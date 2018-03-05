#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <omp.h>

#include "barcodes.h"
#include "preprocess.h"
#include "samdict.h"
#include "samrecord.h"
#include "split.h"
#include "bwabridge.h"
#include "main.h"
#include "align.h"

void init_cloud(Cloud *c)
{
	static int cloud_id = 0;
	c->exp_cov = 0.0;
	c->parent = NULL;
	c->child = NULL;
	c->id = cloud_id++;
	c->bad = 0;
}

int is_pair(SAMRecord *r1, SAMRecord *r2)
{
	if (r1->rev == r2->rev || r1->chrom != r2->chrom)
		return 0;

	if (r2->rev) {
		SAMRecord *rt = r2;
		r2 = r1;
		r1 = rt;
	}

	const int64_t d = r1->pos - r2->pos;
	return INSERT_MIN <= d && d <= INSERT_MAX;
}

int is_pair_relaxed(SAMRecord *r1, SAMRecord *r2)
{
	if (r1->chrom != r2->chrom)
		return 0;

	if (r1->pos < r2->pos) {
		SAMRecord *rt = r2;
		r2 = r1;
		r1 = rt;
	}

	const int64_t d = r1->pos - r2->pos;
	return INSERT_MIN <= d && d <= INSERT_MAX;
}


/* mate 1 is reversed */
static double mate_dist_penalty(const int64_t mate1_pos, const int64_t mate2_pos)
{
	const int64_t d = mate1_pos - mate2_pos;

	if (INSERT_MIN <= d && d <= INSERT_MAX) {
		return 0.0;
	} else {
		return UNPAIRED_PENALTY;
	}
}

/* name comparator for `SAMRecords`s */
static int name_cmp(const void *v1, const void *v2)
{
	SAMRecord **r1 = (SAMRecord **)v1;
	SAMRecord **r2 = (SAMRecord **)v2;
	uint32_t m1 = (*r1)->mate;
	uint32_t m2 = (*r2)->mate;

	const int cmp1 = strcmp((*r1)->ident, (*r2)->ident);
	const int cmp2 = (m1 > m2) - (m1 < m2);

	return (cmp1 != 0) ? cmp1 : cmp2;
}

static void normalize_cloud_probabilities(Cloud *clouds, const size_t nc)
{
	for (size_t i = 0; i < nc; i++) {
		Cloud *c = &clouds[i];

		if (c->parent != NULL)
			continue;

		double total = 0.0;

		for (Cloud *child = c; child != NULL; child = child->child) {
			total += child->weight;
		}

		for (Cloud *child = c; child != NULL; child = child->child) {
			child->weight /= total;
		}
	}
}

static void read_special_fastq(FILE *fq, FASTQRecord **fq1_recs, FASTQRecord **fq2_recs);

static int read_fastq_rec_bc_group(FILE *in,
                                   FASTQRecord *start,
                                   FASTQRecord **out,
                                   size_t *n_recs,
                                   size_t *out_cap);

static int seek_next_barcode_group(FASTQRecord *fq_recs_full,
                                   FASTQRecord **fq_recs,
                                   size_t *n_recs,
                                   FASTQRecord *latest_fqr,
                                   size_t latest_n);

static void append_alignments(bwaidx_t *ref,
                              mem_opt_t *opts,
                              FASTQRecord *m1,
                              FASTQRecord *m2,
                              SAMRecord **out,
                              size_t *n_recs,
                              size_t *out_cap);

static bwaidx_t *ref;
static mem_opt_t *opts;

void bwa_init(const char *ref_path)
{
	fprintf(stderr, "BWA initialization...\n");
	ref = load_reference(ref_path);
	opts = mem_opt_init();
	opts->max_occ = 3000;
}

void bwa_dealloc(void)
{
	bwa_idx_destroy(ref);
}

/* caution: caller must ensure (fq1 != NULL && fq2 != NULL) ^ (fqx != NULL) */
void find_clouds_and_align(FILE *fq1, FILE *fq2, FILE *fqx, FILE *out_file, const int apply_opt)
{
#define STANDARD_FASTQ() (fqx == NULL)

	omp_lock_t in_lock;
	omp_lock_t out_lock;

	omp_init_lock(&in_lock);
	omp_init_lock(&out_lock);

	fprintf(stderr, "Processing reads...\n");

	/* SAM header */
	// HD
	fprintf(out_file, "@HD\tVN:1.3\tSO:unsorted\n");

	// SQ
	for (int32_t i = 0; i < ref->bns->n_seqs; i++) {
		fprintf(out_file, "@SQ\tSN:%s\tLN:%d\n", ref->bns->anns[i].name, ref->bns->anns[i].len);
	}

	// RG
	if (rg != NULL)
		fprintf(out_file, "%s\n", rg);
	const char *rg_id = (rg != NULL) ? (strstr(rg, "ID:") + 3) : NULL;  // pre-validated

	// PG
	fprintf(out_file, "@PG\tID:ema\tPN:ema\tVN:%s\tCL:%s", VERSION, pg_argv[0]);
	for (int i = 1; i < pg_argc; i++)
		fprintf(out_file, " %s", pg_argv[i]);
	fprintf(out_file, "\n");

	/* for our special FASTQs */
	FASTQRecord *fq1_recs_full = NULL;
	FASTQRecord *fq2_recs_full = NULL;
	// for coordinating the distribution of reads between threads:
	int done = 0;
	FASTQRecord *latest_fqr1 = NULL;
	size_t latest_n1 = 0;
	FASTQRecord *latest_fqr2 = NULL;
	size_t latest_n2 = 0;

	if (!STANDARD_FASTQ())
		read_special_fastq(fqx, &fq1_recs_full, &fq2_recs_full);

	#pragma omp parallel num_threads(NUM_THREADS)
	{
		arena_init();
		size_t nc = 0;
		Cloud *clouds = safe_malloc(MAX_CLOUDS_PER_BC * sizeof(*clouds));
		SAMDict *sd = sam_dict_new();

#define INIT_FASTQ_CAP 5000
#define INIT_ALIGN_CAP 1000000

		/*
		 * FASTQ I/O is broken into two cases:
		 *   (i) Standard paired-end bc-sorted FASTQs:
		 *         We read each barcode group as we need to from the file
		 *  (ii) Our special FASTQ format (described below):
		 *         We read the entire file, sort it, then partition it by
		 *         barcode as needed.
		 */

		FASTQRecord *fq1_recs = NULL;       // for standard FASTQs
		FASTQRecord *fq2_recs = NULL;       // for standard FASTQs

		size_t n_fq1_recs = 0;
		size_t fq1_recs_cap = INIT_FASTQ_CAP;

		size_t n_fq2_recs = 0;
		size_t fq2_recs_cap = INIT_FASTQ_CAP;

		if (STANDARD_FASTQ()) {
			fq1_recs = safe_malloc(fq1_recs_cap * sizeof(*fq1_recs));
			fq2_recs = safe_malloc(fq2_recs_cap * sizeof(*fq2_recs));
		}

		size_t n_records = 0;
		size_t recs_cap = INIT_ALIGN_CAP;
		SAMRecord *records = safe_malloc(recs_cap * sizeof(*records));

		FASTQRecord start1;
		FASTQRecord start2;
		start1.id[0] = '\0';
		start2.id[0] = '\0';

		while (1) {
			int good;

			{
				omp_set_lock(&in_lock);
				if (STANDARD_FASTQ()) {
					good = (read_fastq_rec_bc_group(fq1, &start1, &fq1_recs, &n_fq1_recs, &fq1_recs_cap) |
					        read_fastq_rec_bc_group(fq2, &start2, &fq2_recs, &n_fq2_recs, &fq2_recs_cap)) == 0;
				} else {
					if (!done) {
						good = (seek_next_barcode_group(fq1_recs_full, &fq1_recs, &n_fq1_recs, latest_fqr1, latest_n1) |
						        seek_next_barcode_group(fq2_recs_full, &fq2_recs, &n_fq2_recs, latest_fqr2, latest_n2)) == 0;
					} else {
						good = 0;
					}

					if (good) {
						latest_fqr1 = fq1_recs;
						latest_n1 = n_fq1_recs;

						latest_fqr2 = fq2_recs;
						latest_n2 = n_fq2_recs;
					} else {
						done = 1;
					}
				}
				omp_unset_lock(&in_lock);
			}

			if (!good)
				break;

			const bc_t bc = fq1_recs[0].bc;
			assert(n_fq1_recs == n_fq2_recs);
			const int worth_doing_full_em = (n_fq1_recs >= 30);

			for (size_t i = 0; i < n_fq1_recs; i++) {
				append_alignments(ref, opts, &fq1_recs[i], &fq2_recs[i], &records, &n_records, &recs_cap);
			}

			qsort(records, n_records, sizeof(*records), record_cmp);
			SAMRecord *record = &records[0];

			/* find and process clouds */
			while (record->bc == bc) {
				SAMRecord *r = record;
				Cloud *c = &clouds[nc];

				init_cloud(c);
				sam_dict_add(sd, r, c, 0);
				size_t cov = 1;

				int collision_detected = 0;

				while ((r+1)->bc == r->bc &&
				       (r+1)->chrom == r->chrom &&
				       (r+1)->pos - r->pos <= DIST_THRESH) {

					++r;

					if (!collision_detected && sam_dict_add(sd, r, c, 0)) {
						collision_detected = 1;

						for (size_t i = 0; i < cov; i++) {
							sam_dict_del(sd, &record[i]);
						}
					}

					++cov;
				}

				if (collision_detected) {
					c->bad = 1;
					SAMRecord **cloud_to_split = safe_malloc(cov * sizeof(*cloud_to_split));

					for (size_t i = 0; i < cov; i++) {
						cloud_to_split[i] = &record[i];  // note: `record`; NOT `records`
					}

					qsort(cloud_to_split, cov, sizeof(*cloud_to_split), name_cmp);

					if (apply_opt)
						mark_optimal_alignments_in_cloud(cloud_to_split, cov);

					for (size_t i = 0; i < cov; i++) {
						sam_dict_add(sd, cloud_to_split[i], &clouds[nc], 1);
					}

					free(cloud_to_split);
				}

				++nc;
				record = r+1;
			}

			/* initializations */
			for (SAMDictEnt *e = sd->head; e != NULL; e = e->link_next) {
				Cloud **clouds = e->cand_clouds;
				double *gammas = e->gammas;
				const size_t num_cands = e->num_cands;

				normalize_log_probs(gammas, num_cands);

				for (size_t i = 0; i < num_cands; i++) {
					clouds[i]->exp_cov += gammas[i];
				}
			}

			for (size_t i = 0; i < nc; i++) {
				Cloud *c = &clouds[i];
				c->weight = c->exp_cov;
			}
			normalize_cloud_probabilities(clouds, nc);

			/* EM iterations */
			for (int q = 0; q < EM_ITERS; q++) {
				if (!worth_doing_full_em)
					break;

				double old_gammas[MAX_CANDIDATES];
				double max_change = 0;

				for (size_t i = 0; i < nc; i++) {
					clouds[i].exp_cov = 0.0;
				}

				/* recompute gammas */
				for (SAMDictEnt *e = sd->head; e != NULL; e = e->link_next) {
					SAMDictEnt *mate = e->mate;
					SAMRecord **records = e->cand_records;
					Cloud **clouds = e->cand_clouds;
					double *gammas = e->gammas;
					size_t num_cands = e->num_cands;

					SAMRecord **mate_records = NULL;
					Cloud **mate_clouds = NULL;
					double *mate_gammas = NULL;
					size_t mate_num_cands = 0;

					if (mate != NULL) {
						mate_records = mate->cand_records;
						mate_clouds = mate->cand_clouds;
						mate_gammas = mate->gammas;
						mate_num_cands = mate->num_cands;
					}

					int max_old_gamma_index = 0;
					int max_gamma_index = 0;

					for (size_t i = 0; i < num_cands; i++) {
						double best_mate_score = UNPAIRED_PENALTY;

						if (mate != NULL) {
							for (size_t j = 0; j < mate_num_cands; j++) {
								if (mate_records[j]->chrom == records[i]->chrom &&
								    mate_records[j]->rev != records[i]->rev &&
								    mate_clouds[j] == clouds[i] &&
								    mate_gammas[j] != 0.0) {

									const double penalty = records[i]->rev ? mate_dist_penalty(records[i]->pos, mate_records[j]->pos) :
									                                         mate_dist_penalty(mate_records[j]->pos, records[i]->pos);
									const double mate_score = penalty + log(mate_gammas[j]);

									if (mate_score > best_mate_score) {
										best_mate_score = mate_score;
									}
								}
							}
						}

						old_gammas[i] = gammas[i];

						if (old_gammas[i] > old_gammas[max_old_gamma_index])
							max_old_gamma_index = i;

						gammas[i] = records[i]->score + log(clouds[i]->weight) + best_mate_score;
					}

					normalize_log_probs(gammas, num_cands);

					for (size_t i = 0; i < num_cands; i++) {
						const double change  = fabs(old_gammas[i] - gammas[i]);
						if (change > max_change)
							max_change = change;
						if (gammas[i] > gammas[max_gamma_index])
							max_gamma_index = i;
					}
				}

				for (SAMDictEnt *e = sd->head; e != NULL; e = e->link_next) {
					SAMRecord **records = e->cand_records;
					Cloud **clouds = e->cand_clouds;
					double *gammas = e->gammas;
					size_t num_cands = e->num_cands;

					for (size_t i = 0; i < num_cands; i++) {
						if (records[i]->active && !records[i]->duplicate)
							clouds[i]->exp_cov += gammas[i];
					}
				}

				/* recompute cloud scores */
				for (size_t i = 0; i < nc; i++) {
					Cloud *c = &clouds[i];
					c->weight = c->exp_cov;
				}
				normalize_cloud_probabilities(clouds, nc);
			}

			/* find best alignments */
			SAMDictEnt *e = sd->head;

			while (e != NULL) {
				if (!e->visited) {
					struct xa alts1[MAX_ALTS];
					struct xa alts2[MAX_ALTS];
					size_t n_alts1 = 0;
					size_t n_alts2 = 0;
					SAMDictEnt *m = e->mate;
					double best_gamma = 0.0;
					double best_gamma_mate = 0.0;
					Cloud *cloud1 = NULL;
					Cloud *cloud2 = NULL;
					SAMRecord *best = find_best_record(e, &best_gamma, &cloud1, alts1, &n_alts1);
					SAMRecord *best_mate = (m != NULL) ? find_best_record(m, &best_gamma_mate, &cloud2, alts2, &n_alts2) : NULL;

					e->visited = 1;
					if (m != NULL) {
						m->visited = 1;
						m->mate = NULL;
					}

					{
						omp_set_lock(&out_lock);
						print_sam_record(best, best_mate, best_gamma, cloud1, out_file, rg_id, alts1, n_alts1);
						print_sam_record(best_mate, best, best_gamma_mate, cloud2, out_file, rg_id, alts2, n_alts2);
						omp_unset_lock(&out_lock);
					}
				}

				SAMDictEnt *t = e;
				e = e->link_next;
				free(t);
			}

			nc = 0;
			n_records = 0;
			sam_dict_clear(sd);
			arena_clear();
		}

		free(clouds);
		free(sd);
		free(records);

		if (STANDARD_FASTQ()) {
			free(fq1_recs);
			free(fq2_recs);
		}
	}

	if (!STANDARD_FASTQ()) {
		free(fq1_recs_full);
		free(fq2_recs_full);
	}

	arena_destroy();
	omp_destroy_lock(&in_lock);
	omp_destroy_lock(&out_lock);
}

static bc_t get_bc_direct(FASTQRecord *rec)
{
	char *bc_str = strrchr(rec->id, ':');
	assert(bc_str != NULL);
	*bc_str = '\0';
	return encode_bc(&bc_str[1]);
}

static int read_fastq_rec(FILE *in, FASTQRecord *out)
{
	char sep[256];

	if (!fgets(out->id, sizeof(out->id), in))
		return 1;

	assert(fgets(out->read, sizeof(out->read), in));
	assert(fgets(sep, sizeof(sep), in));
	assert(fgets(out->qual, sizeof(out->qual), in));

	out->bc = get_bc_direct(out);
	out->rlen = strlen(out->read) - 1;  // -1 because of newline
	return 0;
}

static int read_fastq_rec_bc_group(FILE *in,
                                   FASTQRecord *start,
                                   FASTQRecord **out,
                                   size_t *n_recs,
                                   size_t *out_cap)
{
	*n_recs = 0;

	if (IS_SENTINEL(*start) && read_fastq_rec(in, start) != 0)
		return 1;

	const bc_t bc = start->bc;
	(*out)[(*n_recs)++] = *start;

	do {
		if (*n_recs == *out_cap) {
			*out_cap = (*out_cap * 3)/2 + 1;
			*out = safe_realloc(*out, *out_cap * sizeof(**out));
		}

		if (read_fastq_rec(in, &(*out)[*n_recs]) != 0) {
			SET_SENTINEL(*start);
			return 0;
		}

		if ((*out)[*n_recs].bc != bc) {
			*start = (*out)[*n_recs];
			return 0;
		}

		++(*n_recs);
	} while (1);
}

/*
 * Special FASTQ file format:
 *   barcode read1 id qual1 read2 qual2
 *
 * `read1` and `qual1` have the barcode bases/quals trimmed.
 */
static int special_fastq_record_cmp(const void *v1, const void *v2)
{
	const char *c1 = *(const char **)v1;
	const char *c2 = *(const char **)v2;
	return strncmp(c1, c2, BC_LEN);
}

static void read_special_fastq(FILE *fq, FASTQRecord **fq1_recs, FASTQRecord **fq2_recs)
{
	char bc_str[BC_LEN + 1];
	char buf[5000];
	const size_t n_records_actual = count_lines(fq) + 1;
	char **records = safe_malloc(n_records_actual * sizeof(*records));
	size_t n_records = 0;

	while (fgets(buf, sizeof(buf), fq)) {
		char *record = safe_malloc(strlen(buf) + 1);
		strcpy(record, buf);
		records[n_records++] = record;
	}

	qsort(records, n_records, sizeof(*records), special_fastq_record_cmp);

	*fq1_recs = safe_malloc((n_records + 1) * sizeof(**fq1_recs));
	*fq2_recs = safe_malloc((n_records + 1) * sizeof(**fq2_recs));

	for (size_t i = 0; i < n_records; i++) {
		char *p = records[i];
		copy_until_space(bc_str, &p);
		const bc_t bc = encode_bc(bc_str);

		FASTQRecord *fq1 = &(*fq1_recs)[i];
		FASTQRecord *fq2 = &(*fq2_recs)[i];

		fq1->bc = bc;
		fq2->bc = bc;

		copy_until_space(fq1->id, &p);
		strcpy(fq2->id, fq1->id);
		copy_until_space(fq1->read, &p);
		copy_until_space(fq1->qual, &p);
		copy_until_space(fq2->read, &p);
		copy_until_space(fq2->qual, &p);

		free(records[i]);
	}

	SET_SENTINEL((*fq1_recs)[n_records]);
	SET_SENTINEL((*fq2_recs)[n_records]);

	free(records);
}

static int seek_next_barcode_group(FASTQRecord *fq_recs_full,
                                   FASTQRecord **fq_recs,
                                   size_t *n_recs,
                                   FASTQRecord *latest_fqr,
                                   size_t latest_n)
{
	if (latest_fqr != NULL) {
		*fq_recs = latest_fqr;
		*n_recs = latest_n;
	}

	if (*fq_recs == NULL) {
		*fq_recs = fq_recs_full;
	} else {
		*fq_recs = &(*fq_recs)[*n_recs];
	}

	if (IS_SENTINEL(**fq_recs)) {
		*n_recs = 0;
		return 1;
	}

	const bc_t bc = (*fq_recs)->bc;
	FASTQRecord *p = *fq_recs;
	*n_recs = 0;

	while (!IS_SENTINEL(*p) && p->bc == bc) {
		p++;
		(*n_recs)++;
	}

	latest_fqr = *fq_recs;
	latest_n = *n_recs;

	return 0;
}


static void score_alignment(SingleReadAlignment *r, SAMRecord *s)
{
	static double LOG_MATCH_SCORE;
	static double LOG_MISMATCH_SCORE;
	static double LOG_INDEL_SCORE;
	static double LOG_CLIP_SCORE;
	static double LOG10_MISMATCH_SCORE;
	static double LOG10_INDEL_SCORE;
	static double LOG10_CLIP_SCORE;
	static int init = 0;

	if (!init) {
		LOG_MATCH_SCORE = log(1 - ERROR_RATE);
		LOG_MISMATCH_SCORE = log(ERROR_RATE);
		LOG_INDEL_SCORE = log(INDEL_RATE);
		LOG_CLIP_SCORE = log(CLIP_RATE);
		LOG10_MISMATCH_SCORE = log10(ERROR_RATE);
		LOG10_INDEL_SCORE = log10(INDEL_RATE);
		LOG10_CLIP_SCORE = log10(CLIP_RATE);
		init = 1;
	}

	int matches = 0;
	int mismatches = 0;
	int indels = 0;
	int indels_to_count = 0;
	int clipping = 0;

	const uint32_t *cigar = r->cigar;
	const int cigar_len = r->n_cigar;
	for (int i = 0; i < cigar_len; i++) {
		const uint32_t op   = cigar[i];
		const uint32_t type = op & 0xf;
		const uint32_t n    = op >> 4;

		switch (type) {
		case 0:  // M
			matches += n;
			break;
		case 1:  // I
		case 2:  // D
			indels += n;
			++indels_to_count;
			break;
		case 3:  // S
		case 4:  // H
			clipping += n;
			break;
		default:
			assert(0);
		}
	}

	mismatches = r->edit_dist - indels;
	matches -= mismatches;

	s->mapq = r->mapq;

	s->score = matches*LOG_MATCH_SCORE +
	           mismatches*LOG_MISMATCH_SCORE +
	           indels_to_count*LOG_INDEL_SCORE +
	           clipping*LOG_CLIP_SCORE;

	s->score_mapq = (int)(60.0 +
	                      mismatches*LOG10_MISMATCH_SCORE +
	                      indels_to_count*LOG10_INDEL_SCORE +
	                      clipping*LOG10_CLIP_SCORE);
}

static void alignment_to_sam_rec(FASTQRecord *fq,
                                 FASTQRecord *fq_mate,
                                 SingleReadAlignment *r,
                                 SAMRecord *s,
                                 unsigned mate,
                                 int clip,
                                 int clip_edit_dist)
{
	s->bc = fq->bc;
	s->chrom = chrom_index(r->chrom);
	s->pos = r->pos + 1;

	char *c = &fq->id[1];  // skip first '@'
	size_t i;
	for (i = 0; *c && *c != '\n'; i++)
		s->ident[i] = *c++;
	s->ident[i] = '\0';

	score_alignment(r, s);

	s->clip = clip;
	s->clip_edit_dist = clip_edit_dist;
	s->hash = 0;
	s->mate_hash = 0;
	s->hashed = 0;
	s->mate_hashed = 0;
	s->unique = 0;
	s->duplicate = 0;
	s->active = 1;
	s->mate = mate;
	s->rev = r->rev;
	s->fq = fq;
	s->fq_mate = fq_mate;
	s->aln = *r;
}

// adapted from BWA
static int mem_approx_mapq_se_insist(const mem_opt_t *opt, const mem_alnreg_t *a)
{
	int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
	double identity;
	sub = a->csub > sub? a->csub : sub;
	if (sub >= a->score) return 0;
	l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
	identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
	if (a->score == 0) {
		mapq = 0;
	} else if (opt->mapQ_coef_len > 0) {
		double tmp;
		tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
		tmp *= identity * identity;
		mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
	} else {
		mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
		mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	}
	if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
	//if (mapq > 60) mapq = 60;
	if (mapq > 254) mapq = 254;
	if (mapq < 0) mapq = 0;
	mapq = (int)(mapq * (1. - a->frac_rep) + .499);
	return mapq;
}

static void append_alignments(bwaidx_t *ref,
                              mem_opt_t *opts,
                              FASTQRecord *m1,
                              FASTQRecord *m2,
                              SAMRecord **out,
                              size_t *n_recs,
                              size_t *out_cap)
{
#define REALLOC_IF_NECESSARY() \
	do { \
		if (*n_recs + 1 == *out_cap) { \
			*out_cap = (*out_cap * 3)/2 + 1; \
			*out = safe_realloc(*out, *out_cap * sizeof(**out)); \
		} \
	} while (0)

	const int mate1_len = m1->rlen;
	const int mate2_len = m2->rlen;

	EasyAlignmentPairs alignments = bwa_mem_mate_sw(ref, opts, m1->read, mate1_len, m2->read, mate2_len, 25);
	int best_dist = -1;
	size_t aligns_added1 = 0;
	size_t aligns_added2 = 0;

	for (size_t i = 0; i < alignments.len1; i++) {
		EasyAlignment *a = &alignments.a1[i];
		SingleReadAlignment r;
		bwa_smith_waterman(ref, opts, m1->read, mate1_len, a->chained_hit, &r);

		const int clip = (mate1_len - (a->read_e - a->read_s));

		if (clip >= mate1_len/2)
			break;

		const int dist = r.edit_dist + clip;
		if (i == 0)
			best_dist = dist;
		else if (dist - best_dist > EXTRA_SEARCH_DEPTH)
			break;

		REALLOC_IF_NECESSARY();
		r.mapq = mem_approx_mapq_se_insist(opts, a->chained_hit);
		alignment_to_sam_rec(m1, m2, &r, &(*out)[(*n_recs)++], 0, clip, dist);
		++aligns_added1;
	}

	if (aligns_added1 == 1)
		(*out)[*n_recs - 1].unique = 1;

	for (size_t i = 0; i < alignments.len2; i++) {
		EasyAlignment *a = &alignments.a2[i];
		SingleReadAlignment r;
		bwa_smith_waterman(ref, opts, m2->read, mate2_len, a->chained_hit, &r);

		const int clip = (mate2_len - (a->read_e - a->read_s));

		if (clip >= mate2_len/2)
			break;

		const int dist = r.edit_dist + clip;
		if (i == 0)
			best_dist = dist;
		else if (dist - best_dist > EXTRA_SEARCH_DEPTH)
			break;

		REALLOC_IF_NECESSARY();
		r.mapq = mem_approx_mapq_se_insist(opts, a->chained_hit);
		alignment_to_sam_rec(m2, m1, &r, &(*out)[(*n_recs)++], 1, clip, dist);
		++aligns_added2;
	}

	if (aligns_added2 == 1)
		(*out)[*n_recs - 1].unique = 1;

	(*out)[*n_recs].bc = 0;
}

