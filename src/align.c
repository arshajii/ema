#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "barcodes.h"
#include "preprocess.h"
#include "samdict.h"
#include "samrecord.h"
#include "bwabridge.h"
#include "main.h"
#include "align.h"

// mate 1 is reversed
static double mate_dist_penalty(const int64_t mate1_pos, const int64_t mate2_pos)
{
	const int64_t d = mate1_pos - mate2_pos;

	if (INSERT_MIN <= d && d <= INSERT_MAX) {
		return 0.0;
	} else {
		return UNPAIRED_PENALTY;
	}
}

static void init_cloud(Cloud *c)
{
	c->lo = 0;
	c->hi = 0;
	c->cov = 0;
	c->exp_cov = 0.0;
	c->parent = NULL;
	c->child = NULL;
}

static void rec2cloud(SAMRecord *r, Cloud *c)
{
	const uint32_t pos = r->pos;
	if (pos > c->hi)
		c->hi = pos;
	if (pos < c->lo)
		c->lo = pos;
}

struct mmap_dist {
	SAMRecord *rec;
	Cloud *cloud;
	size_t rec_idx;
	size_t cloud_idx;
	size_t dist;
};

/* mate-aware read-to-cloud distance */
static size_t dist(SAMRecord *r, SAMRecord *mate, Cloud *c)
{
	const int64_t a = c->lo;
	const int64_t b = c->hi + MATE2_LEN;
	const int64_t x = mate ? ((r->pos + mate->pos)/2) : r->pos;
	return llabs(x - (a+b)/2);
}

/* distance comparator for `struct mmap_dist`s */
static int dist_cmp(const void *v1, const void *v2)
{
	size_t d1 = (*((struct mmap_dist **)v1))->dist;
	size_t d2 = (*((struct mmap_dist **)v2))->dist;
	return (d1 > d2) - (d1 < d2);
}

/* name comparator for `SAMRecords`s */
static int name_cmp(const void *v1, const void *v2)
{
	SAMRecord **r1 = (SAMRecord **)v1;
	SAMRecord **r2 = (SAMRecord **)v2;
	uint32_t m1 = (*r1)->mate;
	uint32_t m2 = (*r2)->mate;

	const int cmp1 = (m1 > m2) - (m1 < m2);

	if (cmp1 == 0)
		return strcmp((*r1)->ident, (*r2)->ident);
	else
		return cmp1;
}

/* caution: `records` should be name-sorted */
static Cloud **split_cloud(SAMRecord **records, const size_t n_records, Cloud *clouds, size_t *n_split_clouds)
{
#define K 5

	struct multimapped_read {
		size_t idx;  // index in records
		int n;       // how many mappings in the to-be-split cloud
	} mmaps[MAX_CLOUDS_PER_BC];
	size_t n_mmaps = 0;

	Cloud **cloud_assignments = safe_calloc(n_records, sizeof(*cloud_assignments));

#if DEBUG
	int found_target_read = 0;
#endif

	/* find the multi-mapped reads, record highest */
	size_t i = 0;
	size_t n_clouds = 1;
	size_t mmap_max_idx = 0;
	while (i < n_records) {

#if DEBUG
		if (strcmp(records[i]->ident, TRACK_READ) == 0)
			found_target_read = 1;
#endif

		size_t j = i+1;
		while (j < n_records && record_eq(records[j], records[i])) {
			++j;
		}

		size_t n = j-i;

		if (n > 1) {
			mmaps[n_mmaps].idx = i;
			mmaps[n_mmaps].n = n;

			if (n > n_clouds) {
				n_clouds = n;
				mmap_max_idx = n_mmaps;
			}

			++n_mmaps;
		}

		i = j;
	}

	assert(n_clouds > 1);

	/* bail-out if too many splits
	     --or--
	   too few multi-mappings to get reliable splits */
	if (n_clouds > 10 || n_mmaps < 5) {
		*n_split_clouds = 1;

		Cloud *c = &clouds[0];
		init_cloud(c);
		c->lo = 0xffffffff;
		c->hi = 0x00000000;
		c->cov = n_records;

		for (size_t i = 0; i < n_records; i++) {
			rec2cloud(records[i], c);
			cloud_assignments[i] = c;
		}

		c->cov = n_records;
		return cloud_assignments;
	}

	*n_split_clouds = n_clouds;

	/* match up most likely mates */
	SAMRecord **mate_map = safe_calloc(n_records, sizeof(*mate_map));
	for (size_t i = 0; i < n_records; i++) {
		if (mate_map[i] != NULL)
			continue;

		SAMRecord *rec = records[i];
		int64_t closest_dist = 0xffffffff;
		SAMRecord *mate = NULL;
		size_t mate_idx = -1;

		for (size_t j = 0; j < n_records; j++) {
			SAMRecord *rec2 = records[j];
			if (strcmp(rec->ident, rec2->ident) == 0 && rec->mate != rec2->mate && rec->rev != rec2->rev) {
				int64_t p1 = rec->pos;
				int64_t p2 = rec2->pos;
				int64_t dist = llabs(llabs(p1 - p2) - INSERT_AVG);

				if (dist < closest_dist && dist < INSERT_MAX) {
					mate = rec2;
					closest_dist = dist;
					mate_idx = j;
				}
			}
		}

		if (mate != NULL) {
			mate_map[i] = mate;
			mate_map[mate_idx] = rec;
		}
	}

#if DEBUG
	if (found_target_read) {
		printf("n_clouds: %lu\n", n_clouds);
		printf("n_records: %lu\n", n_records);
		printf("n_mmaps: %lu\n", n_mmaps);
	}
#endif

	const size_t n_clouds_2 = n_clouds*n_clouds;
	struct mmap_dist *mmap_dists = safe_malloc(n_clouds_2 * sizeof(*mmap_dists));
	struct mmap_dist **mmap_dist_buf = safe_malloc(n_clouds_2 * sizeof(*mmap_dist_buf));
	char *mmap_rec_cache = safe_malloc(n_clouds);
	char *mmap_cloud_cache = safe_malloc(n_clouds);

	/* initialize clouds */
	for (size_t i = 0; i < n_clouds; i++) {
		init_cloud(&clouds[i]);
	}

	/* preliminary cloud assignments */
	for (size_t i = 0; i < n_clouds; i++) {
		const size_t idx = mmaps[mmap_max_idx].idx + i;
		Cloud *cloud = &clouds[i];
		cloud_assignments[idx] = cloud;
		cloud->hi = cloud->lo = records[idx]->pos;
		++(cloud->cov);
	}

#if DEBUG
	if (found_target_read) {
		for (size_t i = 0; i < n_clouds; i++) {
			printf("(%u %u) ", clouds[i].lo, clouds[i].hi);
		}
		printf("\n");
	}
#endif

	/* assign each record to its "nearest" cloud */
	for (size_t k = 0; k < K; k++) {

#if DEBUG
		int move_made = 0;
#endif

		size_t next_mmap = 0;
		size_t i = 0;

		while (i < n_records) {
			if (next_mmap < n_mmaps && i == mmaps[next_mmap].idx) {
				i += mmaps[next_mmap++].n;
				continue;
			}

			SAMRecord *rec = records[i];
			size_t closest_cloud = 0;
			size_t closest_cloud_dist = -1;
			for (size_t j = 0; j < n_clouds; j++) {
				const size_t d = dist(rec, mate_map[i], &clouds[j]);
				if (d < closest_cloud_dist) {
					closest_cloud = j;
					closest_cloud_dist = d;
				}
			}

			Cloud *new = &clouds[closest_cloud];

#if DEBUG
			Cloud *old = cloud_assignments[i];

			if (old != new)
				move_made = 1;
#endif

			cloud_assignments[i] = new;
			i++;
		}

		for (size_t i = 0; i < n_mmaps; i++) {
			const size_t base_idx = mmaps[i].idx;
			const size_t n = mmaps[i].n;
			size_t n_dists = 0;

			for (size_t r = 0; r < n; r++) {
				for (size_t c = 0; c < n_clouds; c++) {
					const size_t idx = base_idx + r;
					SAMRecord *rec = records[idx];
					Cloud *cloud = &clouds[c];
					const size_t dists_idx = r*n_clouds + c;

					mmap_dists[dists_idx].rec = rec;
					mmap_dists[dists_idx].cloud = cloud;
					mmap_dists[dists_idx].rec_idx = r;
					mmap_dists[dists_idx].cloud_idx = c;
					mmap_dists[dists_idx].dist = dist(rec, mate_map[idx], cloud);
					mmap_dist_buf[n_dists++] = &mmap_dists[dists_idx];
				}
			}

			qsort(mmap_dist_buf, n_dists, sizeof(*mmap_dist_buf), dist_cmp);

			for (size_t r = 0; r < n; r++) {
				mmap_rec_cache[r] = 0;
			}

			for (size_t c = 0; c < n_clouds; c++) {
				mmap_cloud_cache[c] = 0;
			}

			for (size_t j = 0; j < n_dists; j++) {
				struct mmap_dist *dist = mmap_dist_buf[j];

				const size_t rec_idx = dist->rec_idx;
				const size_t cloud_idx = dist->cloud_idx;

				if (!mmap_rec_cache[rec_idx] && !mmap_cloud_cache[cloud_idx]) {
					const size_t rec_idx_abs = base_idx + rec_idx;
					Cloud *new = dist->cloud;
					cloud_assignments[rec_idx_abs] = new;
					mmap_rec_cache[rec_idx] = 1;
					mmap_cloud_cache[cloud_idx] = 1;
				}
			}
		}

		for (size_t i = 0; i < n_clouds; i++) {
			Cloud *c = &clouds[i];
			c->lo = 0xffffffff;
			c->hi = 0x00000000;
			c->cov = 0;
		}

		for (size_t i = 0; i < n_records; i++) {
			rec2cloud(records[i], cloud_assignments[i]);
			++(cloud_assignments[i]->cov);
		}

#if DEBUG
		if (found_target_read) {
			for (size_t i = 0; i < n_clouds; i++) {
				printf("(%u %u)", clouds[i].lo, clouds[i].hi);
			}
			printf(" %d\n", move_made);
		}
#endif
	}

#if DEBUG
	if (found_target_read) {
		for (size_t i = 0; i < n_records; i++) {
			SAMRecord *rec = records[i];
			printf("%s/%d %u %f %p\n", rec->ident, 1+rec->mate, rec->pos, rec->score, cloud_assignments[i]);
		}
	}
#endif

	free(mmap_dists);
	free(mmap_dist_buf);
	free(mmap_rec_cache);
	free(mmap_cloud_cache);
	free(mate_map);
	return cloud_assignments;

#undef K
}

static void normalize_cloud_probabilities(Cloud *clouds, const size_t nc)
{
	for (size_t i = 0; i < nc; i++) {
		Cloud *c = &clouds[i];

		if (c->parent != NULL)
			continue;

		double total = 0;

		for (Cloud *child = c; child != NULL; child = child->child) {
			total += child->reads_prob;
		}

		for (Cloud *child = c; child != NULL; child = child->child) {
			child->reads_prob /= total;
		}
	}
}

static int read_fastq_rec_bc_group(FILE *in,
                                   FASTQRecord *start,
                                   FASTQRecord **out,
                                   size_t *n_recs,
                                   size_t *out_cap);

static void append_alignments(bwaidx_t *ref,
                              mem_opt_t *opts,
                              FASTQRecord *m1,
                              FASTQRecord *m2,
                              SAMRecord **out,
                              size_t *n_recs,
                              size_t *out_cap);

void find_clouds_and_align(FILE *fq1, FILE *fq2, const char *ref_path, FILE *out_file, const char *rg)
{
	// BWA
	fprintf(stderr, "BWA initialization...\n");
	arena_init();
	bwaidx_t *ref = load_reference(ref_path);
	mem_opt_t *opts = mem_opt_init();
	opts->max_occ = 3000;

	// SAM header
	fprintf(stderr, "Processing reads...\n");
	fprintf(out_file, "@HD\tVN:1.3\tSO:unsorted\n");
	for (int32_t i = 0; i < ref->bns->n_seqs; i++) {
		fprintf(out_file, "@SQ\tSN:%s\tLN:%d\n", ref->bns->anns[i].name, ref->bns->anns[i].len);
	}
	if (rg != NULL)
		fprintf(out_file, "%s\n", rg);
	const char *rg_id = (rg != NULL) ? (strstr(rg, "ID:") + 3) : NULL;  // pre-validated

	size_t nc = 0;
	Cloud *clouds = safe_malloc(MAX_CLOUDS_PER_BC * sizeof(*clouds));
	SAMDict *sd = sam_dict_new();

#define INIT_FASTQ_CAP 5000
#define INIT_ALIGN_CAP 1000000
	size_t n_fq1_recs = 0;
	size_t fq1_recs_cap = INIT_FASTQ_CAP;
	FASTQRecord *fq1_recs = safe_malloc(fq1_recs_cap * sizeof(*fq1_recs));

	size_t n_fq2_recs = 0;
	size_t fq2_recs_cap = INIT_FASTQ_CAP;
	FASTQRecord *fq2_recs = safe_malloc(fq2_recs_cap * sizeof(*fq2_recs));

	size_t n_records = 0;
	size_t recs_cap = INIT_ALIGN_CAP;
	SAMRecord *records = safe_malloc(recs_cap * sizeof(*records));

	FASTQRecord start1;
	FASTQRecord start2;
	start1.id[0] = '\0';
	start2.id[0] = '\0';

	while ((read_fastq_rec_bc_group(fq1, &start1, &fq1_recs, &n_fq1_recs, &fq1_recs_cap) |
	        read_fastq_rec_bc_group(fq2, &start2, &fq2_recs, &n_fq2_recs, &fq2_recs_cap)) == 0) {
		const bc_t bc = fq1_recs[0].bc;

		assert(n_fq1_recs == n_fq2_recs);
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
			c->lo = r->pos;
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

					/* unlink the cloud with a collision */
					if (c->parent != NULL)
						c->parent->child = c->child;

					if (c->child != NULL)
						c->child->parent = c->parent;

					c->child = NULL;
					c->parent = NULL;
				}

				++cov;
			}

			if (collision_detected) {
				SAMRecord **cloud_to_split = safe_malloc(cov * sizeof(*cloud_to_split));

				for (size_t i = 0; i < cov; i++) {
					cloud_to_split[i] = &record[i];
				}

				qsort(cloud_to_split, cov, sizeof(*cloud_to_split), name_cmp);
				size_t n_split_clouds;
				Cloud **assignments = split_cloud(cloud_to_split, cov, &clouds[nc], &n_split_clouds);

				for (size_t i = 0; i < cov; i++) {
					sam_dict_add(sd, cloud_to_split[i], assignments[i], 1);
				}

				free(cloud_to_split);
				free(assignments);

				nc += n_split_clouds;
			} else {
				c->cov = cov;
				c->hi = r->pos + MATE2_LEN;
				++nc;
			}

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
			c->reads_prob = c->exp_cov;
		}

		normalize_cloud_probabilities(clouds, nc);

		/* EM iterations */
		for (int q = 0; q < 7; q++) {
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

				for (size_t i = 0; i < num_cands; i++) {
					int pair_found = 0;
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

								if (!pair_found || mate_score > best_mate_score) {
									best_mate_score = mate_score;
									pair_found = 1;
								}
							}
						}
					}

					gammas[i] = records[i]->score + log(clouds[i]->reads_prob) + best_mate_score;
				}

				normalize_log_probs(gammas, num_cands);
			}

			for (SAMDictEnt *e = sd->head; e != NULL; e = e->link_next) {
				Cloud **clouds = e->cand_clouds;
				double *gammas = e->gammas;
				size_t num_cands = e->num_cands;

				for (size_t i = 0; i < num_cands; i++) {
					clouds[i]->exp_cov += gammas[i];
				}
			}

			/* recompute cloud scores */
			for (size_t i = 0; i < nc; i++) {
				Cloud *c = &clouds[i];
				c->reads_prob = c->exp_cov;
			}

			normalize_cloud_probabilities(clouds, nc);
		}

		/* find best alignments */
		SAMDictEnt *e = sd->head;

		while (e != NULL) {
			if (!e->visited) {
				SAMDictEnt *m = e->mate;
				double best_gamma = 0.0;
				double best_gamma_mate = 0.0;
				SAMRecord *best = find_best_record(e, &best_gamma);
				SAMRecord *best_mate = (m != NULL) ? find_best_record(m, &best_gamma_mate) : NULL;

				e->visited = 1;
				if (m != NULL) {
					m->visited = 1;
					m->mate = NULL;
				}

				print_sam_record(best, best_mate, best_gamma, out_file, rg_id);
				print_sam_record(best_mate, best, best_gamma_mate, out_file, rg_id);
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
	free(fq1_recs);
	free(fq2_recs);
	arena_destroy();
	bwa_idx_destroy(ref);
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
	char sep[64];

	if (!fgets(out->id, sizeof(out->id), in))
		return 1;

	assert(fgets(out->read, sizeof(out->read), in));
	assert(fgets(sep, sizeof(sep), in));
	assert(fgets(out->qual, sizeof(out->qual), in));

	out->bc = get_bc_direct(out);
	return 0;
}

static int read_fastq_rec_bc_group(FILE *in,
                                   FASTQRecord *start,
                                   FASTQRecord **out,
                                   size_t *n_recs,
                                   size_t *out_cap)
{
	*n_recs = 0;

	if (start->id[0] == '\0' && read_fastq_rec(in, start) != 0)
		return 1;

	const bc_t bc = start->bc;
	(*out)[(*n_recs)++] = *start;

	do {
		if (*n_recs == *out_cap) {
			*out_cap = (*out_cap * 3)/2 + 1;
			*out = safe_realloc(*out, *out_cap * sizeof(**out));
		}

		if (read_fastq_rec(in, &(*out)[*n_recs]) != 0) {
			start->id[0] = '\0';
			return 0;
		}

		if ((*out)[*n_recs].bc != bc) {
			*start = (*out)[*n_recs];
			return 0;
		}

		++(*n_recs);
	} while (1);
}

static double score_alignment(SingleReadAlignment *r)
{
	static double LOG_MATCH_SCORE;
	static double LOG_MISMATCH_SCORE;
	static double LOG_INDEL_SCORE;
	static double LOG_CLIP_SCORE;
	static int init = 0;

	if (!init) {
		LOG_MATCH_SCORE = log(1 - ERROR_RATE);
		LOG_MISMATCH_SCORE = log(ERROR_RATE);
		LOG_INDEL_SCORE = log(INDEL_RATE);
		LOG_CLIP_SCORE = log(CLIP_RATE);
		init = 1;
	}

	int matches = 0;
	int mismatches = 0;
	int indels = 0;
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

	return matches*LOG_MATCH_SCORE +
	       mismatches*LOG_MISMATCH_SCORE +
	       indels*LOG_INDEL_SCORE +
	       clipping*LOG_CLIP_SCORE;
}

static void alignment_to_sam_rec(FASTQRecord *fq, FASTQRecord *fq_mate, SingleReadAlignment *r, SAMRecord *s, unsigned mate)
{
	s->bc = fq->bc;
	s->chrom = chrom_index(r->chrom);
	s->pos = r->pos + 1;

	char *c = &fq->id[1];  // skip first '@'
	size_t i;
	for (i = 0; *c != '\n'; i++)
		s->ident[i] = *c++;
	s->ident[i] = '\0';

	s->score = score_alignment(r);

	s->hash = 0;
	s->mate_hash = 0;
	s->hashed = 0;
	s->mate_hashed = 0;
	s->mate = mate;
	s->rev = r->rev;
	s->fq = fq;
	s->fq_mate = fq_mate;
	s->aln = *r;
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

	EasyAlignmentPairs alignments = bwa_mem_mate_sw(ref, opts, m1->read, MATE1_LEN, m2->read, MATE2_LEN, 25);
	int best_dist = -1;

	for (size_t i = 0; i < alignments.len1; i++) {
		EasyAlignment *a = &alignments.a1[i];
		SingleReadAlignment r;
		bwa_smith_waterman(ref, opts, m1->read, MATE1_LEN, a->chained_hit, &r);

		const int dist = r.edit_dist + (MATE1_LEN - (a->read_e - a->read_s));
		if (i == 0)
			best_dist = dist;
		else if (dist - best_dist > EXTRA_SEARCH_DEPTH)
			break;

		REALLOC_IF_NECESSARY();
		alignment_to_sam_rec(m1, m2, &r, &(*out)[(*n_recs)++], 0);
	}

	for (size_t i = 0; i < alignments.len2; i++) {
		EasyAlignment *a = &alignments.a2[i];
		SingleReadAlignment r;
		bwa_smith_waterman(ref, opts, m2->read, MATE2_LEN, a->chained_hit, &r);

		const int dist = r.edit_dist + (MATE2_LEN - (a->read_e - a->read_s));
		if (i == 0)
			best_dist = dist;
		else if (dist - best_dist > EXTRA_SEARCH_DEPTH)
			break;

		REALLOC_IF_NECESSARY();
		alignment_to_sam_rec(m2, m1, &r, &(*out)[(*n_recs)++], 1);
	}

	(*out)[*n_recs].bc = 0;
}

