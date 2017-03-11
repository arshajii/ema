#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "barcodes.h"
#include "preprocess.h"
#include "samdict.h"
#include "samrecord.h"
#include "main.h"
#include "align.h"

static double mate_dist_penalty(const int64_t mate1_pos, const int64_t mate2_pos)
{
	int64_t d = llabs(mate1_pos - mate2_pos);

	if (INSERT_MIN <= d && d <= INSERT_MAX) {
		return 0.0;
	} else if (d <= INSERT_CAP) {
		return LONG_INS_PENALTY;
	} else {
		return UNPAIRED_PENALTY;
	}
}

static void init_cloud(Cloud *c)
{
	c->lo = 0;
	c->hi = 0;
	c->cov = 0;
	c->snvs = 0;
	c->indels = 0;
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
	uint8_t m1 = (*r1)->mate;
	uint8_t m2 = (*r2)->mate;

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

	Cloud **cloud_assignments = calloc(n_records, sizeof(*cloud_assignments));

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
	SAMRecord **mate_map = calloc(n_records, sizeof(*mate_map));
	for (size_t i = 0; i < n_records; i++) {
		if (mate_map[i] != NULL)
			continue;

		SAMRecord *rec = records[i];
		int64_t closest_dist = 0xffffffff;
		SAMRecord *mate = NULL;
		size_t mate_idx = -1;

		for (size_t j = 0; j < n_records; j++) {
			SAMRecord *rec2 = records[j];
			if (strcmp(rec->ident, rec2->ident) == 0 && rec->mate != rec2->mate) {
				int64_t p1 = rec->pos;
				int64_t p2 = rec2->pos;
				int64_t dist = llabs(llabs(p1 - p2) - INSERT_AVG);

				if (dist < closest_dist && dist < INSERT_CAP) {
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
	struct mmap_dist *mmap_dists = malloc(n_clouds_2 * sizeof(*mmap_dists));
	struct mmap_dist **mmap_dist_buf = malloc(n_clouds_2 * sizeof(*mmap_dist_buf));
	char *mmap_rec_cache = malloc(n_clouds);
	char *mmap_cloud_cache = malloc(n_clouds);

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
	if (x) {
		for (size_t i = 0; i < n_records; i++) {
			SAMRecord *rec = records[i];
			printf("%s/%d %u %f %p\n", rec->ident, rec->mate, rec->pos, rec->score, cloud_assignments[i]);
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

void find_clouds_and_align(SAMRecord *records, const size_t n_records, const char *out)
{
	FILE *out_file = (out == NULL ? stdout : fopen(out, "w"));

	if (!out_file) {
		IOERROR(out);
	}

	qsort(records, n_records, sizeof(*records), record_cmp);
	//records = remove_dups(records, &num_records);

	size_t nc = 0;
	Cloud *clouds = malloc(MAX_CLOUDS_PER_BC * sizeof(*clouds));
	SAMDict *sd = sam_dict_new();
	SAMRecord *record = &records[0];

	while (record->bc != 0) {
		const bc_t bc = record->bc;

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
				SAMRecord **cloud_to_split = malloc(cov * sizeof(*cloud_to_split));

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
			SAMRecord **records = e->cand_records;
			Cloud **clouds = e->cand_clouds;
			double *gammas = e->gammas;
			const size_t num_cands = e->num_cands;

			for (size_t i = 0; i < num_cands; i++) {
				Cloud *c = clouds[i];
				c->snvs += records[i]->snvs;
				c->indels += records[i]->indels;
			}

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

				SAMRecord **mate_records;
				Cloud **mate_clouds;
				double *mate_gammas;
				size_t mate_num_cands;

				if (mate != NULL) {
					mate_records = mate->cand_records;
					mate_clouds = mate->cand_clouds;
					mate_gammas = mate->gammas;
					mate_num_cands = mate->num_cands;
				}

				for (size_t i = 0; i < num_cands; i++) {
					double best_mate_score = UNPAIRED_PENALTY;

					if (mate != NULL) {
						for (size_t j = 0; j < mate_num_cands; j++) {
							if (mate_records[j]->chrom == records[i]->chrom &&
							    mate_clouds[j] == clouds[i]) {
								const double mate_score =
								  mate_dist_penalty(records[i]->pos, mate_records[j]->pos) + log(mate_gammas[j]);

								if (mate_score > best_mate_score) {
									best_mate_score = mate_score;
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
			SAMRecord *best = find_best_record(e);
			fprintf(out_file, "%s %d %s %u\n", best->ident, best->mate, chrom_lookup(best->chrom), best->pos);
			SAMDictEnt *t = e;
			e = e->link_next;
			free(t);
		}

		nc = 0;
		sam_dict_clear(sd);
	}

	free(clouds);
	free(sd);
}

