#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "samdict.h"
#include "samrecord.h"
#include "main.h"
#include "align.h"
#include "util.h"
#include "split.h"

static double log_density_prob(unsigned int density)
{
	static double density_probs_log[SIZE(tech->density_probs)];
	static int init = 0;

	const size_t size = tech->n_density_probs;

	if (!init) {
		for (size_t i = 0; i < size; i++) {
			const double p = tech->density_probs[i];
			density_probs_log[i] = log(p);
		}

		init = 1;
	}

	if (density < size)
		return density_probs_log[density];
	else
		return density_probs_log[size-1] - (density - size + 1)*log(2.0);  // scale down larger density probs exponentially
}

/* caution: `records` should be name-sorted */
void mark_optimal_alignments_in_cloud(SAMRecord **records, size_t n_records)
{
#define MAX_NO_MOVE 500
#define BUF_SIZE    50000

#define BIN_IDX_FOR_POS(pos, lo)   (((pos) - (lo))/BIN_SIZE)
#define BIN_FOR_POS(bins, pos, lo) ((bins)[BIN_IDX_FOR_POS((pos), (lo))])

#define UPDATE_CLOUD_BOUNDARIES(r, lo, hi) \
	do {                                   \
		if ((r)->pos < lo)                 \
			lo = (r)->pos;                 \
		if ((r)->pos > hi)                 \
			hi = (r)->pos;                 \
	} while (0)

	static int rand_init = 0;
	if (!rand_init) {
		time_t t;
		srand((unsigned)time(&t));
		rand_init = 1;
	}

	double log_config_prob = 0;
	unsigned short bins[MAX_BINS] = {0};

	struct uniquemapped_read {
		size_t idx;
	} umaps[BUF_SIZE];

	struct multimapped_read {
		size_t idx;  // index in records
		int n;       // how many mappings
		int mate_umap;  // mate
		int mate_mmap;  // mate
		int active;  // which of these alignments is active
	} mmaps[BUF_SIZE];

	size_t n_umaps = 0;
	size_t n_mmaps = 0;

	uint32_t cloud_lo = 0xffffffff;
	uint32_t cloud_hi = 0x00000000;

	if (n_records >= BUF_SIZE || n_records <= 5)
		return;

	/* remove records that are too far from the lowest edit distance */
	SAMRecord **records_clean = safe_malloc(n_records * sizeof(*records_clean));
	size_t n_records_clean = 0;

	for (size_t i = 0; i < n_records;) {
		size_t j = i+1;

		while (j < n_records && record_eq(records[j], records[i])) {
			++j;
		}

		const size_t n = j-i;

		if (n > 1) {
			size_t min_edit_dist = 0;
			for (size_t k = 0; k < n; k++) {
				if (records[i+k]->clip_edit_dist < records[i+min_edit_dist]->clip_edit_dist)
					min_edit_dist = k;
			}

			SAMRecord *best = records[i+min_edit_dist];
			const int too_much_clipping = 0;//(best->clip > SPLIT_CLIP_THRESH);

			if (!too_much_clipping) {
				const int edit_dist_cutoff = best->clip_edit_dist + SPLIT_EXTRA_SEARCH_DEPTH;

				for (size_t k = 0; k < n; k++) {
					if (records[i+k]->clip_edit_dist <= edit_dist_cutoff)
						records_clean[n_records_clean++] = records[i+k];
					else
						records[i+k]->active = 0;
				}
			} else {
				// don't add the records...
				for (size_t k = 0; k < n; k++) {
					records[i+k]->active = 0;
				}
				best->active = 1;
			}
		} else {
			//if (records[i]->clip <= SPLIT_CLIP_THRESH)
				records_clean[n_records_clean++] = records[i];
		}

		i = j;
	}

	records = records_clean;
	n_records = n_records_clean;

	/* find the multi-mapped reads, record highest */
	for (size_t i = 0; i < n_records;) {
		UPDATE_CLOUD_BOUNDARIES(records[i], cloud_lo, cloud_hi);
		size_t j = i+1;

		while (j < n_records && record_eq(records[j], records[i])) {
			UPDATE_CLOUD_BOUNDARIES(records[j], cloud_lo, cloud_hi);
			++j;
		}

		const size_t n = j-i;

		if (n > 1) {
			size_t max_score = 0;
			for (size_t k = 0; k < n; k++) {
				if (records[i+k]->score > records[i+max_score]->score)
					max_score = k;
			}

			// find mate
			int mate_umap = -1;
			int mate_mmap = -1;

			for (size_t k = 0; k < n_umaps; k++) {
				if (record_eq_mate(records[i], records[umaps[k].idx])) {
					mate_umap = k;
					break;
				}
			}

			if (mate_umap < 0) {
				for (size_t k = 0; k < n_mmaps; k++) {
					if (record_eq_mate(records[i], records[mmaps[k].idx])) {
						mate_mmap = k;
						mmaps[k].mate_mmap = n_mmaps;  // point mate to us
						break;
					}
				}
			}

			mmaps[n_mmaps++] = (struct multimapped_read){.idx = i,
			                                             .n = n,
			                                             .mate_umap = mate_umap,
			                                             .mate_mmap = mate_mmap,
			                                             .active = max_score};

			log_config_prob += records[i+max_score]->score/SCORE_SCALE;  // add up active alignment scores
		} else {
			for (size_t k = 0; k < n_mmaps; k++) {
				if (record_eq_mate(records[i], records[mmaps[k].idx])) {
					mmaps[k].mate_umap = n_umaps;  // point mate to us
					break;
				}
			}

			umaps[n_umaps++] = (struct uniquemapped_read){.idx = i};
			log_config_prob += records[i]->score/SCORE_SCALE;  // add up active alignment scores
		}

		i = j;
	}

	/* get initial configuration probability */
	const size_t n_bins = (cloud_hi - cloud_lo)/BIN_SIZE + 1;

	if (n_bins >= MAX_BINS || n_records <= 5 || n_mmaps == 0) {
		free(records);
		return;
	}

	for (size_t i = 0; i < n_records; i++) {
		records[i]->active = 0;  // we will re-set the active ones later
	}

	for (size_t i = 0; i < n_umaps; i++) {
		SAMRecord *rec = records[umaps[i].idx];
		++BIN_FOR_POS(bins, rec->pos, cloud_lo);
	}

	for (size_t i = 0; i < n_mmaps; i++) {
		SAMRecord *rec = records[mmaps[i].idx + mmaps[i].active];
		++BIN_FOR_POS(bins, rec->pos, cloud_lo);
	}

	for (size_t i = 0; i < n_bins; i++) {
		log_config_prob += log_density_prob(bins[i]);
	}

	/* simulated annealing to minimize distance variance */
	int no_move_count = 0;
	for (size_t k = 0; k < SIM_ANNEAL_ITERS; k++) {
		const double t = pow(10.0, TMAX_LOG - ((TMAX_LOG - TMIN_LOG)*k)/SIM_ANNEAL_ITERS);

		// pick a read to move, and a target alignment
		size_t r = rand() % n_mmaps;
		size_t r_old = mmaps[r].active;
		size_t r_new = rand() % (mmaps[r].n - 1);
		if (r_new >= r_old)
			++r_new;  // make sure r_new != r_old, but still random

		SAMRecord *active_mate = NULL;
		size_t mate_r = 0;
		int mate_is_mmap = 0;

		if (mmaps[r].mate_umap >= 0) {
			mate_r = mmaps[r].mate_umap;
			active_mate = records[umaps[mate_r].idx];
		} else if (mmaps[r].mate_mmap >= 0) {
			mate_r = mmaps[r].mate_mmap;
			active_mate = records[mmaps[mate_r].idx + mmaps[mate_r].active];
			mate_is_mmap = 1;
		}

		SAMRecord *rec_old = records[mmaps[r].idx + r_old];
		SAMRecord *rec_new = records[mmaps[r].idx + r_new];

		double density_prob_change = 0.0;
		double score_prob_change = 0.0;
		int force_move = 0;
		int mate_new_active = -1;
		size_t mate_old_bin = 0;
		size_t mate_new_bin = 0;
		const int old_paired = active_mate != NULL && is_pair(rec_old, active_mate);
		const int new_paired = active_mate != NULL && is_pair(rec_new, active_mate);

		if (!old_paired && new_paired) {
			force_move = 1;
		} else if (old_paired && !new_paired && mate_is_mmap) {
			// try to move mate to preserve paired-end goodness
			for (int i = 0; i < mmaps[mate_r].n; i++) {
				SAMRecord *mate_rec_new = records[mmaps[mate_r].idx + i];
				if (is_pair(rec_new, mate_rec_new)) {
					SAMRecord *mate_rec_old = active_mate;
					mate_new_active = i;

					mate_old_bin = BIN_IDX_FOR_POS(mate_rec_old->pos, cloud_lo);
					mate_new_bin = BIN_IDX_FOR_POS(mate_rec_new->pos, cloud_lo);

					score_prob_change += (mate_rec_new->score - mate_rec_old->score)/SCORE_SCALE;
					break;
				}
			}
		}

		// compute change in probability
		const size_t old_bin = BIN_IDX_FOR_POS(rec_old->pos, cloud_lo);
		const size_t new_bin = BIN_IDX_FOR_POS(rec_new->pos, cloud_lo);
		const int p1 = (mate_new_active >= 0 && old_bin == mate_old_bin) ? 2 : 1;  // check for source/dest overlap
		const int p2 = (mate_new_active >= 0 && new_bin == mate_new_bin) ? 2 : 1;
		const double old_bin_prob_old = log_density_prob(bins[old_bin]);
		const double old_bin_prob_new = log_density_prob(bins[old_bin] - p1);
		const double new_bin_prob_old = log_density_prob(bins[new_bin]);
		const double new_bin_prob_new = log_density_prob(bins[new_bin] + p2);
		density_prob_change += (old_bin_prob_new - old_bin_prob_old) +
		                       (new_bin_prob_new - new_bin_prob_old);

		if (p1 == 1 && mate_new_active >= 0) {
			const double old_bin_prob_old = log_density_prob(bins[mate_old_bin]);
			const double old_bin_prob_new = log_density_prob(bins[mate_old_bin] - 1);
			density_prob_change += (old_bin_prob_new - old_bin_prob_old);
		}

		if (p2 == 1 && mate_new_active >= 0) {
			const double new_bin_prob_old = log_density_prob(bins[mate_new_bin]);
			const double new_bin_prob_new = log_density_prob(bins[mate_new_bin] + 1);
			density_prob_change += (new_bin_prob_new - new_bin_prob_old);
		}

		score_prob_change += (rec_new->score - rec_old->score)/SCORE_SCALE;
		const double prob_change = density_prob_change + score_prob_change;

		if (force_move || prob_change > 0 || exp(prob_change/t) >= ((double)rand())/RAND_MAX) {
			/* accept change */
			log_config_prob += prob_change;
			mmaps[r].active = r_new;
			bins[old_bin] -= 1;
			bins[new_bin] += 1;

			if (mate_new_active >= 0) {
				mmaps[mate_r].active = mate_new_active;
				bins[mate_old_bin] -= 1;
				bins[mate_new_bin] += 1;
			}
		} else {
			/* revert change */
			++no_move_count;
		}

		if (no_move_count >= MAX_NO_MOVE)
			goto exit;
	}
	exit:

	/* set active records */
	for (size_t i = 0; i < n_umaps; i++) {
		records[umaps[i].idx]->active = 1;
	}

	for (size_t i = 0; i < n_mmaps; i++) {
		records[mmaps[i].idx + mmaps[i].active]->active = 1;
	}

	free(records);  // recall: we reassigned `records` to a new, dup-free array
}

