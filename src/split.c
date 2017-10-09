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

/*
 * This is all for resolving these types of alignments:
 *
 *      align1        align2
 *     |------|      |------|
 * ...========================...  <= 50kb
 *           \        /
 *            \      /
 *             \    /
 *              \  /
 *               \/
 *            |------|
 *              read
 *
 * i.e. same read mapping multiple times in the same cloud.
 * The approach is to split the cloud into sub-clouds that
 * are as similar to one another as possible, thereby making
 * it likely that the resulting alignment patterns correspond
 * to the true source fragment. We use simulated annealing
 * to minimize the variance between distances of corresponding
 * reads within each subcloud.
 */

static double var(const double *x, const size_t n)
{
	if (n == 0)
		return 0.0;

	double var = 0.0;
	double avg = 0.0;

	for (size_t i = 0; i < n; i++)
		avg += x[i];

	avg /= n;

	for (size_t i = 0; i < n; i++)
		var += pow(x[i] - avg, 2);

	return var;
}

/* distance-variance objective we are trying to minimize when splitting */
static double dvar_for_rec(size_t record_idx,
                           SAMRecord **records,
                           Cloud **cloud_assignments,
                           const size_t n_records)
{
#define MAX_PAIRS 500
	size_t n_targets = 0;
	double dist_buf[MAX_PAIRS];
	size_t n_pairs = 0;
	double total_var = 0.0;

	SAMRecord *record = records[record_idx];

	/* find first */
	while (record_idx > 0 && record_eq(records[record_idx - 1], record))
		--record_idx;

	/* figure out where our target read appears */
	for (size_t i = record_idx; i < n_records; i++) {
		if (record_eq(records[i], record)) {
			if (cloud_assignments[i] == NULL)
				return 0.0;

			++n_targets;
		} else {
			break;
		}
	}

	/* compute distance variances between each other read */
	for (size_t i = 0; i < n_records;) {
		SAMRecord *r = records[i];

		if (record_eq(r, record)) {
			i++;
			continue;
		}

		size_t k;
		for (k = i; k < n_records; k++) {
			SAMRecord *r2 = records[k];
			if (record_eq(r, r2)) {
				/* check if our two reads actually appear in the same cloud */
				Cloud *c2 = cloud_assignments[k];
				int shared_cloud = -1;

				if (c2 != NULL) {
					for (size_t l = 0; l < n_targets; l++) {
						if (cloud_assignments[record_idx + l] == c2) {
							shared_cloud = l;
							break;
						}
					}

					if (shared_cloud >= 0) {
						assert(n_pairs < MAX_PAIRS-1);
						dist_buf[n_pairs++] = ((double)(records[record_idx + shared_cloud]->pos)) - ((double)(r2->pos));
					}
				}
			} else {
				break;
			}
		}

		i = k;
		total_var += var(dist_buf, n_pairs);
		n_pairs = 0;
	}

	return total_var;
#undef MAX_PAIRS
}

static void init_cloud_endpoints(Cloud *c)
{
	c->hi = 0x00000000;
	c->lo = 0xffffffff;
}

static void rec_to_cloud(SAMRecord *r, Cloud *c)
{
	const uint32_t pos0 = r->pos;
	const uint32_t pos1 = pos0 + (r->mate == 0 ? MATE1_LEN : MATE2_LEN);
	if (pos1 > c->hi)
		c->hi = pos1;
	if (pos0 < c->lo)
		c->lo = pos0;
}

static int clouds_overlap(Cloud *c1, Cloud *c2)
{
	const uint32_t a1 = c1->lo;
	const uint32_t b1 = c1->hi;
	const uint32_t a2 = c2->lo;
	const uint32_t b2 = c2->hi;
	assert(a1 < b1 && a2 < b2);
	return MAX(a1, a2) <= MIN(b1, b2);
}

/* caution: `records` should be name-sorted */
Cloud **split_cloud_sim_anneal(SAMRecord **records, size_t n_records, Cloud *clouds, size_t *n_split_clouds)
{
#define MAX_NO_MOVE 500

	static int rand_init = 0;
	if (!rand_init) {
		time_t t;
		srand((unsigned)time(&t));
		rand_init = 1;
	}

	struct multimapped_read {
		size_t idx;  // index in records
		int n;       // how many mappings in the to-be-split cloud
	} mmaps[MAX_CLOUDS_PER_BC];
	size_t n_mmaps = 0;

	Cloud **cloud_assignments = safe_calloc(n_records, sizeof(*cloud_assignments));
	double *dvar_cache = safe_calloc(n_records, sizeof(*dvar_cache));

	/* find the multi-mapped reads, record highest */
	size_t i = 0;
	size_t n_clouds = 1;
	while (i < n_records) {
		size_t j = i+1;
		while (j < n_records && record_eq(records[j], records[i])) {
			++j;
		}

		const size_t n = j-i;

		if (n > 1) {
			mmaps[n_mmaps].idx = i;
			mmaps[n_mmaps].n = n;

			if (n > n_clouds) {
				n_clouds = n;
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
		goto bailout;
	}

	*n_split_clouds = n_clouds;

	/* initialize clouds */
	for (size_t i = 0; i < n_clouds; i++) {
		init_cloud(&clouds[i]);
		init_cloud_endpoints(&clouds[i]);
	}

	/* preliminary cloud assignments */
	for (size_t i = 0; i < n_mmaps; i++) {
		const size_t base = mmaps[i].idx;
		for (int j = 0; j < mmaps[i].n; j++) {
			cloud_assignments[base + j] = &clouds[j];
		}
	}

	/* simulated annealing to minimize distance variance */
	int no_move_count = 0;
	for (size_t k = 0; k < SIM_ANNEAL_ITERS; k++) {
		const double t = pow(10.0, TMAX_LOG - ((TMAX_LOG - TMIN_LOG)*k)/SIM_ANNEAL_ITERS);

		size_t c1 = rand() % n_clouds;
		size_t c2 = rand() % (n_clouds-1);
		if (c2 >= c1)  // we want unique clouds
			++c2;

		const size_t r = rand() % n_mmaps;
		int found_c1 = -1;
		int found_c2 = -1;
		const size_t base = mmaps[r].idx;

		for (int i = 0; i < mmaps[r].n; i++) {
			if (cloud_assignments[base + i] == &clouds[c1])
				found_c1 = base + i;
			else if (cloud_assignments[base + i] == &clouds[c2])
				found_c2 = base + i;
		}

		if (found_c1 < 0 && found_c2 < 0)
			continue;

		if (found_c1 < 0) {
			const size_t tmp = c1;
			c1 = c2;
			c2 = tmp;
			found_c1 = found_c2;
			found_c2 = -1;
		}

		double s0;
		if (dvar_cache[base] == 0.0) {
			s0 = dvar_for_rec(base, records, cloud_assignments, n_records);
			dvar_cache[base] = s0;
		} else {
			s0 = dvar_cache[base];
		}

		if (found_c2 < 0) {  // only c1 has read
			/* make the move */
			cloud_assignments[found_c1] = &clouds[c2];

			const double s = dvar_for_rec(base, records, cloud_assignments, n_records);

			if (s < s0 || exp(-sqrt(s - s0)/t) >= ((double)rand())/RAND_MAX) {
				/* accept change */
				dvar_cache[base] = s;
				no_move_count = 0;
			} else {
				/* revert change */
				cloud_assignments[found_c1] = &clouds[c1];
				++no_move_count;
			}
		} else {  // both c1 and c2 have read
			/* make the move */
			cloud_assignments[found_c1] = &clouds[c2];
			cloud_assignments[found_c2] = &clouds[c1];

			const double s = dvar_for_rec(base, records, cloud_assignments, n_records);

			if (s < s0 || exp(-sqrt(s - s0)/t) >= ((double)rand())/RAND_MAX) {
				/* accept change */
				dvar_cache[base] = s;
				no_move_count = 0;
			} else {
				/* revert change */
				cloud_assignments[found_c1] = &clouds[c1];
				cloud_assignments[found_c2] = &clouds[c2];
				++no_move_count;
			}
		}

		if (no_move_count >= MAX_NO_MOVE)
			break;
	}

	/* shouldn't have too many clouds, so use naive overlap check */
	for (size_t i = 0; i < n_records; i++) {
		if (cloud_assignments[i] != NULL)
			rec_to_cloud(records[i], cloud_assignments[i]);
	}

	for (size_t i = 0; i < n_records; i++) {
		Cloud *c1 = cloud_assignments[i];
		if (c1 == NULL)
			continue;

		for (size_t j = i+1; j < n_records; j++) {
			Cloud *c2 = cloud_assignments[j];
			if (c2 == NULL)
				continue;

			if (clouds_overlap(c1, c2)) {
				goto bailout;
			}
		}
	}

	free(dvar_cache);
	return cloud_assignments;

	bailout:
	*n_split_clouds = 1;

	Cloud *c = &clouds[0];
	init_cloud(c);

	for (size_t i = 0; i < n_records; i++) {
		cloud_assignments[i] = c;
	}

	return cloud_assignments;

#undef MAX_NO_MOVE
}

static double log_density_prob(unsigned int density)
{
	/* density probabilities per 1k bin */
	static double density_probs[] = {/*0.6*/ 0.01, 0.1, 0.15, 0.05, 0.04};
	static double density_probs_log[SIZE(density_probs)];
	static const size_t size = SIZE(density_probs);
	static double log_remainder;
	static int init = 0;

	if (!init) {
		double total = 0;
		for (size_t i = 0; i < size; i++) {
			const double p = density_probs[i];
			density_probs_log[i] = log(p);
			total += p;
		}

		log_remainder = log(1.0 - total);
		init = 1;
	}

	if (density < size)
		return density_probs_log[density];
	else
		return log_remainder - (density - size + 1)*log(2.0);  // scale down larger density probs exponentially
}

/* caution: `records` should be name-sorted */
void mark_optimal_alignments_in_cloud(SAMRecord **records, size_t n_records)
{
#define MAX_NO_MOVE 500
#define BUF_SIZE    30000

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
		int active;  // which of these alignments is active
	} mmaps[BUF_SIZE];

	size_t n_umaps = 0;
	size_t n_mmaps = 0;

	uint32_t cloud_lo = 0xffffffff;
	uint32_t cloud_hi = 0x00000000;

	if (n_records >= BUF_SIZE)
		return;

	/* get rid of duplicates */
	/*
	SAMRecord **records_no_dups = safe_malloc(n_records * sizeof(*records_no_dups));
	size_t n_records_no_dups = 0;

	for (size_t i = 0; i < n_records; i++) {
		SAMRecord *rec = records[i];
		if (!rec->duplicate)
			records_no_dups[n_records_no_dups++] = rec;
	}

	records = records_no_dups;
	n_records = n_records_no_dups;
	*/

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
			mmaps[n_mmaps++] = (struct multimapped_read){.idx = i, .n = n, .active = 0};
		} else {
			umaps[n_umaps++] = (struct uniquemapped_read){.idx = i};
		}

		log_config_prob += records[i]->score;  // add up active alignment scores
		i = j;
	}

	/* get initial configuration probability */
	const size_t n_bins = (cloud_hi - cloud_lo)/BIN_SIZE + 1;

	if (n_bins >= MAX_BINS) {
		//free(records);
		return;
	}

	for (size_t i = 0; i < n_records; i++) {
		SAMRecord *rec = records[i];
		rec->active = 0;  // we will re-set the active ones later
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

		SAMRecord *rec_old = records[mmaps[r].idx + r_old];
		SAMRecord *rec_new = records[mmaps[r].idx + r_new];

		// compute change in probability
		const size_t old_bin = BIN_IDX_FOR_POS(rec_old->pos, cloud_lo);
		const size_t new_bin = BIN_IDX_FOR_POS(rec_new->pos, cloud_lo);

		// bear with me here...
		const double old_bin_prob_old = log_density_prob(bins[old_bin]);
		const double old_bin_prob_new = log_density_prob(bins[old_bin] - 1);
		const double new_bin_prob_old = log_density_prob(bins[new_bin]);
		const double new_bin_prob_new = log_density_prob(bins[new_bin] + 1);
		const double density_prob_change = (old_bin_prob_new - old_bin_prob_old) +
		                                   (new_bin_prob_new - new_bin_prob_old);

		const double score_prob_change = rec_new->score - rec_old->score;
		const double prob_change = density_prob_change + score_prob_change;

		if (prob_change > 0 || exp(prob_change/t) >= ((double)rand())/RAND_MAX) {
			/* accept change */
			log_config_prob += prob_change;
			mmaps[r].active = r_new;
			bins[old_bin] -= 1;
			bins[new_bin] += 1;
		} else {
			/* revert change */
			++no_move_count;
		}

		if (no_move_count >= MAX_NO_MOVE)
			break;
	}

	/* set active records */
	for (size_t i = 0; i < n_umaps; i++) {
		records[umaps[i].idx]->active = 1;
	}

	for (size_t i = 0; i < n_mmaps; i++) {
		records[mmaps[i].idx + mmaps[i].active]->active = 1;
	}

	//free(records);  // recall: we reassigned `records` to a new, dup-free array
}

