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
		*n_split_clouds = 1;

		Cloud *c = &clouds[0];
		init_cloud(c);

		for (size_t i = 0; i < n_records; i++) {
			cloud_assignments[i] = c;
		}

		return cloud_assignments;
	}

	*n_split_clouds = n_clouds;

	/* initialize clouds */
	for (size_t i = 0; i < n_clouds; i++) {
		init_cloud(&clouds[i]);
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

	free(dvar_cache);
	return cloud_assignments;

#undef MAX_NO_MOVE
}

