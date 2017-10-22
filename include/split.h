#ifndef SPLIT_H
#define SPLIT_H

#include <stdlib.h>
#include "samrecord.h"
#include "align.h"

#define TMAX_LOG 0.0
#define TMIN_LOG (-12.0)
#define SIM_ANNEAL_ITERS 50000

//Cloud **split_cloud_sim_anneal(SAMRecord **records, size_t n_records, Cloud *clouds, size_t *n_split_clouds);

#define BIN_SIZE 1000
#define MAX_FRAG 1000000
#define MAX_BINS (MAX_FRAG/BIN_SIZE)
#define SCORE_SCALE 20

void mark_optimal_alignments_in_cloud(SAMRecord **records, size_t n_records);

#endif /* SPLIT_H */

