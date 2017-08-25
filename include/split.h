#ifndef SPLIT_H
#define SPLIT_H

#include <stdlib.h>
#include "samrecord.h"
#include "align.h"

#define TMAX_LOG 4.0
#define TMIN_LOG (-11.0)
#define SIM_ANNEAL_ITERS 50000

Cloud **split_cloud_sim_anneal(SAMRecord **records, size_t n_records, Cloud *clouds, size_t *n_split_clouds);

#endif /* SPLIT_H */

