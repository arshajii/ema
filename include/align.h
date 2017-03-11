#ifndef ALIGN_H
#define ALIGN_H

#include <stdlib.h>
#include <stdint.h>
#include "samrecord.h"

void find_clouds_and_align(SAMRecord *records, const size_t n_records, const char *out);

typedef struct cloud {
	uint32_t lo;
	uint32_t hi;
	uint32_t cov;

	int snvs;
	int indels;

	double exp_cov;
	double reads_prob;

	/* disjoint set for probability normalization */
	struct cloud *parent;
	struct cloud *child;
} Cloud;

/* reasonable clouds/barcode upper bound */
#define MAX_CLOUDS_PER_BC 100000

/* cloud distance threshold */
#define DIST_THRESH 50000

/* lengths */
#define BC_LEN     16
#define MATE1_TRIM 7
#define READ_LEN   151
#define MATE1_LEN  (READ_LEN - BC_LEN - MATE1_TRIM)
#define MATE2_LEN  READ_LEN

/* fragment properties */
#define MEAN_FRAG_LEN  1e5
#define MEAN_FRAG_COV  300.0
#define FRAG_COV_SCALE 20.0

/* read properties */
#define INSERT_AVG 250
#define INSERT_MIN 15
#define INSERT_MAX 500
#define INSERT_CAP 1500
#define UNPAIRED_PENALTY (-60.0)
#define LONG_INS_PENALTY (-15.0)

#define ERROR_RATE 0.001
#define INDEL_RATE 0.0001

#endif /* ALIGN_H */

