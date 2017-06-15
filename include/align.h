#ifndef ALIGN_H
#define ALIGN_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

void find_clouds_and_align(FILE *fq1, FILE *fq2, const char *ref_path, FILE *out, const char *rg);

typedef struct cloud {
	uint32_t lo;
	uint32_t hi;
	uint32_t cov;

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

/* read properties */
#define INSERT_AVG 250
#define INSERT_MIN (-35)
#define INSERT_MAX 750
#define UNPAIRED_PENALTY (-3.0)

#define ERROR_RATE 0.001
#define INDEL_RATE 0.0001
#define CLIP_RATE  0.01

#define EXTRA_SEARCH_DEPTH 7

#endif /* ALIGN_H */

