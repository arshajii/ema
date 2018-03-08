#ifndef ALIGN_H
#define ALIGN_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

void bwa_init(const char *ref_path);
void bwa_dealloc(void);

void find_clouds_and_align(FILE *fq1, FILE *fq2, FILE *fqx, FILE *out, const int apply_opt);

typedef struct cloud {
	double exp_cov;
	double weight;

	uint32_t hi, lo;  // start+end points

	/* disjoint set for probability normalization */
	struct cloud *parent;
	struct cloud *child;

	int id;
	unsigned bad : 1;
} Cloud;

struct xa {
	char *chrom;
	uint32_t pos;
	int edit_dist;
	int cigar_len;
	int rev;
	uint32_t cigar[64];
};

void init_cloud(Cloud *c);

struct sam_record;
int is_pair(struct sam_record *r1, struct sam_record *r2);
int is_pair_relaxed(struct sam_record *r1, struct sam_record *r2);

#define EM_ITERS 5

/* reasonable clouds/barcode upper bound */
#define MAX_CLOUDS_PER_BC_SMALL   100000
#define MAX_CLOUDS_PER_BC_LARGE 10000000

/* lengths */
#define BC_LEN     16
#define MATE1_TRIM  7

#define MIN_READ_LEN (BC_LEN + MATE1_TRIM + 1)
#define MAX_READ_LEN 200
#define MAX_ID_LEN   100

/* read properties */
#define INSERT_AVG 250
#define INSERT_MIN (-35)
#define INSERT_MAX 750
#define UNPAIRED_PENALTY (-15.0)

#define INDEL_RATE 0.0001
#define CLIP_RATE  0.03

#define EXTRA_SEARCH_DEPTH       12
#define SPLIT_EXTRA_SEARCH_DEPTH 5
#define SPLIT_CLIP_THRESH        15

#define SECONDARY_ALIGN_THRESH 0.9
#define MAX_ALTS 15

#endif /* ALIGN_H */

