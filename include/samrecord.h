#ifndef SAMRECORD_H
#define SAMRECORD_H

#include <stdint.h>
#include "align.h"
#include "bwabridge.h"
#include "util.h"

typedef struct fastq_record {
	bc_t bc;
	unsigned short rlen;
	char id[100];
	char read[MAX_READ_LEN+2];  // +2 for newline and null-term
	char qual[MAX_READ_LEN+2];
} FASTQRecord;

#define IS_SENTINEL(rec)  ((rec).id[0] == '\0')
#define SET_SENTINEL(rec) ((rec).id[0] = '\0')

/* compact representation of one line in a SAM file */
typedef struct sam_record {
	bc_t bc;
	chrom_t chrom;
	uint32_t pos;

	char ident[100];
	double score;
	int mapq;  // BWA MEM's mapq
	int score_mapq;  // our own computed mapq
	int clip;
	int clip_edit_dist;  // edit distance including clipping

	uint32_t hash;
	uint32_t mate_hash;

	uint32_t mate        : 1;
	uint32_t rev         : 1;
	uint32_t hashed      : 1;
	uint32_t mate_hashed : 1;
	uint32_t duplicate   : 1;
	uint32_t unique      : 1;
	uint32_t active      : 1;
	uint32_t visited     : 1;

	struct xa alts[MAX_ALTS];
	size_t n_alts;
	double gamma;
	Cloud *cloud;

	FASTQRecord *fq;
	FASTQRecord *fq_mate;
	struct sam_record *selected_mate;
	SingleReadAlignment aln;
} SAMRecord;

uint32_t record_hash(SAMRecord *record);
uint32_t record_hash_mate(SAMRecord *record);
uint32_t record_eq(SAMRecord *r1, SAMRecord *r2);
uint32_t record_eq_mate(SAMRecord *r1, SAMRecord *r2);
int record_cmp(const void *v1, const void *v2);

struct xa;
void print_sam_record(SAMRecord *rec,
                      SAMRecord *mate,
                      FILE *out,
                      const char *rg_id);

/* SAM specs */
#define SAM_NUM_FIELDS 11

/* SAM flags */
#define SAM_READ_PAIRED   1
#define SAM_READ_PROPER   2
#define SAM_READ_UNMAPPED 4
#define SAM_MATE_UNMAPPED 8
#define SAM_READ_REVERSED 16
#define SAM_MATE_REVERSED 32
#define SAM_1ST_IN_PAIR   64
#define SAM_2ND_IN_PAIR   128
#define SAM_READ_IS_A_DUP 1024

#endif /* SAMRECORD_H */

