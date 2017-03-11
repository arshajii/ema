#ifndef SAMRECORD_H
#define SAMRECORD_H

#include <stdint.h>
#include "util.h"

/* compact representation of one line in a SAM file */
typedef struct {
	bc_t bc;
	chrom_t chrom;
	uint32_t pos;
	uint8_t mate;
	char ident[64];
	double score;
	int snvs;
	int indels;

	uint32_t hash;
	uint32_t mate_hash;

	uint32_t rev : 1;
} SAMRecord;

uint32_t record_hash(SAMRecord *record);
uint32_t record_hash_mate(SAMRecord *record);
uint32_t record_eq(SAMRecord *r1, SAMRecord *r2);
uint32_t record_eq_mate(SAMRecord *r1, SAMRecord *r2);
int record_cmp(const void *v1, const void *v2);
int is_dup(SAMRecord *r1, SAMRecord *r2);
SAMRecord *remove_dups(SAMRecord *records, size_t *n_records_ptr);
int parse_sam_record(char *record, SAMRecord *out);

/* SAM flags */
#define SAM_READ_UNMAPPED 4
#define SAM_MATE_UNMAPPED 8
#define SAM_READ_REVERSED 16
#define SAM_MATE_REVERSED 32
#define SAM_READ_IS_A_DUP 1024

#endif /* SAMRECORD_H */

