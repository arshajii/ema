#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "main.h"
#include "align.h"
#include "util.h"
#include "samrecord.h"

uint32_t record_hash(SAMRecord *record)
{
	if (record->hash != 0)
		return record->hash;

	const uint32_t hash = hash_ident(record->ident)*record->mate;
	record->hash = hash;
	return hash;
}

uint32_t record_hash_mate(SAMRecord *record)
{
	if (record->mate_hash != 0)
		return record->mate_hash;

	const uint32_t hash = hash_ident(record->ident) * ((record->mate == 1) ? 2 : 1);
	record->mate_hash = hash;
	return hash;
}

uint32_t record_eq(SAMRecord *r1, SAMRecord *r2)
{
	return record_hash(r1) == record_hash(r2) &&
	       r1->mate == r2->mate &&
	       strcmp(r1->ident, r2->ident) == 0;
}

uint32_t record_eq_mate(SAMRecord *r1, SAMRecord *r2)
{
	return record_hash_mate(r1) == record_hash(r2) &&
	       r1->mate != r2->mate &&
	       strcmp(r1->ident, r2->ident) == 0;
}

int record_cmp(const void *v1, const void *v2) {
	const SAMRecord *r1 = (SAMRecord *)v1;
	const SAMRecord *r2 = (SAMRecord *)v2;

	const bc_t bc1 = r1->bc;
	const bc_t bc2 = r2->bc;
	int c = (bc1 > bc2) - (bc1 < bc2);

	if (c != 0) return c;

	const uint8_t chrom1 = r1->chrom;
	const uint8_t chrom2 = r2->chrom;
	c = (chrom1 > chrom2) - (chrom1 < chrom2);

	if (c != 0) return c;

	const uint32_t pos1 = r1->pos;
	const uint32_t pos2 = r2->pos;
	c = (pos1 > pos2) - (pos1 < pos2);

	if (c != 0) return c;

	return strcmp(r1->ident, r2->ident);
}

int is_dup(SAMRecord *r1, SAMRecord *r2)
{
	return strcmp(r1->ident, r2->ident) == 0 &&
	       r1->chrom == r2->chrom &&
	       r1->pos == r2->pos &&
	       r1->mate == r2->mate &&
	       r1->rev == r2->rev &&
	       r1->score == r2->score;
}

/* caution: frees given records buffer */
SAMRecord *remove_dups(SAMRecord *records, size_t *n_records_ptr)
{
	const size_t n_records = *n_records_ptr;
	SAMRecord *records_no_dups = malloc((n_records+1) * sizeof(*records_no_dups));
	size_t n_records_no_dups = 0;
	size_t i = 0;

	while (i < n_records) {
		SAMRecord *rec = &records[i++];
		records_no_dups[n_records_no_dups++] = *rec;

		while (i < n_records && is_dup(rec, &records[i])) {
			i++;
		}
	}

	free(records);
	*n_records_ptr = n_records_no_dups;
	records_no_dups[n_records_no_dups].bc = 0;
	return records_no_dups;
}

static double score_cigar(char *cigar, int *snvs, int *indels);

int parse_sam_record(char *record, SAMRecord *out)
{
	char *split_buf[64];
	split_line(record, split_buf);

	if (split_buf[2][0] == '*' || split_buf[3][0] == '0' || split_buf[5][0] == '*') {
		return 0;
	}

	/* find barcode */
	char *bc_str = split_buf[1];
	while (*bc_str != ':') --bc_str;

	const bc_t bc = encode_bc(bc_str + 1);

	const int flags = atoi(split_buf[1]);

	if (flags & SAM_READ_IS_A_DUP)
		return 0;

	const uint8_t chrom = chrom_index(split_buf[2]);
	const uint32_t pos = atol(split_buf[3]);

	out->bc = bc;
	out->chrom = chrom;
	out->pos = pos;

	size_t i = 0;
	for (char *c = record; c != bc_str; c++) {
		out->ident[i++] = *c;
	}
	out->ident[i] = '\0';

	out->score = score_cigar(split_buf[5], &out->snvs, &out->indels);

	size_t rlen = 0;
	while (!isspace(split_buf[9][++rlen]));
	out->mate = ((rlen == MATE1_LEN) ? 1 : 2);
	out->hash = 0;
	out->mate_hash = 0;
	out->rev = ((flags & SAM_READ_REVERSED) != 0);

	return 1;
}

double LOG_MATCH_SCORE;
double LOG_MISMATCH_SCORE;
double LOG_INDEL_SCORE;

static double score_cigar(char *cigar, int *snvs, int *indels)
{
	static int init = 0;

	if (!init) {
		LOG_MATCH_SCORE = log(1 - ERROR_RATE);
		LOG_MISMATCH_SCORE = log(ERROR_RATE);
		LOG_INDEL_SCORE = log(INDEL_RATE);
		init = 1;
	}

	char *p = &cigar[0];
	double score = 0;

	int snvs_ = 0;
	int indels_ = 0;

	while (isdigit(*p)) {
		const int n = atoi(p);
		while (isdigit(*++p));
		const char code = *p++;

		switch (code) {
		case 'M':
			fprintf(stderr, "error: ambiguous cigar character: M\n");
			goto fail;
		case '=':
			score += n*LOG_MATCH_SCORE;
			break;
		case 'X':
		case 'S':
			score += n*LOG_MISMATCH_SCORE;
			snvs_ += n;
			break;
		case 'I':
		case 'N':
		case 'D':
			score += n*LOG_INDEL_SCORE;
			indels_ += n;
			break;
		default:
			fprintf(stderr, "error: invalid cigar character: %c\n", code);
			goto fail;
		}
	}

	if (snvs != NULL)
		*snvs = snvs_;

	if (indels != NULL)
		*indels = indels_;

	return score;

	fail:
	exit(EXIT_FAILURE);
}

