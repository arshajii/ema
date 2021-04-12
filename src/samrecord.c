#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include "main.h"
#include "align.h"
#include "util.h"
#include "samrecord.h"

uint32_t record_hash(SAMRecord *record)
{
	if (record->hashed)
		return record->hash;

	const uint32_t hash = hash_ident(record->ident)*(record->mate + 1);
	record->hash = hash;
	record->hashed = 1;
	return hash;
}

uint32_t record_hash_mate(SAMRecord *record)
{
	if (record->mate_hashed)
		return record->mate_hash;

	const uint32_t hash = hash_ident(record->ident)*(2 - record->mate);
	record->mate_hash = hash;
	record->mate_hashed = 1;
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

int record_cmp(const void *v1, const void *v2)
{
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

// from BWA source
static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; k++) {
		int op = cigar[k] & 0xf;
		if (op == 0 || op == 2)
			l += (cigar[k] >> 4);
	}
	return l;
}

static inline char rc(const char c)
{
	switch (c) {
	case 'A':
		return 'T';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	case 'T':
		return 'A';
	case 'N':
		return 'N';
	}

	assert(0);
}

void print_sam_record(SAMRecord *rec,
                      SAMRecord *mate,
                      FILE *out,
                      const char *rg_id,
		      const int is_haplotag)
{
	assert(rec != NULL || mate != NULL);
	int flag = SAM_READ_PAIRED;
	char *ident;
	char *chrom = "*";
	uint32_t pos = 0;
	int mapq = 0;
	int read_len;
	bc_t bc;
	SingleReadAlignment *r = NULL;
	FASTQRecord *fq;

	struct xa *alts = NULL;
	size_t n_alts = 0;
	double gamma = 0;
	Cloud *cloud = 0;

	if (rec != NULL) {
		ident = rec->ident;
		chrom = chrom_lookup(rec->chrom);
		pos = rec->pos;
		read_len = rec->fq->rlen;
		bc = rec->bc;
		r = &rec->aln;
		fq = rec->fq;

		alts = rec->alts;
		n_alts = rec->n_alts;
		gamma = rec->gamma;
		cloud = rec->cloud;

		assert(!isnan(gamma));

		const int gamma_mapq = ((gamma <= 0.999999) ? (int)(-10*log10(1 - gamma)) : 60);
		const int score_mapq = rec->score_mapq;
		const int bwa_mapq   = rec->mapq;
		mapq = MIN(gamma_mapq, score_mapq);
		mapq = MIN(mapq, bwa_mapq);
		mapq = MAX(mapq, 0);
		mapq = MIN(mapq, 60);

		if (rec->rev)
			flag |= SAM_READ_REVERSED;

		if (rec->duplicate)
			flag |= SAM_READ_IS_A_DUP;

		flag |= ((rec->mate == 0) ? SAM_1ST_IN_PAIR : SAM_2ND_IN_PAIR);
	} else {
		ident = mate->ident;
		read_len = mate->fq_mate->rlen;
		bc = mate->bc;
		fq = mate->fq_mate;
		flag |= SAM_READ_UNMAPPED;
		flag |= ((mate->mate == 0) ? SAM_2ND_IN_PAIR : SAM_1ST_IN_PAIR);
	}

	if (mate != NULL) {
		if (rec != NULL && is_pair(rec, mate))
			flag |= SAM_READ_PROPER;

		if (mate->rev)
			flag |= SAM_MATE_REVERSED;
	} else {
		flag |= SAM_MATE_UNMAPPED;
	}

	// basics
	fprintf(out, "%s\t%d\t%s\t%u\t%d\t", ident, flag, chrom, pos, mapq);

	// cigar
	if (rec != NULL) {
		const uint32_t *cigar = r->cigar;
		const int cigar_len = r->n_cigar;
		for (int i = 0; i < cigar_len; i++) {
			const uint32_t op   = cigar[i];
			const uint32_t type = op & 0xf;
			const uint32_t n    = op >> 4;
			fprintf(out, "%u%c", n, "MIDSS"[type]);  // note: convert hard to soft clipping
		}
	} else {
		fputc('*', out);
	}

	// mate mapping
	if (mate != NULL) {
		const int same_chrom = (rec != NULL) && (mate->chrom == rec->chrom);
		fprintf(out, "\t%s\t%d", same_chrom ? "=" : chrom_lookup(mate->chrom), mate->pos);

		SingleReadAlignment *s = &mate->aln;
		if (same_chrom) {
			int64_t p0 = r->pos + (r->rev ? get_rlen(r->n_cigar, r->cigar) - 1 : 0);
			int64_t p1 = s->pos + (s->rev ? get_rlen(s->n_cigar, s->cigar) - 1 : 0);
			if (s->n_cigar == 0 || r->n_cigar == 0)
				fprintf(out, "\t0");
			else
				fprintf(out, "\t%ld", (long)(-(p0 - p1 + (p0 > p1 ? 1 : p0 < p1 ? -1 : 0))));
		} else {
			fprintf(out, "\t0");
		}
	} else {
		fprintf(out, "\t*\t0\t0");
	}

	// seq and qual
	fputc('\t', out);
	if (rec != NULL && rec->rev) {
		for (int i = read_len - 1; i >= 0; i--) {
			fputc(rc(fq->read[i]), out);
		}

		fputc('\t', out);

		for (int i = read_len - 1; i >= 0; i--) {
			fputc(fq->qual[i], out);
		}
	} else {
		for (int i = 0; i < read_len; i++) {
			fputc(fq->read[i], out);
		}

		fputc('\t', out);

		for (int i = 0; i < read_len; i++) {
			fputc(fq->qual[i], out);
		}
	}

	// tags
	if (is_haplotag)
	{
		char bc_str[13] = {0};
		decode_bc(bc, bc_str, 1);
		if (rec != NULL) {
			fprintf(out, "\tNM:i:%d\tBX:Z:%s\tXG:f:%.5g\tMI:i:%d\tXF:i:%d", r->edit_dist, bc_str, gamma, cloud->id, cloud->bad);
		} else {
			fprintf(out, "\tBX:Z:%s", bc_str);
		}
	}
	else
	{
		char bc_str[BC_LEN + 1] = {0};
		decode_bc(bc, bc_str, 0);
		if (rec != NULL) {
			fprintf(out, "\tNM:i:%d\tBX:Z:%s-%s\tXG:f:%.5g\tMI:i:%d\tXF:i:%d", r->edit_dist, bc_str, bx_index, gamma, cloud->id, cloud->bad);
		} else {
			fprintf(out, "\tBX:Z:%s-1", bc_str);
		}
	}
	
	if (rg_id != NULL) {
		fprintf(out, "\tRG:Z:");
		for (size_t i = 0; rg_id[i] != '\0' && !isspace(rg_id[i]); i++)
			fputc(rg_id[i], out);
	}

	if (n_alts > 0) {
		fprintf(out, "\tXA:Z:");

		for (size_t i = 0; i < n_alts; i++) {
			struct xa *alt = &alts[i];
			fprintf(out, "%s,%s%d,", alt->chrom, alt->rev ? "-" : "+", alt->pos);
			const uint32_t *cigar = alt->cigar;
			const int cigar_len = alt->cigar_len;
			for (int i = 0; i < cigar_len; i++) {
				const uint32_t op   = cigar[i];
				const uint32_t type = op & 0xf;
				const uint32_t n    = op >> 4;
				fprintf(out, "%u%c", n, "MIDSS"[type]);  // note: convert hard to soft clipping
			}
			fprintf(out, ",%d;", alt->edit_dist);
		}
	}

	fputc('\n', out);
}

