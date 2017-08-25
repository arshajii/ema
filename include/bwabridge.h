#ifndef BWABRIDGE_H
#define BWABRIDGE_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "bwa/bwa.h"
#include "bwa/bwt.h"
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"

/***************************************************************************
 * This is all translated from Lariat's BWA bridge.                        *
 * https://github.com/10XGenomics/lariat/blob/master/go/src/gobwa/gobwa.go *
 ***************************************************************************/

void arena_init(void);
void arena_clear(void);
void arena_destroy(void);
void arena_push(void *p);

/* BWA defs */
typedef struct {
		int64_t rbeg;
		int32_t qbeg, len;
		int score;
} mem_seed_t;  // unaligned memory

typedef struct {
		int n, m, first, rid;
		uint32_t w:29, kept:2, is_alt:1;
		float frac_rep;
		int64_t pos;
		mem_seed_t *seeds;
} mem_chain_t;

typedef struct { size_t n, m; mem_chain_t *a;  } mem_chain_v;

/* Convenience defs */
typedef struct {
	char *name;
	int32_t len;
} ContigInfo;

typedef struct {
	int64_t offset;
	int64_t aln_end;
	char *contig;
	int score;
	int read_s;
	int read_e;
	mem_alnreg_t *chained_hit;
	unsigned rev : 1;
	unsigned sec : 1;
} EasyAlignment;

typedef struct {
	EasyAlignment *a1;
	EasyAlignment *a2;
	size_t len1;
	size_t len2;
} EasyAlignmentPairs;

typedef struct {
	int64_t pos;
	char *chrom;
	int flag;
	int alt;
	int mapq;
	int edit_dist;
	uint32_t *cigar;
	int n_cigar;
	char *alt_mappings;
	int score;
	int sub;
	int alt_sc;
	//int read_s;
	//int read_e;
	//mem_aln_t *raw;
	unsigned rev : 1;
} SingleReadAlignment;

typedef struct {
	int64_t offset;
	char *contig;
	mem_chain_t *chain;
	unsigned rev : 1;
} Chain;

bwaidx_t *load_reference(const char *path);
void get_seq(bwaidx_t *ref, char *chrom, int64_t start, int64_t end, int rev, char *out);
EasyAlignment *bwa_align(bwaidx_t *ref, mem_opt_t *opts, char *seq, const size_t len);
Chain *bwa_chain(bwaidx_t *ref, mem_opt_t *opts, char *seq, const size_t len);
EasyAlignmentPairs bwa_mem_mate_sw(bwaidx_t *ref,
                                   mem_opt_t *opts,
                                   char *read1,
                                   const size_t len1,
                                   char *read2,
                                   const size_t len2,
                                   const int score_delta);
void bwa_smith_waterman(bwaidx_t *ref, mem_opt_t *opts, char *seq, const size_t len, mem_alnreg_t *aln, SingleReadAlignment *res);
void interpret_align(bwaidx_t *ref, mem_alnreg_t *caln, EasyAlignment *res);
void interpret_chain(bwaidx_t *ref, mem_chain_t *chn, Chain *res);
void interpret_single_read_alignment(bwaidx_t *ref, mem_aln_t *aln, SingleReadAlignment *res);

#endif /* BWABRIDGE_H */

