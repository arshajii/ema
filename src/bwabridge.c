#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "util.h"
#include "bwabridge.h"

/***************************************************************************
 * This is all translated from Lariat's BWA bridge.                        *
 * https://github.com/10XGenomics/lariat/blob/master/go/src/gobwa/gobwa.go *
 ***************************************************************************/

extern mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf);
extern mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf);
extern mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar);
extern int mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], const mem_alnreg_t *a, int l_ms, const uint8_t *ms, mem_alnreg_v *ma);
extern uint8_t *bns_fetch_seq(const bntseq_t *bns, const uint8_t *pac, int64_t *beg, int64_t mid, int64_t *end, int *rid);


/* Arena */

#define ARENA_INIT_CAP 100

typedef struct {
	void **p;
	size_t len;
	size_t cap;
} Arena;

static Arena arena = {NULL, 0ULL, 0ULL};

#pragma omp threadprivate(arena)
void arena_init(void)
{
	if (arena.p == NULL) {
		arena.p = safe_malloc(ARENA_INIT_CAP * sizeof(void *));
		arena.len = 0;
		arena.cap = ARENA_INIT_CAP;
	}
}

#pragma omp threadprivate(arena)
void arena_clear(void)
{
	for (size_t i = 0; i < arena.len; i++)
		free(arena.p[i]);
	arena.len = 0;
}

#pragma omp threadprivate(arena)
void arena_destroy(void)
{
	arena_clear();
	free(arena.p);
	arena.p = NULL;
	arena.len = 0;
	arena.cap = 0;
}

#pragma omp threadprivate(arena)
void arena_push(void *p)
{
	if (arena.len == arena.cap) {
		const size_t new_cap = (arena.len * 3)/2 + 1;
		arena.p = safe_realloc(arena.p, new_cap * sizeof(void *));
		arena.cap = new_cap;
	}

	arena.p[arena.len++] = p;
}


/* Convenience wrappers */

char **contig_ids;

bwaidx_t *load_reference(const char *path)
{
	bwaidx_t *ref = bwa_idx_load(path, BWA_IDX_ALL);

	if (ref == NULL) {
		fprintf(stderr, "error: could not load reference at %s\n", path);
		exit(EXIT_FAILURE);
	}

	bntseq_t *contigs = ref->bns;
	int32_t n_contigs = contigs->n_seqs;
	contig_ids = safe_malloc((n_contigs + 1) * sizeof(*contig_ids));

	for (int32_t i = 0; i < n_contigs; i++) {
		contig_ids[i] = contigs->anns[i].name;
	}

	contig_ids[n_contigs] = NULL;
	return ref;
}

static int32_t get_contig_id(const char *chrom)
{
	for (int32_t i = 0; contig_ids[i] != NULL; i++) {
		if (strcmp(chrom, contig_ids[i]) == 0)
			return i;
	}

	return -1;
}

void get_seq(bwaidx_t *ref, char *chrom, int64_t start, int64_t end, int rev, char *out)
{
	static const char two_bit_to_seq[]      = {'A', 'C', 'G', 'T'};
	static const char two_bit_to_seq_comp[] = {'T', 'G', 'C', 'A'};

	bntseq_t *contigs = ref->bns;
	int32_t contig_id = get_contig_id(chrom);
	assert(contig_id >= 0);
	bntann1_t *contig = &contigs->anns[contig_id];

	int64_t offset = contig->offset;
	int64_t offstart = start + offset;
	int64_t offend = end + offset;

	uint8_t *result = bns_fetch_seq(ref->bns, ref->pac, &offstart, (offstart + offend) >> 1, &offend, &contig_id);
	size_t len = offend - offstart;

	if (rev) {
		for (size_t i = 0; i < len; i++)
			out[len - i - 1] = two_bit_to_seq_comp[result[i]];
	} else {
		for (size_t i = 0; i < len; i++)
			out[i] = two_bit_to_seq[result[i]];
	}
	out[len] = '\0';

	free(result);
}

ContigInfo *get_ref_contigs_info(bwaidx_t *ref)
{
	bntseq_t *contigs = ref->bns;
	int32_t n_contigs = contigs->n_seqs;
	ContigInfo *cinfo = safe_malloc(n_contigs * sizeof(*cinfo));

	for (int32_t i = 0; i < n_contigs; i++) {
		bntann1_t *contig = &contigs->anns[i];
		cinfo[i] = (ContigInfo){ .name = contig->name, .len = contig->len };
	}

	return cinfo;
}

static char *seq_convert(const char *seq, const size_t len)
{
	char *buf = safe_malloc(len * sizeof(*buf));
	for (size_t i = 0; i < len; i++)
		buf[i] = nst_nt4_table[(int)seq[i]];
	return buf;
}

typedef struct {
	int64_t pos;     // forward strand 5'-end mapping position
	int rid;         // reference sequence index in bntseq_t; <0 for unmapped
	int flag;        // extra flag
	uint32_t flag2;  //uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
	int n_cigar;     // number of CIGAR operations
	uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
	char *XA;        // alternative mappings
	int score, sub, alt_sc;
} mem_aln_compact_flag_t;

EasyAlignment *bwa_align(bwaidx_t *ref, mem_opt_t *opts, char *seq, const size_t len)
{
	char *seq_conv = seq_convert(seq, len);
	mem_alnreg_v results = mem_align1_core(opts, ref->bwt, ref->bns, ref->pac, len, seq_conv, NULL);
	EasyAlignment *aligns = safe_malloc((results.n + 1) * sizeof(*aligns));

	arena_push(results.a);

	for (size_t i = 0; i < results.n; i++) {
		mem_alnreg_t *a = &results.a[i];
		interpret_align(ref, a, &aligns[i]);
	}

	aligns[results.n].offset = -1;

	free(seq_conv);
	return aligns;
}

Chain *bwa_chain(bwaidx_t *ref, mem_opt_t *opts, char *seq, const size_t len)
{
	uint8_t *seq_conv = (uint8_t *)seq_convert(seq, len);
	mem_chain_v results = mem_chain(opts, ref->bwt, ref->bns, len, seq_conv, NULL);
	Chain *chains = safe_malloc(results.n * sizeof(*chains));

	for (size_t i = 0; i < results.n; i++) {
		mem_chain_t *a = &results.a[i];
		interpret_chain(ref, a, &chains[i]);
	}

	free(seq_conv);
	return chains;
}

EasyAlignmentPairs bwa_mem_mate_sw(bwaidx_t *ref,
                                   mem_opt_t *opts,
                                   char *read1,
                                   const size_t len1,
                                   char *read2,
                                   const size_t len2,
                                   const int score_delta)
{
#define N_PES 4
	static int init = 0;
	static mem_pestat_t pes[N_PES];

	if (!init) {
		pes[0].failed = 1;
		pes[1].failed = 0;
		pes[2].failed = 1;
		pes[3].failed = 1;

		for (int i = 0; i < N_PES; i++) {
			pes[i].low = -35;
			pes[i].high = 500;
			pes[i].avg = 200.0;
			pes[i].std = 100.0;
		}

		init = 1;
	}
#undef N_PES

	char *seq1_conv = seq_convert(read1, len1);
	char *seq2_conv = seq_convert(read2, len2);

	mem_alnreg_v results1 = mem_align1_core(opts, ref->bwt, ref->bns, ref->pac, len1, seq1_conv, NULL);
	mem_alnreg_v results2 = mem_align1_core(opts, ref->bwt, ref->bns, ref->pac, len2, seq2_conv, NULL);

	EasyAlignment *aligns1 = safe_malloc(results1.n * sizeof(*aligns1));
	EasyAlignment *aligns2 = safe_malloc(results2.n * sizeof(*aligns2));

	int best_score1 = 0;
	int best_score2 = 0;

	for (size_t i = 0; i < results1.n; i++) {
		mem_alnreg_t *a = &results1.a[i];
		EasyAlignment *p = &aligns1[i];
		interpret_align(ref, a, p);

		if (p->score > best_score1)
			best_score1 = p->score;
	}

	for (size_t i = 0; i < results2.n; i++) {
		mem_alnreg_t *a = &results2.a[i];
		EasyAlignment *p = &aligns2[i];
		interpret_align(ref, a, p);

		if (p->score > best_score2)
			best_score2 = p->score;
	}

	int num = 0;
	for (size_t i = 0; i < results2.n && num < 50; i++) {
		if (aligns2[i].score >= best_score2 - score_delta) {
			++num;
			mem_matesw(opts, ref->bns, ref->pac, &pes[0], aligns2[i].chained_hit, len1, (uint8_t *)seq1_conv, &results1);
		}
	}

	aligns1 = safe_realloc(aligns1, results1.n * sizeof(*aligns1));
	for (size_t i = 0; i < results1.n; i++) {
		mem_alnreg_t *a = &results1.a[i];
		interpret_align(ref, a, &aligns1[i]);
	}

	num = 0;
	for (size_t i = 0; i < results1.n && num < 50; i++) {
		if (aligns1[i].score >= best_score1 - score_delta) {
			++num;
			mem_matesw(opts, ref->bns, ref->pac, &pes[0], aligns1[i].chained_hit, len2, (uint8_t *)seq2_conv, &results2);
		}
	}

	arena_push(results1.a);
	arena_push(results2.a);

	aligns2 = safe_realloc(aligns2, results2.n * sizeof(*aligns2));
	for (size_t i = 0; i < results2.n; i++) {
		mem_alnreg_t *a = &results2.a[i];
		interpret_align(ref, a, &aligns2[i]);
	}

	free(seq1_conv);
	free(seq2_conv);
	arena_push(aligns1);
	arena_push(aligns2);
	return (EasyAlignmentPairs){aligns1, aligns2, results1.n, results2.n};
}

void bwa_smith_waterman(bwaidx_t *ref, mem_opt_t *opts, char *seq, const size_t len, mem_alnreg_t *aln, SingleReadAlignment *res)
{
	char *seq_conv = seq_convert(seq, len);
	mem_aln_t result = mem_reg2aln(opts, ref->bns, ref->pac, len, seq_conv, aln);

	arena_push(result.cigar);
	arena_push(result.XA);

	interpret_single_read_alignment(ref, &result, res);
	free(seq_conv);
}

void interpret_align(bwaidx_t *ref, mem_alnreg_t *caln, EasyAlignment *res)
{
	bntseq_t *contigs = ref->bns;
	const int contig_id = caln->rid;
	bntann1_t *contig = &contigs->anns[contig_id];

	if (caln->rb < contigs->l_pac) {
		res->offset = caln->rb - contig->offset;
		res->rev = 0;
	} else {
		res->offset = contigs->l_pac*2 - 1 - caln->rb - contig->offset;
		res->rev = 1;
	}

	if (caln->re < contigs->l_pac) {
		res->aln_end = caln->re - contig->offset;
	} else {
		res->aln_end = contigs->l_pac*2 - 1 - caln->re - contig->offset;
	}

	res->contig = contig->name;
	res->sec = (caln->secondary >= 0 || caln->secondary_all > 0) ? 1 : 0;
	res->chained_hit = caln;
	res->score = caln->score;
	res->read_s = caln->qb;
	res->read_e = caln->qe;
}

void interpret_chain(bwaidx_t *ref, mem_chain_t *chn, Chain *res)
{
	bntseq_t *contigs = ref->bns;
	const int contig_id = chn->rid;
	bntann1_t *contig = &contigs->anns[contig_id];

	if (chn->pos < contigs->l_pac) {
		res->offset = chn->pos - contig->offset;
		res->rev = 0;
	} else {
		res->offset = contigs->l_pac*2 - 1 - chn->pos - contig->offset;
		res->rev = 1;
	}

	res->contig = contig->name;
	res->chain = chn;
}

void interpret_single_read_alignment(bwaidx_t *ref, mem_aln_t *aln, SingleReadAlignment *res)
{
	const int fixed_flags = ((mem_aln_compact_flag_t *)aln)->flag2;
	bntseq_t *contigs = ref->bns;
	const int contig_id = aln->rid;
	bntann1_t *contig = &contigs->anns[contig_id];

	res->pos = aln->pos;
	res->chrom = contig->name;
	res->cigar = aln->cigar;
	res->n_cigar = aln->n_cigar;

	res->alt = (fixed_flags & 0x2) >> 1;
	res->mapq = (fixed_flags & 0x3fc) >> 2;//(fixed_flags & 0x2c) >> 2;
	res->rev = fixed_flags & 0x1;
	res->score = aln->score;
	res->sub = aln->sub;
	res->edit_dist = fixed_flags >> 10;
	res->alt_sc = aln->alt_sc;
	//res->raw = aln;
}

