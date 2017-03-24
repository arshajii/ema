#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include "align.h"
#include "barcodes.h"
#include "util.h"
#include "preprocess.h"

#define PREPROC_EXT                   ".preproc.fastq"
#define NO_BC_EXT                     ".no_bc.fastq"
#define INTERLEAVED_PREPROC_EXT(mate) ("_" mate PREPROC_EXT)
#define INTERLEAVED_NO_BC_EXT(mate)   ("_" mate NO_BC_EXT)

static void trim_after_space(char *s)
{
	for (size_t i = 0; s[i] != '\0'; i++) {
		if (isspace(s[i])) {
			s[i] = '\0';
			break;
		}
	}
}

static FILE *open_fastq_with_ext(char *fq, const char *s)
{
	char *ext = strrchr(fq, '.');

	if (ext == NULL)
		ext = &fq[strlen(fq)];

	strcpy(ext, s);

	return fopen(fq, "w");
}

/*
 * This is essentially a translation of Long Ranger's barcode correction scheme.
 */
static int correct_barcode(char *barcode, char *barcode_qual, BarcodeDict *wl)
{
#define ILLUMINA_QUAL_OFFSET 33

	uint8_t quals[BC_LEN];

	int n_count = 0;
	for (size_t i = 0; i < BC_LEN; i++) {
		barcode[i] = toupper(barcode[i]);
		quals[i] = barcode_qual[i] - ILLUMINA_QUAL_OFFSET;

		if (!IS_ACGT(barcode[i])) {
			++n_count;
		}
	}

	const bc_t bc0 = ((n_count == 0) ? encode_bc(barcode) : 0);
	BarcodeInfo *bcinfo = ((n_count == 0) ? wl_lookup(wl, bc0) : NULL);

#define CAND_BUF_SIZE (120 * 4 * 4)
	struct bc_string { char bc_str[BC_LEN]; };
	struct bc_string bc_cands[CAND_BUF_SIZE];
	double bc_cand_probs[CAND_BUF_SIZE];
	size_t n_cands = 0;
#undef CAND_BUF_SIZE

	if (bcinfo == NULL) {
		if (n_count > 1) {
			return 0;
		}

		/* examine Hamming-1 neighbors */
		for (size_t i = 0; i < BC_LEN; i++) {
			const char prev = barcode[i];

			if (n_count > 0 && IS_ACGT(prev))
				continue;

			for (size_t j = 0; j < 4; j++) {
				const char new = "ACGT"[j];

				if (new == prev)
					continue;

				barcode[i] = new;
				const bc_t bc = encode_bc(barcode);
				bcinfo = wl_lookup(wl, bc);

				if (bcinfo != NULL) {
					const double prior = bcinfo->prior;
					const double edit_log_prob = MIN(33.0, (double)(quals[i]));
					const double edit_prob = pow(10.0, (-edit_log_prob/10.0));
					const double p = prior*edit_prob;

					memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
					bc_cand_probs[n_cands] = p;
					++n_cands;
				}
			}

			barcode[i] = prev;
		}
	} else {
		assert(n_count == 0);
		memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
		bc_cand_probs[n_cands] = bcinfo->prior;
		++n_cands;

		/* examine Hamming-2 neighbors */
		for (size_t i1 = 0; i1 < BC_LEN; i1++) {
			const char prev1 = barcode[i1];

			for (size_t j1 = 0; j1 < 4; j1++) {
				const char new1 = "ACGT"[j1];

				if (new1 == prev1)
					continue;

				barcode[i1] = new1;

				for (size_t i2 = i1 + 1; i2 < BC_LEN; i2++) {
					const char prev2 = barcode[i2];

					for (size_t j2 = 0; j2 < 4; j2++) {
						const char new2 = "ACGT"[j2];

						if (new2 == prev2)
							continue;

						barcode[i2] = new2;
						const bc_t bc = encode_bc(barcode);
						bcinfo = wl_lookup(wl, bc);

						if (bcinfo != NULL) {
							const double prior = bcinfo->prior;

							const double e1 = MAX(3.0, (double)(quals[i1]) - 1.0);
							const double e2 = MAX(3.0, (double)(quals[i2]) - 1.0);

							const double edit1_log_prob = MIN(33.0, e1);
							const double edit2_log_prob = MIN(33.0, e2);

							const double edit_prob = pow(10.0, (-edit1_log_prob/10.0)) * pow(10.0, (-edit2_log_prob/10.0));
							const double p = prior*edit_prob;

							memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
							bc_cand_probs[n_cands] = p;
							++n_cands;
						}
					}

					barcode[i2] = prev2;
				}

				barcode[i1] = prev1;
			}
		}
	}

	if (n_cands > 0) {
		double total = bc_cand_probs[0];
		size_t max = 0;

		for (size_t i = 1; i < n_cands; i++) {
			total += bc_cand_probs[i];

			if (bc_cand_probs[i] > bc_cand_probs[max])
				max = i;
		}

		if (bc_cand_probs[max]/total > BC_CONF_THRESH) {
			memcpy(barcode, bc_cands[max].bc_str, BC_LEN);
			return 1;
		}
	}

	return 0;

#undef ILLUMINA_QUAL_OFFSET
}

void preprocess_fastqs(const char *fq1, const char *fq2, const char *wl_path)
{
#define BUF_SIZE 1024

	char namebuf[BUF_SIZE];

	assert(strlen(fq1) < (sizeof(namebuf) - strlen(PREPROC_EXT) - strlen(NO_BC_EXT) - 1));
	assert(strlen(fq2) < (sizeof(namebuf) - strlen(PREPROC_EXT) - strlen(NO_BC_EXT) - 1));

	FILE *fq1_file = fopen(fq1, "r");

	if (fq1_file == NULL)
		IOERROR(fq1);

	FILE *fq2_file;

	const int interleaved = (strcmp(fq1, fq2) == 0);

	if (!interleaved) {
		fq2_file = fopen(fq2, "r");

		if (fq2_file == NULL) {
			IOERROR(fq2);
		}
	} else {
		fq2_file = fq1_file;
	}

	strcpy(namebuf, fq1);
	FILE *fq1_preproc_file = open_fastq_with_ext(namebuf, interleaved ? INTERLEAVED_PREPROC_EXT("1") : PREPROC_EXT);
	if (fq1_preproc_file == NULL) {
		IOERROR(namebuf);
	}

	strcpy(namebuf, fq1);
	FILE *fq1_nobc_file = open_fastq_with_ext(namebuf, interleaved ? INTERLEAVED_NO_BC_EXT("1") : NO_BC_EXT);
	if (fq1_nobc_file == NULL) {
		IOERROR(namebuf);
	}

	strcpy(namebuf, fq2);
	FILE *fq2_preproc_file = open_fastq_with_ext(namebuf, interleaved ? INTERLEAVED_PREPROC_EXT("2") : PREPROC_EXT);
	if (fq2_preproc_file == NULL) {
		IOERROR(namebuf);
	}

	strcpy(namebuf, fq2);
	FILE *fq2_nobc_file = open_fastq_with_ext(namebuf, interleaved ? INTERLEAVED_NO_BC_EXT("1") : NO_BC_EXT);
	if (fq2_nobc_file == NULL) {
		IOERROR(namebuf);
	}

	BarcodeDict wl;
	wl_read(&wl, wl_path);

	char barcode[BC_LEN+1];
	char barcode_qual[BC_LEN+1];
	char id1[BUF_SIZE];
	char read1[BUF_SIZE];
	char sep1[BUF_SIZE];
	char qual1[BUF_SIZE];
	char id2[BUF_SIZE];
	char read2[BUF_SIZE];
	char sep2[BUF_SIZE];
	char qual2[BUF_SIZE];

	barcode[BC_LEN] = '\0';
	barcode_qual[BC_LEN] = '\0';

	/* perform initial barcode counting */
	while (fgets(id1, BUF_SIZE, fq1_file)) {
		assert(fgets(read1, BUF_SIZE, fq1_file));
		assert(fgets(sep1, BUF_SIZE, fq1_file));
		assert(fgets(qual1, BUF_SIZE, fq1_file));

		for (size_t i = 0; i < BC_LEN; i++) {
			if (IS_ACGT(read1[i])) {
				barcode[i] = read1[i];
			} else {
				barcode[0] = '\0';
				break;
			}
		}

		if (barcode[0] != '\0') {
			const bc_t bc = encode_bc(barcode);
			wl_increment(&wl, bc);
		}

		if (interleaved) {
			assert(fgets(id1, BUF_SIZE, fq2_file));
			assert(fgets(read1, BUF_SIZE, fq2_file));
			assert(fgets(sep1, BUF_SIZE, fq2_file));
			assert(fgets(qual1, BUF_SIZE, fq2_file));
		}
	}

	wl_compute_priors(&wl);
	rewind(fq1_file);

	while (fgets(id1, BUF_SIZE, fq1_file)) {
		assert(fgets(read1, BUF_SIZE, fq1_file));
		assert(fgets(sep1, BUF_SIZE, fq1_file));
		assert(fgets(qual1, BUF_SIZE, fq1_file));

		assert(fgets(id2, BUF_SIZE, fq2_file));
		assert(fgets(read2, BUF_SIZE, fq2_file));
		assert(fgets(sep2, BUF_SIZE, fq2_file));
		assert(fgets(qual2, BUF_SIZE, fq2_file));

		trim_after_space(id1);
		trim_after_space(id2);

		size_t id1_len = strlen(id1);
		size_t id2_len = strlen(id2);
		size_t read1_len = strlen(read1);
		size_t read2_len = strlen(read2);
		size_t qual1_len = strlen(qual1);
		size_t qual2_len = strlen(qual2);

		assert(id1_len == id2_len);
		assert(read1_len == read2_len && qual1_len == qual2_len && read1_len == qual1_len);
		assert(read1_len > (BC_LEN + MATE1_TRIM));

		/* get rid of any /1 or /2 */
		if (id1[id1_len - 1] == '1' && id2[id2_len - 1] == '2') {
			assert(id1[id1_len - 2] == '/' && id2[id2_len - 2] == '/');

			id1[id1_len - 2] = '\0';
			id2[id2_len - 2] = '\0';

			id1_len -= 2;
			id2_len -= 2;
		}

		assert(strcmp(id1, id2) == 0);

		memcpy(barcode, read1, BC_LEN);
		memcpy(barcode_qual, qual1, BC_LEN);

		/* trim off barcode and first MATE1_TRIM bases */
		memmove(read1, &read1[BC_LEN + MATE1_TRIM], read1_len - (BC_LEN + MATE1_TRIM) + 1);
		memmove(qual1, &qual1[BC_LEN + MATE1_TRIM], qual1_len - (BC_LEN + MATE1_TRIM) + 1);

		const int good_barcode = correct_barcode(barcode, barcode_qual, &wl);
		FILE *dest1 = good_barcode ? fq1_preproc_file : fq1_nobc_file;
		FILE *dest2 = good_barcode ? fq2_preproc_file : fq2_nobc_file;

		fputs(id1, dest1);
		fputs(id2, dest2);

		if (good_barcode) {
			fputs(":", dest1);
			fputs(barcode, dest1);

			fputs(":", dest2);
			fputs(barcode, dest2);
		}

		fputs("\n", dest1);
		fputs("\n", dest2);

		fputs(read1, dest1);
		fputs(read2, dest2);

		fputs(sep1, dest1);
		fputs(sep2, dest2);

		fputs(qual1, dest1);
		fputs(qual2, dest2);
	}

	wl_dealloc(&wl);
	fclose(fq1_file);
	fclose(fq2_file);
	fclose(fq1_preproc_file);
	fclose(fq2_preproc_file);
	fclose(fq1_nobc_file);
	fclose(fq2_nobc_file);

#undef BUF_SIZE
}

