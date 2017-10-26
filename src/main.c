#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include "samrecord.h"
#include "align.h"
#include "barcodes.h"
#include "preprocess.h"
#include "util.h"
#include "main.h"

#define MAX_CHROM_NAME_LEN 64
static struct { char chrom_name[MAX_CHROM_NAME_LEN]; } *chroms;

char *chrom_lookup(const chrom_t chrom)
{
	return chroms[chrom].chrom_name;
}

chrom_t chrom_index(const char *chrom)
{
	size_t len = 0;
	while (!isspace(chrom[len]) && chrom[len] != '\0')
		++len;

	for (size_t i = 0; chroms[i].chrom_name[0] != '\0'; i++) {
		if (strncmp(chrom, chroms[i].chrom_name, len) == 0) {
			return i;
		}
	}

	assert(0);
	return 0;
}

void read_fai(FILE *fai_file)
{
	const size_t lines = count_lines(fai_file);
	chroms = safe_malloc(lines * sizeof(*chroms));
	size_t i = 0;

	while (fgets(chroms[i].chrom_name, MAX_CHROM_NAME_LEN, fai_file)) {
		size_t j = 0;
		while (!isspace(chroms[i].chrom_name[j]))
			++j;
		chroms[i++].chrom_name[j] = '\0';

	}
	chroms[i].chrom_name[0] = '\0';
}

static int validate_read_group(const char *rg)
{
	return strstr(rg, "@RG\t") == &rg[0] && strstr(rg, "ID:") != NULL;
}

static void print_help_and_exit(const char *argv0, int error)
{
#define P(...) fprintf(out, __VA_ARGS__)

	FILE *out = error ? stderr : stdout;
	P("usage: %s <preproc|align|help> [options]\n", argv0);
	P("\n");
	P("preproc: preprocess barcoded FASTQ files\n");
	P("  -1 <fastq1 path>: specify first FASTQ file [required]\n");
	P("  -2 <fastq2 path>: specify second FASTQ file [none]\n");
	P("  -w <whitelist path>: specify whitelist [required]\n");
	P("  -n <num buckets>: number of barcode buckets to make [20]\n");
	P("  -c <counts file>: specify preexisting barcode counts [none]\n");
	P("\n");
	P("sort: sort preprocessed FASTQs by barcode\n");
	P("  -1 <fastq1 path>: specify first FASTQ file [required]\n");
	P("  -2 <fastq2 path>: specify second FASTQ file [required]\n");
	P("\n");
	P("count: performs preliminary barcode count\n");
	P("  -1 <fastq1 path>: specify first FASTQ file [required]\n");
	P("  -w <whitelist path>: specify whitelist [required]\n");
	P("  -i: indicates FASTQ is interleaved\n");
	P("  -o <output file>: specify output file [stdout]\n");
	P("\n");
	P("align: choose best alignments based on barcodes\n");
	P("  -1 <FASTQ1 path>: first (preprocessed and sorted) FASTQ file [required]\n");
	P("  -2 <FASTQ2 path>: second (preprocessed and sorted) FASTQ file [required]\n");
	P("  -r <FASTA path>: indexed reference [required]\n");
	P("  -o <SAM file>: output SAM file [default: stdout]\n");
	P("  -R <RG string>: full read group string (e.g. $'@RG\\tID:foo\\tSM:bar') [default: none]\n");
	P("  -d: apply fragment read density optimization\n");
	P("\n");
	P("help: print this help message\n");
	exit(error ? EXIT_FAILURE : EXIT_SUCCESS);

#undef P
}

int main(const int argc, char *argv[])
{
#define EQ(s1, s2) (strcmp((s1), (s2)) == 0)

	const char *argv0 = argv[0];

	if (argc < 2) {
		fprintf(stderr, "EMA version %s\n", VERSION);
		fprintf(stderr, "note: use '%s help' for usage information.\n", argv0);
		return EXIT_SUCCESS;
	}

	const char *mode = argv[1];

	if (EQ(mode, "preproc")) {
		char *fq1 = NULL;
		char *fq2 = NULL;
		char *wl = NULL;
		int n = 20;
		char *counts = NULL;
		char c;

		while ((c = getopt(argc, argv, "w:1:2:n:c:")) != -1) {
			switch (c) {
			case '1':
				fq1 = strdup(optarg);
				break;
			case '2':
				fq2 = strdup(optarg);
				break;
			case 'w':
				wl = strdup(optarg);
				break;
			case 'n':
				n = atoi(optarg);
				break;
			case 'c':
				counts = strdup(optarg);
				break;
			default:
				print_help_and_exit(argv0, 1);
			}
		}

		if (fq1 == NULL) {
			fprintf(stderr, "error: specify paired-end FASTQs with -1 and -2 or interleaved FASTQ with -1\n");
			exit(EXIT_FAILURE);
		}

		if (fq2 == NULL)
			fq2 = fq1;

		if (EQ(fq1, "-") ^ EQ(fq2, "-")) {
			fprintf(stderr, "error: can only read interleaved FASTQ over stdin\n");
			exit(EXIT_FAILURE);
		}

		if (wl == NULL && counts == NULL) {
			fprintf(stderr, "error: specify barcode whitelist with -w OR specify barcode counts file with -c\n");
			exit(EXIT_FAILURE);
		}

		preprocess_fastqs(fq1, fq2, wl, n, counts);
		return EXIT_SUCCESS;
	}

	if (EQ(mode, "sort")) {
		char *fq1 = NULL;
		char *fq2 = NULL;
		char c;

		while ((c = getopt(argc, argv, "1:2:")) != -1) {
			switch (c) {
			case '1':
				fq1 = strdup(optarg);
				break;
			case '2':
				fq2 = strdup(optarg);
				break;
			default:
				print_help_and_exit(argv0, 1);
			}
		}

		if (fq1 == NULL || fq2 == NULL) {
			fprintf(stderr, "error: specify paired-end FASTQs with -1 and -2\n");
			exit(EXIT_FAILURE);
		}

		sort_fastq(fq1, fq2);
		return EXIT_SUCCESS;
	}

	if (EQ(mode, "count")) {
		char *fq1 = NULL;
		char *wl = NULL;
		int interleaved = 0;
		char *out = NULL;
		char c;

		while ((c = getopt(argc, argv, "1:w:io:")) != -1) {
			switch (c) {
			case '1':
				fq1 = strdup(optarg);
				break;
			case 'w':
				wl = strdup(optarg);
				break;
			case 'i':
				interleaved = 1;
				break;
			case 'o':
				out = strdup(optarg);
				break;
			default:
				print_help_and_exit(argv0, 1);
			}
		}

		if (fq1 == NULL) {
			fprintf(stderr, "error: specify first FASTQs with -1\n");
			exit(EXIT_FAILURE);
		}

		if (wl == NULL) {
			fprintf(stderr, "error: specify barcode whitelist with -w\n");
			exit(EXIT_FAILURE);
		}

		const int read_from_stdin = EQ(fq1, "-");
		FILE *fq1_file = read_from_stdin ? stdin : fopen(fq1, "r");

		if (!read_from_stdin && fq1_file == NULL) {
			IOERROR(fq1);
		}

		FILE *out_file = (out != NULL) ? fopen(out, "wb") : stdout;

		if (out != NULL && out_file == NULL) {
			IOERROR(out);
		}

		BarcodeDict wldict;
		wl_read(&wldict, wl);
		count_barcodes(&wldict, fq1_file, interleaved);
		wl_serialize(&wldict, out_file);

		if (!read_from_stdin)
			fclose(fq1_file);

		if (out != NULL)
			fclose(out_file);

		wl_dealloc(&wldict);

		return EXIT_SUCCESS;
	}

	if (EQ(mode, "align")) {
		char *ref = NULL;
		char *fq1 = NULL;
		char *fq2 = NULL;
		char *fai = NULL;
		char *out = NULL;
		char *rg  = NULL;
		int apply_opt = 0;
		char c;

		while ((c = getopt(argc, argv, "r:1:2:o:R:d")) != -1) {
			switch (c) {
			case 'r':
				ref = strdup(optarg);
				break;
			case '1':
				fq1 = strdup(optarg);
				break;
			case '2':
				fq2 = strdup(optarg);
				break;
			case 'o':
				out = strdup(optarg);
				break;
			case 'R':
				rg  = strdup(optarg);
				break;
			case 'd':
				apply_opt = 1;
				break;
			default:
				print_help_and_exit(argv0, 1);
			}
		}

		if (ref == NULL) {
			fprintf(stderr, "error: specify reference FASTA with -r\n");
			exit(EXIT_FAILURE);
		}

		if (fq1 == NULL || fq2 == NULL) {
			fprintf(stderr, "error: specify input paired-end FASTQs with -1 and -2\n");
			exit(EXIT_FAILURE);
		}

		if (rg != NULL && !validate_read_group(rg)) {
			fprintf(stderr, "error: malformed read group: '%s' (remember to escape the argument string!)\n", rg);
			exit(EXIT_FAILURE);
		}

		FILE *fq1_file = fopen(fq1, "r");

		if (!fq1_file) {
			IOERROR(fq1);
		}

		FILE *fq2_file = fopen(fq2, "r");

		if (!fq2_file) {
			IOERROR(fq2);
		}

#define FAI_EXT ".fai"
		fai = safe_malloc(strlen(ref) + strlen(FAI_EXT) + 1);
		strcpy(fai, ref);
		strcat(fai, FAI_EXT);
#undef FAI_EXT

		FILE *fai_file = fopen(fai, "r");

		if (!fai_file) {
			IOERROR(fai);
		}

		read_fai(fai_file);
		fclose(fai_file);

		FILE *out_file = (out == NULL ? stdout : fopen(out, "w"));

		if (!out_file) {
			IOERROR(out);
		}

		find_clouds_and_align(fq1_file, fq2_file, ref, out_file, rg, apply_opt);

		fclose(fq1_file);
		fclose(fq2_file);
		fclose(out_file);

		free(fai);
		free(chroms);

		return EXIT_SUCCESS;
	}

	if (EQ(mode, "help")) {
		print_help_and_exit(argv0, 0);
	}

	fprintf(stderr, "error: unrecognized mode\n");
	print_help_and_exit(argv0, 1);

#undef EQ
}

