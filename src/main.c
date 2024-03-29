#define __STDC_WANT_LIB_EXT2__ 1
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>

#include <omp.h>

#include "samrecord.h"
#include "align.h"
#include "barcodes.h"
#include "util.h"
#include "main.h"

// preprocess
#include "cpp/main.h"
#include "cpp/count.h"
#include "cpp/correct.h"

int num_threads_per_file = 1;
int num_threads_for_files = 1;
char *rg = "@RG\tID:rg1\tSM:sample1";
char *bx_index = "1";
char **pg_argv = NULL;
int pg_argc = 0;

int BC_LEN;
PlatformProfile *tech;

#define MAX_CHROM_NAME_LEN 256
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
	return strstr(rg, "@RG\t") == &rg[0] && strstr(rg, "\tID:") != NULL;
}

static void print_help_and_exit(const char *argv0, int error)
{
#define P(...) fprintf(out, __VA_ARGS__)

	FILE *out = error ? stderr : stdout;
	P("usage: %s <count|preproc|align|help> [options]\n", argv0);
	P("\n");
	P("count: perform preliminary barcode count (takes interleaved FASTQ via stdin)\n");
	P("  -w <whitelist path>: specify barcode whitelist [required, unless using -p option]\n");
	P("  -o <output prefix>: specify output prefix [required]\n");
	P("  -p: using haplotag barcodes\n");
	P("\n");
	P("preproc: preprocess barcoded FASTQ files (takes interleaved FASTQ via stdin)\n");
	P("  -w <whitelist path>: specify whitelist [required, unless using -p option]\n");
	P("  -n <num buckets>: number of barcode buckets to make [500]\n");
	P("  -h: apply Hamming-2 correction [off]\n");
	P("  -o: <output directory> specify output directory [required]\n");
	P("  -b: output BX:Z-formatted FASTQs [off]\n");
	P("  -t <threads>: set number of threads [1]\n");
	P("  -p: using haplotag barcodes\n");
	P("  all other arguments: list of all output prefixes generated by count stage\n");
	P("\n");
	P("align: choose best alignments based on barcodes\n");
	P("  -1 <FASTQ1 path>: first (preprocessed and sorted) FASTQ file [none]\n");
	P("  -2 <FASTQ2 path>: second (preprocessed and sorted) FASTQ file [none]\n");
	P("  -s <EMA-FASTQ path>: specify special FASTQ path [none]\n");
	P("  -x: multi-input mode; takes input files after flags and spawns a thread for each [off]\n");
	P("  -r <FASTA path>: indexed reference [required]\n");
	P("  -o <SAM file>: output SAM file [stdout]\n");
	P("  -R <RG string>: full read group string (e.g. '@RG\\tID:foo\\tSM:bar') [none]\n");
	P("  -d: apply fragment read density optimization [off]\n");
	P("  -p <platform>: sequencing platform (one of '10x', 'tru', 'haplotag', 'dbs', 'cpt') [10x]\n");
	P("  -i <index>: index to follow 'BX' tag in SAM output [1]");
	P("  -t <threads>: set number of threads [1]\n");
	P("  all other arguments (only for -x): list of all preprocessed inputs\n");
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

	pg_argv = argv;
	pg_argc = argc;

	const char *mode = argv[1];

	if (EQ(mode, "preproc")) {
		cppinit();

		char *wl = NULL;
		int nbuckets = 500;
		int do_h2 = 0;
		char *out = NULL;
		int t = 1;
		char c;
		char do_bx_format = 0;
		int is_haplotag = 0;

		while ((c = getopt(argc-1, &argv[1], "w:n:hbpo:t:")) != -1) {
			switch (c) {
			case 'w':
				wl = strdup(optarg);
				break;
			case 'n':
				nbuckets = atoi(optarg);
				break;
			case 'h':
				do_h2 = 1;
				break;
			case 'o':
				out = strdup(optarg);
				break;
			case 't':
				t = atoi(optarg);
				break;
			case 'b':
				do_bx_format = 1;
				break;
			case 'p':
				is_haplotag = 1;
				break;
			default:
				print_help_and_exit(argv0, 1);
			}
		}

		if (wl == NULL && !is_haplotag) {
			fprintf(stderr, "error: specify barcode whitelist with -w\n");
			exit(EXIT_FAILURE);
		}

		if (out == NULL) {
			fprintf(stderr, "error: specify output directory with -o\n");
			exit(EXIT_FAILURE);
		}

		const size_t n_inputs = argc - optind - 1;

		if (n_inputs == 0) {
			fprintf(stderr, "warning: no input files specified; nothing to do\n");
			exit(EXIT_SUCCESS);
		}

		const char **inputs = safe_malloc(n_inputs * sizeof(*inputs));

		for (int i = optind + 1; i < argc; i++) {
			const int j = i - (optind + 1);  // w.r.t. `inputs`
			inputs[j] = strdup(argv[i]);
		}

		correct(wl, inputs, n_inputs, out, do_h2, 10 * MB, do_bx_format, t, nbuckets, is_haplotag);
		return EXIT_SUCCESS;
	}

	if (EQ(mode, "count")) {
		cppinit();

		char *wl = NULL;
		char *out = NULL;
		char c;
		int is_haplotag = 0;

		while ((c = getopt(argc-1, &argv[1], "w:o:p")) != -1) {
			switch (c) {
			case 'w':
				wl = strdup(optarg);
				break;
			case 'o':
				out = strdup(optarg);
				break;
			case 'p':
				is_haplotag = 1;
				break;
			default:
				print_help_and_exit(argv0, 1);
			}
		}

		if (wl == NULL && !is_haplotag) {
			fprintf(stderr, "error: specify barcode whitelist with -w\n");
			exit(EXIT_FAILURE);
		}

		if (out == NULL) {
			fprintf(stderr, "error: specify output prefix with -o\n");
			exit(EXIT_FAILURE);
		}

		count(wl, out, 1 * GB, is_haplotag);
		return EXIT_SUCCESS;
	}

	if (EQ(mode, "align")) {
		char *ref = NULL;
		char *fq1 = NULL;
		char *fq2 = NULL;
		char *fqx = NULL;
		char *fai = NULL;
		char *out = NULL;
		int apply_opt = 0;
		int multi_input = 0;
		char *platform = "10x";
		int t = 1;
		char c;

		while ((c = getopt(argc-1, &argv[1], "r:1:2:s:xo:R:dp:i:t:")) != -1) {
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
			case 's':
				fqx = strdup(optarg);
				break;
			case 'x':
				multi_input = 1;
				break;
			case 'o':
				out = strdup(optarg);
				break;
			case 'R':
				rg  = escape(strdup(optarg));
				break;
			case 'd':
				apply_opt = 1;
				break;
			case 'p':
				platform = strdup(optarg);
				break;
			case 'i':
				bx_index = strdup(optarg);
				break;
			case 't':
				t = atoi(optarg);
				break;
			default:
				print_help_and_exit(argv0, 1);
			}
		}

		if (multi_input + (fqx != NULL) + (fq1 != NULL || fq2 != NULL) != 1) {
			fprintf(stderr, "error: must specify *exactly one* of -1/-2, -s or -x\n");
			exit(EXIT_FAILURE);
		}

		if (fq1 == NULL && fq2 != NULL) {
			fprintf(stderr, "error: cannot specify -2 without -1\n");
			exit(EXIT_FAILURE);
		}

		if (ref == NULL) {
			fprintf(stderr, "error: specify reference FASTA with -r\n");
			exit(EXIT_FAILURE);
		}

		if (rg != NULL && !validate_read_group(rg)) {
			fprintf(stderr, "error: malformed read group: '%s'\n", rg);
			exit(EXIT_FAILURE);
		}

		if ((tech = get_platform_profile_by_name(platform)) == NULL) {
			fprintf(stderr, "error: invalid platform name: '%s'\n", platform);
			exit(EXIT_FAILURE);
		}

		FILE *fq1_file = NULL;
		FILE *fq2_file = NULL;
		FILE *fqx_file = NULL;
		FILE *out_file = (out == NULL ? stdout : fopen(out, "w"));
		BC_LEN = tech->bc_len;

		if (!out_file) {
			IOERROR(out);
		}

		if (fqx != NULL) {
			fqx_file = fopen(fqx, "r");

			if (!fqx_file) {
				IOERROR(fqx);
			}
		} else if (fq1 != NULL) {
			fq1_file = fopen(fq1, "r");

			if (!fq1_file) {
				IOERROR(fq1);
			}

			if (fq2 != NULL) {
				fq2_file = fopen(fq2, "r");

				if (!fq2_file) {
					IOERROR(fq2);
				}
			} else {
				fq2_file = fq1_file;
			}
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
		bwa_init(ref);
		write_sam_header(out_file);

		if (multi_input) {
			num_threads_for_files = t;
			num_threads_per_file  = 1;

			const size_t n_inputs = argc - optind - 1;

			if (n_inputs == 0) {
				fprintf(stderr, "warning: no input files specified; nothing to do\n");
				exit(EXIT_SUCCESS);
			}

			FILE **inputs  = safe_malloc(n_inputs * sizeof(*inputs));

			for (int i = optind + 1; i < argc; i++) {
				const int j = i - (optind + 1);  // w.r.t. `inputs` and `outputs`

				const char *filename = strdup(argv[i]);
				inputs[j] = fopen(filename, "r");

				if (!inputs[j]) {
					IOERROR(filename);
				}
			}

			omp_set_nested(1);
			omp_lock_t out_lock;
			omp_init_lock(&out_lock);

			#pragma omp parallel for num_threads(num_threads_for_files)
			for (size_t i = 0; i < n_inputs; i++) {
				find_clouds_and_align(NULL, NULL, inputs[i], out_file, apply_opt, NULL, &out_lock);
				fclose(inputs[i]);
			}

			omp_destroy_lock(&out_lock);
			free(inputs);
		} else {
			num_threads_for_files = 1;
			num_threads_per_file  = t;

			find_clouds_and_align(fq1_file, fq2_file, fqx_file, out_file, apply_opt, NULL, NULL);

			if (fq1_file) fclose(fq1_file);
			if (fq2_file && fq2_file != fq1_file) fclose(fq2_file);
			if (fqx_file) fclose(fqx_file);
		}

		fclose(out_file);
		bwa_dealloc();
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

