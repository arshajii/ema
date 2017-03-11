#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define IOERROR(filename) \
	do { \
		fprintf(stderr, "error: file %s could not be opened\n", (filename)); \
		exit(EXIT_FAILURE); \
	} while (0)

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define IS_ACGT(c) ((c) == 'A' || (c) == 'C' || (c) == 'G' || (c) == 'T')

typedef uint32_t bc_t;
typedef uint8_t chrom_t;

bc_t encode_bc(const char *bc);
size_t count_lines(FILE *f);
void split_line(const char *str, char **out);
uint32_t hash_ident(const char *ident);

void normalize_log_probs(double *p, const size_t n);

void *safe_malloc(const size_t n);
void *safe_calloc(size_t num, size_t size);

#endif /* UTIL_H */

