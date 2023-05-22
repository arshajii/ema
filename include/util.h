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

#define SIZE(a) (sizeof(a)/sizeof((a)[0]))

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define IS_ACGT(c) ((c) == 'A' || (c) == 'C' || (c) == 'G' || (c) == 'T')

typedef uint32_t bc_t;
typedef uint32_t chrom_t;

void copy_until_space(char *dest, char **src);
char *escape(char *s);

bc_t encode_bc_default(const char *bc);
bc_t encode_bc_haptag(const char *bc);
bc_t encode_bc(const char *bc, const int is_haplotag);
void decode_bc_default(bc_t bc, char *out);
void decode_bc_haptag(bc_t bc, char *out);
void decode_bc(bc_t bc, char *out, const int is_haplotag);
size_t count_lines(FILE *f);
size_t trim_after_space(char *s);
uint32_t hash_ident(const char *ident);

void normalize_log_probs(double *p, const size_t n);

void *safe_malloc(const size_t n);
void *safe_calloc(size_t num, size_t size);
void *safe_realloc(void *p, size_t n);

void serialize_uint64(FILE *out, const uint64_t x);
void serialize_uint32(FILE *out, const uint32_t x);
void serialize_uint8(FILE *out, const uint8_t x);
uint64_t read_uint64(FILE *in);
uint32_t read_uint32(FILE *in);
uint8_t read_uint8(FILE *in);

#endif /* UTIL_H */

