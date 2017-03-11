#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "util.h"

bc_t encode_bc(const char *bc)
{
#define BC_ADD_BASE(x) (encoded_bc |= (x))

	bc_t encoded_bc = 0UL;
	char *base = (char *)&bc[15];
	for (int i = 0; i < 16; i++) {
		encoded_bc <<= 2;
		switch (*base--) {
		case 'A': case 'a': BC_ADD_BASE(0UL); break;
		case 'C': case 'c': BC_ADD_BASE(1UL); break;
		case 'G': case 'g': BC_ADD_BASE(2UL); break;
		case 'T': case 't': BC_ADD_BASE(3UL); break;
		default: assert(0); break;
		}
	}

	return encoded_bc;

#undef BC_ADD_BASE
}

size_t count_lines(FILE *f)
{
	size_t lines = 0;
	while (!feof(f)) {
		if (fgetc(f) == '\n')
			++lines;
	}
	rewind(f);
	return lines + 1;
}

void split_line(const char *str, char **out)
{
	char *p = (char *)str;
	size_t i = 0;
	while (*p) {
		out[i++] = p;
		while (*p != '\t' && *p != '\n') ++p;
		++p;
	}
	out[i] = NULL;
}

uint32_t hash_ident(const char *ident)
{
	uint32_t h = 0;
	for (size_t i = 0; ident[i] != '\0'; i++) {
		h = 31*h + ident[i];
	}
	return h;
}

void normalize_log_probs(double *p, const size_t n)
{
#define EPSILON 1e-50

	if (n == 1) {
		p[0] = 1.0;
		return;
	}

	const double thresh = log(EPSILON) - log(n);
	double p_max = p[0];
	for (size_t i = 1; i < n; i++) {
		const double p_cur = p[i];
		if (p_cur > p_max)
			p_max = p_cur;
	}

	double total = 0;
	for (size_t i = 0; i < n; i++) {
		p[i] -= p_max;

		if (p[i] < thresh)
			p[i] = 0;
		else
			p[i] = exp(p[i]);

		total += p[i];
	}

	for (size_t i = 0; i < n; i++) {
		p[i] /= total;
	}

#undef EPSILON
}

void *safe_malloc(const size_t n)
{
	void *p = malloc(n);
	assert(n == 0 || p != NULL);
	return p;
}

void *safe_calloc(size_t num, size_t size)
{
	void *p = calloc(num, size);
	assert(p != NULL || num == 0 || size == 0);
	return p;
}

