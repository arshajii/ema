#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "main.h"
#include "util.h"
#include "barcodes.h"

#define POW_2_24 (1UL << 24)
#define HI24(bc) (((bc) & 0xFFFFFF00) >> 8)

static int bcinfo_cmp(const void *v1, const void *v2) {
	const BarcodeInfo *b1 = (BarcodeInfo *)v1;
	const BarcodeInfo *b2 = (BarcodeInfo *)v2;
	const bc_t bc1 = b1->bc;
	const bc_t bc2 = b2->bc;
	return (bc1 > bc2) - (bc1 < bc2);
}

void wl_read(BarcodeDict *bcdict, const char *whitelist_path)
{
	char buf[1024];

	FILE *wl_file = fopen(whitelist_path, "r");

	if (wl_file == NULL) {
		IOERROR(whitelist_path);
	}

	const size_t n_lines = count_lines(wl_file);
	size_t wl_size = 0;

	BarcodeInfo *whitelist = safe_malloc(n_lines * sizeof(*whitelist));

	while (fgets(buf, sizeof(buf), wl_file)) {
		if (strchr(buf, '#'))
			continue;

		whitelist[wl_size].bc = encode_bc(buf);
		whitelist[wl_size].count = 0;
		++wl_size;
	}

	fclose(wl_file);

	qsort(whitelist, wl_size, sizeof(*whitelist), bcinfo_cmp);

	uint32_t *wl_jumpgate = safe_malloc(POW_2_24 * sizeof(*wl_jumpgate));

	wl_jumpgate[0] = 0;
	uint32_t last_hi = 0;
	for (size_t i = 0; i < wl_size; i++) {
		const bc_t bc = whitelist[i].bc;

		const uint32_t hi = HI24(bc);

		if (hi != last_hi) {
			assert(hi > last_hi);

			for (size_t j = (last_hi + 1); j <= hi; j++)
				wl_jumpgate[j] = i;

			last_hi = hi;
		}
	}

	if (last_hi != 0xFFFFFF) {
		for (size_t j = (last_hi + 1); j < POW_2_24; j++)
			wl_jumpgate[j] = wl_size;
	}

	bcdict->jumpgate = wl_jumpgate;
	bcdict->entries = whitelist;
	bcdict->size = wl_size;
	bcdict->unfound = 0;
}

void wl_dealloc(BarcodeDict *bcdict)
{
	free(bcdict->jumpgate);
	free(bcdict->entries);
}

BarcodeInfo *wl_lookup(BarcodeDict *bcdict, bc_t key)
{
	const uint32_t *wl_jumpgate = bcdict->jumpgate;
	const BarcodeInfo *whitelist = bcdict->entries;
	const size_t wl_size = bcdict->size;

	const uint32_t bc_hi = HI24(key);

	const uint32_t lo = wl_jumpgate[bc_hi];

	if (lo == wl_size) {
		return NULL;
	}

	const uint32_t hi = (bc_hi == 0xFFFFFF ? wl_size : wl_jumpgate[bc_hi + 1]);

	if (lo == hi) {
		return NULL;
	}

	assert(hi > lo);
	BarcodeInfo findme = (BarcodeInfo){ .bc = key };
	BarcodeInfo *target = bsearch(&findme, &whitelist[lo], (hi - lo), sizeof(*whitelist), bcinfo_cmp);
	return target;
}

int wl_increment(BarcodeDict *bcdict, bc_t key)
{
	BarcodeInfo *bcinfo = wl_lookup(bcdict, key);

	if (bcinfo != NULL) {
		++(bcinfo->count);
		return 1;
	} else {
		++(bcdict->unfound);
		return 0;
	}
}

void wl_compute_priors(BarcodeDict *bcdict)
{
	uint64_t total = 0;
	BarcodeInfo *whitelist = bcdict->entries;
	const size_t size = bcdict->size;

	for (size_t i = 0; i < size; i++) {
		total += whitelist[i].count + 1;  // +1 as a pseudocount
	}

	for (size_t i = 0; i < size; i++) {
		whitelist[i].prior = (whitelist[i].count + 1.0)/total;
	}
}

int wl_get_bucket(BarcodeDict *bcdict, BarcodeInfo *bc, const int n_buckets)
{
	return ((bc - bcdict->entries)*n_buckets)/(bcdict->size);
}

void wl_serialize(BarcodeDict *bcdict, FILE *out)
{
	uint32_t *jumpgate = bcdict->jumpgate;
	BarcodeInfo *entries = bcdict->entries;
	const size_t size = bcdict->size;

	for (size_t i = 0; i < POW_2_24; i++) {
		serialize_uint32(out, jumpgate[i]);
	}

	serialize_uint64(out, size);

	for (size_t i = 0; i < size; i++) {
		serialize_uint32(out, entries[i].bc);
		serialize_uint32(out, entries[i].count);
	}
}

void wl_deserialize(BarcodeDict *bcdict, FILE *in)
{
	uint32_t *jumpgate = safe_malloc(POW_2_24 * sizeof(*jumpgate));

	for (size_t i = 0; i < POW_2_24; i++) {
		jumpgate[i] = read_uint32(in);
	}

	const size_t size = read_uint64(in);
	BarcodeInfo *entries = safe_malloc(size * sizeof(*entries));

	for (size_t i = 0; i < size; i++) {
		entries[i].bc = read_uint32(in);
		entries[i].count = read_uint32(in);
	}

	bcdict->jumpgate = jumpgate;
	bcdict->entries = entries;
	bcdict->size = size;
	bcdict->unfound = 0;
}

