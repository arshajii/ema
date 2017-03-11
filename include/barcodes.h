#ifndef BARCODES_H
#define BARCODES_H

#include <stdlib.h>
#include <stdint.h>
#include "main.h"
#include "util.h"

/*
 * Useful structures for computing barcode information -
 */

typedef struct {
	bc_t bc;
	uint32_t count;
	double prior;
} BarcodeInfo;

typedef struct {
	uint32_t *jumpgate;
	BarcodeInfo *entries;
	size_t size;
	uint32_t unfound;
} BarcodeDict;

void wl_read(BarcodeDict *bcdict, const char *whitelist_path);
void wl_dealloc(BarcodeDict *bcdict);
BarcodeInfo *wl_lookup(BarcodeDict *bcdict, bc_t key);
int wl_increment(BarcodeDict *bcdict, bc_t key);
void wl_compute_priors(BarcodeDict *bcdict);

#endif /* BARCODES_H */

