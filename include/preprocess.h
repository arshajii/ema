#ifndef PREPROCESS_H
#define PREPROCESS_H

/* corrects barcodes and generate new FASTQs */
void preprocess_fastqs(const char *fq1, const char *fq2, const char *wl_path, const int n_buckets, const char *counts);

/* performs initial barcode count */
void count_barcodes(BarcodeDict *bcdict, FILE *fq, const int interleaved);

/* barcode min confidence */
#define BC_CONF_THRESH 0.975

#endif /* PREPROCESS_H */

