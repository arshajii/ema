#ifndef PREPROCESS_H
#define PREPROCESS_H

/* corrects barcodes and generate new FASTQs */
void preprocess_fastqs(const char *fq1, const char *fq2, const char *wl_path, const int n_buckets, const char *counts);

/* performs initial barcode count */
void count_barcodes(BarcodeDict *bcdict, FILE *fq, const int interleaved);

/* sorts preprocessed FASTQs by barcode */
void sort_fastq(const char *fq1, const char *fq2);

/* barcode min confidence */
#define BC_CONF_THRESH 0.975

#endif /* PREPROCESS_H */

