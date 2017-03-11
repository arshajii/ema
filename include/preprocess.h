#ifndef PREPROCESS_H
#define PREPROCESS_H

/* corrects barcodes and generate new FASTQs */
void preprocess_fastqs(const char *fq1, const char *fq2, const char *wl_path);

/* barcode min confidence */
#define BC_CONF_THRESH 0.975

#endif /* PREPROCESS_H */

