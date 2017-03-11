#ifndef MAIN_H
#define MAIN_H

#include <stdlib.h>
#include <stdint.h>
#include "samrecord.h"
#include "util.h"

const char *chrom_lookup(const chrom_t chrom);
chrom_t chrom_index(const char *chrom);
SAMRecord *read_sam(FILE *sam_file, size_t *num_records);
void read_fai(FILE *fai_file);

/* debugging */
#define DEBUG 0
#define TRACK_READ ""

#define VERSION "0.0.1"

#endif /* MAIN_H */

