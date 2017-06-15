#ifndef MAIN_H
#define MAIN_H

#include <stdlib.h>
#include <stdint.h>
#include "samrecord.h"
#include "util.h"

char *chrom_lookup(const chrom_t chrom);
chrom_t chrom_index(const char *chrom);
void read_fai(FILE *fai_file);

/* debugging */
#define DEBUG 0
#define TRACK_READ ""

#define VERSION "0.0.2"

#endif /* MAIN_H */

