#ifndef MAIN_H
#define MAIN_H

#include <stdlib.h>
#include <stdint.h>
#include "util.h"
#include "techs.h"

char *chrom_lookup(const chrom_t chrom);
chrom_t chrom_index(const char *chrom);
void read_fai(FILE *fai_file);

extern int NUM_THREADS;
extern char *rg;
extern char **pg_argv;
extern int pg_argc;

extern PlatformProfile *tech;

/* debugging */
#define DEBUG 0
#define TRACK_READ ""

#define VERSION "0.1.0"

#endif /* MAIN_H */

