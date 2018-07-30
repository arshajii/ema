#ifndef MAIN_H
#define MAIN_H

#include <stdlib.h>
#include <stdint.h>
#include "util.h"
#include "techs.h"

char *chrom_lookup(const chrom_t chrom);
chrom_t chrom_index(const char *chrom);
void read_fai(FILE *fai_file);

extern int num_threads_per_file;
extern int num_threads_for_files;
extern char *rg;
extern char *bx_index;
extern char **pg_argv;
extern int pg_argc;

extern PlatformProfile *tech;

/* debugging */
#define DEBUG 0
#define TRACK_READ ""

#define VERSION "0.6.2"

#define KB 1024
#define MB (1024 * KB)
#define GB (1024 * MB)

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#endif /* MAIN_H */

