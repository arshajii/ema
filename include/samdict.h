#ifndef SAMDICT_H
#define SAMDICT_H

#include <stdlib.h>
#include <stdint.h>
#include "samrecord.h"
#include "align.h"

#define MAX_CANDIDATES 5000
#define DEFAULT_CANDS  5
#define SAM_DICT_CAP_SMALL (1 << 16)
#define SAM_DICT_CAP_LARGE (1 << 25)

/* internal dictionary entries */
typedef struct sam_dict_ent {
	SAMRecord *key;
	struct sam_dict_ent *clash_next;
	struct sam_dict_ent *link_next;  // for easy clearing
	struct sam_dict_ent *mate;

	size_t num_cands;
	size_t cap_cands;
	SAMRecord **cand_records;
	Cloud **cand_clouds;
	double *gammas;

	unsigned visited : 1;
} SAMDictEnt;

/* SAM record -to- mapping information dictionary */
typedef struct {
	SAMDictEnt *head;
	SAMDictEnt **entries;
	uint32_t count;
	uint32_t cap;
} SAMDict;

SAMDictEnt *sde_new(SAMRecord *key, Cloud *v);
void sde_free(SAMDictEnt *sde);

SAMDict *sam_dict_new(void);
void sam_dict_free(SAMDict *sd);

int sam_dict_add(SAMDict *sd, SAMRecord *k, Cloud *v, const int force);
void sam_dict_del(SAMDict *sd, SAMRecord *k);
void sam_dict_clear(SAMDict *sd);
SAMRecord *find_best_record(SAMDictEnt *e);

#endif /* SAMDICT_H */

