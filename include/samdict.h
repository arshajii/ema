#ifndef SAMDICT_H
#define SAMDICT_H

#include <stdlib.h>
#include <stdint.h>
#include "samrecord.h"
#include "align.h"

#define MAX_CANDIDATES 5000
#define DICT_CAP (1 << 16)

/* internal dictionary entries */
typedef struct sam_dict_ent {
	SAMRecord *key;
	struct sam_dict_ent *clash_next;
	struct sam_dict_ent *link_next;  // for easy clearing
	struct sam_dict_ent *mate;

	uint16_t num_cands;
	SAMRecord *cand_records[MAX_CANDIDATES];
	Cloud *cand_clouds[MAX_CANDIDATES];
	double gammas[MAX_CANDIDATES];
} SAMDictEnt;

/* SAM record -to- mapping information dictionary */
typedef struct {
	SAMDictEnt *head;
	SAMDictEnt *entries[DICT_CAP];
	uint32_t count;
} SAMDict;

SAMDict *sam_dict_new(void);
int sam_dict_add(SAMDict *sd, SAMRecord *k, Cloud *v, const int force);
void sam_dict_del(SAMDict *sd, SAMRecord *k);
void sam_dict_clear(SAMDict *sd);
SAMRecord *find_best_record(SAMDictEnt *e);

#endif /* SAMDICT_H */

