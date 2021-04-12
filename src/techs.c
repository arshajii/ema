#include <stdlib.h>
#include <string.h>
#include "techs.h"

static bc_t extract_bc_haptag(FASTQRecord *rec)
{
	char *bc_str = strrchr(rec->id, ':');
	assert(bc_str != NULL);
	*bc_str = '\0';

	char *space = strchr(rec->id, ' ');  // for Long Ranger formatted FASTQs
	if (space != NULL)
		*space = '\0';

	return encode_bc(&bc_str[1], 1);
}

static bc_t extract_bc_10x(FASTQRecord *rec)
{
	char *bc_str = strrchr(rec->id, ':');
	assert(bc_str != NULL);
	*bc_str = '\0';

	char *space = strchr(rec->id, ' ');  // for Long Ranger formatted FASTQs
	if (space != NULL)
		*space = '\0';

	return encode_bc(&bc_str[1], 0);
}

static bc_t extract_bc_truseq(FASTQRecord *rec)
{
	char *bc_str = rec->id[0] == '@' ? &(rec->id)[1] : rec->id;
	return atoi(bc_str);
}

static bc_t extract_bc_cptseq(FASTQRecord *rec)
{
	char *bc_str = strrchr(rec->id, ':');
	assert(bc_str != NULL);
	*bc_str = '\0';
	return atoi(&bc_str[3]);
}

static PlatformProfile profiles[] = {
	{.name = "haptag",
	 .extract_bc = extract_bc_haptag,
	 .many_clouds = 0,
	 .dist_thresh = 50000,
	 .error_rate = 0.001,
	 .n_density_probs = 4,
	 .density_probs = {0.6, 0.05, 0.2, 0.01}},
	
	{.name = "10x",
	 .extract_bc = extract_bc_10x,
	 .many_clouds = 0,
	 .dist_thresh = 50000,
	 .error_rate = 0.001,
	 .n_density_probs = 4,
	 .density_probs = {0.6, 0.05, 0.2, 0.01}},

	{.name = "tru",
	 .extract_bc = extract_bc_truseq,
	 .many_clouds = 1,
	 .dist_thresh = 15000,
	 .error_rate = 0.001,
	 .n_density_probs = 4,
	 .density_probs = {0.6, 0.05, 0.2, 0.01}},

	{.name = "cpt",
	 .extract_bc = extract_bc_cptseq,
	 .many_clouds = 1,
	 .dist_thresh = 3500,
	 .error_rate = 0.01,
	 .n_density_probs = 9,
	 .density_probs = {0.6, 0.01, 0.15, 0.001, 0.05, 0.001, 0.02, 0.001, 0.01}},

	{.name = NULL},
};

PlatformProfile *get_platform_profile_by_name(const char *name)
{
	for (size_t i = 0; profiles[i].name != NULL; i++) {
		if (strcmp(name, profiles[i].name) == 0)
			return &profiles[i];
	}

	return NULL;
}

