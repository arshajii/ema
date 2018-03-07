#ifndef TECHS_H
#define TECHS_H

#include <stdint.h>
#include "samrecord.h"
#include "util.h"

typedef bc_t (*bc_extract_routine_t)(FASTQRecord *);

typedef struct {
	const char *name;

	bc_extract_routine_t extract_bc;
	int many_clouds;
	uint32_t dist_thresh;
	double error_rate;

	size_t n_density_probs;
	double density_probs[100];
} PlatformProfile;

PlatformProfile *get_platform_profile_by_name(const char *name);

#endif /* TECHS_H */

