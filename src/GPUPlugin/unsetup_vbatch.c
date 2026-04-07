#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "dmrg_vbatch.h"
#include "setup_vbatch.h"

void unsetup_vbatch(IntegerType** p_left_patch_start,
                    IntegerType** p_right_patch_start,
                    IntegerType** p_xy_patch_start,

                    FpType** p_Abatch_,
                    FpType** p_Bbatch_)
{
	const IntegerType idebug = 1;
	if (idebug >= 1) {
		printf("unsetup_vbatch called\n");
	};

	assert(p_Abatch_ != NULL);
	assert(p_Bbatch_ != NULL);

	assert(p_left_patch_start != NULL);
	assert(p_right_patch_start != NULL);
	assert(p_xy_patch_start != NULL);

	if (*p_Abatch_ != NULL) {
		dmrg_free((void*)*p_Abatch_);
		*p_Abatch_ = NULL;
	};

	if (*p_Bbatch_ != NULL) {
		dmrg_free((void*)*p_Bbatch_);
		*p_Bbatch_ = NULL;
	};

	if (*p_left_patch_start != NULL) {
		free(*p_left_patch_start);
		*p_left_patch_start = NULL;
	};

	if (*p_right_patch_start != NULL) {
		free(*p_right_patch_start);
		*p_right_patch_start = NULL;
	};

	if (*p_xy_patch_start != NULL) {
		free(*p_xy_patch_start);
		*p_xy_patch_start = NULL;
	};
}
