#include "analysis.h"
#include "dmrg_vbatch.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const SizeType idebug = 0;

SizeType MOD(const SizeType& x, const SizeType& y) { return ((x) % (y)); }

SizeType indx2f2(const SizeType& i, const SizeType& j, const SizeType& lda)
{
	return (((i)-1) + ((j)-1) * (lda));
}

double nchoosek(SizeType n, SizeType k)
{
	k = MIN(k, n - k);
	if (k == 0) {
		return ((double)1.0);
	};

	double   dval = 1.0;
	SizeType i    = 0;

	for (i = 1; i <= k; i++) {
		double num = (double)(n - i + 1);
		double den = (double)i;
		dval *= (num / den);
	};
	return (dval);
}

double bincoeff(SizeType n, SizeType k) { return (nchoosek(n, k)); }

SizeType gen_patches_comb(SizeType        left_size,
                          SizeType        right_size,
                          SizeType        target_up,
                          SizeType        target_down,
                          SizeType        keep_left_states,
                          SizeType        keep_right_states,
                          VectorSizeType& left_patch_size_,
                          VectorSizeType& right_patch_size_,
                          VectorSizeType& left_patch_up_,
                          VectorSizeType& left_patch_down_,
                          VectorSizeType& right_patch_up_,
                          VectorSizeType& right_patch_down_)
{

	assert((bincoeff(10, 3) == 120) && (bincoeff(10, 7) == bincoeff(10, 3)));

	SizeType nleft_up    = 0;
	SizeType nleft_down  = 0;
	SizeType nright_up   = 0;
	SizeType nright_down = 0;

	SizeType npatch = 0;

	/*
	function [left_patch_size,right_patch_size] =  ...
	          gen_patches_comb( left_size, right_size, ...
	                       target_up, target_down, ...
	                       keep_left_states, keep_right_states)
	% ------------------------------------------------
	%  [left_patch_size,right_patch_size] =  ...
	%          gen_patches_comb( left_size, right_size, ...
	%		       target_up, target_down, ...
	%                      keep_left_states, keep_right_states)
	%
	%  estimate the size of patches based on
	%  combinatorial arguments
	% ------------------------------------------------
	*/

	SizeType max_left_up    = MIN(left_size, target_up);
	SizeType max_left_down  = MIN(left_size, target_down);
	SizeType max_right_up   = MIN(right_size, target_up);
	SizeType max_right_down = MIN(right_size, target_down);

	/*
	% ---------------------------------------------------------
	% setup logical mask to identify valid (nleft_up,nleft_down)
	% and valid matching (nright_up,nright_down)
	% ---------------------------------------------------------
	isvalid_left = zeros(1+max_left_up,1+max_left_down);
	isvalid_right = zeros(1+max_right_up,1+max_right_down);
	*/

	const SizeType ld_isvalid_left  = (1 + max_left_up);
	const SizeType ld_isvalid_right = (1 + max_right_up);

	VectorSizeType isvalid_left_(ld_isvalid_left * (1 + max_left_down), 0);
	VectorSizeType isvalid_right_(ld_isvalid_right * (1 + max_right_down), 0);

	{
		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
				isvalid_left_[indx2f2(
				    1 + nleft_up, 1 + nleft_down, ld_isvalid_left)]
				    = 0;
			};
		};

		for (nright_down = 0; nright_down <= max_right_down; nright_down++) {
			for (nright_up = 0; nright_up <= max_right_up; nright_up++) {
				isvalid_right_[indx2f2(
				    1 + nright_up, 1 + nright_down, ld_isvalid_right)]
				    = 0;
			};
		};
	}

	/*
	for nleft_up=0:max_left_up,
	for nleft_down=0:max_left_down,
	*/
	for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
		for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
			SizeType nright_up   = target_up - nleft_up;
			SizeType nright_down = target_down - nleft_down;

			SizeType isvalid = (nright_up <= max_right_up) && (nright_down <= max_right_down);
			isvalid_left_[indx2f2(1 + nleft_up, 1 + nleft_down, ld_isvalid_left)]
			    = isvalid ? 1 : 0;
			isvalid_right_[indx2f2(1 + nright_up, 1 + nright_down, ld_isvalid_right)]
			    = isvalid ? 1 : 0;
		};
	};

	/*
	--------------------------------
	npatch = sum( isvalid_left(:) );
	--------------------------------
	*/

	{
		SizeType isum = 0;
		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
				SizeType isvalid = isvalid_left_[indx2f2(
				    1 + nleft_up, 1 + nleft_down, ld_isvalid_left)];
				if (isvalid) {
					isum = isum + 1;
				};
			};
		};
		npatch = isum;

		if (idebug >= 1) {
			printf("gen_patches_comb: isum=%lu\n", isum);
		};
	}

	/*
	left_patch_size  = zeros(npatch,1);
	right_patch_size = zeros(npatch,1);
	*/

	/*
	left_ways  = zeros( max_left_up+1, max_left_down+1);
	right_ways = zeros( max_right_up+1, max_right_down+1);
	*/

	const SizeType ld_left_ways  = max_left_up + 1;
	const SizeType ld_right_ways = max_right_up + 1;

	std::vector<double> left_ways_(ld_left_ways * (max_left_down + 1), 0.0);
	std::vector<double> right_ways_(ld_right_ways * (max_right_down + 1), 0.0);

	{
		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
				left_ways_[indx2f2(1 + nleft_up, 1 + nleft_down, ld_left_ways)] = 0;
			};
		};

		for (nright_down = 0; nright_down <= max_right_down; nright_down++) {
			for (nright_up = 0; nright_up <= max_right_up; nright_up++) {
				right_ways_[indx2f2(1 + nright_up, 1 + nright_down, ld_right_ways)]
				    = 0;
			};
		};
	}

	/*
	% ---------------------------------------------
	% vectorized way to compute
	% Cmat_left(1+nleft_up,1+nleft_down) = ...
	%               nchoosek(left_size,nleft_up) *  ...
	%               nchoosek(left_size,nleft_down)
	% ---------------------------------------------
	*/

	/*
	nleft_up     = 0:max_left_up;
	nleft_down   = 0:max_left_down;
	C_nleft_up   = bincoeff( left_size, nleft_up);
	C_nleft_down = bincoeff( left_size, nleft_down);
	*/

	std::vector<double> C_nleft_up_(max_left_up + 1, 0.0);
	std::vector<double> C_nleft_down_(max_left_down + 1, 0.0);

	{
		SizeType nleft_up   = 0;
		SizeType nleft_down = 0;

		for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
			C_nleft_up_[nleft_up] = bincoeff(left_size, nleft_up);
		};

		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			C_nleft_down_[nleft_down] = bincoeff(left_size, nleft_down);
		};
	}

	/*
	%   ------------------------------------------------------
	%   Cmat_left = C_nleft_up(:) * transpose(C_nleft_down(:));
	%   ------------------------------------------------------
	*/

	const SizeType ld_Cmat_left = (1 + max_left_up);

	std::vector<double> Cmat_left_(ld_Cmat_left * (1 + max_left_down), 0.0);

	{
		SizeType nleft_up   = 0;
		SizeType nleft_down = 0;

		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
				Cmat_left_[indx2f2(1 + nleft_up, 1 + nleft_down, ld_Cmat_left)]
				    = C_nleft_up_[nleft_up] * C_nleft_down_[nleft_down];
			};
		};
	}

	/*
	for nleft_down=0:max_left_down,
	for nleft_up=0:max_left_up,
	*/

	{
		SizeType nleft_up   = 0;
		SizeType nleft_down = 0;

		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {

				if (isvalid_left_[indx2f2(
				        1 + nleft_up, 1 + nleft_down, ld_isvalid_left)]) {
					/*
					%    ---------------------------------
					%    note, beware of possible overflow
					%    or loss of accuracy
					%    since nchoosek(n,k) can quickly become
					%    very large numbers
					%
					%    for example, nchoosek(144,72) = 1.4802*10^42
					%    ---------------------------------
					*/
					left_ways_[indx2f2(
					    1 + nleft_up, 1 + nleft_down, ld_left_ways)]
					    = C_nleft_up_[nleft_up] * C_nleft_down_[nleft_down];
				}
			};
		};
	}

	/*
	% ---------------------------------------------
	% vectorized way to compute
	% Cmat_right(1+nright_up,1+nright_down) = ...
	%               nchoosek(right_size,nright_up) *  ...
	%               nchoosek(right_size,nright_down)
	% ---------------------------------------------
	*/

	/*
	nright_up     = 0:max_right_up;
	nright_down   = 0:max_right_down;
	C_nright_up   = bincoeff( right_size, nright_up);
	C_nright_down = bincoeff( right_size, nright_down);
	*/

	std::vector<double> C_nright_up_(max_right_up + 1, 0.0);
	std::vector<double> C_nright_down_(max_right_down + 1, 0.0);

	{
		SizeType nright_up   = 0;
		SizeType nright_down = 0;

		for (nright_up = 0; nright_up <= max_right_up; nright_up++) {
			C_nright_up_[nright_up] = bincoeff(right_size, nright_up);
		};

		for (nright_down = 0; nright_down <= max_right_down; nright_down++) {
			C_nright_down_[nright_down] = bincoeff(right_size, nright_down);
		};
	};

	/*
	%   ----------------------------------------------------------
	%   Cmat_right = C_nright_up(:) * transpose(C_nright_down(:));
	%   ----------------------------------------------------------
	*/

	const SizeType ld_Cmat_right = (1 + max_right_up);

	std::vector<double> Cmat_right_(ld_Cmat_right * (1 + max_right_down), 0.0);

	/*
	for nright_down=0:max_right_down,
	for nright_up = 0:max_right_up;
	*/

	for (nright_down = 0; nright_down <= max_right_down; nright_down++) {
		for (nright_up = 0; nright_up <= max_right_up; nright_up++) {
			Cmat_right_[indx2f2(1 + nright_up, 1 + nright_down, ld_Cmat_right)]
			    = C_nright_up_[nright_up] * C_nright_down_[nright_down];
		};
	};

	{
		SizeType nright_up   = 0;
		SizeType nright_down = 0;

		for (nright_down = 0; nright_down <= max_right_down; nright_down++) {
			for (nright_up = 0; nright_up <= max_right_up; nright_up++) {

				if (isvalid_right_[indx2f2(
				        1 + nright_up, 1 + nright_down, ld_isvalid_right)]) {
					/*
					%    ---------------------------------
					%    note, beware of possible overflow
					%    or loss of accuracy
					%    since nchoosek(n,k) can quickly become
					%    very large numbers
					%
					%    for example, nchoosek(144,72) = 1.4802*10^42
					%    ---------------------------------
					*/
					right_ways_[indx2f2(
					    1 + nright_up, 1 + nright_down, ld_right_ways)]
					    = C_nright_up_[nright_up] * C_nright_down_[nright_down];
				}
			};
		};
	}

	/*
	total_left_ways = sum( left_ways(:) );
	total_right_ways = sum( right_ways(:) );
	*/

	double total_left_ways  = 0;
	double total_right_ways = 0;
	{

		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
				total_left_ways += left_ways_[indx2f2(
				    1 + nleft_up, 1 + nleft_down, ld_left_ways)];
			};
		};

		for (nright_down = 0; nright_down <= max_right_down; nright_down++) {
			for (nright_up = 0; nright_up <= max_right_up; nright_up++) {
				total_right_ways += right_ways_[indx2f2(
				    1 + nright_up, 1 + nright_down, ld_right_ways)];
			};
		};

		if (idebug >= 1) {
			printf("gen_patches_comb: total_left_ways=%lf, total_right_ways=%lf\n",
			       total_left_ways,
			       total_right_ways);
		};
	}

	/*
	% --------------------------------------------------------------------------------------
	%
	% set use_ceil = 1  to generate more patches
	% set use_ceil = 0  to generate fewer patches
	%
	% --------------------------------------------------------------------------------------
	% use ceiling() function  would generate at least 1 row per valid (nleft_up,nleft_down)
	% this may generate more patches
	%
	% use round() function  would remove many unlikely or small states
	% this may generate fewer patches
	%
	% using more patches would spread the rows across more patches and
	% may perform less work overall
	% --------------------------------------------------------------------------------------
	*/

	SizeType use_ceil = 0;
	/*
	-------------------------------
	if (use_ceil){
	  left_ways = ceil((left_ways/total_left_ways) * keep_left_states );
	  right_ways = ceil((right_ways/total_right_ways) * keep_right_states);

	  }
	else {
	  left_ways = round((left_ways/total_left_ways) * keep_left_states );
	  right_ways = round((right_ways/total_right_ways) * keep_right_states);
	};
	-------------------------------
	*/

	for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
		for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
			double frac
			    = (left_ways_[indx2f2(1 + nleft_up, 1 + nleft_down, ld_left_ways)]
			       / total_left_ways * keep_left_states);

			left_ways_[indx2f2(1 + nleft_up, 1 + nleft_down, ld_left_ways)]
			    = (use_ceil) ? ceil(frac) : round(frac);
		};
	};

	for (nright_up = 0; nright_up <= max_right_up; nright_up++) {
		for (nright_down = 0; nright_down <= max_right_down; nright_down++) {
			double frac
			    = (right_ways_[indx2f2(1 + nright_up, 1 + nright_down, ld_right_ways)]
			       / total_right_ways * keep_right_states);

			right_ways_[indx2f2(1 + nright_up, 1 + nright_down, ld_right_ways)]
			    = (use_ceil) ? ceil(frac) : round(frac);
		};
	};

	npatch = 0;

	/*
	for nleft_down=0:max_left_down,
	for nleft_up=0:max_left_up,
	*/

	for (nleft_down = 0; nleft_down <= max_left_down; nleft_down++) {
		for (nleft_up = 0; nleft_up <= max_left_up; nleft_up++) {
			nright_up   = target_up - nleft_up;
			nright_down = target_down - nleft_down;

			if (isvalid_left_[indx2f2(1 + nleft_up, 1 + nleft_down, ld_isvalid_left)]) {
				SizeType lsize = left_ways_[indx2f2(
				    1 + nleft_up, 1 + nleft_down, ld_left_ways)];
				SizeType rsize = right_ways_[indx2f2(
				    1 + nright_up, 1 + nright_down, ld_right_ways)];

				SizeType has_work = (lsize >= 1) && (rsize >= 1);
				if (has_work) {
					npatch                       = npatch + 1;
					left_patch_size_[npatch - 1] = lsize;
					left_patch_up_[npatch - 1]   = nleft_up;
					left_patch_down_[npatch - 1] = nleft_down;

					right_patch_size_[npatch - 1] = rsize;
					right_patch_up_[npatch - 1]   = nright_up;
					right_patch_down_[npatch - 1] = nright_down;
				};
			};
		};
	};

	/*
	% -------------------------
	% keep only non-zero patches
	% -------------------------
	*/

	/*
	idx_nonzero = find( (left_patch_size > 0) & (right_patch_size > 0) );
	left_patch_size = left_patch_size_[idx_nonzero-1];
	right_patch_size = right_patch_size_[idx_nonzero-1];
	*/
	{
		SizeType ipatch = 0;
		for (ipatch = 1; ipatch <= npatch; ipatch++) {
			SizeType isvalid = (left_patch_size_[ipatch - 1] > 0)
			    && (right_patch_size_[ipatch - 1] > 0);
			assert(isvalid);
		};
	}

	return (npatch);
}
