
#include "dmrg_types.h"

#include "dmrg_vbatch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include <cassert>
#include <vector>

void apply_Htarget_vbatch(SizeType            noperator,
                          SizeType            npatches,
                          VectorSizeType      left_patch_start_,
                          VectorSizeType      right_patch_start_,
                          VectorSizeType      xy_patch_start_,
                          std::vector<FpType> Abatch_,
                          SizeType            ld_Abatch,
                          std::vector<FpType> Bbatch_,
                          SizeType            ld_Bbatch,
                          FpType*             X_,
                          FpType*             Y_)
{
	const IntegerType idebug = 1;
	const IntegerType ialign = 32;
	const double      giga   = 1000.0 * 1000.0 * 1000.0;

	double gflops1         = 0.0;
	double gflops2         = 0.0;
	double time_1st_vbatch = 0.0;
	double time_2nd_vbatch = 0.0;

	size_t nbytes_BX = 0;

	/*
	 ------------------
	 compute  Y = H * X
	 ------------------
	*/
	IntegerType left_max_states  = left_patch_start_[npatches] - 1;
	IntegerType right_max_states = right_patch_start_[npatches] - 1;

	std::vector<IntegerType> left_patch_size_(npatches);
	std::vector<IntegerType> right_patch_size_(npatches);

	for (SizeType ipatch = 1; ipatch <= npatches; ipatch++) {
		IntegerType L1 = left_patch_start_[ipatch - 1];
		IntegerType L2 = left_patch_start_[ipatch] - 1;

		left_patch_size_[ipatch - 1] = L2 - L1 + 1;
	};
	for (SizeType ipatch = 1; ipatch <= npatches; ipatch++) {
		IntegerType R1                = right_patch_start_[ipatch - 1];
		IntegerType R2                = right_patch_start_[ipatch] - 1;
		right_patch_size_[ipatch - 1] = R2 - R1 + 1;
	};

	IntegerType ngroups        = npatches;
	IntegerType ngroups_dim    = ialign * ICEIL(ngroups, ialign);
	IntegerType batch_size     = ngroups * noperator;
	IntegerType batch_size_dim = ialign * ICEIL(batch_size, ialign);

	std::vector<FpType>  alpha_array_(ngroups_dim);
	std::vector<FpType>  beta_array_(ngroups_dim);
	std::vector<FpType*> a_array_(batch_size_dim);
	std::vector<FpType*> b_array_(batch_size_dim);
	std::vector<FpType*> c_array_(batch_size_dim);

	std::vector<IntegerType> m_array_(ngroups_dim);
	std::vector<IntegerType> n_array_(ngroups_dim);
	std::vector<IntegerType> k_array_(ngroups_dim);
	std::vector<IntegerType> group_size_(ngroups_dim);
	std::vector<IntegerType> lda_array_(batch_size_dim);
	std::vector<IntegerType> ldb_array_(batch_size_dim);
	std::vector<IntegerType> ldc_array_(batch_size_dim);

	std::vector<char> transa_array_(ngroups_dim);
	std::vector<char> transb_array_(ngroups_dim);

	IntegerType nrowA = left_max_states;
	IntegerType ncolA = nrowA;
	IntegerType nrowB = right_max_states;
	IntegerType ncolB = nrowB;

	IntegerType nrowBX = nrowB;
	IntegerType ncolBX = (ncolA * noperator);
	IntegerType ld_BX  = ialign * ICEIL(nrowBX, ialign);

	nbytes_BX = ((sizeof(FpType) * ld_BX) * (ncolA * noperator));
	// FpType *BX_ = (FpType *) dmrg_malloc( nbytes_BX );
	FpType* BX_ = new FpType[ld_BX * ncolA * noperator];
	assert(BX_ != NULL);

	IntegerType idx = 1;
	for (SizeType jpatch = 1; jpatch <= npatches; jpatch++) {
		IntegerType igroup = jpatch;
		IntegerType j1     = xy_patch_start_[jpatch - 1];
		IntegerType j2     = xy_patch_start_[jpatch] - 1;
		IntegerType nrowX  = right_patch_size_[jpatch - 1];
		IntegerType ncolX  = left_patch_size_[jpatch - 1];
		assert((j2 - j1 + 1) == (nrowX * ncolX));
		(void)j2;
		(void)nrowX;
		(void)ncolX;

		/*
		 --------------------------------------
		 XJ = reshape( X(j1:j2), nrowX, ncolX )
		 --------------------------------------
		 */
		FpType*     XJ    = &(X_[j1 - 1]);
		IntegerType ld_XJ = nrowX;

		IntegerType R1   = right_patch_start_[jpatch - 1];
		IntegerType R2   = right_patch_start_[jpatch] - 1;
		IntegerType L1   = left_patch_start_[jpatch - 1];
		IntegerType L2   = left_patch_start_[jpatch] - 1;
		IntegerType kmax = noperator;
		IntegerType k    = 0;

		/*
		 -------------------------------
		 independent DGEMM in same group
		 -------------------------------
		 */
		group_size_[igroup - 1] = kmax;
		for (k = 1; k <= kmax; k++) {
			IntegerType offsetB  = (k - 1) * ncolB;
			IntegerType offsetBX = (k - 1) * ncolA;

			/*
			------------------------------------------------------------------------
			BX(1:nrowBX, offsetBX + (L1:L2)) = Bbatch(1:nrowBX, offsetB + (R1:R2) ) *
			                                     XJ( 1:(R2-R1+1), 1:(L2-L1+1));
			------------------------------------------------------------------------
			*/
			transa_array_[igroup - 1] = 'N';
			transb_array_[igroup - 1] = 'N';
			IntegerType mm            = nrowBX;
			IntegerType nn            = L2 - L1 + 1;
			IntegerType kk            = R2 - R1 + 1;
			m_array_[igroup - 1]      = mm;
			n_array_[igroup - 1]      = nn;
			k_array_[igroup - 1]      = kk;

			gflops1 += ((2.0 * mm) * nn) * kk;

			alpha_array_[igroup - 1] = (FpType)1;
			beta_array_[igroup - 1]  = (FpType)0;

			c_array_[idx - 1]      = &(BX_[indx2f(1, offsetBX + L1, ld_BX)]);
			ldc_array_[igroup - 1] = ld_BX;

			a_array_[idx - 1]      = &(Bbatch_[indx2f(1, offsetB + R1, ld_Bbatch)]);
			lda_array_[igroup - 1] = ld_Bbatch;

			b_array_[idx - 1]      = XJ;
			ldb_array_[igroup - 1] = ld_XJ;
			idx                    = idx + 1;
		};
	};
	/*
	 ------------------
	 first vbatch DGEMM
	 ------------------
	 */

	time_1st_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch(transa_array_.data(),
	                  transb_array_.data(),
	                  m_array_.data(),
	                  n_array_.data(),
	                  k_array_.data(),
	                  alpha_array_.data(),
	                  a_array_.data(),
	                  lda_array_.data(),
	                  b_array_.data(),
	                  ldb_array_.data(),
	                  beta_array_.data(),
	                  c_array_.data(),
	                  ldc_array_.data(),
	                  ngroups,
	                  group_size_.data());
	time_1st_vbatch += dmrg_get_wtime();
	gflops1 = gflops1 / (giga);

	/*
	 -------------------------------------------------
	 perform computations with  Y += (BX)*transpose(A)
	 -------------------------------------------------
	*/
	for (SizeType ipatch = 1; ipatch <= npatches; ipatch++) {
		IntegerType igroup = ipatch;

		IntegerType i1 = xy_patch_start_[ipatch - 1];
		IntegerType i2 = xy_patch_start_[ipatch] - 1;

		IntegerType R1 = right_patch_start_[ipatch - 1];
		IntegerType R2 = right_patch_start_[ipatch] - 1;

		IntegerType L1 = left_patch_start_[ipatch - 1];
		IntegerType L2 = left_patch_start_[ipatch] - 1;

		IntegerType isok = ((R2 - R1 + 1) == right_patch_size_[ipatch - 1])
		    && ((L2 - L1 + 1) == left_patch_size_[ipatch - 1]);
		assert(isok);
		(void)i2;
		(void)isok;

		FpType*     YI     = &(Y_[i1 - 1]);
		IntegerType nrowYI = R2 - R1 + 1;
		IntegerType ld_YI  = nrowYI;
		IntegerType ncolYI = L2 - L1 + 1;
		assert((i2 - i1 + 1) == (nrowYI * ncolYI));

		/*
		   --------------------------------------------------------------------
		   YI(1:(R2-R1+1),1:(L2-L1+1)) = BX( R1:R2,1:ncolBX) *
		                                    transpose( Abatch( L1:L2,1:ncolBX) );
		   --------------------------------------------------------------------
		 */
		group_size_[igroup - 1]   = 1;
		transa_array_[igroup - 1] = 'N';
		transb_array_[igroup - 1] = 'T';
		IntegerType mm            = nrowYI;
		IntegerType nn            = ncolYI;
		IntegerType kk            = ncolBX;
		m_array_[igroup - 1]      = mm;
		n_array_[igroup - 1]      = nn;
		k_array_[igroup - 1]      = kk;
		gflops2 += ((2.0 * mm) * nn) * kk;
		alpha_array_[igroup - 1] = (FpType)1;
		beta_array_[igroup - 1]  = (FpType)0;
		a_array_[igroup - 1]     = &(BX_[indx2f(R1, 1, ld_BX)]);
		lda_array_[igroup - 1]   = ld_BX;
		b_array_[igroup - 1]     = &(Abatch_[indx2f(L1, 1, ld_Abatch)]);
		ldb_array_[igroup - 1]   = ld_Abatch;
		c_array_[igroup - 1]     = YI;
		ldc_array_[igroup - 1]   = ld_YI;
	};
	ngroups = npatches;

	/*
	 ------------------
	 second vbatch DGEMM
	 ------------------
	 */
	time_2nd_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch(transa_array_.data(),
	                  transb_array_.data(),
	                  m_array_.data(),
	                  n_array_.data(),
	                  k_array_.data(),
	                  alpha_array_.data(),
	                  a_array_.data(),
	                  lda_array_.data(),
	                  b_array_.data(),
	                  ldb_array_.data(),
	                  beta_array_.data(),
	                  c_array_.data(),
	                  ldc_array_.data(),
	                  ngroups,
	                  group_size_.data());
	time_2nd_vbatch += dmrg_get_wtime();
	gflops2 = gflops2 / (giga);

	if (idebug >= 1) {
		printf("1st vbatch %lf gflops/sec (gflops1=%lf,time=%lf)\n",
		       gflops1 / time_1st_vbatch,
		       gflops1,
		       time_1st_vbatch);
		printf("2nd vbatch %lf gflops/sec (gflops2=%lf,time=%lf)\n",
		       gflops2 / time_2nd_vbatch,
		       gflops2,
		       time_2nd_vbatch);

		printf("overall %lf gflops/sec\n",
		       (gflops1 + gflops2) / (time_1st_vbatch + time_2nd_vbatch));
		printf("memory BX(%lf GBytes)\n", (double)nbytes_BX / (giga));
	};

	// dmrg_free( (void *) BX_ );
	delete[] BX_;
}
