#include "dmrg_types.h"
#include "dmrg_vbatch.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void apply_Htarget_sparse(IntegerType noperator,
                          IntegerType npatches,
                          IntegerType left_patch_start_[],
                          IntegerType right_patch_start_[],
                          IntegerType xy_patch_start_[],

                          IntegerType nC_[],
                          FpType*     gAbatch_[],
                          IntegerType ld_gAbatch_[],
                          FpType*     gBbatch_[],
                          IntegerType ld_gBbatch_[],

                          FpType Xin_[],
                          FpType Yout_[])

// #define Abatch(i,j) Abatch_[ indx2f(i,j,ld_Abatch) ]
// #define Bbatch(i,j) Bbatch_[ indx2f(i,j,ld_Bbatch) ]
{

	const double      giga   = 1000.0 * 1000.0 * 1000.0;
	const IntegerType idebug = 1;
	const IntegerType ialign = 32;

	double total_time = -dmrg_get_wtime();

	size_t nbytes_X  = 0;
	size_t nbytes_Y  = 0;
	size_t nbytes_BX = 0;

	double gflops1         = 0.0;
	double gflops2         = 0.0;
	double time_1st_vbatch = 0.0;
	double time_2nd_vbatch = 0.0;

	FpType* X_ = &(Xin_[0]);
	FpType* Y_ = &(Yout_[0]);
#ifdef USE_MAGMA
	const IntegerType ltrue           = (1 == 1);
	IntegerType       need_allocate_X = ltrue;
	IntegerType       need_allocate_Y = ltrue;
	/*
	IntegerType need_allocate_X = !dmrg_is_managed( Xin_ );
	IntegerType need_allocate_Y = !dmrg_is_managed( Yout_ );
	*/

#endif

	/*
	 ------------------
	 compute  Y = H * X
	 ------------------
	*/
	IntegerType ipatch = 0;
	IntegerType jpatch = 0;

	IntegerType left_patch_size_[npatches];
	IntegerType right_patch_size_[npatches];

	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		IntegerType L1 = left_patch_start_[ipatch - 1];
		IntegerType L2 = left_patch_start_[ipatch] - 1;

		left_patch_size_[ipatch - 1] = L2 - L1 + 1;
	};
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		IntegerType R1                = right_patch_start_[ipatch - 1];
		IntegerType R2                = right_patch_start_[ipatch] - 1;
		right_patch_size_[ipatch - 1] = R2 - R1 + 1;
	};

#ifdef USE_MAGMA
	long xy_size = 0;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		xy_size += left_patch_size_[ipatch - 1] * right_patch_size_[ipatch - 1];
	};
	long xy_size_dim = ialign * ICEIL(xy_size, ialign);

	/*
	 * -----------------------
	 * allocate unified memory
	 * -----------------------
	 */
	if (need_allocate_X) {
		nbytes_X   = sizeof(FpType) * xy_size_dim;
		FpType* X_ = new FpType[xy_size_dim];

		assert(X_ != NULL);
#ifdef USE_MAGMA
		{
			IntegerType ngpu         = dmrg_get_ngpu();
			IntegerType igpu         = 0;
			size_t      X_inc        = xy_size_dim / ngpu;
			size_t      X_inc_nbytes = sizeof(FpType) * X_inc;
			for (igpu = 0; igpu < ngpu; igpu++) {
				size_t  offset = igpu * X_inc;
				FpType* pX     = &(X_[offset]);
				dmrg_set_readonly((void*)pX, X_inc_nbytes, igpu);
				dmrg_prefetch_to_device((void*)pX, X_inc_nbytes, igpu);
			};
		}
#endif
		void*  dest  = (void*)&(X_[0]);
		void*  src   = (void*)&(Xin_[0]);
		size_t count = sizeof(FpType) * xy_size;
		dmrg_memcpy(dest, src, count);
	};

	if (need_allocate_Y) {
		nbytes_Y   = sizeof(FpType) * xy_size_dim;
		FpType* Y_ = new FpType[xy_size_dim];
		assert(Y_ != NULL);
	};

#endif

	IntegerType nnz_nC = 0;
	IntegerType sum_nC = 0;
	IntegerType max_nC = 0;

	long BX_sizes_[npatches];

	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		BX_sizes_[ipatch - 1] = 0;

		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			IntegerType nrowA = left_patch_size_[ipatch - 1];
			IntegerType ncolA = left_patch_size_[jpatch - 1];

			IntegerType nrowB = right_patch_size_[ipatch - 1];
			IntegerType ncolB = right_patch_size_[jpatch - 1];

			IntegerType nrowX = ncolB;
			IntegerType ncolX = ncolA;

			IntegerType ncolBX = ncolX;
			IntegerType ld_B   = ld_gBbatch_[ipatch - 1];

			assert((1 <= nrowB) && (1 <= ncolB));
			assert((1 <= nrowA) && (1 <= ncolA));
			assert((1 <= nrowB) && (nrowB <= ld_B));
			assert((1 <= nrowX) && (1 <= ncolX));

			IntegerType nconnector = nC_[indx2f(ipatch, jpatch, npatches)];
			if (nconnector >= 1) {
				IntegerType ld_BX = ld_B;

				sum_nC += nconnector;
				max_nC = MAX(max_nC, nconnector);
				nnz_nC++;
				BX_sizes_[ipatch - 1] += ld_BX * (ncolBX * nconnector);
			};
		};
	};
	long sum_BX_sizes = 0;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		sum_BX_sizes += BX_sizes_[ipatch - 1];
	};

	/*
	 ---------------
	 setup gBXbatch
	 ---------------
	 */
	FpType* gBXbatch_[npatches];

	nbytes_BX      = sizeof(FpType) * sum_BX_sizes;
	FpType* pBXmem = new FpType[sum_BX_sizes];

	if (pBXmem == NULL) {
		printf("apply_Htarget_sparse: sum_BX_sizes=%le\n", (double)sum_BX_sizes);
		printf("max_nC=%d, sum_nC=%d, nnz_nC=%d\n", max_nC, sum_nC, nnz_nC);

		IntegerType ipatch = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			printf("ld_gBbatch(%d)=%d, ld_gAbatch(%d)=%d\n",
			       ipatch,
			       ld_gBbatch_[ipatch - 1],
			       ipatch,
			       ld_gAbatch_[ipatch - 1]);
		};
	};
	assert(pBXmem != NULL);

	{
		FpType* BXbatch = pBXmem;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			gBXbatch_[ipatch - 1] = BXbatch;
			BXbatch += BX_sizes_[ipatch - 1];
		};
	};

	/*
	 --------------------------
	 for simplicity and follow magma vbatch
	 assume each group has size only 1
	 --------------------------
	 */
	IntegerType ngroups        = MAX(npatches, sum_nC);
	IntegerType ngroups_dim    = ialign * ICEIL(ngroups, ialign);
	IntegerType batch_size     = ngroups;
	IntegerType batch_size_dim = ialign * ICEIL(batch_size, ialign);

	FpType alpha_array_[ngroups_dim];
	FpType beta_array_[ngroups_dim];

	FpType* a_array_[batch_size_dim];
	FpType* b_array_[batch_size_dim];
	FpType* c_array_[batch_size_dim];

	IntegerType m_array_[ngroups_dim];
	IntegerType n_array_[ngroups_dim];
	IntegerType k_array_[ngroups_dim];
	IntegerType group_size_[ngroups_dim];
	IntegerType lda_array_[batch_size_dim];
	IntegerType ldb_array_[batch_size_dim];
	IntegerType ldc_array_[batch_size_dim];

	char transa_array_[ngroups_dim];
	char transb_array_[ngroups_dim];

#define transa_array(i) transa_array_[(i) - 1]
#define transb_array(i) transb_array_[(i) - 1]
#define m_array(i) m_array_[(i) - 1]
#define n_array(i) n_array_[(i) - 1]
#define k_array(i) k_array_[(i) - 1]
#define alpha_array(i) alpha_array_[(i) - 1]
#define beta_array(i) beta_array_[(i) - 1]
#define lda_array(i) lda_array_[(i) - 1]
#define ldb_array(i) ldb_array_[(i) - 1]
#define ldc_array(i) ldc_array_[(i) - 1]
#define a_array(i) a_array_[(i) - 1]
#define b_array(i) b_array_[(i) - 1]
#define c_array(i) c_array_[(i) - 1]
#define group_size(i) group_size_[(i) - 1]

	IntegerType ibatch = 1;
	gflops1            = 0;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		FpType* BXbatch = gBXbatch_[ipatch - 1];
		FpType* Bbatch  = gBbatch_[ipatch - 1];

		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			IntegerType nconnection = nC_[indx2f(ipatch, jpatch, npatches)];
			if (nconnection >= 1) {

				IntegerType nrowA = left_patch_size_[ipatch - 1];
				IntegerType ncolA = left_patch_size_[jpatch - 1];

				IntegerType nrowB = right_patch_size_[ipatch - 1];
				IntegerType ncolB = right_patch_size_[jpatch - 1];

				IntegerType nrowX = ncolB;
				IntegerType ncolX = ncolA;

				IntegerType nrowY = nrowB;
				IntegerType ncolY = nrowA;

				IntegerType nrowBX = nrowB;
				IntegerType ncolBX = ncolX;

				IntegerType ld_B  = ld_gBbatch_[ipatch - 1];
				IntegerType ld_BX = ld_B;

				IntegerType ix1   = xy_patch_start_[jpatch - 1];
				FpType*     XJ    = &(X_[ix1 - 1]);
				IntegerType ld_XJ = nrowX;

				assert((1 <= nrowY) && (1 <= ncolY));
				assert((1 <= nrowX) && (1 <= ncolX));

				IntegerType iconnection = 0;
				for (iconnection = 1; iconnection <= nconnection; iconnection++) {
					/*
					 ------------------------------------------------------------------
					 BX(1:nrowBX,1:ncolBX) = B(1:nrowB, 1:ncolB) * XJ(1:nrowX,
					 1:ncolX)
					 ------------------------------------------------------------------
					 */
					char transA = 'N';
					char transB = 'N';

					IntegerType mm = nrowBX;
					IntegerType nn = ncolBX;
					IntegerType kk = ncolB;

					FpType alpha = 1;
					FpType beta  = 0;

					FpType* pA = Bbatch;
					FpType* pB = XJ;
					FpType* pC = BXbatch;

					IntegerType ld1 = ld_B;
					IntegerType ld2 = ld_XJ;
					IntegerType ld3 = ld_BX;

					transa_array(ibatch) = transA;
					transb_array(ibatch) = transB;

					m_array(ibatch) = mm;
					n_array(ibatch) = nn;
					k_array(ibatch) = kk;

					alpha_array(ibatch) = alpha;
					beta_array(ibatch)  = beta;

					lda_array(ibatch) = ld1;
					ldb_array(ibatch) = ld2;
					ldc_array(ibatch) = ld3;

					a_array(ibatch) = pA;
					b_array(ibatch) = pB;
					c_array(ibatch) = pC;

					gflops1 += ((2.0 * mm) * nn) * kk;

					Bbatch += ld_B * ncolB;
					BXbatch += ld_BX * ncolBX;
					ibatch++;
				}; /* end for iconnection */
			};

		}; /* end for jpatch */
	}; /* end for ipatch */

	ngroups = (ibatch - 1);
	assert(ngroups == sum_nC);
	{
		IntegerType i = 0;
		for (i = 1; i <= ngroups; i++) {
			group_size(i) = 1;
		};
	};

	/*
	 ------------------
	 first vbatch DGEMM
	 ------------------
	 */

	time_1st_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch(transa_array_,
	                  transb_array_,
	                  m_array_,
	                  n_array_,
	                  k_array_,
	                  alpha_array_,
	                  a_array_,
	                  lda_array_,
	                  b_array_,
	                  ldb_array_,
	                  beta_array_,
	                  c_array_,
	                  ldc_array_,
	                  ngroups,
	                  group_size_);
	time_1st_vbatch += dmrg_get_wtime();
	gflops1 = gflops1 / giga;

#ifdef USE_MAGMA
	{
		FpType*     BX_           = pBXmem;
		IntegerType ngpu          = dmrg_get_ngpu();
		IntegerType igpu          = 0;
		size_t      BX_inc        = sum_BX_sizes / ngpu;
		size_t      BX_inc_nbytes = sizeof(FpType) * BX_inc;
		for (igpu = 0; igpu < ngpu; igpu++) {
			size_t  offset = igpu * BX_inc;
			FpType* pBX    = &(BX_[offset]);
			dmrg_set_readonly((void*)pBX, BX_inc_nbytes, igpu);
		};
	}
#endif

	/*
	   --------------------------------
	   perform  Y = (BX) * transpose(A)
	   --------------------------------
	*/
	gflops2 = 0;
	ibatch  = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		FpType* Abatch  = gAbatch_[ipatch - 1];
		FpType* BXbatch = gBXbatch_[ipatch - 1];

		IntegerType iy = xy_patch_start_[ipatch - 1];
		FpType*     YI = &(Y_[iy - 1]);

		IntegerType nrowA = left_patch_size_[ipatch - 1];
		IntegerType nrowB = right_patch_size_[ipatch - 1];
		IntegerType nrowY = nrowB;
		IntegerType ncolY = nrowA;

		IntegerType ld_A  = ld_gAbatch_[ipatch - 1];
		IntegerType ld_BX = ld_gBbatch_[ipatch - 1];
		IntegerType ld_YI = nrowY;

		/*
		   -------------------------------
		   compute total number of columns
		   -------------------------------
		*/
		IntegerType total_columns = 0;
		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			IntegerType nconnector = nC_[indx2f(ipatch, jpatch, npatches)];
			if (nconnector >= 1) {

				IntegerType ncolA = left_patch_size_[jpatch - 1];

				IntegerType ncolB = right_patch_size_[jpatch - 1];

				IntegerType nrowX = ncolB;
				IntegerType ncolX = ncolA;

				assert((1 <= nrowX) && (1 <= ncolX));

				IntegerType nrowBX = nrowB;
				IntegerType ncolBX = ncolX;

				assert((1 <= nrowBX) && (1 <= ncolBX));

				total_columns += (nconnector * ncolBX);
			};
		}; /* end for jpatch */

		/*
		 --------------------------------------------------------------------------------
		 YI(nrowY,ncolY) = BX(1:nrowBX, 1:ncolumns) * transpose( A(1:nrowA, 1:ncolumns) )
		 --------------------------------------------------------------------------------
		 */
		IntegerType ncolumns = total_columns;
		IntegerType mm       = nrowY;
		IntegerType nn       = ncolY;
		IntegerType kk       = ncolumns;

		char transA = 'N';
		char transB = 'T';

		FpType alpha = 1;
		FpType beta  = 0;

		FpType* pA = BXbatch;
		FpType* pB = Abatch;
		FpType* pC = YI;

		IntegerType ld1 = ld_BX;
		IntegerType ld2 = ld_A;
		IntegerType ld3 = ld_YI;

		transa_array(ibatch) = transA;
		transb_array(ibatch) = transB;
		m_array(ibatch)      = mm;
		n_array(ibatch)      = nn;
		k_array(ibatch)      = kk;

		alpha_array(ibatch) = alpha;
		beta_array(ibatch)  = beta;

		a_array(ibatch) = pA;
		b_array(ibatch) = pB;
		c_array(ibatch) = pC;

		lda_array(ibatch) = ld1;
		ldb_array(ibatch) = ld2;
		ldc_array(ibatch) = ld3;

		gflops2 += ((2.0 * mm) * nn) * kk;
		ibatch++;
	}; /* end for ipatch */

	ngroups = (ibatch - 1);
	assert(ngroups == npatches);
	{
		IntegerType i = 0;
		for (i = 1; i <= ngroups; i++) {
			group_size(i) = 1;
		};
	}

	/*
	 ------------------
	 second vbatch DGEMM
	 ------------------
	 */
	time_2nd_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch(transa_array_,
	                  transb_array_,
	                  m_array_,
	                  n_array_,
	                  k_array_,
	                  alpha_array_,
	                  a_array_,
	                  lda_array_,
	                  b_array_,
	                  ldb_array_,
	                  beta_array_,
	                  c_array_,
	                  ldc_array_,
	                  ngroups,
	                  group_size_);
	time_2nd_vbatch += dmrg_get_wtime();
	gflops2 = gflops2 / giga;

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
	};

	assert(pBXmem != NULL);

	delete[] pBXmem;
	// dmrg_free( (void *) pBXmem );

#ifdef USE_MAGMA
	/*
	 * -------------------
	 * free unified memory
	 * -------------------
	 */

	if (need_allocate_X) {
		//    dmrg_free( X_ );
		delete[] X_;
		X_ = NULL;
	};
	if (need_allocate_Y) {
		void*  dest  = &(Yout_[0]);
		void*  src   = &(Y_[0]);
		size_t count = sizeof(FpType) * xy_size;
		dmrg_memcpy(dest, src, count);
		//   dmrg_free( Y_ );
		delete[] Y_;
		Y_ = NULL;
	};
#endif

	total_time += dmrg_get_wtime();
	if (idebug >= 1) {
		printf("apply_Htarget_sparse: total_time=%lf \n", total_time);
		printf("apply_Htarget_sparse:memory BX (%f GBytes) X (%f GBytes) Y (%f GBytes) \n",
		       (double)nbytes_BX / (giga),
		       (double)nbytes_X / (giga),
		       (double)nbytes_Y / (giga));
	};
}
