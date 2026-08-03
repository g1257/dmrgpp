#include "DMRGConfig.h"
#include "dmrg_types.h"
#include "dmrg_vbatch.h"
#include <algorithm>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T>
void apply_Htarget_sparse(SizeType                 noperator,
                          SizeType                 npatches,
                          VectorSizeType           left_patch_start_,
                          VectorSizeType           right_patch_start_,
                          VectorSizeType           xy_patch_start_,
                          VectorSizeType           nC_,
                          T**                      gAbatch_,
                          const VectorIntegerType& ld_gAbatch_,
                          T**                      gBbatch_,
                          const VectorIntegerType& ld_gBbatch_,
                          T*                       Xin_,
                          T*                       Yout_)
{

	typedef std::vector<T*> VectorPointerType;

	const double      giga   = 1000.0 * 1000.0 * 1000.0;
	const IntegerType idebug = 1;
	const SizeType    ialign = 32;

	double total_time = -dmrg_get_wtime();

	size_t nbytes_X  = 0;
	size_t nbytes_Y  = 0;
	size_t nbytes_BX = 0;

	double gflops1         = 0.0;
	double gflops2         = 0.0;
	double time_1st_vbatch = 0.0;
	double time_2nd_vbatch = 0.0;

	T* X_ = &(Xin_[0]);
	T* Y_ = &(Yout_[0]);
#ifdef USE_MAGMA
	IntegerType need_allocate_X = !dmrg_is_managed(Xin_);
	IntegerType need_allocate_Y = !dmrg_is_managed(Yout_);
#endif
	assert(X_ != nullptr);
	assert(Y_ != nullptr);

	/*
	 ------------------
	 compute  Y = H * X
	 ------------------
	*/
	SizeType ipatch = 0;
	SizeType jpatch = 0;

	VectorSizeType left_patch_size_(npatches, 0);
	VectorSizeType right_patch_size_(npatches, 0);
	// IntegerType left_patch_size_[npatches];
	// IntegerType right_patch_size_[npatches];

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
	SizeType xy_size = 0;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		xy_size += left_patch_size_[ipatch - 1] * right_patch_size_[ipatch - 1];
	};
	SizeType xy_size_dim = ialign * ICEIL(xy_size, ialign);

	/*
	 * -----------------------
	 * allocate unified memory
	 * -----------------------
	 */
	if (need_allocate_X) {
		nbytes_X = sizeof(T) * xy_size_dim;
		X_       = dmrg_malloc<T>(nbytes_X, nbytes_X);
		assert(X_ != nullptr);

#ifdef USE_MAGMA
		{
			IntegerType ngpu         = dmrg_get_ngpu();
			IntegerType igpu         = 0;
			size_t      X_inc        = xy_size_dim / ngpu;
			size_t      X_inc_nbytes = sizeof(T) * X_inc;
			for (igpu = 0; igpu < ngpu; igpu++) {
				size_t offset = igpu * X_inc;
				T*     pX     = &(X_[offset]);
				dmrg_set_readonly((void*)pX, X_inc_nbytes, igpu);
				dmrg_prefetch_to_device((void*)pX, X_inc_nbytes, igpu);
			};
		}
#endif
		void*  dest  = (void*)&(X_[0]);
		void*  src   = (void*)&(Xin_[0]);
		size_t count = sizeof(T) * xy_size;
		dmrg_memcpy(dest, src, count);
	}

	if (need_allocate_Y) {
		nbytes_Y = sizeof(T) * xy_size_dim;
		Y_       = dmrg_malloc<T>(nbytes_Y, nbytes_Y);
		assert(Y_ != nullptr);
	}

#endif

	IntegerType nnz_nC = 0;
	SizeType    sum_nC = 0;
	IntegerType max_nC = 0;

	VectorSizeType BX_sizes_(npatches, 0);
	// long BX_sizes_[npatches];

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
			(void)nrowA;
			(void)nrowB;
			(void)nrowX;

			IntegerType nconnector = nC_[indx2f(ipatch, jpatch, npatches)];
			if (nconnector >= 1) {
				IntegerType ld_BX = ld_B;

				sum_nC += nconnector;
				max_nC = std::max(max_nC, nconnector);
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
	// T* gBXbatch_[npatches];
	VectorPointerType gBXbatch_(npatches, nullptr);

	nbytes_BX = sizeof(T) * sum_BX_sizes;
	T* pBXmem = dmrg_malloc<T>(nbytes_BX, nbytes_BX);

	if (pBXmem == NULL) {
		printf("apply_Htarget_sparse: sum_BX_sizes=%le\n", (double)sum_BX_sizes);
		printf("max_nC=%d, sum_nC=%lu, nnz_nC=%d\n", max_nC, sum_nC, nnz_nC);

		SizeType ipatch = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			printf("ld_gBbatch(%lu)=%d, ld_gAbatch(%lu)=%d\n",
			       ipatch,
			       ld_gBbatch_[ipatch - 1],
			       ipatch,
			       ld_gAbatch_[ipatch - 1]);
		};
	};
	assert(pBXmem != NULL);

	{
		T* BXbatch = pBXmem;
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
	SizeType ngroups        = std::max(npatches, sum_nC);
	SizeType ngroups_dim    = ialign * ICEIL(ngroups, ialign);
	SizeType batch_size     = ngroups;
	SizeType batch_size_dim = ialign * ICEIL(batch_size, ialign);

	T* alpha_array_ = dmrg_malloc<T>(sizeof(T) * ngroups_dim, sizeof(T) * ngroups_dim);
	T* beta_array_  = dmrg_malloc<T>(sizeof(T) * ngroups_dim, sizeof(T) * ngroups_dim);
	assert(alpha_array_ != nullptr);
	assert(beta_array_ != nullptr);

	T** a_array_ = dmrg_malloc<T*>(sizeof(T*) * batch_size_dim, sizeof(T*) * batch_size_dim);
	T** b_array_ = dmrg_malloc<T*>(sizeof(T*) * batch_size_dim, sizeof(T*) * batch_size_dim);
	T** c_array_ = dmrg_malloc<T*>(sizeof(T*) * batch_size_dim, sizeof(T*) * batch_size_dim);
	assert(a_array_ != nullptr);
	assert(b_array_ != nullptr);
	assert(c_array_ != nullptr);

	IntegerType* m_array_    = dmrg_malloc<IntegerType>(sizeof(IntegerType) * ngroups_dim,
                                                         sizeof(IntegerType) * ngroups_dim);
	IntegerType* n_array_    = dmrg_malloc<IntegerType>(sizeof(IntegerType) * ngroups_dim,
                                                         sizeof(IntegerType) * ngroups_dim);
	IntegerType* k_array_    = dmrg_malloc<IntegerType>(sizeof(IntegerType) * ngroups_dim,
                                                         sizeof(IntegerType) * ngroups_dim);
	IntegerType* group_size_ = dmrg_malloc<IntegerType>(sizeof(IntegerType) * ngroups_dim,
	                                                    sizeof(IntegerType) * ngroups_dim);
	IntegerType* lda_array_  = dmrg_malloc<IntegerType>(sizeof(IntegerType) * batch_size_dim,
                                                           sizeof(IntegerType) * batch_size_dim);
	IntegerType* ldb_array_  = dmrg_malloc<IntegerType>(sizeof(IntegerType) * batch_size_dim,
                                                           sizeof(IntegerType) * batch_size_dim);
	IntegerType* ldc_array_  = dmrg_malloc<IntegerType>(sizeof(IntegerType) * batch_size_dim,
                                                           sizeof(IntegerType) * batch_size_dim);
	assert(m_array_ != nullptr);
	assert(n_array_ != nullptr);
	assert(k_array_ != nullptr);
	assert(group_size_ != nullptr);
	assert(lda_array_ != nullptr);
	assert(ldb_array_ != nullptr);
	assert(ldc_array_ != nullptr);

	SizeType nbytes_trans  = (sizeof(char) * ngroups_dim);
	char*    transa_array_ = dmrg_malloc<char>(nbytes_trans, nbytes_trans);
	char*    transb_array_ = dmrg_malloc<char>(nbytes_trans, nbytes_trans);
	assert(transa_array_ != nullptr);
	assert(transb_array_ != nullptr);

	IntegerType ibatch = 1;
	gflops1            = 0;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		T* BXbatch = gBXbatch_[ipatch - 1];
		T* Bbatch  = gBbatch_[ipatch - 1];

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
				T*          XJ    = &(X_[ix1 - 1]);
				IntegerType ld_XJ = nrowX;

				assert((1 <= nrowY) && (1 <= ncolY));
				assert((1 <= nrowX) && (1 <= ncolX));
				(void)nrowY;
				(void)ncolY;

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

					T alpha = 1;
					T beta  = 0;

					T* pA = Bbatch;
					T* pB = XJ;
					T* pC = BXbatch;

					IntegerType ld1 = ld_B;
					IntegerType ld2 = ld_XJ;
					IntegerType ld3 = ld_BX;

					transa_array_[ibatch - 1] = transA;
					transb_array_[ibatch - 1] = transB;

					m_array_[ibatch - 1] = mm;
					n_array_[ibatch - 1] = nn;
					k_array_[ibatch - 1] = kk;

					alpha_array_[ibatch - 1] = alpha;
					beta_array_[ibatch - 1]  = beta;

					lda_array_[ibatch - 1] = ld1;
					ldb_array_[ibatch - 1] = ld2;
					ldc_array_[ibatch - 1] = ld3;

					a_array_[ibatch - 1] = pA;
					b_array_[ibatch - 1] = pB;
					c_array_[ibatch - 1] = pC;

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
		SizeType i = 0;
		for (i = 1; i <= ngroups; i++) {
			group_size_[i - 1] = 1;
		};
	};

	/*
	 ------------------
	 first vbatch DGEMM
	 ------------------
	 */

	time_1st_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch<T>(transa_array_,
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
		T*          BX_           = pBXmem;
		IntegerType ngpu          = dmrg_get_ngpu();
		IntegerType igpu          = 0;
		size_t      BX_inc        = sum_BX_sizes / ngpu;
		size_t      BX_inc_nbytes = sizeof(T) * BX_inc;
		for (igpu = 0; igpu < ngpu; igpu++) {
			size_t offset = igpu * BX_inc;
			T*     pBX    = &(BX_[offset]);
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
		T* Abatch  = gAbatch_[ipatch - 1];
		T* BXbatch = gBXbatch_[ipatch - 1];

		IntegerType iy = xy_patch_start_[ipatch - 1];
		T*          YI = &(Y_[iy - 1]);

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
				(void)nrowX;

				IntegerType nrowBX = nrowB;
				IntegerType ncolBX = ncolX;

				assert((1 <= nrowBX) && (1 <= ncolBX));
				(void)nrowBX;

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

		T alpha = 1;
		T beta  = 0;

		T* pA = BXbatch;
		T* pB = Abatch;
		T* pC = YI;

		IntegerType ld1 = ld_BX;
		IntegerType ld2 = ld_A;
		IntegerType ld3 = ld_YI;

		transa_array_[ibatch - 1] = transA;
		transb_array_[ibatch - 1] = transB;
		m_array_[ibatch - 1]      = mm;
		n_array_[ibatch - 1]      = nn;
		k_array_[ibatch - 1]      = kk;

		alpha_array_[ibatch - 1] = alpha;
		beta_array_[ibatch - 1]  = beta;

		a_array_[ibatch - 1] = pA;
		b_array_[ibatch - 1] = pB;
		c_array_[ibatch - 1] = pC;

		lda_array_[ibatch - 1] = ld1;
		ldb_array_[ibatch - 1] = ld2;
		ldc_array_[ibatch - 1] = ld3;

		gflops2 += ((2.0 * mm) * nn) * kk;
		ibatch++;
	}; /* end for ipatch */

	ngroups = (ibatch - 1);
	assert(ngroups == npatches);
	{
		SizeType i = 0;
		for (i = 1; i <= ngroups; i++) {
			group_size_[i - 1] = 1;
		};
	}

	/*
	 ------------------
	 second vbatch DGEMM
	 ------------------
	 */
	time_2nd_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch<T>(transa_array_,
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
	dmrg_free(pBXmem);
	pBXmem = nullptr;

#ifdef USE_MAGMA
	/*
	 * -------------------
	 * free unified memory
	 * -------------------
	 */

	if (need_allocate_X) {
		dmrg_free(X_);
		X_ = nullptr;
	};
	if (need_allocate_Y) {
		void*  dest  = &(Yout_[0]);
		void*  src   = &(Y_[0]);
		size_t count = sizeof(T) * xy_size;
		dmrg_memcpy(dest, src, count);
		dmrg_free(Y_);
		Y_ = nullptr;
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

template void apply_Htarget_sparse<MYTYPE>(SizeType,
                                           SizeType,
                                           VectorSizeType,
                                           VectorSizeType,
                                           VectorSizeType,
                                           VectorSizeType,
                                           MYTYPE**,
                                           const VectorIntegerType&,
                                           MYTYPE**,
                                           const VectorIntegerType&,
                                           MYTYPE*,
                                           MYTYPE*);

#if defined(USE_COMPLEX_Z)

template void apply_Htarget_sparse<double>(SizeType,
                                           SizeType,
                                           VectorSizeType,
                                           VectorSizeType,
                                           VectorSizeType,
                                           VectorSizeType,
                                           double**,
                                           const VectorIntegerType&,
                                           double**,
                                           const VectorIntegerType&,
                                           double*,
                                           double*);
#else

template void apply_Htarget_sparse<std::complex<double>>(SizeType,
                                                         SizeType,
                                                         VectorSizeType,
                                                         VectorSizeType,
                                                         VectorSizeType,
                                                         VectorSizeType,
                                                         std::complex<double>**,
                                                         const VectorIntegerType&,
                                                         std::complex<double>**,
                                                         const VectorIntegerType&,
                                                         std::complex<double>*,
                                                         std::complex<double>*);
#endif
