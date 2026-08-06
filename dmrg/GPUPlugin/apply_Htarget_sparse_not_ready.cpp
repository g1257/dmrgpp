#include "GPUPluginConfig.h"
#include "dmrg_types.h"
#include "dmrg_vbatch.h"

#include <Kokkos_Core.hpp>
#include <PsimagLite/KokkosType.h>

#include <algorithm>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T, typename KokkosScalar>
void apply_Htarget_sparse(SizeType                                               noperator,
                          SizeType                                               npatches,
                          VectorSizeType                                         left_patch_start_,
                          VectorSizeType                                         right_patch_start_,
                          VectorSizeType                                         xy_patch_start_,
                          VectorSizeType                                         nC_,
                          T**                                                    gAbatch_,
                          const VectorIntegerType&                               ld_gAbatch_,
                          T**                                                    gBbatch_,
                          const VectorIntegerType&                               ld_gBbatch_,
                          Kokkos::View<const KokkosScalar*, Kokkos::SharedSpace> Xin_,
                          Kokkos::View<KokkosScalar*, Kokkos::SharedSpace>       Yout_)
{
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

	Kokkos::DefaultExecutionSpace exec;
	Kokkos::SharedSpace           mem;

	auto X_ = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Xin_);
	auto Y_ = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Yout_);

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

	{
		IntegerType ngpu         = dmrg_get_ngpu();
		IntegerType igpu         = 0;
		size_t      X_inc        = xy_size_dim / ngpu;
		size_t      X_inc_nbytes = sizeof(T) * X_inc;
		for (igpu = 0; igpu < ngpu; igpu++) {
			size_t              offset = igpu * X_inc;
			const KokkosScalar* pX     = &(X_[offset]);
			dmrg_set_readonly((void*)pX, X_inc_nbytes, igpu);
			dmrg_prefetch_to_device((void*)pX, X_inc_nbytes, igpu);
		};
	}
	Kokkos::deep_copy(exec, Y_, Xin_);
#endif

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
	KokkosScalar** gBXbatch_ = reinterpret_cast<KokkosScalar**>(
	    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(npatches * sizeof(KokkosScalar*)));

	nbytes_BX = sizeof(T) * sum_BX_sizes;
	Kokkos::View<KokkosScalar*, Kokkos::SharedSpace> pBXmem("pBXmem", sum_BX_sizes);

	{
		KokkosScalar* BXbatch = pBXmem.data();
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

	T alpha = 1.;
	T beta  = 0.;

	KokkosScalar const** a_array_ = reinterpret_cast<KokkosScalar const**>(
	    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(sizeof(KokkosScalar*) * batch_size_dim));
	KokkosScalar const** b_array_ = reinterpret_cast<KokkosScalar const**>(
	    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(sizeof(KokkosScalar*) * batch_size_dim));
	KokkosScalar** c_array_ = reinterpret_cast<KokkosScalar**>(
	    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(sizeof(KokkosScalar*) * batch_size_dim));
	assert(a_array_ != nullptr);
	assert(b_array_ != nullptr);
	assert(c_array_ != nullptr);

	Kokkos::View<IntegerType*, Kokkos::SharedSpace> m_array_("m_array", ngroups_dim);
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> n_array_("n_array", ngroups_dim);
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> k_array_("k_array", ngroups_dim);
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> group_size_("group_size", ngroups_dim);

	Kokkos::View<IntegerType*, Kokkos::SharedSpace> lda_array_("lda_array", batch_size_dim);
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> ldb_array_("ldb_array", batch_size_dim);
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> ldc_array_("ldc_array", batch_size_dim);

	IntegerType ibatch = 1;
	gflops1            = 0;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		KokkosScalar* BXbatch = gBXbatch_[ipatch - 1];
		T*            Bbatch  = gBbatch_[ipatch - 1];

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

				IntegerType         ix1   = xy_patch_start_[jpatch - 1];
				const KokkosScalar* XJ    = &(X_[ix1 - 1]);
				IntegerType         ld_XJ = nrowX;

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

					IntegerType mm = nrowBX;
					IntegerType nn = ncolBX;
					IntegerType kk = ncolB;

					T*                  pA = Bbatch;
					const KokkosScalar* pB = XJ;
					KokkosScalar*       pC = BXbatch;

					IntegerType ld1 = ld_B;
					IntegerType ld2 = ld_XJ;
					IntegerType ld3 = ld_BX;

					m_array_[ibatch - 1] = mm;
					n_array_[ibatch - 1] = nn;
					k_array_[ibatch - 1] = kk;

					lda_array_[ibatch - 1] = ld1;
					ldb_array_[ibatch - 1] = ld2;
					ldc_array_[ibatch - 1] = ld3;

					a_array_[ibatch - 1]
					    = reinterpret_cast<const KokkosScalar*>(pA);
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

	char transA     = 'N';
	char transB     = 'N';
	time_1st_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch<T>(transA,
	                     transB,
	                     m_array_.data(),
	                     n_array_.data(),
	                     k_array_.data(),
	                     alpha,
	                     reinterpret_cast<const T**>(a_array_),
	                     lda_array_.data(),
	                     reinterpret_cast<const T**>(b_array_),
	                     ldb_array_.data(),
	                     beta,
	                     reinterpret_cast<T**>(c_array_),
	                     ldc_array_.data(),
	                     ngroups,
	                     group_size_.data());
	time_1st_vbatch += dmrg_get_wtime();
	gflops1 = gflops1 / giga;

#ifdef USE_MAGMA
	{
		KokkosScalar* BX_           = pBXmem.data();
		IntegerType   ngpu          = dmrg_get_ngpu();
		IntegerType   igpu          = 0;
		size_t        BX_inc        = sum_BX_sizes / ngpu;
		size_t        BX_inc_nbytes = sizeof(T) * BX_inc;
		for (igpu = 0; igpu < ngpu; igpu++) {
			size_t        offset = igpu * BX_inc;
			KokkosScalar* pBX    = BX_ + offset;
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
		T*            Abatch  = gAbatch_[ipatch - 1];
		KokkosScalar* BXbatch = gBXbatch_[ipatch - 1];

		IntegerType   iy = xy_patch_start_[ipatch - 1];
		KokkosScalar* YI = &(Y_[iy - 1]);

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

		KokkosScalar* pA = BXbatch;
		T*            pB = Abatch;
		KokkosScalar* pC = YI;

		IntegerType ld1 = ld_BX;
		IntegerType ld2 = ld_A;
		IntegerType ld3 = ld_YI;

		m_array_[ibatch - 1] = mm;
		n_array_[ibatch - 1] = nn;
		k_array_[ibatch - 1] = kk;

		a_array_[ibatch - 1] = reinterpret_cast<KokkosScalar*>(pA);
		b_array_[ibatch - 1] = reinterpret_cast<KokkosScalar*>(pB);
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
	transA          = 'N';
	transB          = 'T';
	time_2nd_vbatch = -dmrg_get_wtime();
	dmrg_Xgemm_vbatch<T>(transA,
	                     transB,
	                     m_array_.data(),
	                     n_array_.data(),
	                     k_array_.data(),
	                     alpha,
	                     reinterpret_cast<const T**>(a_array_),
	                     lda_array_.data(),
	                     reinterpret_cast<const T**>(b_array_),
	                     ldb_array_.data(),
	                     beta,
	                     reinterpret_cast<T**>(c_array_),
	                     ldc_array_.data(),
	                     ngroups,
	                     group_size_.data());
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

	Kokkos::kokkos_free<Kokkos::SharedSpace>(a_array_);
	Kokkos::kokkos_free<Kokkos::SharedSpace>(b_array_);
	Kokkos::kokkos_free<Kokkos::SharedSpace>(c_array_);
	Kokkos::kokkos_free<Kokkos::SharedSpace>(gBXbatch_);
	Kokkos::deep_copy(exec, Yout_, Y_);

	total_time += dmrg_get_wtime();
	if (idebug >= 1) {
		printf("apply_Htarget_sparse: total_time=%lf \n", total_time);
		printf("apply_Htarget_sparse:memory BX (%f GBytes) X (%f GBytes) Y (%f GBytes) \n",
		       (double)nbytes_BX / (giga),
		       (double)nbytes_X / (giga),
		       (double)nbytes_Y / (giga));
	};
}

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
                                           Kokkos::View<const double*, Kokkos::SharedSpace>,
                                           Kokkos::View<double*, Kokkos::SharedSpace>);

template void apply_Htarget_sparse<std::complex<double>>(
    SizeType,
    SizeType,
    VectorSizeType,
    VectorSizeType,
    VectorSizeType,
    VectorSizeType,
    std::complex<double>**,
    const VectorIntegerType&,
    std::complex<double>**,
    const VectorIntegerType&,
    Kokkos::View<const Kokkos::complex<double>*, Kokkos::SharedSpace>,
    Kokkos::View<Kokkos::complex<double>*, Kokkos::SharedSpace>);
