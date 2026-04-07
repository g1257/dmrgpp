#include "setup_sparse_batch.h"
#include "dmrg_lapack.h"
#include "setup_nC.h"

#include <cassert>
#include <cmath>
#include <iostream>

template <typename T>
void setup_sparse_batch(
    SizeType noperator,
    SizeType npatches,
    /*
                                            ----------
                                            IntegerTypeent(in)
                                            ----------
                                            */
    const VectorSizeType&  left_patch_size_, /* sizes of left patches (INPUT) */
    const VectorSizeType&  right_patch_size_, /* sizes of right patches (INPUT) */
    const std::vector<T*>& Amatrix_, /* array of pointers to small Amat matrices */
    const VectorSizeType&  ld_Amatrix_, /* array of leading dimension of the small Amat matrices */
    const std::vector<T*>& Bmatrix_, /* array of pointers to small Bmat matrices */
    const VectorSizeType&  ld_Bmatrix_, /* array of leading dimension of the small Bmat matrices */
    /*
                                            -----------
                                            IntegerTypeent(out)
                                            -----------
                                            */
    VectorSizeType&    left_patch_start_, /* cumulative sum  of left_patch_size (OUTPUT) */
    VectorSizeType&    right_patch_start_, /* cumulative sum of right_patch_size (OUTPUT) */
    VectorSizeType&    xy_patch_start_, /* cumulative sum of product patch size (OUTPUT) */
    VectorSizeType&    nC_,
    T***               pgAbatch,
    VectorIntegerType& ld_gAbatch_,
    T***               pgBbatch,
    VectorIntegerType& ld_gBbatch_)
/*
 ---------------------------------------------------------
 returns the T *pointer to the Amat or Bmat matrices
 ---------------------------------------------------------
*/
{
	const IntegerType idebug     = 1;
	const IntegerType ialign     = 32;
	const IntegerType lfalse     = 0;
	const IntegerType ltrue      = !lfalse;
	const IntegerType use_Xlacpy = ltrue;
	const double      giga       = 1000.0 * 1000.0 * 1000.0;

	double total_time = -dmrg_get_wtime();

	size_t nbytes_Abatch = 0;
	size_t nbytes_Bbatch = 0;

	SizeType ipatch = 0;

	VectorSizeType gnnz_A_(npatches * npatches * noperator, 0);
	VectorSizeType gnnz_B_(npatches * npatches * noperator, 0);

	SizeType left_sum_state  = 0;
	SizeType right_sum_state = 0;
	{
		SizeType ipatch = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			left_sum_state += left_patch_size_[ipatch - 1];
			right_sum_state += right_patch_size_[ipatch - 1];

			assert(left_patch_size_[ipatch - 1] >= 1);
			assert(right_patch_size_[ipatch - 1] >= 1);
		};

		assert(left_sum_state >= 1);
		assert(right_sum_state >= 1);
	}
	if (idebug >= 1) {
		printf("setup_sparse_batch:npatches=%lu,noperator=%lu,", npatches, noperator);
		printf("left_sum_state=%lu,right_sum_state=%lu\n", left_sum_state, right_sum_state);
	};

	left_patch_start_.resize(npatches + 1);
	right_patch_start_.resize(npatches + 1);
	xy_patch_start_.resize(npatches + 1);
	nC_.resize(npatches * npatches);

	left_patch_start_[0] = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		left_patch_start_[ipatch]
		    = left_patch_start_[ipatch - 1] + left_patch_size_[ipatch - 1];
	};

	right_patch_start_[0] = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		right_patch_start_[ipatch]
		    = right_patch_start_[ipatch - 1] + right_patch_size_[ipatch - 1];
	};

	xy_patch_start_[0] = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		IntegerType nrowX       = right_patch_size_[ipatch - 1];
		IntegerType ncolX       = left_patch_size_[ipatch - 1];
		xy_patch_start_[ipatch] = xy_patch_start_[ipatch - 1] + (nrowX * ncolX);
	};

	/*
	 --------------------------------
	 setup nC(1:npatches, 1:npatches)
	 number of non-zero pairs of Amat(k), Bmat(k)
	 in block (ipatch,jpatch)
	 --------------------------------
*/

	/*
   IntegerType *nC_ = new IntegerType[npatches*npatches];
   assert( nC_ != NULL );
*/
	setup_nC(noperator,
	         npatches,
	         Amatrix_,
	         ld_Amatrix_,
	         Bmatrix_,
	         ld_Bmatrix_,
	         left_patch_size_,
	         right_patch_size_,
	         /*
	          * ------------
	          * IntegerTypeent(inout)
	          * ------------
	          */
	         nC_,
	         gnnz_A_,
	         gnnz_B_);
	/*
	                 nC_[indx2f(1,1,npatches)],
	                 gnnz_A_[index3(1,1,1,npatches)],
	                 gnnz_B_[index3(1,1,1,npatches)]
	                  );
*/

	/*
   --------------------------------------------------
   setup storage for gAbatch{ipatch}, gBbatch{ipatch}
   --------------------------------------------------
   */
	{
		SizeType ipatch = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			IntegerType nrowA = left_patch_size_[ipatch - 1];
			IntegerType nrowB = right_patch_size_[ipatch - 1];

			ld_gAbatch_[ipatch - 1] = ialign * ICEIL(nrowA, ialign);
			ld_gBbatch_[ipatch - 1] = ialign * ICEIL(nrowB, ialign);

			assert((1 <= nrowA) && (nrowA <= ld_gAbatch_[ipatch - 1]));
			assert((1 <= nrowB) && (nrowB <= ld_gBbatch_[ipatch - 1]));
		};
	};

	/*
   --------------------------------------------
   calculate the amount of storage needed per
   block row  gAbatch{ipatch},  gBbatch{ipatch}
   --------------------------------------------
  */

	VectorSizeType Abatch_sizes_(npatches, 0);
	VectorSizeType Bbatch_sizes_(npatches, 0);
	{
		SizeType ipatch = 0;
		SizeType jpatch = 0;

		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			Abatch_sizes_[ipatch - 1] = 0;
			Bbatch_sizes_[ipatch - 1] = 0;

			for (jpatch = 1; jpatch <= npatches; jpatch++) {
				SizeType nconnection = nC_[indx2f(ipatch, jpatch, npatches)];
				if (nconnection >= 1) {
					/*
			----------------------------
			accumulate number of columns
			----------------------------
			*/
					SizeType ncolA = left_patch_size_[jpatch - 1];
					SizeType ncolB = right_patch_size_[jpatch - 1];
					Abatch_sizes_[ipatch - 1] += (nconnection * ncolA);
					Bbatch_sizes_[ipatch - 1] += (nconnection * ncolB);
				};
			};
			/*
	  -----------------------
	  total number of entries
	  -----------------------
	  */
			Abatch_sizes_[ipatch - 1] *= ld_gAbatch_[ipatch - 1];
			Bbatch_sizes_[ipatch - 1] *= ld_gBbatch_[ipatch - 1];

		}; /* end for ipatch */
	};
	/*
   ------------------------------
   calculate total storage needed
   ------------------------------
   */
	SizeType sum_Abatch_sizes = 0;
	SizeType sum_Bbatch_sizes = 0;
	{
		SizeType ipatch = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			sum_Abatch_sizes += Abatch_sizes_[ipatch - 1];
			sum_Bbatch_sizes += Bbatch_sizes_[ipatch - 1];
		};
	};
	/*
  -----------------------------------
  allocate storage for Abatch, Bbatch
  Each Abatch{ipatch} is a contiguous block of memory
  Each Bbatch{ipatch} is a contiguous block of memory

  Moreover, all from a very large contiguous block of memory
  -----------------------------------
  */

	T** gAbatch_ = new T*[npatches];
	T** gBbatch_ = new T*[npatches];

	// printf("sizeof(FbType)=%lf \n",double(sizeof(T)));

	nbytes_Abatch = sizeof(T) * sum_Abatch_sizes;
	nbytes_Bbatch = sizeof(T) * sum_Bbatch_sizes;

	T* pAmem = new T[sum_Abatch_sizes];
	T* pBmem = new T[sum_Bbatch_sizes];
	assert(pAmem != NULL);
	assert(pBmem != NULL);

	{
		SizeType ipatch = 0;
		SizeType ipA    = 0;
		SizeType ipB    = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			// gAbatch_[ipatch-1] = &(pAmem[ipA]);
			gAbatch_[ipatch - 1] = pAmem + ipA;
			ipA += Abatch_sizes_[ipatch - 1];
			// gBbatch_[ipatch-1] = &(pBmem[ipB]);
			gBbatch_[ipatch - 1] = pBmem + ipB;
			ipB += Bbatch_sizes_[ipatch - 1];
		};
	};

	/*
  -------------------------
  fill in Abatch and Bbatch
  -------------------------
  */

	{

		SizeType ipatch    = 0;
		SizeType jpatch    = 0;
		SizeType ioperator = 0;

		/*
	--------------------------------
	Note the order of loops is important:

	for ipatch is  outer most loop
	for ioperator is inner most loop
	--------------------------------

	*/

		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			T* Adest = gAbatch_[ipatch - 1];
			T* Bdest = gBbatch_[ipatch - 1];

			for (jpatch = 1; jpatch <= npatches; jpatch++) {
				for (ioperator = 1; ioperator <= noperator; ioperator++) {

					T* Asrc
					    = Amatrix_[index3(ipatch, jpatch, ioperator, npatches)];
					T* Bsrc
					    = Bmatrix_[index3(ipatch, jpatch, ioperator, npatches)];

					SizeType nrowA = left_patch_size_[ipatch - 1];
					SizeType ncolA = left_patch_size_[jpatch - 1];

					SizeType nrowB = right_patch_size_[ipatch - 1];
					SizeType ncolB = right_patch_size_[jpatch - 1];

					SizeType nnz_Amat
					    = gnnz_A_[index3(ipatch, jpatch, ioperator, npatches)];
					SizeType nnz_Bmat
					    = gnnz_B_[index3(ipatch, jpatch, ioperator, npatches)];

					IntegerType has_work = (nnz_Amat >= 1) && (nnz_Bmat >= 1);
					if (has_work) {
						const char* uplo = " ";
						IntegerType m    = nrowA;
						IntegerType n    = ncolA;
						IntegerType ld1  = ld_Amatrix_[index3(
                                                    ipatch, jpatch, ioperator, npatches)];
						IntegerType ld2  = ld_gAbatch_[ipatch - 1];

						IntegerType isok = (1 <= m) && (1 <= n)
						    && (m <= ld1) && (m <= ld2);
						assert(isok);

						if (use_Xlacpy) {
							Xlacpy_(
							    uplo, &m, &n, Asrc, &ld1, Adest, &ld2);
						} else {
							dmrg_lacpy(
							    uplo, m, n, Asrc, ld1, Adest, ld2);
						};

						Adest += (ld2 * n);
					};

					if (has_work) {
						const char* uplo = " ";
						IntegerType m    = nrowB;
						IntegerType n    = ncolB;
						IntegerType ld1  = ld_Bmatrix_[index3(
                                                    ipatch, jpatch, ioperator, npatches)];
						IntegerType ld2  = ld_gBbatch_[ipatch - 1];

						IntegerType isok = (1 <= m) && (1 <= n)
						    && (m <= ld1) && (m <= ld2);
						assert(isok);

						if (use_Xlacpy) {
							Xlacpy_(
							    uplo, &m, &n, Bsrc, &ld1, Bdest, &ld2);
						} else {
							dmrg_lacpy(
							    uplo, m, n, Bsrc, ld1, Bdest, ld2);
						};

						Bdest += (ld2 * n);
					};
				}; /* end for ioperator */
			}; /* end for jpatch */

		}; /* end for ipatch */

#ifdef USE_MAGMA

		{
			IntegerType ngpu              = dmrg_get_ngpu();
			size_t      Abatch_inc        = (sum_Abatch_sizes) / ngpu;
			size_t      Bbatch_inc        = (sum_Bbatch_sizes) / ngpu;
			size_t      Abatch_inc_nbytes = sizeof(T) * Abatch_inc;
			size_t      Bbatch_inc_nbytes = sizeof(T) * Bbatch_inc;

			IntegerType igpu = 0;
			for (igpu = 0; igpu < ngpu; igpu++) {
				size_t Aoffset = igpu * Abatch_inc;
				size_t Boffset = igpu * Bbatch_inc;
				T*     pA      = &(pAmem[Aoffset]);
				T*     pB      = &(pBmem[Boffset]);

				dmrg_set_readonly((void*)pA, Abatch_inc_nbytes, igpu);
				dmrg_prefetch_to_device((void*)pA, Abatch_inc_nbytes, igpu);

				dmrg_set_readonly((void*)pB, Bbatch_inc_nbytes, igpu);
				dmrg_prefetch_to_device((void*)pB, Bbatch_inc_nbytes, igpu);
			};
		};
#endif
	};

	*pgAbatch = gAbatch_;
	*pgBbatch = gBbatch_;

	total_time += dmrg_get_wtime();
	if (idebug >= 1) {
		printf("setup_sparse_batch: total_time = %lf \n", total_time);
		printf("setup_sparse_batch:memory Abatch (%lf GBytes) Bbatch (%lf GBytes)\n",
		       (double)nbytes_Abatch / (giga),
		       (double)nbytes_Bbatch / (giga));
	};
}

template void setup_sparse_batch<MYTYPE>(
    SizeType noperator,
    SizeType npatches,
    /*
                                            ----------
                                            IntegerTypeent(in)
                                            ----------
                                            */
    const VectorSizeType&, /* sizes of left patches (INPUT) */
    const VectorSizeType&, /* sizes of right patches (INPUT) */
    const std::vector<MYTYPE*>&, /* array of pointers to small Amat matrices */
    const VectorSizeType&, /* array of leading dimension of the small Amat matrices */
    const std::vector<MYTYPE*>&, /* array of pointers to small Bmat matrices */
    const VectorSizeType&, /* array of leading dimension of the small Bmat matrices */
    /*
                                            -----------
                                            IntegerTypeent(out)
                                            -----------
                                            */
    VectorSizeType&, /* cumulative sum  of left_patch_size (OUTPUT) */
    VectorSizeType&, /* cumulative sum of right_patch_size (OUTPUT) */
    VectorSizeType&, /* cumulative sum of product patch size (OUTPUT) */
    VectorSizeType&,
    MYTYPE***,
    VectorIntegerType&,
    MYTYPE***,
    VectorIntegerType&);

template <typename T> void unsetup_sparse_batch(T*** pgAbatch, T*** pgBbatch)
{
	assert(pgAbatch != NULL);
	assert(pgBbatch != NULL);
	T** gAbatch = *pgAbatch;
	T** gBbatch = *pgBbatch;
	T*  pAmem   = gAbatch[0];
	T*  pBmem   = gBbatch[0];
	assert(pAmem != NULL);
	assert(pBmem != NULL);

	delete[] pAmem;
	delete[] pBmem;

	pAmem = nullptr;
	pBmem = nullptr;

	delete[] gAbatch;
	delete[] gBbatch;

	gAbatch = nullptr;
	gBbatch = nullptr;
}
template void unsetup_sparse_batch<MYTYPE>(MYTYPE***, MYTYPE***);
#if defined(USE_COMPLEX_Z)
template void unsetup_sparse_batch<double>(double***, double***);
#else
template void unsetup_sparse_batch<std::complex<double>>(std::complex<double>***,
                                                         std::complex<double>***);
#endif
