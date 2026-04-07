#include "setup_vbatch.h"
#include "dmrg_lapack.h"
#include "dmrg_types.h"
#include "dmrg_vbatch.h"
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#if defined(USE_COMPLEX_Z)
double ABS(std::complex<double> x) { return std::abs(x); }
#else
double ABS(double x) { return fabs(x); }
#endif

SizeType icall = 0;
void     print_nnz(SizeType              noperator,
                   SizeType              npatches,
                   SizeType              ncall,
                   const VectorSizeType& left_patch_size_,
                   const VectorSizeType& right_patch_size_,
                   VectorSizeType        gnnz_A_,
                   VectorSizeType        gnnz_B_)
{
	SizeType const MAX_STRLEN = 2 * 1024;

	char filename_A[MAX_STRLEN + 1];
	char filename_B[MAX_STRLEN + 1];

	snprintf(&(filename_A[0]), MAX_STRLEN, "Abatch.%05lu.m", ncall);
	snprintf(&(filename_B[0]), MAX_STRLEN, "Bbatch.%05lu.m", ncall);

	FILE* fid_A = fopen(filename_A, "w");
	FILE* fid_B = fopen(filename_B, "w");

	assert(fid_A != NULL);
	assert(fid_B != NULL);

	SizeType ipatch    = 0;
	SizeType jpatch    = 0;
	SizeType ioperator = 0;

	fprintf(fid_A, "left_patch_size = zeros(%lu,1);\n", npatches);
	fprintf(fid_B, "right_patch_size = zeros(%lu,1);\n", npatches);

	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		fprintf(
		    fid_A, "left_patch_size(%lu) = %lu;\n", ipatch, left_patch_size_[(ipatch)-1]);
		fprintf(
		    fid_B, "right_patch_size(%lu) = %lu;\n", ipatch, right_patch_size_[(ipatch)-1]);
	};

	fprintf(fid_A, "gnnz_A = zeros([%lu,%lu,%lu]);\n", npatches, npatches, noperator);
	fprintf(fid_B, "gnnz_B = zeros([%lu,%lu,%lu]);\n", npatches, npatches, noperator);

	for (ioperator = 1; ioperator <= noperator; ioperator++) {
		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			for (ipatch = 1; ipatch <= npatches; ipatch++) {
				if (gnnz_A_[index3(ipatch, jpatch, ioperator, npatches)] > 0) {
					fprintf(
					    fid_A,
					    "gnnz_A(%lu,%lu,%lu) = %lu;\n",
					    ipatch,
					    jpatch,
					    ioperator,
					    gnnz_A_[index3(ipatch, jpatch, ioperator, npatches)]);
				};

				if (gnnz_B_[index3(ipatch, jpatch, ioperator, npatches)] > 0) {
					fprintf(
					    fid_B,
					    "gnnz_B(%lu,%lu,%lu) = %lu;\n",
					    ipatch,
					    jpatch,
					    ioperator,
					    gnnz_B_[index3(ipatch, jpatch, ioperator, npatches)]);
				};
			};
		};
	};

	SizeType istat_A = fclose(fid_A);
	SizeType istat_B = fclose(fid_B);

	assert(istat_A == 0);
	assert(istat_B == 0);
}

template <typename T>
void setup_vbatch(
    SizeType               noperator, /* number of connections (INPUT) */
    SizeType               npatches, /* number of patches (INPUT) */
    const VectorSizeType&  left_patch_size_, /* sizes of left patches (INPUT) */
    const VectorSizeType&  right_patch_size_, /* sizes of right patches (INPUT) */
    VectorSizeType&        left_patch_start_, /* cumulative sum  of left_patch_size (OUTPUT) */
    VectorSizeType&        right_patch_start_, /* cumulative sum of right_patch_size (OUTPUT) */
    VectorSizeType&        xy_patch_start_, /* cumulative sum of product patch size (OUTPUT) */
    std::vector<T>&        pAbatch, /* Abatch matrix (OUTPUT) */
    SizeType&              pld_Abatch, /* leading dimension of Abatch matrix (OUTPUT) */
    std::vector<T>&        pBbatch, /* Bbatch matrix (OUTPUT) */
    SizeType&              pld_Bbatch, /*  leading dimension of Bbatch matrix (OUTPUT) */
    const std::vector<T*>& Amatrix_, /* array of pointers to small Amat matrices */
    const VectorSizeType&  ld_Amatrix_, /* array of leading dimension of the small Amat matrices */
    const std::vector<T*>& Bmatrix_, /* array of pointers to small Bmat matrices */
    const VectorSizeType&  ld_Bmatrix_ /* array of leading dimension of the small Bmat matrices */
)
/*
 ---------------------------------------------------------
 returns the T *pointer to the Amat or Bmat matrices
 ---------------------------------------------------------
*/
{
	const SizeType idebug = 2;

	const double giga          = 1000.0 * 1000.0 * 1000.0;
	SizeType     nbytes_Abatch = 0;
	SizeType     nbytes_Bbatch = 0;

	std::vector<double> nnz_A(noperator, 0.0);
	std::vector<double> nnz_B(noperator, 0.0);
	std::vector<double> nnz_Adiag(noperator, 0.0);
	std::vector<double> nnz_Bdiag(noperator, 0.0);
	VectorSizeType      gnnz_A_(npatches * npatches * noperator, 0);
	VectorSizeType      gnnz_B_(npatches * npatches * noperator, 0);

	if (idebug >= 1) {
		SizeType i = 0;
		for (i = 0; i < noperator; i++) {
			nnz_A[i]     = 0;
			nnz_B[i]     = 0;
			nnz_Adiag[i] = 0;
			nnz_Bdiag[i] = 0;
		};
	};

	SizeType left_max_state  = 0;
	SizeType right_max_state = 0;
	{
		SizeType ipatch = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			left_max_state += left_patch_size_[(ipatch)-1];
			right_max_state += right_patch_size_[(ipatch)-1];

			assert(left_patch_size_[(ipatch)-1] >= 1);
			assert(right_patch_size_[(ipatch)-1] >= 1);
		};
	}

	SizeType nrow_Abatch = left_max_state;
	SizeType ncol_Abatch = left_max_state * noperator;

	SizeType nrow_Bbatch = right_max_state;
	SizeType ncol_Bbatch = right_max_state * noperator;

	const SizeType ialign    = 32;
	SizeType       ld_Abatch = ialign * ICEIL(nrow_Abatch, ialign);
	SizeType       ld_Bbatch = ialign * ICEIL(nrow_Bbatch, ialign);

	nbytes_Abatch = (sizeof(T) * ld_Abatch * ncol_Abatch);
	nbytes_Bbatch = (sizeof(T) * ld_Bbatch * ncol_Bbatch);
	/*
	  T *Abatch_ = new T[ld_Abatch * ncol_Abatch];
	  T *Bbatch_ = new T[ld_Bbatch * ncol_Bbatch];

	  assert( Abatch_ != 0 );
	  assert( Bbatch_ != 0 );
	*/
	pAbatch.resize(ld_Abatch * ncol_Abatch, T(0.0));
	pBbatch.resize(ld_Abatch * ncol_Abatch, T(0.0));

	//  std::vector<T> Abatch_(ld_Abatch*ncol_Abatch,T(0.0));
	//  std::vector<T> Bbatch_(ld_Abatch*ncol_Abatch,T(0.0));

	SizeType ipatch = 0;

	left_patch_start_[0] = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		left_patch_start_[ipatch]
		    = left_patch_start_[(ipatch)-1] + left_patch_size_[(ipatch)-1];
	};

	right_patch_start_[0] = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		right_patch_start_[ipatch]
		    = right_patch_start_[(ipatch)-1] + right_patch_size_[(ipatch)-1];
	};

	xy_patch_start_[0] = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		SizeType nrowX          = right_patch_size_[(ipatch)-1];
		SizeType ncolX          = left_patch_size_[(ipatch)-1];
		xy_patch_start_[ipatch] = xy_patch_start_[(ipatch)-1] + (nrowX * ncolX);
	};

	/*
	 ---------------------------
	 set seed to be reproducible
	 ---------------------------
	 */
	{
		unsigned int iseed = 13;
		srand(iseed);
	}

	/*
	 -------------------------
	 fill in Abatch and Bbatch
	 -------------------------
	 */
	{
		const size_t size_Abatch = ((size_t)ld_Abatch) * left_max_state * noperator;
		const size_t size_Bbatch = ((size_t)ld_Bbatch) * right_max_state * noperator;

		assert(size_Abatch >= 1);
		assert(size_Bbatch >= 1);

		{
#ifdef USE_GETSET
			T* hAbatch_ = new T[size_Abatch];
#else
			std::vector<T*> hAbatch_(size_Abatch, nullptr);
			for (SizeType i = 0; i < size_Abatch; i++) {
				hAbatch_[i] = &(pAbatch[i]);
			}
#endif

			SizeType ipatch    = 0;
			SizeType jpatch    = 0;
			SizeType ioperator = 0;

			for (ioperator = 1; ioperator <= noperator; ioperator++) {
				for (jpatch = 1; jpatch <= npatches; jpatch++) {
					for (ipatch = 1; ipatch <= npatches; ipatch++) {
						T*       Asrc      = Amatrix_[index3(
                                                    ipatch, jpatch, ioperator, npatches)];
						SizeType ld_Asrc   = ld_Amatrix_[index3(
                                                    ipatch, jpatch, ioperator, npatches)];
						SizeType nrow_Asrc = left_patch_size_[(ipatch)-1];
						SizeType ncol_Asrc = left_patch_size_[(jpatch)-1];

						SizeType ia = left_patch_start_[(ipatch)-1];
						SizeType ja = left_patch_start_[(jpatch)-1];

						T*       Adest    = hAbatch_[indx2f(
                                                    ia,
                                                    ja + (ioperator - 1) * left_max_state,
                                                    ld_Abatch)];
						SizeType ld_Adest = ld_Abatch;

						{
							char        uplo = ' ';
							IntegerType m    = nrow_Asrc;
							IntegerType n    = ncol_Asrc;
							IntegerType ld1  = ld_Asrc;
							IntegerType ld2  = ld_Adest;

							SizeType isok = (1 <= m) && (1 <= n)
							    && (m <= ld1) && (m <= ld2);
							if (!isok) {
								printf("setup_vbatch:ipatch=%lu,"
								       "jpatch=%lu,ioperator=%lu\n",
								       ipatch,
								       jpatch,
								       ioperator);
								printf("npatches=%lu,noperator=%lu,"
								       " left_max_state=%lu\n",
								       (SizeType)npatches,
								       (SizeType)noperator,
								       (SizeType)left_max_state);
								printf("m=%d,n=%d,ld1=%d,ld2=%d\n",
								       m,
								       n,
								       ld1,
								       ld2);
								printf("ia=%lu, ja=%lu\n", ia, ja);
							};

							assert(m >= 1);
							assert(n >= 1);
							assert(m <= ld1);
							assert(m <= ld2);

							//  ------------------------------------------
							//  copy m by n submatrix from Asrc  to
							//  Adest
							//  ------------------------------------------

							Xlacpy_(
							    &uplo, &m, &n, Asrc, &ld1, Adest, &ld2);

							if (idebug >= 1) {
								SizeType     i      = 0;
								SizeType     j      = 0;
								double       lnnz_A = 0;
								const double tiny
								    = 1.0 / pow(2, 50);
								for (j = 0; j < SizeType(n); j++) {
									for (i = 0; i < SizeType(m);
									     i++) {
										double abs_aij = ABS(
										    Asrc[i
										         + j * ld1]);
										lnnz_A += (abs_aij
										           < tiny)
										    ? 0
										    : 1;
									};
								};
								nnz_A[ioperator - 1] += lnnz_A;

								if (ipatch == jpatch) {
									nnz_Adiag[ioperator - 1]
									    += lnnz_A;
								};

								gnnz_A_[index3(ipatch,
								               jpatch,
								               ioperator,
								               npatches)]
								    = lnnz_A;
							};
						};
					};
				};
			};

#ifdef USE_GETSET
			{
				const SizeType incx = 1;
				const SizeType incy = 1;

				SizeType       n      = size_Abatch;
				SizeType       istart = 1;
				const SizeType nb     = 1024 * 1024 * 1024;

				for (istart = 1; istart <= n; istart += nb) {
					SizeType iend  = MIN(n, istart + nb - 1);
					SizeType isize = (iend - istart + 1);
					T*       src   = &(hAbatch_[istart - 1]);
					T*       dest  = &(pAbatch[istart - 1]);

					dmrg_Xsetvector(isize, src, incx, dest, incy);
				};
			}
			delete[] hAbatch_;
#endif
		};

		{
#ifdef USE_GETSET

			T* hBbatch_ = new T[size_Bbatch];
#else
			std::vector<T*> hBbatch_(size_Bbatch, nullptr);
			for (SizeType i = 0; i < size_Bbatch; i++) {
				hBbatch_[i] = &(pBbatch[i]);
			}
#endif

			SizeType ioperator = 0;
			for (ioperator = 1; ioperator <= noperator; ioperator++) {
				SizeType jpatch = 0;
				for (jpatch = 1; jpatch <= npatches; jpatch++) {
					SizeType ipatch = 0;
					for (ipatch = 1; ipatch <= npatches; ipatch++) {

						T*       Bsrc      = Bmatrix_[index3(
                                                    ipatch, jpatch, ioperator, npatches)];
						SizeType ld_Bsrc   = ld_Bmatrix_[index3(
                                                    ipatch, jpatch, ioperator, npatches)];
						SizeType nrow_Bsrc = right_patch_size_[(ipatch)-1];
						SizeType ncol_Bsrc = right_patch_size_[(jpatch)-1];

						SizeType ib = right_patch_start_[(ipatch)-1];
						SizeType jb = right_patch_start_[(jpatch)-1];

						T*       Bdest    = hBbatch_[indx2f(
                                                    ib,
                                                    jb + (ioperator - 1) * right_max_state,
                                                    ld_Bbatch)];
						SizeType ld_Bdest = ld_Bbatch;

						{
							char        uplo = ' ';
							IntegerType m    = nrow_Bsrc;
							IntegerType n    = ncol_Bsrc;
							IntegerType ld1  = ld_Bsrc;
							IntegerType ld2  = ld_Bdest;

							assert(1 <= m);
							assert(1 <= n);
							assert(m <= ld1);
							assert(m <= ld2);

							Xlacpy_(
							    &uplo, &m, &n, Bsrc, &ld1, Bdest, &ld2);

							if (idebug >= 1) {
								SizeType     i      = 0;
								SizeType     j      = 0;
								double       lnnz_B = 0;
								const double tiny
								    = 1.0 / pow(2.0, 50);
								for (j = 0; j < SizeType(n); j++) {
									for (i = 0; i < SizeType(m);
									     i++) {
										double abs_bij = ABS(
										    Bsrc[i
										         + j * ld1]);
										lnnz_B += (abs_bij
										           < tiny)
										    ? 0
										    : 1;
									};
								};
								nnz_B[ioperator - 1] += lnnz_B;
								if (ipatch == jpatch) {
									nnz_Bdiag[ioperator - 1]
									    += lnnz_B;
								};
								gnnz_B_[index3(ipatch,
								               jpatch,
								               ioperator,
								               npatches)]
								    = lnnz_B;
							};
						};
					};
				};
			};

#ifdef USE_GETSET
			{
				const SizeType incx = 1;
				const SizeType incy = 1;

				SizeType       n      = size_Bbatch;
				SizeType       istart = 0;
				const SizeType nb     = 1024 * 1024 * 1024;

				for (istart = 1; istart <= n; istart += nb) {
					SizeType iend  = MIN(n, istart + nb - 1);
					SizeType isize = (iend - istart + 1);
					T*       src   = &(hBbatch_[istart - 1]);
					T*       dest  = &(pBbatch[istart - 1]);

					dmrg_Xsetvector(isize, src, incx, dest, incy);
				};
			};

			delete[] hBbatch_;
#endif
		};
	}

	if (idebug >= 1) {
		/*
		 * -------------
		 * estimate work for Y = B * X * transpose(A);
		 * -------------
		 */

		double total_flops_method_1 = 0;
		double total_flops_method_2 = 0;
		double total_flops_min      = 0;

		SizeType ipatch = 0;
		SizeType jpatch = 0;

		SizeType left_sum  = 0;
		SizeType right_sum = 0;
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			SizeType left_size  = left_patch_size_[(ipatch)-1];
			SizeType right_size = right_patch_size_[(ipatch)-1];

			left_sum += left_size;
			right_sum += right_size;

			assert(1 <= left_size);
			assert(1 <= right_size);
		};
		assert(left_sum == left_max_state);
		assert(right_sum == right_max_state);

		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			for (ipatch = 1; ipatch <= npatches; ipatch++) {

				SizeType nrowA = left_patch_size_[(ipatch)-1];
				SizeType ncolA = left_patch_size_[(jpatch)-1];
				SizeType nrowB = right_patch_size_[(ipatch)-1];
				SizeType ncolB = right_patch_size_[(jpatch)-1];

				SizeType nrowX = ncolB;
				SizeType ncolX = ncolA;
				SizeType nrowY = nrowB;
				SizeType ncolY = nrowA;

				SizeType nrowBX = nrowB;
				SizeType ncolBX = ncolX;

				/*
				 ---------------------------------------------------------
				 Method 1:  (i) BX = B * X,    (ii)  Y = BX * transpose(A)
				 ---------------------------------------------------------
				 */
				double flops_bx       = (2.0 * (nrowBX * ncolBX)) * ncolB;
				double flops_bx_a     = (2.0 * (nrowY * ncolY)) * ncolBX;
				double flops_method_1 = flops_bx + flops_bx_a;

				/*
				 ---------------------------------------------------------
				 Method 2:  (i) XAt = X * transpose(A),    (ii)  Y = B * XAt
				 ---------------------------------------------------------
				 */
				double flops_xat      = (2.0 * (nrowX * nrowA)) * ncolX;
				double flops_b_xat    = (2.0 * (nrowY * ncolY)) * ncolB;
				double flops_method_2 = flops_xat + flops_b_xat;

				total_flops_method_1 += flops_method_1;
				total_flops_method_2 += flops_method_2;
				total_flops_min += MIN(flops_method_1, flops_method_2);
			};
		};
		total_flops_method_1 *= ((double)noperator);
		total_flops_method_2 *= ((double)noperator);
		total_flops_min *= ((double)noperator);

		printf("npatches=%lu, noperator=%lu left_max_state=%lu right_max_state=%lu\n",
		       npatches,
		       noperator,
		       left_max_state,
		       right_max_state);

		{
			SizeType ioperator = 0;
			for (ioperator = 0; ioperator < noperator; ioperator++) {

				double sparsity_ratio_A
				    = ((nnz_A[ioperator] / left_max_state) / left_max_state);
				double sparsity_ratio_B
				    = ((nnz_B[ioperator] / right_max_state) / right_max_state);

				double fraction_in_diag_A = nnz_Adiag[ioperator] / nnz_A[ioperator];
				double fraction_in_diag_B = nnz_Bdiag[ioperator] / nnz_B[ioperator];

				printf("ioperator=%lu ", ioperator);
				printf("sparsity ratio A=%6.4f, fraction in diag blocks=%6.4f ",
				       sparsity_ratio_A,
				       fraction_in_diag_A);
				printf("sparsity ratio B=%6.4f, fraction in diag blocks=%6.4f ",
				       sparsity_ratio_B,
				       fraction_in_diag_B);
				printf(" \n");
			};
		};

		printf("total_flops_method_1=%le (%6.2lf), ",
		       total_flops_method_1,
		       total_flops_method_1 / total_flops_min);
		printf("total_flops_method_2=%le (%6.2lf), ",
		       total_flops_method_2,
		       total_flops_method_2 / total_flops_min);
		printf("total_flops_min=%le \n", total_flops_min);

		/*
		 * ----------------------------
		 * write out number of nonzeros
		 * ----------------------------
		 */
		if (idebug >= 2) {
			print_nnz(noperator,
			          npatches,
			          icall,
			          left_patch_size_,
			          right_patch_size_,
			          gnnz_A_,
			          gnnz_B_);
			icall++;
		};
	};

	if (idebug >= 1) {
		printf("setup_vbatch:memory Abatch (%f GBytes) Bbatch (%f GBytes)\n",
		       (double)nbytes_Abatch / (giga),
		       (double)nbytes_Bbatch / (giga));
	};

	pld_Abatch = ld_Abatch;
	pld_Bbatch = ld_Bbatch;
}

int ICEIL2(const SizeType& x, const SizeType& n) { return (((x) + (n)-1) / (n)); }
/*
int indx2f2(const SizeType& i,const SizeType& j,const SizeType& lda)
{
        return (((i)-1) + ((j)-1)*(lda));
}
*/
template void setup_vbatch<MYTYPE>(
    SizeType              noperator, /* number of connections (INPUT) */
    SizeType              npatches, /* number of patches (INPUT) */
    const VectorSizeType& left_patch_size_, /* sizes of left patches (INPUT) */
    const VectorSizeType& right_patch_size_, /* sizes of right patches (INPUT) */
    VectorSizeType&       left_patch_start_, /* cumulative sum  of left_patch_size (OUTPUT) */
    VectorSizeType&       right_patch_start_, /* cumulative sum of right_patch_size (OUTPUT) */
    VectorSizeType&       xy_patch_start_, /* cumulative sum of product patch size (OUTPUT) */
    std::vector<MYTYPE>&  pAbatch, /* Abatch matrix (OUTPUT) */
    SizeType&             ld_Abatch, /* leading dimension of Abatch matrix (OUTPUT) */
    std::vector<MYTYPE>&  pBbatch, /* Bbatch matrix (OUTPUT) */
    SizeType&             ld_Bbatch, /*  leading dimension of Bbatch matrix (OUTPUT) */
    const std::vector<MYTYPE*>& Amatrix_, /* array of pointers to small Amat matrices */
    const VectorSizeType& ld_Amatrix_, /* array of leading dimension of the small Amat matrices */
    const std::vector<MYTYPE*>& Bmatrix_, /* array of pointers to small Bmat matrices */
    const VectorSizeType& ld_Bmatrix_ /* array of leading dimension of the small Bmat matrices */
);
