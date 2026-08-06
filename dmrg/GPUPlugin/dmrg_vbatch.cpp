#include "dmrg_vbatch.h"
#include "dmrg_lapack.h"
#include "dmrg_types.h"

#include "PsimagLite/KokkosType.h"
#include <Kokkos_Core.hpp>

static IntegerType is_initialized = 0;

#ifdef USE_MAGMA

#include "dmrg_magma.h"

#ifndef MAXGPUS
#define MAXGPUS 8

#endif

static IntegerType    max_gpus = MAXGPUS;
static magma_queue_t  queue_array[MAXGPUS];
static magma_device_t device_array[MAXGPUS];
static IntegerType    ngpu = 0;

static magma_queue_t queue  = 0;
static IntegerType   device = 0;

#endif

#ifdef _OPENMP
#include <omp.h>
double dmrg_get_wtime() { return omp_get_wtime(); }
#else
#include <time.h>
double dmrg_get_wtime() { return ((double)clock()) / CLOCKS_PER_SEC; }
#endif

IntegerType ICEIL(const IntegerType& x, const IntegerType& n)
{
	IntegerType value = (((x) + (n)-1) / (n));
	return value;
}

SizeType ICEIL(const SizeType& x, const SizeType& n) { return (((x) + (n)-1) / (n)); }

SizeType indx2f(const SizeType& i, const SizeType& j, const SizeType& lda)
{
	return (((i)-1) + ((j)-1) * (lda));
}

SizeType index3(const SizeType& ipatch,
                const SizeType& jpatch,
                const SizeType& ioperator,
                const SizeType& npatches)
{
	SizeType myint = ((((ipatch)-1) + (((jpatch)-1) * (npatches)))
	                  + (((ioperator)-1) * ((npatches) * (npatches))));
	return myint;
}

void dmrg_finalize()
{
#ifdef USE_MAGMA
	IntegerType idev = 0;
	for (idev = 0; idev < MAXGPUS; idev++) {
		if (queue_array[idev] != 0) {
			magma_queue_destroy(queue_array[idev]);
		};
	};
	magma_finalize();
#endif
}

#ifdef USE_MAGMA
IntegerType dmrg_get_ngpu() { return ngpu; };

IntegerType dmrg_set_max_gpus(const IntegerType max_gpus_in)
{
	const IntegerType prev_max_gpus = max_gpus;
	if (max_gpus > max_gpus_in) {
		max_gpus = max_gpus_in;
	};
	return (prev_max_gpus);
}
#endif

void dmrg_init()
{
	if (!is_initialized) {
		is_initialized = 1;

#ifdef USE_MAGMA
		const IntegerType idebug = 1;

		device = 0;
		magma_init();

		magma_getdevices(device_array, MAXGPUS, &ngpu);
		assert(ngpu >= 1);

		if (idebug >= 1) {
			printf("dmrg_init: ngpu = %d detected \n", ngpu);
		};

		max_gpus = (max_gpus >= MAXGPUS) ? MAXGPUS : max_gpus;
		ngpu     = (ngpu >= max_gpus) ? max_gpus : ngpu;
		if (idebug >= 1) {
			printf("dmrg_init: ngpu = %d used \n", ngpu);
		};

		IntegerType idev = 0;
		for (idev = 0; idev < ngpu; idev++) {
			device = device_array[idev];

			magma_setdevice(device);
			magma_queue_create(device, &queue);
			assert(queue != 0);
			queue_array[idev] = queue;
		};

		idev   = 0;
		queue  = queue_array[idev];
		device = device_array[idev];

#endif
	};
}

template <typename T>
void dmrg_Xgetvector(const IntegerType n,
                     T*                dx_src,
                     const IntegerType incx,
                     T*                hy_dst,
                     const IntegerType incy)
{
#ifdef USE_MAGMA
	magma_Xgetvector(n, dx_src, incx, hy_dst, incy, queue);
#else
	Xcopy_(&n, dx_src, &incx, hy_dst, &incy);
#endif
}

template <typename T>
void dmrg_Xsetvector(const IntegerType n,
                     const T*          hx_src,
                     const IntegerType incx,
                     T*                dy_dst,
                     const IntegerType incy)
{
#ifdef USE_MAGMA
	magma_Xsetvector(n, hx_src, incx, dy_dst, incy, queue);
#else
	Xcopy_(&n, hx_src, &incx, dy_dst, &incy);
#endif
}

template <typename T>
void dmrg_Xgetmatrix(const IntegerType m,
                     const IntegerType n,
                     T*                dA_src,
                     const IntegerType ldda,
                     T*                hB_dst,
                     const IntegerType ldb)
{
#ifdef USE_MAGMA
	magma_Xgetmatrix(m, n, dA_src, ldda, hB_dst, ldb, queue);
#else
	const char* uplo = "A";
	Xlacpy_(uplo, &m, &n, dA_src, &ldda, hB_dst, &ldb);
#endif
}

template <typename T>
void dmrg_Xsetmatrix(const IntegerType m,
                     const IntegerType n,
                     const T*          hA_src,
                     const IntegerType lda,
                     T*                dB_dst,
                     const IntegerType lddb)
{
#ifdef USE_MAGMA
	magma_Xsetmatrix(m, n, hA_src, lda, dB_dst, lddb, queue);
#else
	const char* uplo = "A";
	Xlacpy_(uplo, &m, &n, hA_src, &lda, dB_dst, &lddb);
#endif
}

template <typename T>
void dmrg_Xgemm_vbatch(char         ctransa,
                       char         ctransb,
                       IntegerType* m_array,
                       IntegerType* n_array,
                       IntegerType* k_array,
                       T            alpha,
                       const T**    a_array,
                       IntegerType* lda_array,
                       const T**    b_array,
                       IntegerType* ldb_array,
                       T            beta,
                       T**          c_array,
                       IntegerType* ldc_array,
                       SizeType     group_count,
                       IntegerType* group_size)
{
	using KokkosScalar = PsimagLite::KokkosType<T>::type;

	const IntegerType idebug       = 1;
	const double      giga         = 1000.0 * 1000.0 * 1000.0;
	double            gflops       = 0;
	double            elapsed_time = 0;

	size_t nbytes       = 0;
	size_t nbytes_total = 0;

	if (idebug >= 1) {
		elapsed_time    = -dmrg_get_wtime();
		SizeType igroup = 0;

		for (igroup = 0; igroup < group_count; igroup++) {

			gflops += ((double)m_array[igroup]) * ((double)n_array[igroup])
			    * ((double)k_array[igroup]) * ((double)group_size[igroup]) * 2.0;
		};
		gflops = gflops / (giga);
	};

	/*
	---------------------
	expand out the groups
	---------------------
	*/
	IntegerType batch_size = 0;
	for (SizeType igroup = 0; igroup < group_count; igroup++) {
		batch_size += group_size[igroup];
	};

	/*
   ----------------------------------------------------------
   Note magma_Xgemm_vbatched need arrays of size batch_size+1
   Set vbatch_dim to be multiple of 32 for better alignment in memory
   ----------------------------------------------------------
   */
	const IntegerType ialign     = 32;
	IntegerType       vbatch_dim = ialign * ICEIL((batch_size + 1), ialign);

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> m_vbatch("m_vbatch", vbatch_dim);

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> n_vbatch("n_vbatch", vbatch_dim);

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> k_vbatch("k_vbatch", vbatch_dim);

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> lda_vbatch("lda_vbatch", vbatch_dim);

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> ldb_vbatch("ldb_vbatch", vbatch_dim);

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	Kokkos::View<IntegerType*, Kokkos::SharedSpace> ldc_vbatch("ldc_vbatch", vbatch_dim);

	nbytes = sizeof(KokkosScalar*) * (vbatch_dim);
	nbytes_total += nbytes;
	const KokkosScalar** a_vbatch = reinterpret_cast<const KokkosScalar**>(
	    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(nbytes));

	nbytes = sizeof(KokkosScalar*) * (vbatch_dim);
	nbytes_total += nbytes;
	const KokkosScalar** b_vbatch = reinterpret_cast<const KokkosScalar**>(
	    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(nbytes));

	nbytes = sizeof(KokkosScalar*) * (vbatch_dim);
	nbytes_total += nbytes;
	KokkosScalar** c_vbatch
	    = reinterpret_cast<KokkosScalar**>(Kokkos::kokkos_malloc<Kokkos::SharedSpace>(nbytes));

	IntegerType idx = 0;
	for (SizeType igroup = 0; igroup < group_count; igroup++) {
		for (IntegerType i = 0; i < group_size[igroup]; i++) {
			m_vbatch[idx] = m_array[igroup];
			n_vbatch[idx] = n_array[igroup];
			k_vbatch[idx] = k_array[igroup];

			lda_vbatch[idx] = lda_array[igroup];
			ldb_vbatch[idx] = ldb_array[igroup];
			ldc_vbatch[idx] = ldc_array[igroup];

			a_vbatch[idx] = reinterpret_cast<const KokkosScalar*>(a_array[idx]);
			b_vbatch[idx] = reinterpret_cast<const KokkosScalar*>(b_array[idx]);
			c_vbatch[idx] = reinterpret_cast<KokkosScalar*>(c_array[idx]);

			idx = idx + 1;
		};
	};

#ifdef USE_MAGMA
	/*
	---------------------------------------------------------------
	Note magma_Xgemm_vbatched assumes all transa, transb, alpha, beta are the same
	---------------------------------------------------------------
	*/
	IntegerType isTransA = ((ctransa == 'T') || (ctransa == 't'));
	IntegerType isTransB = ((ctransb == 'T') || (ctransb == 't'));

	magma_trans_t transA = isTransA ? MagmaTrans : MagmaNoTrans;
	magma_trans_t transB = isTransB ? MagmaTrans : MagmaNoTrans;
	{

		/*
		 * ------------
		 * extra checks
		 * ------------
		 */
		if (idebug >= 1) {
			IntegerType i = 0;
			for (i = 0; i < batch_size; i++) {
				IntegerType mm = m_vbatch[i];
				IntegerType nn = n_vbatch[i];
				IntegerType kk = k_vbatch[i];

				IntegerType lda = lda_vbatch[i];
				IntegerType ldb = ldb_vbatch[i];
				IntegerType ldc = ldc_vbatch[i];

				assert(a_vbatch[i] != nullptr);
				assert(b_vbatch[i] != nullptr);
				assert(c_vbatch[i] != nullptr);

				IntegerType is_ok_mm = (mm >= 1);
				IntegerType is_ok_nn = (nn >= 1);
				IntegerType is_ok_kk = (kk >= 1);

				IntegerType is_ok_mnk = is_ok_mm && is_ok_nn && is_ok_kk;
				if (!is_ok_mnk) {
					fprintf(
					    stderr,
					    "dmrg_vbatch:batch_size=%d,i=%d,mm=%d,nn=%d,kk=%d\n",
					    batch_size,
					    i,
					    mm,
					    nn,
					    kk);
					fflush(stderr);
				};

				IntegerType is_ok_lda = (lda >= 1);
				IntegerType is_ok_ldb = (ldb >= 1);
				IntegerType is_ok_ldc = (ldc >= 1);

				IntegerType is_ok_ldabc = is_ok_lda && is_ok_ldb && is_ok_ldc;

				if (!is_ok_ldabc) {
					fprintf(
					    stderr,
					    "dmrg_vbatch:batch_size=%d,i=%d,lda=%d,ldb=%d,ldc=%d\n",
					    batch_size,
					    i,
					    lda,
					    ldb,
					    ldc);
					fflush(stderr);
				};

				assert(mm >= 1);
				assert(nn >= 1);
				assert(kk >= 1);

				assert(lda >= 1);
				assert(ldb >= 1);
				assert(ldc >= 1);
			};
		};

		if (ngpu == 1) {

			{
				IntegerType max_m = 0;
				IntegerType max_n = 0;
				IntegerType max_k = 0;

				for (SizeType i = 0; i < group_count; i++) {
					const IntegerType m = m_array[i];
					const IntegerType n = n_array[i];
					const IntegerType k = k_array[i];

					max_m = std::max(m, max_m);
					max_n = std::max(n, max_n);
					max_k = std::max(k, max_k);
				};
				if (idebug >= 1) {
					printf(
					    "max_m=%d, max_n=%d, max_k=%d\n", max_m, max_n, max_k);
				};

				magmablas_Xgemm_vbatched_max_nocheck(transA,
				                                     transB,
				                                     m_vbatch.data(),
				                                     n_vbatch.data(),
				                                     k_vbatch.data(),
				                                     alpha,
				                                     a_vbatch,
				                                     lda_vbatch.data(),
				                                     b_vbatch,
				                                     ldb_vbatch.data(),
				                                     beta,
				                                     c_vbatch,
				                                     ldc_vbatch.data(),
				                                     batch_size,
				                                     max_m,
				                                     max_n,
				                                     max_k,
				                                     queue);
			}
		} else {
			/*
			 * --------------------------------------------
			 * simple partitioning of work to multiple GPUs
			 * --------------------------------------------
			 */

			IntegerType* gpm_vbatch[MAXGPUS];
			IntegerType* gpn_vbatch[MAXGPUS];
			IntegerType* gpk_vbatch[MAXGPUS];

			KokkosScalar** gpa_vbatch[MAXGPUS];
			KokkosScalar** gpb_vbatch[MAXGPUS];
			KokkosScalar** gpc_vbatch[MAXGPUS];

			IntegerType* gplda_vbatch[MAXGPUS];
			IntegerType* gpldb_vbatch[MAXGPUS];
			IntegerType* gpldc_vbatch[MAXGPUS];

			IntegerType gmax_m[MAXGPUS];
			IntegerType gmax_n[MAXGPUS];
			IntegerType gmax_k[MAXGPUS];

			IntegerType gBatchCount[MAXGPUS];

			IntegerType inc  = (batch_size + (ngpu - 1)) / ngpu;
			IntegerType idev = 0;
			for (idev = 0; idev < ngpu; idev++) {
				IntegerType istart = idev * inc;
				IntegerType iend   = istart + inc - 1;
				if (iend >= (batch_size - 1)) {
					iend = batch_size - 1;
				};
				IntegerType isize = (iend - istart + 1);
				if (idebug >= 1) {
					printf("idev=%d, istart=%d, iend=%d, isize=%d\n",
					       idev,
					       istart,
					       iend,
					       isize);
				};

				device = device_array[idev];
				queue  = queue_array[idev];
				magma_setdevice(device);

				/*
				 * --------------
				 * copy arguments
				 * --------------
				 */

				size_t nbytes = (1 + isize) * sizeof(IntegerType);
				Kokkos::View<IntegerType*, Kokkos::SharedSpace> pm_vbatch(
				    "pm_vbatch", isize + 1);
				Kokkos::View<IntegerType*, Kokkos::SharedSpace> pn_vbatch(
				    "pn_vbatch", isize + 1);
				Kokkos::View<IntegerType*, Kokkos::SharedSpace> pk_vbatch(
				    "pk_vbatch", isize + 1);

				nbytes                   = (1 + isize) * sizeof(T*);
				KokkosScalar** pa_vbatch = reinterpret_cast<KokkosScalar**>(
				    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(nbytes));
				KokkosScalar** pb_vbatch = reinterpret_cast<KokkosScalar**>(
				    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(nbytes));
				KokkosScalar** pc_vbatch = reinterpret_cast<KokkosScalar**>(
				    Kokkos::kokkos_malloc<Kokkos::SharedSpace>(nbytes));

				nbytes = (1 + isize) * sizeof(IntegerType);
				Kokkos::View<IntegerType*, Kokkos::SharedSpace> plda_vbatch(
				    "plda_vbatch", isize + 1);
				Kokkos::View<IntegerType*, Kokkos::SharedSpace> pldb_vbatch(
				    "pldb_vbatch", isize + 1);
				Kokkos::View<IntegerType*, Kokkos::SharedSpace> pldc_vbatch(
				    "pldc_vbatch", isize + 1);

				{
					IntegerType i = 0;
					for (i = 0; i < isize; i++) {
						pm_vbatch[i] = m_vbatch[istart + i];
						pn_vbatch[i] = n_vbatch[istart + i];
						pk_vbatch[i] = k_vbatch[istart + i];

						plda_vbatch[i] = lda_vbatch[istart + i];
						pldb_vbatch[i] = ldb_vbatch[istart + i];
						pldc_vbatch[i] = ldc_vbatch[istart + i];
					};
				};

				double      gmem   = 0;
				double      gflops = 0;
				IntegerType max_m  = 0;
				IntegerType max_n  = 0;
				IntegerType max_k  = 0;
				{
					IntegerType i = 0;
					for (i = 0; i < isize; i++) {
						IntegerType m = pm_vbatch[i];
						IntegerType n = pn_vbatch[i];
						IntegerType k = pk_vbatch[i];

						max_m = std::max(max_m, m);
						max_n = std::max(max_n, n);
						max_k = std::max(max_k, k);

						gmem += m * n + m * k + k * n;
						gflops += 2.0 * m * n * k;
					};

					double gmem_in_bytes = gmem * sizeof(T);

					if (idebug >= 1) {
						printf("device=%d need %le bytes %le flops \n",
						       idev,
						       gmem_in_bytes,
						       gflops);
					};
				}

				/*
				 * -------------
				 * save poIntegerTypeers
				 * -------------
				 */
				gpm_vbatch[idev] = pm_vbatch.data();
				gpn_vbatch[idev] = pn_vbatch.data();
				gpk_vbatch[idev] = pk_vbatch.data();

				gpa_vbatch[idev] = pa_vbatch;
				gpb_vbatch[idev] = pb_vbatch;
				gpc_vbatch[idev] = pc_vbatch;

				gplda_vbatch[idev] = plda_vbatch.data();
				gpldb_vbatch[idev] = pldb_vbatch.data();
				gpldc_vbatch[idev] = pldc_vbatch.data();

				gmax_m[idev] = max_m;
				gmax_n[idev] = max_n;
				gmax_k[idev] = max_k;

				gBatchCount[idev] = isize;

			}; /* end for idev */

			/*
			 * -------------------------
			 * non-blocking computations
			 * -------------------------
			 */
			for (idev = 0; idev < ngpu; idev++) {
				device = device_array[idev];
				queue  = queue_array[idev];
				magma_setdevice(device);

				IntegerType* pm_vbatch = gpm_vbatch[idev];
				IntegerType* pn_vbatch = gpn_vbatch[idev];
				IntegerType* pk_vbatch = gpk_vbatch[idev];

				assert(pm_vbatch != nullptr);
				assert(pn_vbatch != nullptr);
				assert(pk_vbatch != nullptr);

				KokkosScalar** pa_vbatch = gpa_vbatch[idev];
				KokkosScalar** pb_vbatch = gpb_vbatch[idev];
				KokkosScalar** pc_vbatch = gpc_vbatch[idev];

				assert(pa_vbatch != nullptr);
				assert(pb_vbatch != nullptr);
				assert(pc_vbatch != nullptr);

				IntegerType* plda_vbatch = gplda_vbatch[idev];
				IntegerType* pldb_vbatch = gpldb_vbatch[idev];
				IntegerType* pldc_vbatch = gpldc_vbatch[idev];

				assert(plda_vbatch != nullptr);
				assert(pldb_vbatch != nullptr);
				assert(pldc_vbatch != nullptr);

				IntegerType max_m = gmax_m[idev];
				IntegerType max_n = gmax_n[idev];
				IntegerType max_k = gmax_k[idev];

				assert(max_m >= 1);
				assert(max_n >= 1);
				assert(max_k >= 1);

				IntegerType pbatch_size = gBatchCount[idev];

				if (pbatch_size >= 1) {
					magmablas_Xgemm_vbatched_max_nocheck(transA,
					                                     transB,
					                                     pm_vbatch,
					                                     pn_vbatch,
					                                     pk_vbatch,
					                                     alpha,
					                                     pa_vbatch,
					                                     plda_vbatch,
					                                     pb_vbatch,
					                                     pldb_vbatch,
					                                     beta,
					                                     pc_vbatch,
					                                     pldc_vbatch,
					                                     pbatch_size,
					                                     max_m,
					                                     max_n,
					                                     max_k,
					                                     queue);
				};

			}; /* end for idev */

			/*
			 * -------------------------------
			 * wait for computations to finish
			 * -------------------------------
			 */
			for (idev = 0; idev < ngpu; idev++) {
				device = device_array[idev];
				queue  = queue_array[idev];
				magma_setdevice(device);
				magma_queue_sync(queue);
			};

			for (idev = 0; idev < ngpu; idev++) {
				KokkosScalar** pa_vbatch = gpa_vbatch[idev];
				KokkosScalar** pb_vbatch = gpb_vbatch[idev];
				KokkosScalar** pc_vbatch = gpc_vbatch[idev];

				Kokkos::kokkos_free<Kokkos::SharedSpace>(pa_vbatch);
				Kokkos::kokkos_free<Kokkos::SharedSpace>(pb_vbatch);
				Kokkos::kokkos_free<Kokkos::SharedSpace>(pc_vbatch);
			} /* end for idev */

			idev   = 0;
			device = device_array[idev];
			queue  = queue_array[idev];

		}; /* end if (ngpu > 1) */
	};
#else
	{
		IntegerType i = 0;
		// #pragma omp parallel for private(i) schedule(dynamic)
		for (i = 0; i < batch_size; i++) {
			const IntegerType mm  = m_vbatch[i];
			const IntegerType nn  = n_vbatch[i];
			const IntegerType kk  = k_vbatch[i];
			const IntegerType lda = lda_vbatch[i];
			const IntegerType ldb = ldb_vbatch[i];
			const IntegerType ldc = ldc_vbatch[i];
			const T*          pa  = reinterpret_cast<const T*>(a_vbatch[i]);
			const T*          pb  = reinterpret_cast<const T*>(b_vbatch[i]);
			T*                pc  = reinterpret_cast<T*>(c_vbatch[i]);

			Xgemm_(&ctransa,
			       &ctransb,
			       &mm,
			       &nn,
			       &kk,
			       &alpha,
			       pa,
			       &lda,
			       pb,
			       &ldb,
			       &beta,
			       pc,
			       &ldc);
		};
	};
#endif

	Kokkos::kokkos_free<Kokkos::SharedSpace>(a_vbatch);
	Kokkos::kokkos_free<Kokkos::SharedSpace>(b_vbatch);
	Kokkos::kokkos_free<Kokkos::SharedSpace>(c_vbatch);

	if (idebug >= 1) {
		elapsed_time += dmrg_get_wtime();
		double gflops_per_sec = 0;
		if (elapsed_time > 0) {
			gflops_per_sec = gflops / elapsed_time;
		};

		printf("dmrg_vbatch: gflops=%lf, elapsed_time=%lf, gflops/sec=%lf\n",
		       gflops,
		       elapsed_time,
		       gflops_per_sec);
		printf("dmrg_vbatch need %lf GBytes\n", (double)nbytes_total / (giga));
	};
}

template void dmrg_Xgemm_vbatch<double>(char,
                                        char,
                                        IntegerType*,
                                        IntegerType*,
                                        IntegerType*,
                                        double,
                                        const double**,
                                        IntegerType*,
                                        const double**,
                                        IntegerType*,
                                        double,
                                        double**,
                                        IntegerType*,
                                        SizeType,
                                        IntegerType*);

template void dmrg_Xgemm_vbatch<std::complex<double>>(char                         ctransa_array,
                                                      char                         ctransb_array,
                                                      IntegerType*                 m_array,
                                                      IntegerType*                 n_array,
                                                      IntegerType*                 k_array,
                                                      std::complex<double>         alpha,
                                                      const std::complex<double>** a_array,
                                                      IntegerType*                 lda_array,
                                                      const std::complex<double>** b_array,
                                                      IntegerType*                 ldb_array,
                                                      std::complex<double>         beta,
                                                      std::complex<double>**       c_array,
                                                      IntegerType*                 ldc_array,
                                                      SizeType                     group_count,
                                                      IntegerType*                 group_size);

template void
dmrg_Xgetvector<double>(const IntegerType, double*, const IntegerType, double*, const IntegerType);

template void dmrg_Xgetvector<std::complex<double>>(const IntegerType,
                                                    std::complex<double>*,
                                                    const IntegerType,
                                                    std::complex<double>*,
                                                    const IntegerType);

template void dmrg_Xsetvector<double>(const IntegerType,
                                      const double*,
                                      const IntegerType,
                                      double*,
                                      const IntegerType);

template void dmrg_Xsetvector<std::complex<double>>(const IntegerType,
                                                    const std::complex<double>*,
                                                    const IntegerType,
                                                    std::complex<double>*,
                                                    const IntegerType);

template void dmrg_Xgetmatrix<double>(const IntegerType,
                                      const IntegerType,
                                      double*,
                                      const IntegerType,
                                      double*,
                                      const IntegerType);

template void dmrg_Xgetmatrix<std::complex<double>>(const IntegerType,
                                                    const IntegerType,
                                                    std::complex<double>*,
                                                    const IntegerType,
                                                    std::complex<double>*,
                                                    const IntegerType);

template void dmrg_Xsetmatrix<double>(const IntegerType,
                                      const IntegerType,
                                      const double*,
                                      const IntegerType,
                                      double*,
                                      const IntegerType);

template void dmrg_Xsetmatrix<std::complex<double>>(const IntegerType,
                                                    const IntegerType,
                                                    const std::complex<double>*,
                                                    const IntegerType,
                                                    std::complex<double>*,
                                                    const IntegerType);
