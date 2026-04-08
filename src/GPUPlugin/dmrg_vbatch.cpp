#include "dmrg_vbatch.h"
#include "dmrg_lapack.h"
#include "dmrg_types.h"

#define USE_MALLOC

#ifdef USE_INTEL_MKL
#include "mkl.h"
#endif

static IntegerType is_initialized = 0;

#ifdef USE_MAGMA
#undef USE_MALLOC
#define USE_MALLOC

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

SizeType MAX(const SizeType& x, const SizeType& y) { return (((x) > (y)) ? (x) : (y)); }

SizeType MIN(const SizeType& x, const SizeType& y) { return (((x) < (y)) ? (x) : (y)); }

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
	magma_Xgetvector(n, (MAGMA_T*)dx_src, incx, (MAGMA_T*)hy_dst, incy, queue);
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
	magma_Xsetvector(n, (MAGMA_T*)hx_src, incx, (MAGMA_T*)dy_dst, incy, queue);
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
	magma_Xgetmatrix(m, n, (MAGMA_T*)dA_src, ldda, (MAGMA_T*)hB_dst, ldb, queue);
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
	magma_Xsetmatrix(m, n, (MAGMA_T*)hA_src, lda, (MAGMA_T*)dB_dst, lddb, queue);
#else
	const char* uplo = "A";
	Xlacpy_(uplo, &m, &n, hA_src, &lda, dB_dst, &lddb);
#endif
}

template <typename T>
void dmrg_Xgemm_vbatch(char*        ctransa_array,
                       char*        ctransb_array,
                       IntegerType* m_array,
                       IntegerType* n_array,
                       IntegerType* k_array,
                       T*           alpha_array,
                       T**          a_array,
                       IntegerType* lda_array,
                       T**          b_array,
                       IntegerType* ldb_array,
                       T*           beta_array,
                       T**          c_array,
                       IntegerType* ldc_array,
                       SizeType     group_count,
                       IntegerType* group_size)
{

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

#ifdef USE_INTEL_MKL
	{

		MKL_INT m_array_mkl[group_count];
		MKL_INT n_array_mkl[group_count];
		MKL_INT k_array_mkl[group_count];

		MKL_INT lda_array_mkl[group_count];
		MKL_INT ldb_array_mkl[group_count];
		MKL_INT ldc_array_mkl[group_count];

		MKL_INT group_size_mkl[group_count];

		CBLAS_TRANSPOSE transa_array[group_count];
		CBLAS_TRANSPOSE transb_array[group_count];

		MKL_INT group_count_mkl = group_count;

		IntegerType igroup = 0;
		for (igroup = 0; igroup < group_count; igroup++) {
			char        transa    = ctransa_array[igroup];
			char        transb    = ctransb_array[igroup];
			IntegerType is_transa = (transa == 'T') || (transa == 't');
			IntegerType is_transb = (transb == 'T') || (transb == 't');

			transa_array[igroup] = is_transa ? CblasTrans : CblasNoTrans;
			transb_array[igroup] = is_transb ? CblasTrans : CblasNoTrans;

			m_array_mkl[igroup] = m_array[igroup];
			n_array_mkl[igroup] = n_array[igroup];
			k_array_mkl[igroup] = k_array[igroup];

			lda_array_mkl[igroup] = lda_array[igroup];
			ldb_array_mkl[igroup] = ldb_array[igroup];
			ldc_array_mkl[igroup] = ldc_array[igroup];

			group_size_mkl[igroup] = group_size[igroup];
		};

		cblas_Xgemm_batch(CblasColMajor,
		                  transa_array,
		                  transb_array,
		                  m_array_mkl,
		                  n_array_mkl,
		                  k_array_mkl,
		                  alpha_array,
		                  (const T**)&(a_array[0]),
		                  lda_array_mkl,
		                  (const T**)&(b_array[0]),
		                  ldb_array_mkl,
		                  beta_array,
		                  c_array,
		                  ldc_array_mkl,
		                  group_count_mkl,
		                  group_size_mkl);
	};
#else
	/*
	---------------------
	expand out the groups
	---------------------
	*/
	IntegerType batch_size = 0;
	SizeType    igroup     = 0;
	for (igroup = 0; igroup < group_count; igroup++) {
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
#ifdef USE_MALLOC

	nbytes = sizeof(T) * (vbatch_dim);
	nbytes_total += nbytes;
	T* alpha_vbatch = new T[vbatch_dim];

	nbytes = sizeof(T) * (vbatch_dim);
	nbytes_total += nbytes;
	T* beta_vbatch = new T[vbatch_dim];

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	IntegerType* m_vbatch = new IntegerType[vbatch_dim];

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	IntegerType* n_vbatch = new IntegerType[vbatch_dim];

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	IntegerType* k_vbatch = new IntegerType[vbatch_dim];

	nbytes = sizeof(char) * (vbatch_dim);
	nbytes_total += nbytes;
	char* transa_vbatch = new char[vbatch_dim];

	nbytes = sizeof(char) * (vbatch_dim);
	nbytes_total += nbytes;
	char* transb_vbatch = new char[vbatch_dim];

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	IntegerType* lda_vbatch = new IntegerType[vbatch_dim];

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	IntegerType* ldb_vbatch = new IntegerType[vbatch_dim];

	nbytes = sizeof(IntegerType) * (vbatch_dim);
	nbytes_total += nbytes;
	IntegerType* ldc_vbatch = new IntegerType[vbatch_dim];

	nbytes = sizeof(T*) * (vbatch_dim);
	nbytes_total += nbytes;
	T** a_vbatch = new T*[vbatch_dim];

	nbytes = sizeof(T*) * (vbatch_dim);
	nbytes_total += nbytes;
	T** b_vbatch = new T*[vbatch_dim];

	nbytes = sizeof(T*) * (vbatch_dim);
	nbytes_total += nbytes;
	T** c_vbatch = new T*[vbatch_dim];

#else

	std::vector<T> alpha_vbatch(vbatch_dim, 0.0);
	std::vector<T> beta_vbatch(vbatch_dim, 0.0);
	VectorSizeType m_vbatch(vbatch_dim, 0);
	VectorSizeType n_vbatch(vbatch_dim, 0);
	VectorSizeType k_vbatch(vbatch_dim, 0);
	VectorCharType transa_vbatch(vbatch_dim, ' ');
	VectorCharType transb_vbatch(vbatch_dim, ' ');
	VectorSizeType lda_vbatch(vbatch_dim, 0);
	VectorSizeType ldb_vbatch(vbatch_dim, 0);
	VectorSizeType ldc_vbatch(vbatch_dim, 0);

	std::vector<T*> a_vbatch = a_array;
	std::vector<T*> b_vbatch = b_array;
	std::vector<T*> c_vbatch = c_array;

#endif

	IntegerType idx = 0;
	for (igroup = 0; igroup < group_count; igroup++) {
		IntegerType i = 0;
		for (i = 0; i < group_size[igroup]; i++) {
			transa_vbatch[idx] = ctransa_array[igroup];
			transb_vbatch[idx] = ctransb_array[igroup];

			m_vbatch[idx] = m_array[igroup];
			n_vbatch[idx] = n_array[igroup];
			k_vbatch[idx] = k_array[igroup];

			alpha_vbatch[idx] = alpha_array[igroup];
			beta_vbatch[idx]  = beta_array[igroup];

			lda_vbatch[idx] = lda_array[igroup];
			ldb_vbatch[idx] = ldb_array[igroup];
			ldc_vbatch[idx] = ldc_array[igroup];

			a_vbatch[idx] = a_array[idx];
			b_vbatch[idx] = b_array[idx];
			c_vbatch[idx] = c_array[idx];

			idx = idx + 1;
		};
	};

#ifdef USE_MAGMA
	{
		/*
	---------------------------------------------------------------
	Note magma_Xgemm_vbatched assumes all transa, transb, alpha, beta are the same
	---------------------------------------------------------------
	*/
		IntegerType is_same_transa = 1;
		IntegerType is_same_transb = 1;
		IntegerType is_same_alpha  = 1;
		IntegerType is_same_beta   = 1;

		for (idx = 0; idx < batch_size; idx++) {
			is_same_transa = (transa_vbatch[0] == transa_vbatch[idx]);
			if (!is_same_transa) {
				break;
			};
		};
		for (idx = 0; idx < batch_size; idx++) {
			is_same_transb = (transb_vbatch[0] == transb_vbatch[idx]);
			if (!is_same_transb) {
				break;
			};
		};
		for (idx = 0; idx < batch_size; idx++) {
			is_same_alpha = (alpha_vbatch[0] == alpha_vbatch[idx]);
			if (!is_same_alpha) {
				break;
			};
		};
		for (idx = 0; idx < batch_size; idx++) {
			is_same_beta = (beta_vbatch[0] == beta_vbatch[idx]);
			if (!is_same_beta) {
				break;
			};
		};

		IntegerType is_same_all
		    = is_same_transa && is_same_transb && is_same_alpha && is_same_beta;

		if (!is_same_all) {
			if (!is_same_transa) {
				fprintf(
				    stderr, "dmrg_vbatch: is_same_transa = %d\n", is_same_transa);
			};
			if (!is_same_transb) {
				fprintf(
				    stderr, "dmrg_vbatch: is_same_transb = %d\n", is_same_transb);
			};
			if (!is_same_alpha) {
				fprintf(stderr, "dmrg_vbatch: is_same_alpha = %d\n", is_same_alpha);
			};
			if (!is_same_beta) {
				fprintf(stderr, "dmrg_vbatch: is_same_beta = %d\n", is_same_beta);
			};
		};
		assert(is_same_all);

		IntegerType isTransA = ((transa_vbatch[0] == 'T') || (transa_vbatch[0] == 't'));
		IntegerType isTransB = ((transb_vbatch[0] == 'T') || (transb_vbatch[0] == 't'));

		magma_trans_t transA = isTransA ? MagmaTrans : MagmaNoTrans;
		magma_trans_t transB = isTransB ? MagmaTrans : MagmaNoTrans;

		MAGMA_T alpha = *((MAGMA_T*)&(alpha_vbatch[0]));
		MAGMA_T beta  = *((MAGMA_T*)&(beta_vbatch[0]));

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

				T* Amat = a_vbatch[i];
				T* Bmat = b_vbatch[i];
				T* Cmat = c_vbatch[i];

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

				assert(Amat != 0);
				assert(Bmat != 0);
				assert(Cmat != 0);
			};
		};

		if (ngpu == 1) {

			{
				IntegerType max_m = 0;
				IntegerType max_n = 0;
				IntegerType max_k = 0;

				IntegerType i = 0;
				for (i = 0; i < group_count; i++) {
					const IntegerType m = m_array[i];
					const IntegerType n = n_array[i];
					const IntegerType k = k_array[i];

					max_m = MAX(m, max_m);
					max_n = MAX(n, max_n);
					max_k = MAX(k, max_k);
				};
				if (idebug >= 1) {
					printf(
					    "max_m=%d, max_n=%d, max_k=%d\n", max_m, max_n, max_k);
				};

				magmablas_Xgemm_vbatched_max_nocheck(
				    transA,
				    transB,
				    m_vbatch,
				    n_vbatch,
				    k_vbatch,
				    (MAGMA_T)alpha,
				    (MAGMA_T const* const*)a_vbatch,
				    lda_vbatch,
				    (MAGMA_T const* const*)b_vbatch,
				    ldb_vbatch,
				    (MAGMA_T)beta,
				    (MAGMA_T**)c_vbatch,
				    ldc_vbatch,
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

			T** gpa_vbatch[MAXGPUS];
			T** gpb_vbatch[MAXGPUS];
			T** gpc_vbatch[MAXGPUS];

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

				IntegerType* pm_vbatch = new IntegerType[1 + isize];
				IntegerType* pn_vbatch = new IntegerType[1 + isize];
				IntegerType* pk_vbatch = new IntegerType[1 + isize];

				assert(pm_vbatch != NULL);
				assert(pn_vbatch != NULL);
				assert(pk_vbatch != NULL);

				T** pa_vbatch = new T*[1 + isize];
				T** pb_vbatch = new T*[1 + isize];
				T** pc_vbatch = new T*[1 + isize];

				assert(pa_vbatch != NULL);
				assert(pb_vbatch != NULL);
				assert(pc_vbatch != NULL);

				IntegerType* plda_vbatch = new IntegerType[1 + isize];
				IntegerType* pldb_vbatch = new IntegerType[1 + isize];
				IntegerType* pldc_vbatch = new IntegerType[1 + isize];

				assert(plda_vbatch != 0);
				assert(pldb_vbatch != 0);
				assert(pldc_vbatch != 0);

				{
					IntegerType i = 0;
					for (i = 0; i < isize; i++) {
						pm_vbatch[i] = m_vbatch[istart + i];
						pn_vbatch[i] = n_vbatch[istart + i];
						pk_vbatch[i] = k_vbatch[istart + i];

						pa_vbatch[i] = a_vbatch[istart + i];
						pb_vbatch[i] = b_vbatch[istart + i];
						pc_vbatch[i] = c_vbatch[istart + i];

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

						max_m = MAX(max_m, m);
						max_n = MAX(max_n, n);
						max_k = MAX(max_k, k);

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
				gpm_vbatch[idev] = pm_vbatch;
				gpn_vbatch[idev] = pn_vbatch;
				gpk_vbatch[idev] = pk_vbatch;

				gpa_vbatch[idev] = pa_vbatch;
				gpb_vbatch[idev] = pb_vbatch;
				gpc_vbatch[idev] = pc_vbatch;

				gplda_vbatch[idev] = plda_vbatch;
				gpldb_vbatch[idev] = pldb_vbatch;
				gpldc_vbatch[idev] = pldc_vbatch;

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

				assert(pm_vbatch != NULL);
				assert(pn_vbatch != NULL);
				assert(pk_vbatch != NULL);

				T** pa_vbatch = gpa_vbatch[idev];
				T** pb_vbatch = gpb_vbatch[idev];
				T** pc_vbatch = gpc_vbatch[idev];

				assert(pa_vbatch != NULL);
				assert(pb_vbatch != NULL);
				assert(pc_vbatch != NULL);

				IntegerType* plda_vbatch = gplda_vbatch[idev];
				IntegerType* pldb_vbatch = gpldb_vbatch[idev];
				IntegerType* pldc_vbatch = gpldc_vbatch[idev];

				assert(plda_vbatch != NULL);
				assert(pldb_vbatch != NULL);
				assert(pldc_vbatch != NULL);

				IntegerType max_m = gmax_m[idev];
				IntegerType max_n = gmax_n[idev];
				IntegerType max_k = gmax_k[idev];

				assert(max_m >= 1);
				assert(max_n >= 1);
				assert(max_k >= 1);

				IntegerType pbatch_size = gBatchCount[idev];

				if (pbatch_size >= 1) {
					magmablas_Xgemm_vbatched_max_nocheck(
					    transA,
					    transB,
					    pm_vbatch,
					    pn_vbatch,
					    pk_vbatch,
					    (MAGMA_T)alpha,
					    (MAGMA_T const* const*)pa_vbatch,
					    plda_vbatch,
					    (MAGMA_T const* const*)pb_vbatch,
					    pldb_vbatch,
					    (MAGMA_T)beta,
					    (MAGMA_T**)pc_vbatch,
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
				IntegerType* pm_vbatch = gpm_vbatch[idev];
				IntegerType* pn_vbatch = gpn_vbatch[idev];
				IntegerType* pk_vbatch = gpk_vbatch[idev];

				T** pa_vbatch = gpa_vbatch[idev];
				T** pb_vbatch = gpb_vbatch[idev];
				T** pc_vbatch = gpc_vbatch[idev];

				IntegerType* plda_vbatch = gplda_vbatch[idev];
				IntegerType* pldb_vbatch = gpldb_vbatch[idev];
				IntegerType* pldc_vbatch = gpldc_vbatch[idev];

				delete[] pm_vbatch;
				delete[] pn_vbatch;
				delete[] pk_vbatch;

				delete[] pa_vbatch;
				delete[] pb_vbatch;
				delete[] pc_vbatch;

				delete[] plda_vbatch;
				delete[] pldb_vbatch;
				delete[] pldc_vbatch;

			}; /* end for idev */

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
			const IntegerType mm     = m_vbatch[i];
			const IntegerType nn     = n_vbatch[i];
			const IntegerType kk     = k_vbatch[i];
			const char        transa = transa_vbatch[i];
			const char        transb = transb_vbatch[i];
			const T           alpha  = alpha_vbatch[i];
			const T           beta   = beta_vbatch[i];
			const IntegerType lda    = lda_vbatch[i];
			const IntegerType ldb    = ldb_vbatch[i];
			const IntegerType ldc    = ldc_vbatch[i];
			const T*          pa     = a_vbatch[i];
			const T*          pb     = b_vbatch[i];
			T*                pc     = c_vbatch[i];

			Xgemm_(&transa,
			       &transb,
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

#ifdef USE_MALLOC

	delete[] alpha_vbatch;
	delete[] beta_vbatch;

	delete[] m_vbatch;
	delete[] n_vbatch;
	delete[] k_vbatch;

	delete[] transa_vbatch;
	delete[] transb_vbatch;

	delete[] lda_vbatch;
	delete[] ldb_vbatch;
	delete[] ldc_vbatch;

	delete[] a_vbatch;
	delete[] b_vbatch;
	delete[] c_vbatch;

#endif

#endif

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

template void dmrg_Xgemm_vbatch<MYTYPE>(char*,
                                        char*,
                                        IntegerType*,
                                        IntegerType*,
                                        IntegerType*,
                                        MYTYPE*,
                                        MYTYPE**,
                                        IntegerType*,
                                        MYTYPE**,
                                        IntegerType*,
                                        MYTYPE*,
                                        MYTYPE**,
                                        IntegerType*,
                                        SizeType,
                                        IntegerType*);

#if defined(USE_COMPLEX_Z)
template void dmrg_Xgemm_vbatch<double>(char*,
                                        char*,
                                        IntegerType*,
                                        IntegerType*,
                                        IntegerType*,
                                        double*,
                                        double**,
                                        IntegerType*,
                                        double**,
                                        IntegerType*,
                                        double*,
                                        double**,
                                        IntegerType*,
                                        SizeType,
                                        IntegerType*);

#else
template void dmrg_Xgemm_vbatch<std::complex<double>>(char*                  ctransa_array,
                                                      char*                  ctransb_array,
                                                      IntegerType*           m_array,
                                                      IntegerType*           n_array,
                                                      IntegerType*           k_array,
                                                      std::complex<double>*  alpha_array,
                                                      std::complex<double>** a_array,
                                                      IntegerType*           lda_array,
                                                      std::complex<double>** b_array,
                                                      IntegerType*           ldb_array,
                                                      std::complex<double>*  beta_array,
                                                      std::complex<double>** c_array,
                                                      IntegerType*           ldc_array,
                                                      SizeType               group_count,
                                                      IntegerType*           group_size);
#endif

template void
dmrg_Xgetvector<MYTYPE>(const IntegerType, MYTYPE*, const IntegerType, MYTYPE*, const IntegerType);
template void dmrg_Xsetvector<MYTYPE>(const IntegerType,
                                      const MYTYPE*,
                                      const IntegerType,
                                      MYTYPE*,
                                      const IntegerType);
template void dmrg_Xgetmatrix<MYTYPE>(const IntegerType,
                                      const IntegerType,
                                      MYTYPE*,
                                      const IntegerType,
                                      MYTYPE*,
                                      const IntegerType);
template void dmrg_Xsetmatrix<MYTYPE>(const IntegerType,
                                      const IntegerType,
                                      const MYTYPE*,
                                      const IntegerType,
                                      MYTYPE*,
                                      const IntegerType);
