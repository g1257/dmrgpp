#include "dmrg_lapack.h"
#include "dmrg_vbatch.h"

template <typename T>
void dmrg_lacpy(const char*       uplo,
                const IntegerType m,
                const IntegerType n,
                const T*          src,
                const IntegerType ld_src,
                T*                dest,
                const IntegerType ld_dest)
{
#ifdef USE_MAGMA
	const IntegerType is_upper = (*uplo == 'U') || (*uplo == 'u');
	const IntegerType is_lower = (*uplo == 'L') || (*uplo == 'l');
	const IntegerType is_full  = (!is_upper) && (!is_lower);

	IntegerType is_block_copy = is_full && (m == ld_src) && (m == ld_dest);
	if (is_block_copy) {
		const size_t nbytes = sizeof(T) * m * n;
		dmrg_memcpy(dest, src, nbytes);
	} else {
		const IntegerType min_mn = (m <= n) ? m : n;
		const IntegerType ncol   = is_full ? n : min_mn;
		IntegerType       jcol   = 0;

		for (jcol = 0; jcol < ncol; jcol++) {
			const IntegerType irow   = jcol;
			IntegerType       istart = is_upper ? 0 : is_lower ? irow : 0;
			IntegerType       iend   = is_upper ? irow : is_lower ? m - 1 : m - 1;
			IntegerType       count  = iend - istart + 1;
			if (count >= 1) {

				const T* psrc  = src + jcol * ld_src + istart;
				T*       pdest = dest + jcol * ld_dest + istart;

				const size_t nbytes = count * sizeof(T);
				dmrg_memcpy(pdest, psrc, nbytes);
			}
		}
	}
#else
	Xlacpy_(uplo, &m, &n, src, &ld_src, dest, &ld_dest);
#endif
}

template void dmrg_lacpy<double>(const char*       uplo,
                                 const IntegerType m,
                                 const IntegerType n,
                                 const double*     src,
                                 const IntegerType ld_src,
                                 double*           dest,
                                 const IntegerType ld_dest);

template void dmrg_lacpy<std::complex<double>>(const char*                 uplo,
                                               const IntegerType           m,
                                               const IntegerType           n,
                                               const std::complex<double>* src,
                                               const IntegerType           ld_src,
                                               std::complex<double>*       dest,
                                               const IntegerType           ld_dest);
