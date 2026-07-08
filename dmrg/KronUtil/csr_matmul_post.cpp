#include "util.h"

#include "PsimagLiteConfig.h"
#ifdef PSIMAGLITE_USE_KOKKOS
#include <KokkosBatched_Axpy.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBlas2_serial_gemv.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_spmv.hpp>
#include <Kokkos_Complex.hpp>
#include <Kokkos_Core.hpp>
#include <type_traits>
#include <vector>

// helper to pick Kokkos scalar type only when appropriate (file-scope)
template <typename T, bool IsComplex> struct KokkosScalarTypePost {
	using type = T;
};
template <typename T> struct KokkosScalarTypePost<T, true> {
	using type = Kokkos::complex<typename T::value_type>;
};
#endif

template <typename ComplexOrRealType>
void csr_matmul_post(char                                                       trans_A,
                     const PsimagLite::CrsMatrix<ComplexOrRealType>&            a,
                     const int                                                  nrow_Y,
                     const int                                                  ncol_Y,
                     const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                     const int                                                  nrow_X,
                     const int                                                  ncol_X,
                     PsimagLite::MatrixNonOwned<ComplexOrRealType>&             xout)
{
#ifdef PSIMAGLITE_USE_KOKKOS
	Kokkos::Profiling::ScopedRegion region("PsimagLite::csr_matmul_post");
#endif
	/*
	 * -------------------------------------------------------
	 * A in compressed sparse ROW format
	 *
	 * compute   X +=  Y * op(A)
	 * where op(A) is transpose(A)   if trans_A = 'T' or 't'
	 *       op(A) is A              otherwise
	 *
	 * if need transpose(A) then
	 *   X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * tranpose(A(nrow_A,ncol_A))
	 *   requires (nrow_X == nrow_Y) && (ncol_Y == ncol_A) && (ncol_X == nrow_A)
	 *
	 * if need A then
	 *  X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * A(nrow_A,ncol_A)
	 *  requires  (nrow_X == nrow_Y) && ( ncol_Y == nrow_A) && (ncol_X == ncol_A)
	 * -------------------------------------------------------
	 */
	const bool is_complex      = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int  nrow_A          = a.rows();
	int        isTranspose     = (trans_A == 'T') || (trans_A == 't');
	int        isConjTranspose = (trans_A == 'C') || (trans_A == 'c');
	int        isConj          = (trans_A == 'Z') || (trans_A == 'z');

#ifdef PSIMAGLITE_USE_KOKKOS
	// Use KokkosSparse::spmv by transposing operations
	using HostExec = Kokkos::DefaultExecutionSpace;
	HostExec exec;
	using Ordinal = int;

	using KokkosScalar = typename KokkosScalarTypePost<
	    ComplexOrRealType,
	    PsimagLite::IsComplexNumber<ComplexOrRealType>::True>::type;

	const int nnz = a.nonZeros();

	// build host arrays
	std::vector<int> rowptr(nrow_A + 1);
	for (int i = 0; i <= nrow_A; ++i)
		rowptr[i] = a.getRowPtr(i);

	std::vector<int>          cols(nnz);
	std::vector<KokkosScalar> vals(nnz);
	for (int k = 0; k < nnz; ++k) {
		cols[k] = a.getCol(k);
		if constexpr (!PsimagLite::IsComplexNumber<ComplexOrRealType>::True) {
			vals[k] = a.getValue(k);
		} else {
			ComplexOrRealType v = a.getValue(k);
			if (is_complex && (isConj || isConjTranspose))
				v = PsimagLite::conj(v);
			vals[k] = Kokkos::complex<typename ComplexOrRealType::value_type>(v.real(),
			                                                                  v.imag());
		}
	}

	// build CrsMatrix from raw host arrays; constructor will deep-copy to device
	{
		double maxAbs = 0.0;
		for (int i = 0; i < nnz; ++i) {
			double v = 0.0;
			if constexpr (!PsimagLite::IsComplexNumber<ComplexOrRealType>::True)
				v = std::abs((double)vals[i]);
			else
				v = std::abs((double)vals[i].real())
				    + std::abs((double)vals[i].imag());
			if (v > maxAbs)
				maxAbs = v;
		}
		(void)maxAbs;
	}
	KokkosSparse::CrsMatrix<KokkosScalar, Ordinal, HostExec> A_crs(
	    "A_crs", nrow_A, (int)a.cols(), nnz, vals.data(), rowptr.data(), cols.data());

	// For each column of Y perform spmv
	for (int iy = 0; iy < nrow_Y; ++iy) {
		std::vector<KokkosScalar> yhost((size_t)ncol_Y);
		for (int j = 0; j < ncol_Y; ++j) {
			if constexpr (!PsimagLite::IsComplexNumber<ComplexOrRealType>::True) {
				yhost[j] = yin(iy, j);
			} else {
				auto vv  = yin(iy, j);
				yhost[j] = Kokkos::complex<typename ComplexOrRealType::value_type>(
				    vv.real(), vv.imag());
			}
		}

		auto y_dev = Kokkos::create_mirror_view_and_copy(
		    HostExec(),
		    Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
		        yhost.data(), ncol_Y));
		auto x_dev_out = Kokkos::View<KokkosScalar*>("x_dev_out", ncol_X);

		if (isTranspose || isConjTranspose) {
			KokkosSparse::spmv(
			    "N", (KokkosScalar)1.0, A_crs, y_dev, (KokkosScalar)0.0, x_dev_out);
		} else {
			KokkosSparse::spmv(
			    "T", (KokkosScalar)1.0, A_crs, y_dev, (KokkosScalar)0.0, x_dev_out);
		}

		std::vector<KokkosScalar>                      xhost((size_t)ncol_X);
		Kokkos::View<KokkosScalar*, Kokkos::HostSpace> h_xhost(xhost.data(), ncol_X);
		Kokkos::deep_copy(h_xhost, x_dev_out);
		exec.fence();

		for (int jx = 0; jx < ncol_X; ++jx) {
			if constexpr (!PsimagLite::IsComplexNumber<ComplexOrRealType>::True) {
				xout(iy, jx) += xhost[jx];
			} else {
				Kokkos::complex<typename ComplexOrRealType::value_type> c
				    = xhost[jx];
				xout(iy, jx) += ComplexOrRealType(
				    static_cast<typename ComplexOrRealType::value_type>(c.real()),
				    static_cast<typename ComplexOrRealType::value_type>(c.imag()));
			}
		}
	}

	return;
#endif
	if (isTranspose || isConjTranspose) {
		/*
		 *   ----------------------------------------------------------
		 *   X(nrow_X,ncol_X) +=  Y(nrow_Y,ncol_Y) * transpose(A(nrow_A,ncol_A))
		 *   X(ix,jx) +=  Y(iy,jy) * transpose( A(ia,ja) )
		 *   X(ix,jx) += sum( Y(iy, ja) * At(ja,ia), over ja )
		 *   ----------------------------------------------------------
		 */

		assert(nrow_X == nrow_Y);
		assert(static_cast<SizeType>(ncol_Y) == a.cols() && (ncol_X == nrow_A));

		int ia = 0;
		for (ia = 0; ia < nrow_A; ia++) {
			int istart = a.getRowPtr(ia);
			int iend   = a.getRowPtr(ia + 1);
			int k      = 0;
			for (k = istart; k < iend; k++) {
				int               ja   = a.getCol(k);
				ComplexOrRealType aij  = a.getValue(k);
				ComplexOrRealType atji = aij;
				if (is_complex && isConjTranspose) {
					atji = PsimagLite::conj(atji);
				};

				int iy = 0;
				for (iy = 0; iy < nrow_Y; iy++) {
					int ix = iy;
					int jx = ia;

					xout(ix, jx) += (yin(iy, ja) * atji);
				}
			}
		}
	} else {
		/*
		 * ---------------------------------------------
		 * X(nrow_X,ncol_X) += Y(nrow_Y,ncol_Y) * A(nrow_A,ncol_A)
		 * X(ix,jx) += sum( Y(iy,ia) * A(ia,ja), over ia )
		 * ---------------------------------------------
		 */

		assert(nrow_X == nrow_Y);
		assert(ncol_Y == nrow_A && static_cast<SizeType>(ncol_X) == a.cols());

		int ia = 0;
		for (ia = 0; ia < nrow_A; ia++) {
			int istart = a.getRowPtr(ia);
			int iend   = a.getRowPtr(ia + 1);
			int k      = 0;
			for (k = istart; k < iend; k++) {
				int               ja  = a.getCol(k);
				ComplexOrRealType aij = a.getValue(k);
				if (is_complex && isConj) {
					aij = PsimagLite::conj(aij);
				};

				int iy = 0;
				for (iy = 0; iy < nrow_Y; iy++) {
					int ix = iy;
					int jx = ja;

					xout(ix, jx) += (yin(iy, ia) * aij);
				}
			}
		}
	}
}
