#ifndef PSI_GEMMR_H
#define PSI_GEMMR_H
#include "Vector.h"
#include <cstdlib>
#include "BLAS.h"
#include "Parallelizer2.h"

// Originally written by Ed

namespace PsimagLite {

template<typename T>
class GemmR {

public:

	GemmR(bool idebug, SizeType nb, SizeType nthreads)
	    :  idebug_(idebug), nb_(nb), nthreads_(nthreads)
	{
		if (nb_ == 0)
			throw RuntimeError("GemmR: nb cannot be zero\n");
	}

	void operator()(char const transA,
	                char const transB,
	                int const m,
	                int const n,
	                int const k,
	                T const & alpha,
	                T const * const A,
	                int const ldA,
	                T const * const B,
	                int const ldB,
	                T const & beta,
	                T       * const C,
	                int const ldC)
	{
		if (isSmall(m, n)) {
			psimag::BLAS::GEMM( transA,
			                    transB,
			                    m,
			                    n,
			                    k,
			                    alpha,
			                    A,
			                    ldA,
			                    B,
			                    ldB,
			                    beta,
			                    C,
			                    ldC );
			return;
		}

		bigMatrices(transA,
		            transB,
		            m,
		            n,
		            k,
		            alpha,
		            A,
		            ldA,
		            B,
		            ldB,
		            beta,
		            C,
		            ldC);
	}

private:

	void bigMatrices(char const transA,
	                 char const transB,
	                 int const m,
	                 int const n,
	                 int const k,
	                 T const & alpha,
	                 T const * const A,
	                 int const ldA,
	                 T const * const B,
	                 int const ldB,
	                 T const & beta,
	                 T       * const C,
	                 int const ldC)
	{
		assert(!isSmall(m, n));

		//  -------------------------
		//  if n = 101, and nb = 100,
		//  wish to avoid  one call with size 100
		//  and another call with size 1
		//  a more balanced work load is
		//  one call with size 51
		//  and one call with size 50
		//  -------------------------
		const SizeType nblocks_i = (m + (nb_ - 1))/nb_;
		const SizeType nblocks_j = (n + (nb_ - 1))/nb_;
		int const nb_i = ((m % nblocks_i) == 0) ?
		            (m/nblocks_i) :
		            ((m/nblocks_i) + 1);
		int const nb_j = ((n % nblocks_j) == 0) ?
		            (n/nblocks_j) :
		            ((n/nblocks_j) + 1);

		bool const is_transA = (transA == 'T') || (transA == 't');
		bool const is_conjA  = (transA == 'C') || (transA == 'c');
		bool const is_notransA = (!is_transA) && (!is_conjA);

		bool const is_transB = (transB == 'T') || (transB == 't');
		bool const is_conjB  = (transB == 'C') || (transB == 'c');
		bool const is_notransB = (!is_transB) && (!is_conjB);


		if (idebug_ >= 1) {
			std::cout << " GEMMR: "
			          << " m " << m
			          << " n " << n
			          << " k " << k
			          << " nb " << nb_
			          << " nb_i " << nb_i
			          << " nb_j " << nb_j
			          << " nblocks_i " << nblocks_i
			          << " nblocks_j " << nblocks_j
			          << "\n";
		}

		auto lambda = [transA, transB, m, n, k, alpha, &A, ldA, &B, ldB, beta, &C, ldC,
		        nblocks_i, nblocks_j, nb_i, nb_j, is_notransA, is_notransB]
		        (SizeType ij_block, SizeType)
		{
			const SizeType i_block = (ij_block % nblocks_i);
			const SizeType j_block = (ij_block - i_block)/nblocks_i;
			assert( (i_block + j_block*nblocks_i) == ij_block );

			int const ic_start = 1 + (i_block)*nb_i;
			int const ic_end   = std::min( m, ic_start + nb_i - 1);
			int const jc_start = 1 + (j_block)*nb_j;
			int const jc_end   = std::min( n, jc_start + nb_j - 1);

			// --------------------------------------------
			// compute   C(ic_start:ic_end, jc_start:jc_end) =
			//              A(ic_start:ic_end, 1:k) * B( 1:k, jc_start:jc_end)
			// --------------------------------------------

			int const mm = (ic_end - ic_start + 1);
			int const nn = (jc_end - jc_start + 1);

			// -----------------
			// no splitting in k
			// -----------------
			int const kk = k;

			int const ia = (is_notransA) ? ic_start : 1;
			int const ja = (is_notransA) ? 1        : ic_start;

			int const ib = (is_notransB) ? 1        : jc_start;
			int const jb = (is_notransB) ? jc_start : 1;

			int const ic = ic_start;
			int const jc = jc_start;

			T const * const pA = &(A[ ia-1 + (ja-1)*ldA ]);
			T const * const pB = &(B[ ib-1 + (jb-1)*ldB ]);
			T       * const pC = &(C[ ic-1 + (jc-1)*ldC ]);

			psimag::BLAS::GEMM( transA, transB,
			                    mm,nn,kk,
			                    alpha,
			                    pA, ldA,
			                    pB, ldB,
			                    beta,
			                    pC, ldC);
		};

		const SizeType tasks = nblocks_i * nblocks_j;
		CodeSectionParams csp = Concurrency::codeSectionParams;
		csp.npthreads = std::min(tasks, nthreads_);
		Parallelizer2<> parallelizer2(csp);
		parallelizer2.parallelFor(0, tasks, lambda);
	}

	bool isSmall(SizeType m, SizeType n) const
	{
		return ( (m <= nb_) && (n <= nb_) );
	}

	bool idebug_;
	SizeType nb_;
	SizeType nthreads_;
}; // class GemmR

}

#endif
