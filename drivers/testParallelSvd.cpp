#include "BLAS.h"
#include "GemmR.h"
#include "Random48.h"
#include "Matrix.h"
#include "Svd.h"

typedef double RealType;

template<typename SomeRngType>
void fillRandom(PsimagLite::Matrix<RealType>& m,
                RealType min,
                RealType max,
                SomeRngType& rng)
{
	const SizeType rows = m.rows();
	const SizeType cols = m.cols();
	for (SizeType i = 0; i < rows; ++i)
		for (SizeType j = 0; j < cols; ++j)
			m(i, j) = min + rng()*max;
}

bool equalMatrices(PsimagLite::Matrix<RealType>& a,
                   PsimagLite::Matrix<RealType>& b,
                   RealType tolerance)
{
	const SizeType rows = a.rows();
	const SizeType cols = b.cols();
	if (rows != b.rows() || cols != b.cols()) return false;
	for (SizeType i = 0; i < rows; ++i)
		for (SizeType j = 0; j < cols; ++j)
			if (fabs(a(i, j) - b(i, j)) > tolerance)
				return false;

	return true;
}

int main(int argc, char ** argv)
{
	if (argc < 2)
		throw PsimagLite::RuntimeError("USAGE: " + PsimagLite::String(argv[0])
	        + " total nthreadsOuter nthreadsInner\n");

	int total = atoi(argv[1]);
	int nthreadsOuter = atoi(argv[2]);
	int nthreadsInner = atoi(argv[3]);

	PsimagLite::Concurrency concurrency(&argc, &argv, nthreadsInner);

	PsimagLite::Random48<RealType> rng(1234);

	auto lambda = [&rng, nthreadsInner](SizeType, SizeType) {
		SizeType lda = static_cast<SizeType>(rng()*500) + 10;
		SizeType cda = lda;
		PsimagLite::Matrix<RealType> A(lda, cda);

		fillRandom(A, -10, 10, rng);

		PsimagLite::Svd<RealType> svd;

		PsimagLite::Vector<RealType>::Type s(lda);
		PsimagLite::Matrix<RealType> vt;

		svd('A', A, s, vt);
	};

	PsimagLite::CodeSectionParams csp = PsimagLite::Concurrency::codeSectionParams;
	csp.npthreads = nthreadsOuter;

	PsimagLite::Parallelizer2<> parallelizer2(csp);
	parallelizer2.parallelFor(0, total, lambda);

}
