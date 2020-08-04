#include "Matrix.h"
#include <cstdlib> // for atoi
#include "BLAS.h"
#include "Random48.h"
#define USE_PTHREADS_OR_NOT_NG
#include "Parallelizer.h"

typedef double MyRealType;
typedef PsimagLite::Matrix<MyRealType> MatrixType;

class MyBlasWrapper {

public:

	typedef PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;

	MyBlasWrapper(SizeType m, SizeType n, SizeType k, SizeType total)
	    : a_(total), b_(total), c_(total), myrandom48_(1234)
	{
		for (SizeType i = 0; i < total; ++i) {
			SizeType factor = (i & 1) ? 2 : 1;
			a_[i] = new MatrixType(m*factor, k*factor);
			b_[i] = new MatrixType(k*factor, n*factor);
			c_[i] = new MatrixType(m*factor, n*factor);
			fillMatrix(*(a_[i]));
			fillMatrix(*(b_[i]));
		}
	}

	~MyBlasWrapper()
	{
		const SizeType total = a_.size();
		for (SizeType i = 0; i < total; ++i) {
			delete a_[i];
			delete b_[i];
			delete c_[i];
			a_[i] = b_[i] = c_[i] = 0;
		}
	}

	SizeType tasks() const { return a_.size(); }

	void doTask(SizeType ind, SizeType)
	{
		assert(a_[ind] && b_[ind] && c_[ind]);
		SizeType mm = a_[ind]->rows();
		SizeType nn = b_[ind]->cols();
		SizeType kk = a_[ind]->cols();
		assert(kk == b_[ind]->rows());
		assert(c_[ind]->rows() == mm && c_[ind]->cols() == nn);
		SizeType lda = a_[ind]->rows();
		SizeType ldb = b_[ind]->rows();
		SizeType ldc = c_[ind]->rows();
		MyRealType* aptr = &(a_[ind]->operator()(0,0));
		MyRealType* bptr = &(b_[ind]->operator()(0,0));
		MyRealType* cptr = &(c_[ind]->operator()(0,0));

		psimag::BLAS::GEMM( 'N',
		                    'N',
		                    mm, // rows of op(A)
		                    nn, // columns of op(B)
		                    kk, // columns of op(A)
		                    1.0,
		                    aptr,
		                    lda, // first dimension of A
		                    bptr,
		                    ldb, // first dimension of B
		                    0.0,
		                    cptr,
		                    ldc); // first dimension of C

		*(a_[ind]) = *(b_[ind]);
		*(b_[ind]) = *(c_[ind]);
	}

private:

	void fillMatrix(MatrixType& m)
	{
		for (SizeType j = 0; j < m.cols(); ++j)
			for (SizeType i = 0; i < m.rows(); ++i)
				m(i, j) = myrandom48_();
	}

	VectorMatrixType a_;
	VectorMatrixType b_;
	VectorMatrixType c_;
	PsimagLite::Random48<MyRealType> myrandom48_;
};

int main(int argc, char**argv)
{
	if (argc < 6) {
		std::cerr<<"USAGE: "<<argv[0]<<" m k n total threads\n";
		return 1;
	}

	SizeType m = atoi(argv[1]);
	SizeType k = atoi(argv[2]);
	SizeType n = atoi(argv[3]);
	SizeType total = atoi(argv[4]);
	SizeType threads = atoi(argv[5]);

	PsimagLite::CodeSectionParams codeSections(threads);
	PsimagLite::Parallelizer<MyBlasWrapper> parallel(codeSections);
	MyBlasWrapper myblasWrapper(m, n, k, total);

	parallel.loopCreate(myblasWrapper);
}
