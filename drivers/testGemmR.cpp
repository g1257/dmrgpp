#include <iostream>
#include "omp.h"
#include "Matrix.h"
#include "GemmR.h"

typedef std::complex<double> zcomplex;

template<typename T>
T make_val(double const x, double const y) {
	return(x);
}

template<>
double make_val<double>(double const x, double const y) {
	return( x );
}

template<>
zcomplex make_val<zcomplex>(double const x, double const y) {
	zcomplex z(x,y);
	return( z );
}

template<typename T>
int test_GEMMR(int const Mmax,
               int const Nmax,
               int const Kmax,
               int const nb,
               bool needsPrinting)
{
	int const idebug = (needsPrinting) ? 1 : 0;
	int nerrors = 0;

	char const trans_table[3] = {'N', 'T', 'C'};

	T const alpha = make_val<T>(1.1, 2.1);
	T const beta  = make_val<T>(3.1, 4.1);

	PsimagLite::GemmR<T> gemmR(needsPrinting,
	                           nb,
	                           PsimagLite::Concurrency::codeSectionParams.npthreads);

	for(int k=1; k <= Kmax; k += nb) {
		for(int n=1; n <= Nmax; n += nb) {
			for(int m=1; m <= Mmax; m += nb) {
				for(int itransB=0; itransB < 3; itransB++) {
					for(int itransA=0; itransA < 3; itransA++) {

						char const transA = trans_table[itransA];
						char const transB = trans_table[itransB];

						bool const is_transA = (transA == 'T') || (transA == 't');
						bool const is_transB = (transB == 'T') || (transB == 't');
						bool const is_conjA  = (transA == 'C') || (transA == 'c');
						bool const is_conjB  = (transB == 'C') || (transB == 'c');
						bool const is_notransA = (!is_transA) && (!is_conjA);
						bool const is_notransB = (!is_transB) && (!is_conjB);

						int const mC = m;
						int const nC = n;
						int const mA = (is_notransA) ? mC : k;
						int const nA = (is_notransA) ? k  : mC;
						int const mB = (is_notransB) ? k  : nC;
						int const nB = (is_notransB) ? nC : k;

						PsimagLite::Matrix<T> C(mC,nC);
						PsimagLite::Matrix<T> C_gemmr(mC,nC);
						PsimagLite::Matrix<T> A(mA,nA);
						PsimagLite::Matrix<T> B(mB,nB);

						int const ldA = mA;
						int const ldB = mB;
						int const ldC = mC;

						for(int j=0; j < nC; j++) {
							for(int i=0; i < mC; i++) {
								T cij = make_val<T>( 1.0*(i+j)/(mC+nC), 1.0*i*j/(mC*nC) );
								C(i,j) = cij;
								C_gemmr(i,j) = cij;

							}
						}

						for(int j=0; j < nA; j++) {
							for(int i=0; i < mA; i++) {
								T aij = make_val<T>(-1.0*(i+j+1), 1.0*(j-i+1));
								A(i,j) = aij;
							}
						}

						for(int j=0; j < nB; j++) {
							for(int i=0; i < mB; i++) {
								T bij = make_val<T>(1.0*(i+j+1)/(mB*nB),
								                    -1.0*(j-i+1)/(mB*nB));
								B(i,j) = bij;
							}
						}

						gemmR(transA, transB,
						      m,n,k,
						      alpha, &(A(0,0)), ldA,
						      &(B(0,0)), ldB,
						      beta,  &(C_gemmr(0,0)), ldC);

						psimag::BLAS::GEMM(
						            transA, transB,
						            m,n,k,
						            alpha, &(A(0,0)), ldA,
						            &(B(0,0)), ldB,
						            beta,  &(C(0,0)), ldC);

						double max_err = 0;
						double c_norm = 0;
						for(int j=0; j < nC; j++) {
							for(int i=0; i < mC; i++) {
								double const err = std::abs( C(i,j) - C_gemmr(i,j) );
								max_err = std::max( max_err, err );
								c_norm += std::abs( C(i,j) ) ;
							}
						}

						double const tol = 0.0000001;
						bool const isok = (max_err < tol );
						if (!isok) {
							nerrors++;
						}
						if ((!isok) || (idebug >= 1)) {
							std::cout << " transA " << transA
							          << " transB " << transB
							          << " m " << m
							          << " n " << n
							          << " k " << k
							          << " max_err " << max_err
							          << " c_norm " << c_norm
							          << "\n";
						}
					}
				}
			}
		}
	}

	return(nerrors);
}

int main(int argc, char ** argv)
{
	int const Nmax = 300;
	int const Mmax = 301;
	int const Kmax = 302;
	int nerr_zcomplex = 0;

	if (argc < 2)
		throw PsimagLite::RuntimeError("USAGE: " + PsimagLite::String(argv[0])
	        + " nthreads [nb] [debug]\n");

	int nthreads = atoi(argv[1]);

	int const nb = (argc >= 3) ? atoi(argv[2]) : 99;
	const bool needsPrinting = (argc == 4) ? atoi(argv[3]) > 0 : false;

	PsimagLite::Concurrency concurrency(&argc, &argv, nthreads);

	int nerr_double = test_GEMMR<double>(Mmax, Nmax, Kmax, nb, needsPrinting);
	if (nerr_double == 0) {
		nerr_zcomplex = test_GEMMR<zcomplex>(Mmax, Nmax, Kmax, nb, needsPrinting);
	}

	bool const all_passed = (nerr_double == 0) &&
	        (nerr_zcomplex == 0);
	if (all_passed) {
		std::cout << "ALL PASSED "
		          << "\n";
	} else {
		std::cout << " nerr_double = " << nerr_double
		          << " nerr_zcomplex = " << nerr_zcomplex
		          << "\n";
	}
}


