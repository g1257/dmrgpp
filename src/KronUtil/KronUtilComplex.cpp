#include "KronUtil.h"
#include <string>

void csr_kron_mult(const char transA,
                   const char transB,
                   const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& arowptr,
                   const PsimagLite::Vector<int>::Type& acol,
                   const PsimagLite::Vector<std::complex<double> >::Type& aval,

                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& browptr,
                   const PsimagLite::Vector<int>::Type& bcol,
                   const PsimagLite::Vector<std::complex<double> >::Type& bval,

                   const std::complex<double>* yin,
                   std::complex<double>* xout )
{
	std::string msg("csr_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}

void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<std::complex<double> >& a_,

                       const int nrow_B,
                       const int ncol_B,
                       const PsimagLite::Vector<int>::Type& browptr,
                       const PsimagLite::Vector<int>::Type& bcol,
                       const PsimagLite::Vector<std::complex<double> >::Type& bval,

                       const std::complex<double>* yin,
                       std::complex<double>* xout)
{
	std::string msg("den_csr_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}


void den_kron_mult(const char transA,
                   const char transB,


                   const PsimagLite::Matrix<std::complex<double> >& a_,
                   const PsimagLite::Matrix<std::complex<double> >& b_,

                   const std::complex<double>* yin,
                   std::complex<double>* xout)
{
	std::string msg("den_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}

 void csr_den_kron_mult(const char transA,
                               const char transB,

                               const int nrow_A,
                               const int ncol_A,
                               const PsimagLite::Vector<int>::Type& arowptr,
                               const PsimagLite::Vector<int>::Type& acol,
                               const PsimagLite::Vector<std::complex<double> >::Type& aval,
                               const PsimagLite::Matrix<std::complex<double> >& b_,

                               const std::complex<double>* yin,
                               std::complex<double>* xout)
 {
	 std::string msg("csr_den_kron_mult: Unimplemented for complex");
	 throw std::runtime_error(msg);
 }
