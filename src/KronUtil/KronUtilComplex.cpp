#include "KronUtil.h"
#include <string>

void csr_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::CrsMatrix<std::complex<double> >&,
                   const PsimagLite::CrsMatrix<std::complex<double> >&,
                   const std::complex<double>* yin,
                   std::complex<double>* xout )
{
	std::string msg("csr_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}

void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<std::complex<double> >& a_,
                       const PsimagLite::CrsMatrix<std::complex<double> >&,
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
                       const PsimagLite::CrsMatrix<std::complex<double> >&,
                       const PsimagLite::Matrix<std::complex<double> >& b_,
                       const std::complex<double>* yin,
                       std::complex<double>* xout)
{
	std::string msg("csr_den_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}
