#include "KronUtil.h"
#include <string>

void csr_kron_mult(const char transA,
                   const char transB,
                   const int nrow_A,
                   const int ncol_A,
                   const int arowptr[],
                   const int acol[],
                   const std::complex<double> aval[],

                   const int nrow_B,
                   const int ncol_B,
                   const int browptr[],
                   const int bcol[],
                   const std::complex<double> bval[],

                   const std::complex<double> yin[],
                   std::complex<double> xout[] )
{
	std::string msg("csr_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}

void den_csr_kron_mult(const char transA,
                       const char transB,
                       const int nrow_A,
                       const int ncol_A,
                       const std::complex<double> a_[],

                       const int nrow_B,
                       const int ncol_B,
                       const int browptr[],
                       const int bcol[],
                       const std::complex<double> bval[],

                       const std::complex<double> yin[],
                       std::complex<double> xout[])
{
	std::string msg("den_csr_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}


void den_kron_mult(const char transA,
                   const char transB,

                   const int nrow_A,
                   const int ncol_A,
                   const std::complex<double> a_[],

                   const int nrow_B,
                   const int ncol_B,
                   const std::complex<double> b_[],

                   const std::complex<double> yin[],
                   std::complex<double> xout[])
{
	std::string msg("den_kron_mult: Unimplemented for complex");
	throw std::runtime_error(msg);
}

 void csr_den_kron_mult( const char transA,
                               const char transB,

                               const int nrow_A,
                               const int ncol_A,
                               const int arowptr[],
                               const int acol[],
                               const std::complex<double> aval[],

                               const int nrow_B,
                               const int ncol_B,
                               const std::complex<double> b_[],

                               const std::complex<double> yin[],
                               std::complex<double> xout[])
 {
	 std::string msg("csr_den_kron_mult: Unimplemented for complex");
	 throw std::runtime_error(msg);
 }
