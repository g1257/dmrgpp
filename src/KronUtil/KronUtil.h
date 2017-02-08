#ifndef KRON_UTIL_HEADER_H
#define KRON_UTIL_HEADER_H
#include <complex>

extern
void csr_kron_mult(const char transA,
                   const char transB,
                   const int nrow_A,
                   const int ncol_A,
                   const int arowptr[],
                   const int acol[],
                   const double aval[],

                   const int nrow_B,
                   const int ncol_B,
                   const int browptr[],
                   const int bcol[],
                   const double bval[],

                   const double yin[],
                   double xout[] );


extern
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const int nrow_A,
                       const int ncol_A,
                       const double a_[],

                       const int nrow_B,
                       const int ncol_B,
                       const int browptr[],
                       const int bcol[],
                       const double bval[],

                       const double yin[],
                       double xout[] );

extern
void den_kron_mult(const char transA,
                   const char transB,

                   const int nrow_A,
                   const int ncol_A,
                   const double a_[],

                   const int nrow_B,
                   const int ncol_B,
                   const double b_[],

                   const double yin[],
                   double xout[] );

extern void csr_den_kron_mult( const char transA,
                               const char transB,

                               const int nrow_A,
                               const int ncol_A,
                               const int arowptr[],
                               const int acol[],
                               const double aval[],

                               const int nrow_B,
                               const int ncol_B,
                               const double b_[],

                               const double yin[],
                               double xout[] );
// complex functions below

extern
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
                   std::complex<double> xout[] );


extern
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
                       std::complex<double> xout[] );

extern
void den_kron_mult(const char transA,
                   const char transB,

                   const int nrow_A,
                   const int ncol_A,
                   const std::complex<double> a_[],

                   const int nrow_B,
                   const int ncol_B,
                   const std::complex<double> b_[],

                   const std::complex<double> yin[],
                   std::complex<double> xout[] );

extern void csr_den_kron_mult( const char transA,
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
                               std::complex<double> xout[] );
#endif

