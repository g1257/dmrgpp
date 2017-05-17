#include "csr_kron_mult.cpp"
#include "den_csr_kron_mult.cpp"
#include "den_kron_mult.cpp"
#include "csr_den_kron_mult.cpp"

template
void csr_kron_mult<double>(const char transA,
                           const char transB,
                           const PsimagLite::CrsMatrix<double>& a,
                           const PsimagLite::CrsMatrix<double>& b,
                           const PsimagLite::Vector<double>::Type& yin,
                           SizeType offsetY,
                           PsimagLite::Vector<double>::Type& xout,
                           SizeType offsetX);

template
void csr_kron_mult
<std::complex<double> >(const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<std::complex<double> >&,
                        const PsimagLite::CrsMatrix<std::complex<double> >&,
                        const PsimagLite::Vector<std::complex<double> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<double> >::Type& xout,
                        SizeType offsetX);

//-----------------------------------------------------------------------------------

template
void den_csr_kron_mult<double>(const char transA,
                               const char transB,
                               const PsimagLite::Matrix<double>& a_,
                               const PsimagLite::CrsMatrix<double>&,
                               const PsimagLite::Vector<double>::Type& yin,
                               SizeType offsetY,
                               PsimagLite::Vector<double>::Type& xout,
                               SizeType offsetX);
template
void den_csr_kron_mult
<std::complex<double> >(const char transA,
                        const char transB,
                        const PsimagLite::Matrix<std::complex<double> >& a_,
                        const PsimagLite::CrsMatrix<std::complex<double> >&,
                        const PsimagLite::Vector<std::complex<double> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<double> >::Type& xout,
                        SizeType offsetX);


//-----------------------------------------------------------------------------------

template
void den_kron_mult<double>(const char transA,
                           const char transB,
                           const PsimagLite::Matrix<double>& a_,
                           const PsimagLite::Matrix<double>& b_,
                           const PsimagLite::Vector<double>::Type& yin,
                           SizeType offsetY,
                           PsimagLite::Vector<double>::Type& xout,
                           SizeType offsetX);

template
void den_kron_mult
<std::complex<double> >(const char transA,
                        const char transB,
                        const PsimagLite::Matrix<std::complex<double> >& a_,
                        const PsimagLite::Matrix<std::complex<double> >& b_,
                        const PsimagLite::Vector<std::complex<double> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<double> >::Type& xout,
                        SizeType offsetX);


//-----------------------------------------------------------------------------------

template
void csr_den_kron_mult<double>(const char transA,
                                const char transB,
                                const PsimagLite::CrsMatrix<double>&,
                                const PsimagLite::Matrix<double>& b_,
                                const PsimagLite::Vector<double>::Type& yin,
                                SizeType offsetY,
                                PsimagLite::Vector<double>::Type& xout,
                                SizeType offsetX);

template
void csr_den_kron_mult
<std::complex<double> >(const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<std::complex<double> >&,
                        const PsimagLite::Matrix<std::complex<double> >& b_,
                        const PsimagLite::Vector<std::complex<double> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<double> >::Type& xout,
                        SizeType offsetX);


