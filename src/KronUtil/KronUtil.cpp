#include "csr_kron_mult.cpp"
#include "den_csr_kron_mult.cpp"
#include "den_kron_mult.cpp"
#include "csr_den_kron_mult.cpp"
#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

template
void csr_kron_mult<RealType>(const char transA,
                           const char transB,
                           const PsimagLite::CrsMatrix<RealType>& a,
                           const PsimagLite::CrsMatrix<RealType>& b,
                           const PsimagLite::Vector<RealType>::Type& yin,
                           SizeType offsetY,
                           PsimagLite::Vector<RealType>::Type& xout,
                           SizeType offsetX);

template
void csr_kron_mult
<std::complex<RealType> >(const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<std::complex<RealType> >&,
                        const PsimagLite::CrsMatrix<std::complex<RealType> >&,
                        const PsimagLite::Vector<std::complex<RealType> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<RealType> >::Type& xout,
                        SizeType offsetX);

//-----------------------------------------------------------------------------------

template
void den_csr_kron_mult<RealType>(const char transA,
                               const char transB,
                               const PsimagLite::Matrix<RealType>& a_,
                               const PsimagLite::CrsMatrix<RealType>&,
                               const PsimagLite::Vector<RealType>::Type& yin,
                               SizeType offsetY,
                               PsimagLite::Vector<RealType>::Type& xout,
                               SizeType offsetX);
template
void den_csr_kron_mult
<std::complex<RealType> >(const char transA,
                        const char transB,
                        const PsimagLite::Matrix<std::complex<RealType> >& a_,
                        const PsimagLite::CrsMatrix<std::complex<RealType> >&,
                        const PsimagLite::Vector<std::complex<RealType> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<RealType> >::Type& xout,
                        SizeType offsetX);


//-----------------------------------------------------------------------------------

template
void den_kron_mult<RealType>(const char transA,
                           const char transB,
                           const PsimagLite::Matrix<RealType>& a_,
                           const PsimagLite::Matrix<RealType>& b_,
                           const PsimagLite::Vector<RealType>::Type& yin,
                           SizeType offsetY,
                           PsimagLite::Vector<RealType>::Type& xout,
                           SizeType offsetX);

template
void den_kron_mult
<std::complex<RealType> >(const char transA,
                        const char transB,
                        const PsimagLite::Matrix<std::complex<RealType> >& a_,
                        const PsimagLite::Matrix<std::complex<RealType> >& b_,
                        const PsimagLite::Vector<std::complex<RealType> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<RealType> >::Type& xout,
                        SizeType offsetX);


//-----------------------------------------------------------------------------------

template
void csr_den_kron_mult<RealType>(const char transA,
                                const char transB,
                                const PsimagLite::CrsMatrix<RealType>&,
                                const PsimagLite::Matrix<RealType>& b_,
                                const PsimagLite::Vector<RealType>::Type& yin,
                                SizeType offsetY,
                                PsimagLite::Vector<RealType>::Type& xout,
                                SizeType offsetX);

template
void csr_den_kron_mult
<std::complex<RealType> >(const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<std::complex<RealType> >&,
                        const PsimagLite::Matrix<std::complex<RealType> >& b_,
                        const PsimagLite::Vector<std::complex<RealType> >::Type& yin,
                        SizeType offsetY,
                        PsimagLite::Vector<std::complex<RealType> >::Type& xout,
                        SizeType offsetX);


