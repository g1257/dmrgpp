#include "AdditionalOnSiteHamiltonianExt.h"

template<typename T>
void additionalOnSiteHamiltonianExt(PsimagLite::CrsMatrix<T>&,
                                    SizeType,
                                    typename PsimagLite::Real<T>::Type)
{
}

template
void additionalOnSiteHamiltonianExt<double>(PsimagLite::CrsMatrix<double>&,
                                    SizeType,
                                    double);

template
void additionalOnSiteHamiltonianExt<std::complex<double> >
(PsimagLite::CrsMatrix<std::complex<double> >&,
                                    SizeType,
                                    double);
