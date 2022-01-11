#ifndef ADDITIONALONSITEHAMILTONIANEXT_H
#define ADDITIONALONSITEHAMILTONIANEXT_H
#include "CrsMatrix.h"

template<typename T>
void additionalOnSiteHamiltonianExt(PsimagLite::CrsMatrix<T>&,
                                    SizeType,
                                    typename PsimagLite::Real<T>::Type);

#endif // ADDITIONALONSITEHAMILTONIANEXT_H
