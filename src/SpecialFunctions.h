#ifndef PSI_SPECIAL_FUNCTIONS_H
#define PSI_SPECIAL_FUNCTIONS_H
#include "GslWrapper.h"

namespace PsimagLite
{

std::complex<double> LnGammaFunction(const std::complex<double>&);

double Ci(const double& x);

} // namespace PsimagLite
#endif // PSI_SPECIAL_FUNCTIONS_H
