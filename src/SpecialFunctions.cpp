#include "SpecialFunctions.h"

namespace PsimagLite {
std::complex<double> LnGammaFunction(const std::complex<double>& z)
{
	GslWrapper gslWrapper;
	GslWrapper::gsl_sf_result lnr;
	GslWrapper::gsl_sf_result arg;
	gslWrapper.gsl_sf_lngamma_complex_e(PsimagLite::real(z),PsimagLite::imag(z),&lnr,&arg);
	return std::complex<double>(lnr.val,arg.val);
}

double Ci(const double& x)
{
	GslWrapper gslWrapper;
	GslWrapper::gsl_sf_result result;
	gslWrapper.gsl_sf_Ci_e(x,&result);
	return result.val;
}

}

