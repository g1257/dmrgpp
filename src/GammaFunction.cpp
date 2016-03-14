#include "GammaFunction.h"

namespace PsimagLite {
std::complex<double> LnGammaFunction(const std::complex<double>& z)
{
	GslWrapper gslWrapper;
	double lnr = 0.0;
	double arg = 0.0;
	gslWrapper.gsl_sf_lngamma_complex_e(std::real(z),std::imag(z),&lnr,&arg);
	return std::complex<double>(lnr,arg);
}

}

