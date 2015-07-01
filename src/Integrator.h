#ifndef PSI_INTEGRATOR_H
#define PSI_INTEGRATOR_H

#include "GslWrapper.h"

namespace PsimagLite {

template<typename FunctionType>
class Integrator {

public:

	typedef GslWrapper GslWrapperType;
	typedef typename FunctionType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Integrator(FunctionType& function)
	    : pts_(2),
	      epsabs_(1e-9),
	      epsrel_(1e-9),
	      limit_(1000000),
	      result_(0),
	      abserr_(0),
	      workspace_(gslWrapper_.gsl_integration_workspace_alloc(limit_+2))
	{
		f_.function= &FunctionType::function;

		f_.params = &function.params();
	}

	~Integrator()
	{
		gslWrapper_.gsl_integration_workspace_free (workspace_);
	}

	RealType operator()(const VectorRealType& pts)
	{
		pts_ = pts;
		int status = gslWrapper_.gsl_integration_qagp(&f_,
		                                 &(pts_[0]),
		        pts_.size(),
		        epsabs_,
		        epsrel_,
		        limit_,
		        workspace_,
		        &result_,
		        &abserr_);

		if (status)
			gslWrapper_.printError(status);

		return result_;
	}

private:

	GslWrapperType gslWrapper_;
	VectorRealType pts_;
	RealType epsabs_;
	RealType epsrel_;
	SizeType limit_;
	RealType result_;
	RealType abserr_;
	GslWrapperType::gsl_integration_workspace *workspace_;
	GslWrapperType::gsl_function f_;

};
} // namespace PsimagLite
#endif // PSI_INTEGRATOR_H

