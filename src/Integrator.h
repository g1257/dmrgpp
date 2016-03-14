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

	enum IntegrationEnum {
		INTEG_QAG, INTEG_QAGP
	};

	Integrator(FunctionType& function,
	           RealType epsabs = 1e-9,
	           RealType epsrel = 1e-9,
	           SizeType limit = 1000000)
	    : epsabs_(epsabs),
	      epsrel_(epsrel),
	      limit_(limit),
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

	RealType operator()(VectorRealType& pts,
	                    IntegrationEnum integ = INTEG_QAG,
	                    int key = 4)
	{
		switch (integ) {
		default:
			return qag(pts,key);
		case INTEG_QAGP:
			return qagp(pts);
		}
	}

	RealType toInfinity(RealType a,
	                    IntegrationEnum integ = INTEG_QAG,
	                    int key = 4)
	{
		int status = gslWrapper_.gsl_integration_qagiu(&f_,
		                                               a,
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

	RealType operator()()
	{
		int status = gslWrapper_.gsl_integration_qagi(&f_,
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

	RealType qagp(VectorRealType& pts)
	{
		int status = gslWrapper_.gsl_integration_qagp(&f_,
		                                              &(pts[0]),
		        pts.size(),
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

	RealType qag(const VectorRealType& pts, int key)
	{
		int status = gslWrapper_.gsl_integration_qag(&f_,
		                                             pts[0],
		        pts[1],
		        epsabs_,
		        epsrel_,
		        limit_,
		        key,
		        workspace_,
		        &result_,
		        &abserr_);

		if (status)
			gslWrapper_.printError(status);

		return result_;
	}

	GslWrapperType gslWrapper_;
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

