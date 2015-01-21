
/*
*/

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <iostream>
#include "Vector.h"
#include <stdexcept>

#ifdef USE_GSL
extern "C" {
#include <gsl/gsl_multimin.h>
}

namespace PsimagLite {

template<typename FunctionType>
typename FunctionType::FieldType myFunction(const gsl_vector *v, void *params)
{
	FunctionType* ft = (FunctionType *)params;
	return ft->operator()(v->data,v->size);
}

template<typename FunctionType>
void myDfunction(const gsl_vector *v,
                 void *params,
                 gsl_vector* df)
{
	FunctionType* ft = (FunctionType *)params;
	ft->df(v->data,v->size,df->data,df->size);
}

template<typename FunctionType>
void myFdFunction(const gsl_vector *v,
                  void *params,
                  double *f,
                  gsl_vector *df)
{
	*f = myFunction<FunctionType>(v, params);
	myDfunction<FunctionType>(v, params, df);
}

template<typename RealType,typename FunctionType>
class Minimizer {

	typedef typename FunctionType::FieldType FieldType;
	typedef typename Vector<FieldType>::Type VectorType;
	typedef Minimizer<RealType,FunctionType> ThisType;

public:

	enum {GSL_SUCCESS=::GSL_SUCCESS, GSL_CONTINUE=::GSL_CONTINUE};

	Minimizer(FunctionType& function,SizeType maxIter, bool verbose = false)
	    : function_(function),
	      maxIter_(maxIter),
	      verbose_(verbose),
	      status_(100),
	      gslT_(gsl_multimin_fminimizer_nmsimplex2),
	      gslS_(gsl_multimin_fminimizer_alloc(gslT_,function_.size())),
	      gslDt_(gsl_multimin_fdfminimizer_conjugate_fr),
	      gslDs_(gsl_multimin_fdfminimizer_alloc (gslDt_, function_.size()))
	{}

	~Minimizer()
	{
		gsl_multimin_fminimizer_free(gslS_);
		gsl_multimin_fdfminimizer_free(gslDs_);
	}

	int simplex(VectorType& minVector,RealType delta=1e-3,RealType tolerance=1e-3)
	{
		gsl_vector *x;
		/* Starting point,  */
		x = gsl_vector_alloc (function_.size());
		for (SizeType i=0;i<minVector.size();i++)
			gsl_vector_set (x, i, minVector[i]);

		gsl_vector *xs;
		xs = gsl_vector_alloc (function_.size());
		for (SizeType i=0;i<minVector.size();i++)
			gsl_vector_set (xs, i, delta);

		gsl_multimin_function func;
		func.f= myFunction<FunctionType>;
		func.n = function_.size();
		func.params = &function_;
		gsl_multimin_fminimizer_set (gslS_, &func, x, xs);

		SizeType iter = 0;
		RealType prevValue = 0;

		for (;iter<maxIter_;iter++) {
			status_ = gsl_multimin_fminimizer_iterate (gslS_);

			if (status_) {
				PsimagLite::String gslError(gsl_strerror(status_));
				PsimagLite::String msg("Minimizer::simplex(...): GSL Error: ");
				msg += gslError + "\n";
				throw RuntimeError(msg);
			}

			RealType size = gsl_multimin_fminimizer_size(gslS_);
			status_ = gsl_multimin_test_size(size, tolerance);

			if (verbose_) {
				RealType thisValue = function_(gslS_->x->data,func.n);
				RealType diff = fabs(thisValue - prevValue);
				std::cerr<<"simplex: "<<iter<<" "<<thisValue<<" diff= "<<diff;
				std::cerr<<" status= "<<status_<<" size="<<size<<"\n";
				prevValue = thisValue;
			}

			if (status_ == GSL_SUCCESS) break;
		}

		found(minVector,gslS_->x,iter);
		gsl_vector_free (x);
		gsl_vector_free (xs);
		return iter;
	}

	int conjugateGradient(VectorType& minVector,
	                      RealType delta=1e-3,
	                      RealType delta2=1e-3,
	                      RealType tolerance=1e-3)
	{
		gsl_vector *x;
		/* Starting point,  */
		x = gsl_vector_alloc (function_.size());
		for (SizeType i=0;i<minVector.size();i++)
			gsl_vector_set (x, i, minVector[i]);

		gsl_multimin_function_fdf func;
		func.f= myFunction<FunctionType>;
		func.df = myDfunction<FunctionType>;
		func.fdf = myFdFunction<FunctionType>;
		func.n = function_.size();
		func.params = &function_;

		gsl_multimin_fdfminimizer_set(gslDs_, &func, x, delta, delta2);

		RealType prevValue = 0;
		SizeType iter = 0;
		for (;iter<maxIter_;iter++) {
			status_ = gsl_multimin_fdfminimizer_iterate(gslDs_);

			if (status_) {
				PsimagLite::String gslError(gsl_strerror(status_));
				PsimagLite::String msg("Minimizer::conjugateGradient(...): GSL Error: ");
				msg += gslError + "\n";
				throw RuntimeError(msg);
			}

			status_ = gsl_multimin_test_gradient (gslDs_->gradient, tolerance);

			if (verbose_) {
				RealType thisValue = function_(gslDs_->x->data,func.n);
				RealType diff = fabs(thisValue - prevValue);
				std::cerr<<"conjugateGradient: "<<iter<<" "<<thisValue;
				std::cerr<<" diff= "<<diff;
				std::cerr<<" status= "<<status_<<"\n";
				prevValue = thisValue;
			}

			if (status_ == GSL_SUCCESS) break;

		}

		found(minVector,gslDs_->x,iter);
		gsl_vector_free (x);
		return iter;
	}

	int status() const { return status_; }

	PsimagLite::String statusString() const
	{
		switch (status_) {
		case GSL_SUCCESS:
			return "GSL_SUCCESS";
			break;
		case GSL_CONTINUE:
			return "GSL_CONTINUE";
			break;
		default:
			return "UNKNOWN";
			break;
		}
	}

private:

	void found(VectorType& minVector,gsl_vector* x,SizeType iter)
	{
		for (SizeType i=0;i<minVector.size();i++)
			minVector[i] = gsl_vector_get(x,i);
	}

	FunctionType& function_;
	SizeType maxIter_;
	bool verbose_;
	int status_;
	const gsl_multimin_fminimizer_type *gslT_;
	gsl_multimin_fminimizer *gslS_;
	const gsl_multimin_fdfminimizer_type *gslDt_;
	gsl_multimin_fdfminimizer *gslDs_;

}; // class Minimizer
}

#else

namespace PsimagLite {

template<typename RealType,typename FunctionType>
class Minimizer {

	typedef typename FunctionType::FieldType FieldType;
	typedef typename Vector<FieldType>::Type VectorType;

public:

	enum {GSL_SUCCESS=0, GSL_CONTINUE=1};

	Minimizer(FunctionType&,SizeType, bool = false)
	{
		PsimagLite::String str("Minimizer needs the gsl\n");
		throw PsimagLite::RuntimeError(str);
	}

	int simplex(VectorType&,RealType=1e-3,RealType=1e-3)
	{
		PsimagLite::String str("Minimizer needs the gsl\n");
		throw PsimagLite::RuntimeError(str);
	}

	int status() const { return 1; }

	PsimagLite::String statusString() const
	{
		return "Minimizer needs the gsl";
	}
};

} // namespace PsimagLite
#endif

#endif // MINIMIZER_H

