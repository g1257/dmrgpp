
/** \ingroup PsimagLite */
/*@{*/

/*! \file RootFindingBisection.h
 *
 *  RootFindingBisection such as chemical potential
 *
 */
#ifndef ROOT_FIND_BISECT_H
#define ROOT_FIND_BISECT_H
#include "ProgressIndicator.h"
#include "Matrix.h"
#include "Fermi.h"

namespace PsimagLite {
template<typename FunctionType>
class RootFindingBisection {
	typedef typename FunctionType::RealType RealType;

public:
	RootFindingBisection(const FunctionType& function,
	                     RealType a = -100.,
	                     RealType b = 100.,
	                     SizeType maxIter=100,
	                     RealType tolerance=1.0e-3)
	    : function_(function),maxIter_(maxIter),tolerance_(tolerance),
	      a_(a),b_(b)
	{
		RealType f1 = function_(a_) * function_(b_);
		if (f1 >= 0)
			throw PsimagLite::RuntimeError("RootFindingBisection: condition not met\n");
	}

	void operator()(RealType& mu)
	{
		// INPUT: Function f, endpoint values a, b, tolerance TOL,
		//                maximum iterations NMAX
		// CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or
		//                           f(a) > 0 and f(b) < 0
		// OUTPUT: value which differs from a root of f(x)=0 by less than TOL

		SizeType n = 0;
		while (n<maxIter_) {
			RealType c = (a_ + b_)/2; // new midpoint

			if (function_(c) == 0 || (b_ - a_)/2 < tolerance_) { //solution found
				mu = c;
				return;
			}
			n++;
			if (sign(function_(c)) == sign(function_(a_)))
				a_ = c;
			else
				b_ = c;
		}

		throw std::runtime_error("RootFindBisection failed\n");
	}

private:

	int sign(const RealType& x) const
	{
		return (x>=0) ?  1  : -1;
	}

	const FunctionType& function_;
	SizeType maxIter_;
	RealType tolerance_;
	RealType a_,b_;
}; // RootFindingBisection
} // namespace PsimagLite

/*@}*/
#endif// ROOT_FIND_BISECT_H

