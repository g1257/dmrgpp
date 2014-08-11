
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

	RealType operator()() const
	{
		// INPUT: Function f, endpoint values a, b, tolerance TOL,
		//                maximum iterations NMAX
		// CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or
		//                           f(a) > 0 and f(b) < 0
		// OUTPUT: value which differs from a root of f(x)=0 by less than TOL

		RealType a = a_;
		RealType b = b_;
		SizeType n = 0;
		RealType functionAtA = function_(a);
		while (n<maxIter_) {
			RealType c = (a + b)/2; // new midpoint
			RealType tmp = function_(c);
			if (fabs(tmp) < tolerance_) { //solution found
				return c;
			}
			n++;
			if (sign(tmp) == sign(functionAtA)) {
				a = c;
				functionAtA = tmp;
			} else {
				b = c;
			}
		}

		PsimagLite::String msg("RootFindBisection, too many iterations ");
		msg += "maximum = " + ttos(maxIter_) + " tolerance = ";
		msg += ttos(tolerance_) + " a = " + ttos(a_) + " b = " + ttos(b_) + "\n";
		throw std::runtime_error(msg);
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

