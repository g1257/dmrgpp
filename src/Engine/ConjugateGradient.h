/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/

/** \ingroup DMRG */
/*@{*/

/*! \file ConjugateGradient.h
 *
 *  impl. of the conjugate gradient method
 *
 */
#ifndef CONJ_GRAD_H
#define CONJ_GRAD_H

#include "Matrix.h"
#include "Vector.h"
#include "ProgressIndicator.h"

namespace Dmrg {

template<typename MatrixType>
class	ConjugateGradient {
	typedef typename MatrixType::value_type FieldType;
	typedef typename PsimagLite::Vector<FieldType>::Type VectorType;
	typedef typename PsimagLite::Real<FieldType>::Type RealType;

public:
	ConjugateGradient(SizeType max,RealType eps)
	    : progress_("ConjugateGradient"), max_(max), eps_(eps) {}

	//! A and b, the result x, and also the initial solution x0
	void operator()(VectorType& x,
	                const MatrixType& A,
	                const VectorType& b) const
	{
		VectorType v = multiply(A,x);
		VectorType p(b.size());
		VectorType rprev(b.size());
		VectorType rnext;
		for (SizeType i=0;i<rprev.size();i++) {
			rprev[i] = b[i] - v[i];
			p[i] = rprev[i];
		}

		SizeType k = 0;
		while (k<max_) {
			VectorType tmp = multiply(A,p);
			FieldType scalarrprev = scalarProduct(rprev,rprev);
			FieldType val = scalarrprev/scalarProduct(p,tmp);
			v <= x + val * p;
			x = v;
			v <= rprev - val * tmp;
			rnext = v;
			if (PsimagLite::norm(rnext)<eps_) break;
			val = scalarProduct(rnext,rnext)/scalarrprev;
			v <= rnext - val*p;
			p = v;
			rprev = rnext;
			k++;
		}

		PsimagLite::OstringStream msg;
		msg<<"Finished after "<<k<<" steps out of "<<max_;
		msg<<" requested eps= "<<eps_;
		RealType finalEps = PsimagLite::norm(rnext);
		msg<<" actual eps= "<<finalEps;
		progress_.printline(msg,std::cout);

		if (finalEps <= eps_) return;

		PsimagLite::OstringStream msg2;
		msg2<<"WARNING: actual eps "<<finalEps<<" greater than requested eps= "<<eps_;
		progress_.printline(msg2,std::cout);
	}

private:

	FieldType scalarProduct(const VectorType& v1,const VectorType& v2) const
	{
		FieldType sum = 0;
		for (SizeType i=0;i<v1.size();i++) sum += PsimagLite::conj(v1[i])*v2[i];
		return sum;
	}

	VectorType multiply(const MatrixType& A,const VectorType& v) const
	{
		VectorType y(A.rows(),0);
		A.matrixVectorProduct(y,v);
		return y;
	}

	PsimagLite::ProgressIndicator progress_;
	SizeType max_;
	RealType eps_;
}; // class ConjugateGradient

} // namespace Dmrg

/*@}*/
#endif // CONJ_GRAD_H

