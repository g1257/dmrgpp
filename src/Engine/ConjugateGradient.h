/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
		ConjugateGradient(SizeType max=1000,const RealType& eps = 1e-6)
		: progress_("ConjugateGradient"), max_(max), eps_(eps) {}

		//! A and b, the result x, and also the initial solution x0
		void operator()(typename PsimagLite::Vector<VectorType>::Type& x,
		                  const MatrixType& A,
		                  const typename PsimagLite::Vector<FieldType>::Type& b) const
		{
			VectorType v = multiply(A,x[0]);
			typename PsimagLite::Vector<VectorType>::Type r,p;
			r.push_back(b);
			p.push_back(b);
			for (SizeType i=0;i<r[0].size();i++) {
				r[0][i] = b[i] - v[i];
				p[0][i] = r[0][i];
			}
			SizeType k = 0;
			typename PsimagLite::Vector<FieldType>::Type alpha,beta;
			while(k<max_) {
				VectorType tmp = multiply(A,p[k]);
				FieldType val = scalarProduct(r[k],r[k])/
				           scalarProduct(p[k],tmp);
				alpha.push_back(val);
				v = x[k] + alpha[k] * p[k];
				x.push_back(v);
				v = r[k] - alpha[k] * tmp;
				r.push_back(v);
				if (PsimagLite::norm(r[k+1])<eps_) break;
				val = scalarProduct(r[k+1],r[k+1])/scalarProduct(r[k],r[k]);
				beta.push_back(val);
				v = r[k+1] - beta[k]*p[k];
				p.push_back(v);
				k++;
			}

			PsimagLite::OstringStream msg;
			msg<<"Finished after "<<k<<" steps out of "<<max_;
			progress_.printline(msg,std::cout);
		}

	private:

		FieldType scalarProduct(const VectorType& v1,const VectorType& v2) const
		{
			FieldType sum = 0;
			for (SizeType i=0;i<v1.size();i++) sum += std::conj(v1[i])*v2[i];
			return sum;
		}

		VectorType multiply(const MatrixType& A,const VectorType& v) const
		{
			VectorType y(A.rank(),0);
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

