/*
Copyright (c) 2009-2011, 2013, UT-Battelle, LLC
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

/*! \file CorrectionVectorFunction.h
 *
 *  This is an implementation of PRB 60, 335, Eq. (24)
 *
 */
#ifndef CORRECTION_V_FUNCTION_H
#define CORRECTION_V_FUNCTION_H
#include "ConjugateGradient.h"

namespace Dmrg {
template<typename MatrixType,typename InfoType>
class	CorrectionVectorFunction {

	typedef typename MatrixType::value_type FieldType;
	typedef typename PsimagLite::Vector<FieldType>::Type VectorType;
	typedef typename PsimagLite::Real<FieldType>::Type RealType;

	class InternalMatrix {

	public:

		typedef FieldType value_type ;
		InternalMatrix(const MatrixType& m,const InfoType& info,RealType E0)
		    : m_(m),info_(info),E0_(E0)
		{
			if (info_.omega().first != PsimagLite::FREQ_REAL)
				throw PsimagLite::RuntimeError("Matsubara only with KRYLOV\n");
		}

		SizeType rows() const { return m_.rows(); }

		void matrixVectorProduct(VectorType& x,const VectorType& y) const
		{
			RealType eta = info_.eta();
			RealType omegaMinusE0 = info_.omega().second + E0_;
			VectorType xTmp(x.size(),0);
			m_.matrixVectorProduct(xTmp,y); // xTmp = Hy
			VectorType x2(x.size(),0);
			m_.matrixVectorProduct(x2,xTmp); // x2 = H^2 y
			const RealType f1 = (-2.0);
			// this needs fixing
			// preferred:
			// x <= x2 + f1*omegaMinusE0*xTmp + (omegaMinusE0*omegaMinusE0 + eta*eta)*y;
			// equivalent
			for (SizeType i = 0; i < x.size(); ++i)
				x[i] = x2[i] + f1*omegaMinusE0*xTmp[i] +
				        (omegaMinusE0*omegaMinusE0 + eta*eta)*y[i];

			x /= (-eta);
		}

	private:

		const MatrixType& m_;
		const InfoType& info_;
		RealType E0_;
	};

	typedef ConjugateGradient<InternalMatrix> ConjugateGradientType;

public:

	CorrectionVectorFunction(const MatrixType& m,const InfoType& info,RealType E0)
	    : im_(m,info,E0),cg_(info.cgSteps(),info.cgEps())
	{}

	void getXi(VectorType& result,const VectorType& sv) const
	{
		VectorType x0(result.size(),0.0);

		result = x0; // initial ansatz
		cg_(result,im_,sv);
	}

private:

	InternalMatrix im_;
	ConjugateGradientType cg_;
}; // class CorrectionVectorFunction
} // namespace Dmrg

/*@}*/
#endif // CORRECTION_V_FUNCTION_H

