// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009-2011, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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
// END LICENSE BLOCK
/** \ingroup PsimagLite */
/*@{*/

/*! \file KernelPolynomial.h
 *
 *  Given a series of moments, it recreates the function
 *  see http://dx.doi.org/10.1103/RevModPhys.78.275
 *
 */

#ifndef KERNEL_POLYNOMIAL_H
#define KERNEL_POLYNOMIAL_H
// #include "ProgressIndicator.h"
#include "Vector.h"
#include "TypeToString.h"

namespace PsimagLite {

	//! MatrixType must have the following interface:
	//! 	RealType type to indicate the matrix type
	//! 	rank() member function to indicate the rank of the matrix
	//! 	matrixVectorProduct(std::vector<RealType>& x,const std::vector<RealType>& const y) 
	//!    	   member function that implements the operation x += Hy

	template<typename RealType>
	class KernelPolynomial {

	public:

		KernelPolynomial(std::vector<RealType>& moments,
		                 const RealType& a,
		                 const RealType& b)
		: progress_("KernelPolynomial",0),
		  moments_(moments),
		  oneOverA_(1.0/a),
		  b_(b)
		{}

		void jackson(std::vector<RealType>& res,
		             const RealType& start,
		             const RealType& end,
		             const RealType& step) const
		{
			std::vector<RealType> gn(moments_.size());
			initKernelJackson(gn,moments_.size());
			internalMain(res,start,end,gn);
		}
		
	private:
		
		void internalMain(std::vector<RealType>& res,
		                  const RealType& start,
		                  const RealType& end,
		                  const RealType& step,
		                  const std::vector<RealType>& gn)
		{
			std::vector<RealType> gnmun(gn.size());
			computeGnMuN(gnmun,gn);
			
			res.clear();

			RealType omega = start;
			while(omega<end) {
				RealType x = (omega-b_)*oneOverA_;
				RealType den = sqrt(1.0 - x*x);
				res.push_back(calcF(x,gnmun)/den);
				omega += step;
			}
		}
		
		RealType calcF(const RealType& x,const std::vector<RealType>& gnmn) const
		{
			RealType sum = 0.5*gnmn[0];
			for (size_t i=1;i<gnmn.size();i++) sum += gnmn[i]*chebyshev_(i,x);
			return 2.0*sum;
		}
		
		void computeGnMuN( std::vector<RealType>& gnmn,std::vector<RealType>& gn) const
		{
			for (size_t i=0;i<gnmn.size();i++) gnmn[i] = moments_[i] * gn[i];
		}

		void initKernelJackson(std::vector<RealType>& gn) const
		{
			size_t nPlus1 = gn.size()+1;
			RealType cot1 = 1.0/tan(M_PI/nPlus1);
			for (size_t i=0;i<gn.size();i++) {
				gn[i] = (nPlus1-i)*cos(M_PI*i/nPlus1)+sin(M_PI*i/nPlus1)*cot1;
				gn[i] /= nPlus1;
			}
		}

		std::vector<RealType> moments_;
		RealType oneOverA_,b_;
// 		ProgressIndicator progress_;
	}; // class KernelPolynomial
} // namespace PsimagLite

/*@}*/
#endif // KERNEL_POLYNOMIAL_H
