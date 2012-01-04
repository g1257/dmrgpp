
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 0.0.1]
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

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <iostream>
#include <vector>
#include <stdexcept>
// FIXME: DON'T CALL OR INCLUDE GSL directly, use GslWrapper class
extern "C" {
#include <gsl/gsl_multimin.h>
}

namespace PsimagLite {
	
	template<typename FieldType>
	class MockVector {
	public:
		MockVector(const gsl_vector *v) : v_(v)
		{
		}
		const FieldType& operator[](size_t i) const
		{
			return v_->data[i];
		}
		size_t size() const { return v_->size; }
	private:
		const gsl_vector *v_;
	}; // class MockVector
	
	
	template<typename FunctionType>
	typename FunctionType::FieldType myFunction(const gsl_vector *v, void *params)
	{
		MockVector<typename FunctionType::FieldType> mv(v);
		FunctionType* ft = (FunctionType *)params;
		return ft->operator()(mv);
	}
	
	
	template<typename RealType,typename FunctionType>
	class Minimizer {
		
		typedef typename FunctionType::FieldType FieldType;
		typedef std::vector<FieldType> VectorType;
		typedef Minimizer<RealType,FunctionType> ThisType;

		
	public:
		
		Minimizer(FunctionType& function,size_t maxIter)
				: function_(function),
				  maxIter_(maxIter),
				  gslT_(gsl_multimin_fminimizer_nmsimplex2),
				  gslS_(gsl_multimin_fminimizer_alloc(gslT_,function_.size()))
		{
		}
		
		~Minimizer()
		{
			gsl_multimin_fminimizer_free (gslS_);
		}

		int simplex(VectorType& minVector,RealType delta=1e-3,RealType tolerance=1e-3)
		{
			gsl_vector *x;
			/* Starting point,  */
			x = gsl_vector_alloc (function_.size());
			for (size_t i=0;i<minVector.size();i++)
				gsl_vector_set (x, i, minVector[i]);

			gsl_vector *xs;
			xs = gsl_vector_alloc (function_.size());
			for (size_t i=0;i<minVector.size();i++)
				gsl_vector_set (xs, i, delta);

			gsl_multimin_function func;
			func.f= myFunction<FunctionType>;
			func.n = function_.size();
			func.params = &function_;
			gsl_multimin_fminimizer_set (gslS_, &func, x, xs);

			for (size_t iter=0;iter<maxIter_;iter++) {
				int status = gsl_multimin_fminimizer_iterate (gslS_);

				if (status) throw std::runtime_error("Minimizer::simplex(...): Error encountered\n");

				RealType size = gsl_multimin_fminimizer_size(gslS_);
				status = gsl_multimin_test_size(size, tolerance);

				if (status == GSL_SUCCESS) {
					found(minVector,gslS_->x,iter);
					gsl_vector_free (x);
					gsl_vector_free (xs);
					return iter;
				}
			}
			gsl_vector_free (x);
			gsl_vector_free (xs);
			return -1;
		}
		
	private:

		void found(VectorType& minVector,gsl_vector* x,size_t iter)
		{
			for (size_t i=0;i<minVector.size();i++)
				minVector[i] = gsl_vector_get(x,i);
		}
		
		FunctionType& function_;
		size_t maxIter_;
		const gsl_multimin_fminimizer_type *gslT_;
		gsl_multimin_fminimizer *gslS_;
		
	}; // class Minimizer
	
} // namespace PsimagLite
#endif // MINIMIZER_H

