// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file InternalProductStored.h
 *
 *  A class to encapsulate the product x+=Hy, where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef InternalProductStored_HEADER_H
#define InternalProductStored_HEADER_H

#include <vector>
#include "ProgressIndicator.h"

namespace Dmrg {
	template<typename T,typename ModelType>
	class InternalProductStored {
	public:	
		typedef T HamiltonianElementType;
		typedef HamiltonianElementType value_type; 
		typedef typename ModelType::ModelHelperType ModelHelperType;
		typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
		//typedef typename SparseMatrixType::value_type SparseElementType;

		InternalProductStored(ModelType const *model,
				      ModelHelperType const *modelHelper,
				      const ReflectionSymmetryType* rs=0)
		: matrixStored_(2),pointer_(0),progress_("InternalProductStored",0)
		{
			model_ = model;
			modelHelper_=modelHelper;
			std::ostringstream msg;

			if (!rs) {
				matrixStored_[0].clear();
				model->fullHamiltonian(matrixStored_[0],*modelHelper);
				assert(isHermitian(matrixStored_[0],true));
				msg<<"fullHamiltonian has rank="<<matrixStored_[0].rank()<<" nonzeros="<<matrixStored_[0].nonZero()<<"\n";
				progress_.printline(msg,std::cout);
				return;
			}
			SparseMatrixType matrix2;
			model->fullHamiltonian(matrix2,*modelHelper);
			rs->transform(matrixStored_[0],matrixStored_[1],matrix2);
//			msg<<"fullHamiltonian with rank="<<matrixStored_[0].rank()<<" has reflection\n";
//			progress_.printline(msg,std::cout);
		}

		size_t rank() const { return matrixStored_[pointer_].rank(); }

		template<typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType &x, SomeVectorType const &y) const
		{
			 matrixStored_[pointer_].matrixVectorProduct(x,y);
		}

		HamiltonianElementType operator()(size_t i,size_t j) const
		{
			return matrixStored_[pointer_](i,j);
		}

		size_t reflectionSector() const { return pointer_; }

		void reflectionSector(size_t p) { pointer_=p; }

	private:
		ModelType const *model_;
		ModelHelperType const *modelHelper_;
		std::vector<SparseMatrixType> matrixStored_;
		size_t pointer_;
		PsimagLite::ProgressIndicator progress_;
	}; // class InternalProductStored
} // namespace Dmrg

/*@}*/
#endif
