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

/*! \file LeftRightSuper.h
 *
 *  A class that contains the left block or system, the
 *  right block or environ, and the superblock
 */
#ifndef LEFT_RIGHT_SUPER_H
#define LEFT_RIGHT_SUPER_H

#include "ProgressIndicator.h"

namespace Dmrg {
	
	template<typename BasisWithOperatorsType_>
	class LeftRightSuper {
		public:
			typedef BasisWithOperatorsType_ BasisWithOperatorsType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef typename BasisType::BasisDataType BasisDataType;
			typedef typename BasisWithOperatorsType::SparseMatrixType
					SparseMatrixType;
			typedef typename BasisWithOperatorsType::OperatorsType
					OperatorsType;
			typedef typename OperatorsType::OperatorType
					OperatorType;
			typedef typename BasisType::BlockType BlockType;
			typedef PsimagLite::ProgressIndicator ProgressIndicatorType;

			enum {GROW_TO_THE_RIGHT = BasisWithOperatorsType::GROW_RIGHT,
				GROW_TO_THE_LEFT= BasisWithOperatorsType::GROW_LEFT};

			template<typename IoInputter>
			LeftRightSuper(IoInputter& io,bool bogus=false)
			: progress_("LeftRightSuper",0),
			  left_(io,"",bogus),right_(io,"",bogus),super_(io,"",bogus)
			{
			}

			LeftRightSuper(
					const std::string& slabel,
					const std::string& elabel,
					const std::string& selabel)
			: progress_("LeftRightSuper",0),
			  left_(slabel),right_(elabel),super_(selabel)
			  {}
			
			template<typename SomeModelType>
			void growLeftBlock(
					const SomeModelType& model,
					const BasisWithOperatorsType &pS,
					BlockType const &X)
			{
				grow(model,pS,X,GROW_TO_THE_RIGHT);
			}

			template<typename SomeModelType>
			void growRightBlock(
					const SomeModelType& model,
					const BasisWithOperatorsType &pE,
					BlockType const &X)
			{
				grow(right_,model,pE,X,GROW_TO_THE_LEFT);
			}

			void printSizes(const std::string& label,std::ostream& os) const
			{
				std::ostringstream msg;
				msg<<label<<": left-block basis="<<left_.size();
				msg<<", right-block basis="<<right_.size();
				msg<<" sites="<<left_.block().size()<<"+";
				msg<<right_.block().size();
				progress_.printline(msg,os);
			}

			size_t sites() const
			{
				return left_.block().size() + right_.block().size();
			}

			void setToProduct(size_t quantumSector)
			{
				super_.setToProduct(left_,right_,quantumSector);
			}

			const BasisType& left() const { return left_; }

			const BasisType& right() const { return right_; }

			const BasisType& super() const { return super_; }

		private:

			//! add block X to basis pS and put the result in left_:
			template<typename SomeModelType>
			void grow(
					BasisWithOperatorsType& leftOrRight,
					const SomeModelType& model,
					const BasisWithOperatorsType &pS,
					BlockType const &X,
					size_t dir)
			{
				SparseMatrixType hmatrix;
				BasisDataType q;
				std::vector<OperatorType> creationMatrix;
				model.setNaturalBasis(creationMatrix,hmatrix,q,X);
				BasisWithOperatorsType Xbasis("Xbasis");

				Xbasis.setVarious(X,hmatrix,q,creationMatrix);
				leftOrRight.setToProduct(pS,Xbasis,dir);

				SparseMatrixType matrix=leftOrRight.hamiltonian();

				if (dir==GROW_TO_THE_RIGHT)
					model.addHamiltonianConnection(
							matrix,leftOrRight,pS,Xbasis,model.orbitals());
				else
					model.addHamiltonianConnection(
							matrix,leftOrRight,Xbasis,pS,model.orbitals());

				leftOrRight.setHamiltonian(matrix);
			}

			ProgressIndicatorType progress_;
			BasisWithOperatorsType left_;
			BasisWithOperatorsType right_;
			BasisType super_;
			
	}; // class LeftRightSuper

} // namespace Dmrg 

/*@}*/
#endif // LEFT_RIGHT_SUPER_H
