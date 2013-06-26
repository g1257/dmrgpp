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
	
	template<typename BasisWithOperatorsType_,typename SuperBlockType>
	class LeftRightSuper {
		public:
			typedef typename SuperBlockType::RealType RealType;
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
			typedef  LeftRightSuper<
					BasisWithOperatorsType_,SuperBlockType> ThisType;

			enum {GROW_TO_THE_RIGHT = BasisWithOperatorsType::GROW_RIGHT,
				GROW_TO_THE_LEFT= BasisWithOperatorsType::GROW_LEFT};

			template<typename IoInputter>
			LeftRightSuper(IoInputter& io)
			: progress_("LeftRightSuper"),
			  left_(0),right_(0),super_(0),refCounter_(0)
			{
				// watch out: same order as save here:
				super_ = new SuperBlockType(io,"");
				left_ = new BasisWithOperatorsType(io,"");
				right_ = new BasisWithOperatorsType(io,"");
			}

			LeftRightSuper(
					const PsimagLite::String& slabel,
					const PsimagLite::String& elabel,
					const PsimagLite::String& selabel)
			: progress_("LeftRightSuper"),
			  left_(0),right_(0),super_(0),refCounter_(0)
			{
				left_ = new BasisWithOperatorsType(slabel);
				right_ = new BasisWithOperatorsType(elabel);
				super_ = new SuperBlockType(selabel);
			}

			~LeftRightSuper()
			{
				if (refCounter_>0) {
					refCounter_--;
					return;
				}
				delete left_;
				delete right_;
				delete super_;
			}

			LeftRightSuper(
					BasisWithOperatorsType& left,
					BasisWithOperatorsType& right,
					SuperBlockType& super)
			: progress_("LeftRightSuper"),
			  left_(&left),right_(&right),super_(&super),refCounter_(1)
			{
			}

			LeftRightSuper(const ThisType& rls)
			: progress_("LeftRightSuper"),refCounter_(1)
			{
				left_=rls.left_;
				right_=rls.right_;
				super_=rls.super_;
			}

			// deep copy
			ThisType& operator=(const ThisType& lrs)
			{
				deepCopy(lrs);
				return *this;
			}

			template<typename SomeModelType>
			void growLeftBlock(const SomeModelType& model,
					   BasisWithOperatorsType &pS,
					   BlockType const &X,
					   RealType time)
			{
				grow(*left_,model,pS,X,GROW_TO_THE_RIGHT,time);
			}

			template<typename SomeModelType>
			void growRightBlock(const SomeModelType& model,
					    BasisWithOperatorsType &pE,
					    BlockType const &X,
					    RealType time)
			{
				grow(*right_,model,pE,X,GROW_TO_THE_LEFT,time);
			}

			void printSizes(const PsimagLite::String& label,std::ostream& os) const
			{
				PsimagLite::OstringStream msg;
				msg<<label<<": left-block basis="<<left_->size();
				msg<<", right-block basis="<<right_->size();
				msg<<" sites="<<left_->block().size()<<"+";
				msg<<right_->block().size();
				progress_.printline(msg,os);
			}

			SizeType sites() const
			{
				return left_->block().size() + right_->block().size();
			}

			/*!PTEX_LABEL{setToProductLrs} */
			void setToProduct(SizeType quantumSector)
			{
				super_->setToProduct(*left_,*right_,quantumSector);
			}

			template<typename IoOutputType>
			void save(IoOutputType& io) const
			{
				super_->save(io);
				left_->save(io);
				right_->save(io);
			}

			const BasisWithOperatorsType& left()  const { return *left_; }

			const BasisWithOperatorsType& right() const { return *right_; }

			BasisWithOperatorsType& leftNonConst()  { return *left_; }

			BasisWithOperatorsType& rightNonConst() { return *right_; }


			const SuperBlockType& super() const { return *super_; }

			void left(const BasisWithOperatorsType& left)
			{
				if (refCounter_>0) throw PsimagLite::RuntimeError
						("LeftRightSuper::left(...): not the owner\n");
				*left_=left; // deep copy
			}

			void right(const BasisWithOperatorsType& right)
			{
				if (refCounter_>0) throw PsimagLite::RuntimeError
						("LeftRightSuper::right(...): not the owner\n");
				*right_=right; // deep copy
			}


			template<typename IoInputType>
			void load(IoInputType& io)
			{
				super_->load(io);
				left_->load(io);
				right_->load(io);
			}

		private:
			LeftRightSuper(ThisType& rls);

			void deepCopy(const ThisType& rls)
			{
				*left_=*rls.left_;
				*right_=*rls.right_;
				*super_=*rls.super_;
				if (refCounter_>0) refCounter_--;
			}

			//! add block X to basis pS and put the result in left_:
			/**
			The function grow is located in !PTEX\\_REF{HERE}.
			Local operators are set for the basis in question with a call to 
			\\cppClass{BasisWithOperators}'s member function \\cppFunction{setOperators()}.  
			When adding sites to the system or environment the program does a 
			full outer product, i.e., it increases the size of all local operators. 
			This is performed by the call to \\verb!setToProduct(pSprime,pS,Xbasis,dir,option)!
			in the grow function, which actually calls \\verb!pSprime.setToProduct(pS,xBasis,dir)!
			This function also recalculates the Hamiltonian in the outer product 
			of (i) the previous system basis $pS$, and (ii) the basis $Xbasis$ 
			corresponding to the site(s) that is (are) being added.
			To do this, the Hamiltonian connection between the two parts 
			needs to be calculated and added, and this is done in the call to 
			\\cppFunction{addHamiltonianConnection}, see !PTEX\\_REF{295}. 
			Finally, the resulting dmrgBasis object for the outer product, 
			pSprime, is set to contain this full Hamiltonian with the call
			to  \\cppFunction{pSprime.setHamiltonian(matrix)}. 
			*/
			template<typename SomeModelType>
			void grow(BasisWithOperatorsType& leftOrRight,
			          const SomeModelType& model,
				  BasisWithOperatorsType &pS,
				  const BlockType& X,
			          SizeType dir,
				  RealType time)
			{
				SparseMatrixType hmatrix;
				BasisDataType q;
				typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
				model.setNaturalBasis(creationMatrix,hmatrix,q,X,time);
				BasisWithOperatorsType Xbasis("Xbasis");

				Xbasis.setVarious(X,hmatrix,q,creationMatrix);
				leftOrRight.setToProduct(pS,Xbasis,dir);

				SparseMatrixType matrix=leftOrRight.hamiltonian();

				ThisType* lrs;
				BasisType* leftOrRightL =  &leftOrRight;
				if (dir==GROW_TO_THE_RIGHT) {
					lrs = new ThisType(pS,Xbasis,*leftOrRightL);
				} else {
					lrs = new  ThisType(Xbasis,pS,*leftOrRightL);
				}
				//!PTEX_LABEL{295}
				model.addHamiltonianConnection(matrix,*lrs);
				delete lrs;
				leftOrRight.setHamiltonian(matrix);
			}

			ProgressIndicatorType progress_;
			BasisWithOperatorsType* left_;
			BasisWithOperatorsType* right_;
			SuperBlockType* super_;
			SizeType refCounter_;
			
	}; // class LeftRightSuper

} // namespace Dmrg 

/*@}*/
#endif // LEFT_RIGHT_SUPER_H
