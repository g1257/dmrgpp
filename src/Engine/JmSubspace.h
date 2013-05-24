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

/*! \file JmSubspace.h
 *
 *  Each object of this class contains a subspace of the outer product.
 *  States (a,b) in this subspace give rise to a state c in the outer product
 *  with given quantum numbers (j,m,q) [See paper for more info]
 *
 */
#ifndef JM_SUBSPACE_H
#define JM_SUBSPACE_H

#include "Sort.h" // in PsimagLite
#include "PackIndices.h" // in PsimagLite
#include "Su2SymmetryGlobals.h"

namespace Dmrg {
	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	class JmSubspace {
			typedef std::pair<SizeType,SizeType> PairType;
			typedef std::pair<PairType,PairType> TwoPairsType;
			typedef Su2SymmetryGlobals<FieldType> Su2SymmetryGlobalsType;
			typedef typename Su2SymmetryGlobalsType::ClebschGordanType ClebschGordanType;
			typedef PsimagLite::PackIndices PackIndicesType;

		public:
			typedef std::pair<PairType,TwoPairsType> FlavorType;

			JmSubspace(const PairType& jm,SizeType index,const PairType& jm1,const PairType& jm2,SizeType nelectrons,int heavy=1)
			:	jm_(jm),nelectrons_(nelectrons),heavy_(heavy),cgObject_(&(Su2SymmetryGlobalsType::clebschGordanObject))
			{
				
				push(index,jm1,jm2,nelectrons);
				
			}

			static void setToProduct(const SymmetryRelatedType* symm1,const SymmetryRelatedType* symm2,
						const PsimagLite::Vector<SizeType>::Type& ne1,const PsimagLite::Vector<SizeType>::Type& ne2)
			{
				symm1_=symm1;
				symm2_=symm2;
				ne1_=&ne1;
				ne2_=&ne2;
			}

			void push (SizeType index,const PairType& jm1,const PairType& jm2,SizeType nelectrons)
			{
				if (nelectrons!=nelectrons_) throw PsimagLite::RuntimeError("JmSubspace::push(): nelectrons changed!!\n");
				indices_.push_back(index);
				setFlavorsIndex(index,jm1,jm2);
			}

			bool operator==(const std::pair<PairType,SizeType>& nejm) const 
			{
				std::pair<PairType,SizeType> nejmStored=std::pair<PairType,SizeType>(jm_,nelectrons_);
				if (nejm==nejmStored) return true;
				return false;
			}

			SizeType heavy() const { return heavy_; }

			//! This function is performance critical
			SizeType createFactors(SparseMatrixType& factors,SizeType offset)
			{
				flavors_.clear();
				PsimagLite::Vector<SizeType>::Type perm(indices_.size());
				PsimagLite::Sort<PsimagLite::Vector<SizeType>::Type > sort;
				sort.sort(flavorIndices_,perm);
				SizeType flavorSaved=flavorIndices_[0];
				flavors_.push_back(flavorIndices_[0]);
				SizeType counter=0;
				for (SizeType k=0;k<indices_.size();k++) {
					if (flavorIndices_[k]!=flavorSaved) {
						flavors_.push_back(flavorIndices_[k]);
						counter++;
						flavorSaved = flavorIndices_[k];
					}
					if (heavy_) factors.set(indices_[perm[k]],offset + counter, values_[perm[k]]);
				}
				return flavors_.size();
			}

			PairType getJmValue() const 
			{
				return jm_;	
			}

			SizeType numberOfFlavors() const 
			{
				return flavors_.size(); 
			}

			SizeType getFlavor(SizeType i) const
			{
				return flavors_[i];
			}

			SizeType getNe() const { return nelectrons_; }

			SizeType numberOfIndices() const { return indices_.size(); }

// 			void print(std::ostream &os,SizeType ns=0,bool printIndices=true) const
// 			{
// 				os<<"===============================\n";
// 				os<<"jm="<<jm_<<"\n";
// 				os<<"electrons="<<nelectrons_<<"\n";
// 				os<<"number_of_indices="<<indices_.size()<<"\n";
// 				if (printIndices) { 	
// 					for (SizeType i=0;i<indices_.size();i++) {
// 						if (ns>0) {
// 							SizeType y = SizeType(indices_[i]/ns);
// 							SizeType x = indices_[i] % ns;
// 							os<<"indices["<<i<<"]="<<x<<" "<<y<<"\n";
// 						} else {
// 							os<<"indices["<<i<<"]="<<indices_[i]<<"\n";
// 						}
// 					}
// 				}
// 				os<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5\n";
// 			}

			void clear()
			{
				indices_.clear();
				cg_.clear();
				flavors_.clear();
				flavorIndices_.clear();
				values_.clear();
				
			}

			SizeType qn() const
			{
				return SymmetryRelatedType::neJmToIndex(nelectrons_,jm_); 
			}

			SizeType findFreeRow(SizeType counter,SizeType total)
			{
				for (SizeType i=counter;i<total;i++) {
					int x=PsimagLite::isInVector(indices_,i);
					if (x<0) return i;
				}
				throw PsimagLite::RuntimeError("findfreerow: no free rows\n");
			}	

			static SizeType flavor(SizeType f1,SizeType f2,SizeType ne1,SizeType ne2,SizeType j1,SizeType j2)
			{
				SizeType x = f1 + f2*symm1_->flavorsMax();
				SizeType y = ne1 + ne2*symm1_->electronsMax();
				SizeType z = j1 + j2*symm1_->jMax();
				SizeType xmax = symm1_->flavorsMax()*symm2_->flavorsMax();
				SizeType ret = x + y*xmax;
				SizeType ymax =symm1_->electronsMax()*symm2_->electronsMax();
				return ret + z*xmax*ymax;
			}

		private:
			PairType jm_;
			SizeType nelectrons_;
			int heavy_;
			ClebschGordanType* cgObject_;
			PsimagLite::Vector<SizeType>::Type indices_;
			typename PsimagLite::Vector<FieldType>::Type cg_,values_;
			typename PsimagLite::Vector<SizeType>::Type flavors_,flavorIndices_;
			
			static const PsimagLite::Vector<SizeType>::Type* ne1_;
			static const PsimagLite::Vector<SizeType>::Type* ne2_;
			static const SymmetryRelatedType* symm1_;
			static const SymmetryRelatedType* symm2_;
			static SizeType flavorsMaxPerSite__,numberOfSites__,electronsMaxPerSite__,jMaxPerSite__;

			void setFlavorsIndex(SizeType i,const PairType& jm1,const PairType& jm2)
			{
				SizeType alpha=0,beta=0;
				PackIndicesType pack(symm1_->size());
				pack.unpack(alpha,beta,i);
				
				SizeType ne1 = (*ne1_)[alpha];
				SizeType ne2 = (*ne2_)[beta];
				SizeType flavor1 = symm1_->getFlavor(alpha);
				SizeType flavor2 = symm2_->getFlavor(beta);
				PairType flavorPair = PairType(flavor1,flavor2);
				PairType nePair = PairType(ne1,ne2);
				PairType j1j2 = PairType(jm1.first,jm2.first);
				FlavorType flavorPair2 = FlavorType(flavorPair,TwoPairsType(nePair,j1j2));
				flavorIndices_.push_back(calculateFlavor(flavorPair2));

				FieldType value = 0;
				if (heavy_) value=cgObject_->operator()(jm_,jm1,jm2);
				values_.push_back(value);
			}

			static SizeType calculateFlavor(const FlavorType& g)
			{
				// order is : f1, f2, ne1, ne2, j1, j2
				SizeType f1 = (g.first).first;
				SizeType f2 = (g.first).second;
				TwoPairsType tp=g.second;
				SizeType ne1 = (tp.first).first;
				SizeType ne2 = (tp.first).second;
				SizeType j1 = (tp.second).first;
				SizeType j2 = (tp.second).second;
				return flavor(f1,f2,ne1,ne2,j1,j2);
			}
	}; // class JmSubspace
	
	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const SymmetryRelatedType* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::symm1_=0;

	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const SymmetryRelatedType* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::symm2_=0;

	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const PsimagLite::Vector<SizeType>::Type* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::ne1_=0;

	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const PsimagLite::Vector<SizeType>::Type* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::ne2_=0;

} // namespace Dmrg

/*@}*/
#endif
