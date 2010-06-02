// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

#include "Utils.h"
#include "ClebschGordan.h"

namespace Dmrg {
	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	class JmSubspace {
			typedef std::pair<size_t,size_t> PairType;
			typedef std::pair<PairType,PairType> TwoPairsType;
			typedef ClebschGordan<FieldType> ClebschGordanType;
		public:
			typedef std::pair<PairType,TwoPairsType> FlavorType;

			JmSubspace(const PairType& jm,size_t index,const PairType& jm1,const PairType& jm2,size_t nelectrons,int heavy=1)
			:	jm_(jm),nelectrons_(nelectrons),heavy_(heavy)
			{
				
				push(index,jm1,jm2,nelectrons);
				
			}

			static void setToProduct(const SymmetryRelatedType* symm1,const SymmetryRelatedType* symm2,
						const std::vector<size_t>& ne1,const std::vector<size_t>& ne2)
			{
				symm1_=symm1;
				symm2_=symm2;
				ne1_=&ne1;
				ne2_=&ne2;
			}

			void push (size_t index,const PairType& jm1,const PairType& jm2,size_t nelectrons)
			{
				if (nelectrons!=nelectrons_) throw std::runtime_error("JmSubspace::push(): nelectrons changed!!\n");
				indices_.push_back(index);
				setFlavorsIndex(index,jm1,jm2);
			}

			bool operator==(const std::pair<PairType,size_t>& nejm) const 
			{
				std::pair<PairType,size_t> nejmStored=std::pair<PairType,size_t>(jm_,nelectrons_);
				if (nejm==nejmStored) return true;
				return false;
			}

			size_t heavy() const { return heavy_; }

			//! This function is performance critical
			size_t createFactors(SparseMatrixType& factors,size_t offset)
			{
				flavors_.clear();
				std::vector<size_t> perm(indices_.size());
				utils::sort(flavorIndices_,perm);
				size_t flavorSaved=flavorIndices_[0];
				flavors_.push_back(flavorIndices_[0]);
				size_t counter=0;
				for (size_t k=0;k<indices_.size();k++) {
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

			size_t numberOfFlavors() const 
			{
				return flavors_.size(); 
			}

			size_t getFlavor(size_t i) const
			{
				return flavors_[i];
			}

			size_t getNe() const { return nelectrons_; }

			size_t numberOfIndices() const { return indices_.size(); }

// 			void print(std::ostream &os,size_t ns=0,bool printIndices=true) const
// 			{
// 				os<<"===============================\n";
// 				os<<"jm="<<jm_<<"\n";
// 				os<<"electrons="<<nelectrons_<<"\n";
// 				os<<"number_of_indices="<<indices_.size()<<"\n";
// 				if (printIndices) { 	
// 					for (size_t i=0;i<indices_.size();i++) {
// 						if (ns>0) {
// 							size_t y = size_t(indices_[i]/ns);
// 							size_t x = indices_[i] % ns;
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

			size_t qn() const
			{
				return SymmetryRelatedType::neJmToIndex(nelectrons_,jm_); 
			}

			size_t findFreeRow(size_t counter,size_t total)
			{
				for (size_t i=counter;i<total;i++) {
					int x=utils::isInVector(indices_,i);
					if (x<0) return i;
				}
				throw std::runtime_error("findfreerow: no free rows\n");
			}	

			static size_t flavor(size_t f1,size_t f2,size_t ne1,size_t ne2,size_t j1,size_t j2)
			{
				size_t x = f1 + f2*symm1_->flavorsMax();
				size_t y = ne1 + ne2*symm1_->electronsMax();
				size_t z = j1 + j2*symm1_->jMax();
				size_t xmax = symm1_->flavorsMax()*symm2_->flavorsMax();
				size_t ret = x + y*xmax;
				size_t ymax =symm1_->electronsMax()*symm2_->electronsMax();
				return ret + z*xmax*ymax;
			}

		private:
			PairType jm_;
			size_t nelectrons_;
			int heavy_;
			std::vector<size_t> indices_;
			std::vector<FieldType> cg_,values_;
			std::vector<size_t> flavors_,flavorIndices_;
			ClebschGordanType cgObject_;
			
			static const std::vector<size_t>* ne1_;
			static const std::vector<size_t>* ne2_;
			static const SymmetryRelatedType* symm1_;
			static const SymmetryRelatedType* symm2_;
			static size_t flavorsMaxPerSite__,numberOfSites__,electronsMaxPerSite__,jMaxPerSite__;

			void setFlavorsIndex(size_t i,const PairType& jm1,const PairType& jm2)
			{
				size_t alpha=0,beta=0;
				utils::getCoordinates(alpha,beta,i,symm1_->size());
				
				size_t ne1 = (*ne1_)[alpha];
				size_t ne2 = (*ne2_)[beta];
				size_t flavor1 = symm1_->getFlavor(alpha);
				size_t flavor2 = symm2_->getFlavor(beta);
				PairType flavorPair = PairType(flavor1,flavor2);
				PairType nePair = PairType(ne1,ne2);
				PairType j1j2 = PairType(jm1.first,jm2.first);
				FlavorType flavorPair2 = FlavorType(flavorPair,TwoPairsType(nePair,j1j2));
				flavorIndices_.push_back(calculateFlavor(flavorPair2));

				FieldType value = 0;
				if (heavy_) value=cgObject_(jm_,jm1,jm2);
				values_.push_back(value);
			}

			static size_t calculateFlavor(const FlavorType& g)
			{
				// order is : f1, f2, ne1, ne2, j1, j2
				size_t f1 = (g.first).first;
				size_t f2 = (g.first).second;
				TwoPairsType tp=g.second;
				size_t ne1 = (tp.first).first;
				size_t ne2 = (tp.first).second;
				size_t j1 = (tp.second).first;
				size_t j2 = (tp.second).second;
				return flavor(f1,f2,ne1,ne2,j1,j2);
			}
	}; // class JmSubspace
	
	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const SymmetryRelatedType* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::symm1_=0;

	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const SymmetryRelatedType* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::symm2_=0;

	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const std::vector<size_t>* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::ne1_=0;

	template<typename FieldType,typename SparseMatrixType,typename SymmetryRelatedType>
	const std::vector<size_t>* JmSubspace<FieldType,SparseMatrixType,SymmetryRelatedType>::ne2_=0;

} // namespace Dmrg

/*@}*/
#endif
