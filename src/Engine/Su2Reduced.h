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
#ifndef SU2_REDUCED_HEADER_H
#define SU2_REDUCED_HEADER_H

#include "Map.h"
#include "Su2SymmetryGlobals.h"

/** \ingroup DMRG */
/*@{*/

/*! \file Su2Reduced.h
 *
 *  Using WignerEckart Theorem to speed up SU(2) algorithm
 *
 */
namespace Dmrg {
	template<typename LeftRightSuperType>
	class Su2Reduced {
		typedef typename LeftRightSuperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename LeftRightSuperType::BasisWithOperatorsType
		                                     BasisWithOperatorsType;
		typedef typename LeftRightSuperType::BasisType BasisType;
		typedef typename BasisType::RealType RealType;
		typedef std::pair<SizeType,SizeType> PairType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef Su2SymmetryGlobals<RealType> Su2SymmetryGlobalsType;
		typedef typename Su2SymmetryGlobalsType::ClebschGordanType ClebschGordanType;

		static const SizeType System=0,Environ=1;
		//static const BasisType* basis1Ptr_;

	public:
		Su2Reduced(int m,
		           const LeftRightSuperType& lrs,
		           bool = false)
		:m_(m),
		 lrs_(lrs),
		 cgObject_(Su2SymmetryGlobalsType::clebschGordanObject)
		{
			typename PsimagLite::Vector<PairType>::Type jsEffective;
			typename PsimagLite::Vector<SizeType>::Type jvalues;
			// find all possible j values
			SizeType counter=0;
			for (SizeType i=0;i<lrs.left().numberOfOperators();i++) {
				SizeType j = lrs.left().getReducedOperatorByIndex(i).jm.first;
				int x = PsimagLite::isInVector(jvalues,j);
				if (x<0) jvalues.push_back(j);
				counter += (j+1);
			}

			// build all lfactors
			lfactor_.resize(counter);
			counter=0;
			for (SizeType i=0;i<jvalues.size();i++) {
				for (SizeType m=0;m<=jvalues[i];m++) {
					buildAdditional(lfactor_[counter],jvalues[i],m,jvalues[i]-m,jsEffective);
					counter++;
				}
			}

			calcHamFactor(lfactorHamiltonian_);
			calcEffectiveStates(jsEffective);

			//basis1Ptr_ = &basis1;
			createReducedHamiltonian(hamiltonian2_,lrs_.left());

			createReducedHamiltonian(hamiltonian3_,lrs_.right());
		}

		const PairType& reducedEffective(SizeType i) const
		{
			return reducedEffective_[i];
		}

		SizeType reducedEffectiveSize() const
		{
			return 	reducedEffective_.size();
		}

		SizeType flavorMapping(SizeType i1prime,SizeType i2prime) const
		{
			return flavorsOldInverse_[reducedInverse_(i1prime,i2prime)];
		}

		SizeType flavorMapping(SizeType i) const
		{
			return flavorsOldInverse_[i];
		}

		SparseElementType reducedFactor(SizeType angularMomentum,SizeType category,bool flip,SizeType lf1,SizeType lf2) const
		{
			SizeType category2 = category;
			if (flip) category2 = angularMomentum-category;
			SparseElementType lfactor=lfactor_[category2](lf1,lf2);
			return lfactor;
		}

		SparseElementType reducedHamiltonianFactor(SizeType j1,SizeType j2) const
		{
			return lfactorHamiltonian_(j1,j2);
		}

		const SparseMatrixType& hamiltonianLeft() const
		{
			return hamiltonian2_;
		}

		const SparseMatrixType& hamiltonianRight() const
		{
			return hamiltonian3_;
		}

	private:
		void createReducedHamiltonian(SparseMatrixType& hamReduced,const BasisWithOperatorsType& basis)
		{
			SizeType xx = basis.numberOfOperators();
			hamReduced.resize(xx,xx);
			typename PsimagLite::Vector<SizeType>::Type basisrinverse(basis.size());
			for (SizeType i=0;i<basis.size();i++) {
				SizeType f = basis.getFlavor(i);
				SizeType j = basis.jmValue(i).first;
				basisrinverse[i]=findJf(basis,j,f);
			}

			SizeType angularMomentum =0;
			typename PsimagLite::Vector<const OperatorType*>::Type opSrc(angularMomentum+1);

			OperatorType myOp;
			myOp.data = basis.hamiltonian();
			myOp.fermionSign=1;
			myOp.jm=typename OperatorType::PairType(0,0);
			myOp.angularFactor = 1.0;
			opSrc[0]=&myOp;
			createReducedOperator(hamReduced,opSrc,basis,basisrinverse,basis.reducedSize(),0);
		}

		SizeType findJf(const BasisWithOperatorsType& basis,SizeType j,SizeType f)
		{
			for (SizeType i=0;i<basis.reducedSize();i++) {
				PairType jm=basis.jmValue(basis.reducedIndex(i));
				if (jm.first!=j) continue;
				if (basis.getFlavor(basis.reducedIndex(i))!=f) continue;
				return i;
			}
			return 0; // bogus
		}

		void createReducedConj(SizeType k1,
		                       SparseMatrixType& opDest,
		                       const SparseMatrixType& opSrc,
		                       const BasisWithOperatorsType& basis,
		                       SizeType counter)
		{
			SizeType n=opSrc.row();
			opDest=transposeConjugate(opSrc);
			for (SizeType i=0;i<n;i++) {
				PairType jm = basis.jmValue(basis.reducedIndex(i));
				for (int k=opDest.getRowPtr(i);k<opDest.getRowPtr(i+1);k++) {
					SizeType j=opDest.getCol(k);
					PairType jmPrime = basis.jmValue(basis.reducedIndex(j));

					SparseElementType factor=SparseElementType(jm.first+1)/SparseElementType(jmPrime.first+1);
					factor = sqrt(factor);
					int x = k1+jm.first-jmPrime.first;
					x = int(x/2);
					if (x%2!=0) factor = -1.0*factor;
					SparseElementType val=opDest.getValue(k)*factor;
					opDest.setValues(k,val);
				}
			}
		}

		void createReducedOperator(SparseMatrixType& opDest,
		                           const typename PsimagLite::Vector<const OperatorType*>::Type& opSrc,
		                           const BasisWithOperatorsType& basis,
		                           const typename PsimagLite::Vector<SizeType>::Type& basisrInverse,
		                           SizeType n,
		                           SizeType)
		{
			PsimagLite::Matrix<SparseElementType> opDest1(n,n);
			for (SizeType i=0;i<opSrc.size();i++)
				createReducedOperator(opDest1,*opSrc[i],basis,basisrInverse); //,PairType(k,mu1),mysign);
			fullMatrixToCrsMatrix(opDest,opDest1);
		}

		void createReducedOperator(PsimagLite::Matrix<SparseElementType>& opDest1,
		                           const OperatorType& opSrc,
		                           const BasisWithOperatorsType& basis,
		                           const typename PsimagLite::Vector<SizeType>::Type& basisrInverse)
		{
			for (SizeType i=0;i<opSrc.data.row();i++) {
				PairType jm = basis.jmValue(i);
				for (int l=opSrc.data.getRowPtr(i);l<opSrc.data.getRowPtr(i+1);l++) {
					SizeType iprime = opSrc.data.getCol(l);
					PairType jmPrime = basis.jmValue(iprime);

					RealType divisor = opSrc.angularFactor*(jmPrime.first+1);
					opDest1(basisrInverse[i],basisrInverse[iprime]) +=
					                   opSrc.data.element(i,iprime)*cgObject_(jmPrime,jm,opSrc.jm)/divisor;

				}
			}
		}

		void buildAdditional(PsimagLite::Matrix<SparseElementType>& lfactor,
		                     SizeType k,
		                     SizeType mu1,
		                     SizeType mu2,
		                     typename PsimagLite::Vector<PairType>::Type& jsEffective)
		{
			PairType kmu1(k,mu1);
			PairType kmu2(k,mu2);
			SizeType counter=0;
			int offset = lrs_.super().partition(m_);
			PairType jm=lrs_.super().jmValue(offset);
			lfactor.resize(lrs_.left().jMax()*lrs_.right().jMax(),lrs_.left().jMax()*lrs_.right().jMax());
			for (SizeType i1=0;i1<lrs_.left().jVals();i1++) {
				for (SizeType i2=0;i2<lrs_.right().jVals();i2++) {
					for (SizeType i1prime=0;i1prime<lrs_.left().jVals();i1prime++) {
						for (SizeType i2prime=0;i2prime<lrs_.right().jVals();i2prime++) {
							SparseElementType sum=calcLfactor(lrs_.left().jVals(i1),lrs_.right().jVals(i2),
							    lrs_.left().jVals(i1prime),lrs_.right().jVals(i2prime),jm,kmu1,kmu2);
							if (sum!=static_cast<SparseElementType>(0)) {
								counter++;
								PairType jj(PairType(lrs_.left().jVals(i1),lrs_.right().jVals(i2)));
								int x3=PsimagLite::isInVector(jsEffective,jj);
								if (x3<0) jsEffective.push_back(jj);
							}
							lfactor(lrs_.left().jVals(i1)+lrs_.right().jVals(i2)*lrs_.left().jMax(),
							    lrs_.left().jVals(i1prime)+lrs_.right().jVals(i2prime)*lrs_.left().jMax())=sum;

						}
					}
				}
			}
		}

		void calcHamFactor(PsimagLite::Matrix<SparseElementType>& lfactor)
		{
			int offset = lrs_.super().partition(m_);
			PairType jm=lrs_.super().jmValue(offset);
			lfactor.resize(lrs_.left().jMax(),lrs_.right().jMax());
			PairType kmu(0,0);
			for (SizeType i1=0;i1<lrs_.left().jVals();i1++) {
				for (SizeType i2=0;i2<lrs_.right().jVals();i2++) {
					lfactor(lrs_.left().jVals(i1),lrs_.right().jVals(i2))=
					       calcLfactor(lrs_.left().jVals(i1),lrs_.right().jVals(i2),
					       lrs_.left().jVals(i1),lrs_.right().jVals(i2),jm,kmu,kmu);
				}
			}
		}

		void calcEffectiveStates(typename PsimagLite::Vector<PairType>::Type jsEffective)
		{
			int offset = lrs_.super().partition(m_);
			PairType jm=lrs_.super().jmValue(offset);
			reducedEffective_.clear();
			reducedInverse_.resize(lrs_.left().reducedSize(),lrs_.right().reducedSize());
			SizeType electrons = lrs_.super().electrons(offset);
			for (SizeType i=0;i<jsEffective.size();i++) {
				for (SizeType i1=0;i1<lrs_.left().reducedSize();i1++) {
					for (SizeType i2=0;i2<lrs_.right().reducedSize();i2++) {
						if (electrons!=lrs_.left().electrons(lrs_.left().reducedIndex(i1))+
							lrs_.right().electrons(lrs_.right().reducedIndex(i2)))
							continue;
						PairType jj(lrs_.left().jmValue(lrs_.left().reducedIndex(i1)).first,
						     lrs_.right().jmValue(lrs_.right().reducedIndex(i2)).first);
						if (jj!=jsEffective[i]) continue;
						reducedEffective_.push_back(PairType(i1,i2));
						reducedInverse_(i1,i2)=reducedEffective_.size()-1;
					}
				}
			}

			getFlavorMap(jm);

			reorderMap();
		}

		SparseElementType calcLfactor(SizeType j1,
		                              SizeType j2,
		                              SizeType j1prime,
		                              SizeType j2prime,
		                              const PairType& jm,
		                              const PairType& kmu1,
		                              const PairType& kmu2) const
		{
			SparseElementType sum=0;
			SizeType k = kmu1.first;
			SizeType mu1=kmu1.second;
			SizeType mu2=kmu2.second;

			for (SizeType m1=0;m1<=j1;m1++) {
				PairType jm1(j1,m1);
				int x = j1+j2-jm.first;
				if (x%2!=0) continue;
				if (int(int(x/2)-m1 + jm.second)<0) continue;
				SizeType m2 = SizeType(x/2)-m1 + jm.second;
				PairType jm2(j2,m2);
				if (m2>j2) continue;

				x = j1prime-k-j1;
				if (x%2!=0) continue;
				if (static_cast<int>(m1+mu1)+static_cast<int>(x/2)<0) continue;
				SizeType m1prime = m1+mu1+SizeType(x/2);
				if (m1prime>j1prime) continue;
				PairType jm1prime(j1prime,m1prime);

				x = j2prime-k-j2;
				if (x%2!=0) continue;
				if (static_cast<int>(m2+mu2)+static_cast<int>(x/2)<0) continue;
				SizeType m2prime = m2+mu2+SizeType(x/2);
				if (m2prime>j2prime) continue;
				PairType jm2prime(j2prime,m2prime);

				x=j1prime+j2prime-jm.first;
				if (x%2!=0) continue;
				x= x/2;
				x += jm.second -m1prime-m2prime;
				if (x!=0) continue;

				sum +=  cgObject_(jm1prime,jm1,kmu1)*
				        cgObject_(jm2prime,jm2,kmu2)*
				        cgObject_(jm,jm1,jm2)*cgObject_(jm,jm1prime,jm2prime);

			}
			return sum;
		}

		void getFlavorMap(const PairType& jm)
		{
			PsimagLite::Map<SizeType,SizeType>::Type flavorsOldInverse;
			lrs_.super().flavor2Index(flavorsOldInverse,jm);
			flavorsOldInverse_.resize(reducedEffective_.size());

			for (SizeType i=0;i<reducedEffective_.size();i++) {
				SizeType i1= lrs_.left().reducedIndex(reducedEffective_[i].first);
				SizeType f1= lrs_.left().getFlavor(i1);
				SizeType ne1= lrs_.left().electrons(i1);
				SizeType j1= lrs_.left().jmValue(i1).first;

				SizeType i2= lrs_.right().reducedIndex(reducedEffective_[i].second);
				SizeType f2= lrs_.right().getFlavor(i2);
				SizeType ne2= lrs_.right().electrons(i2);
				SizeType j2= lrs_.right().jmValue(i2).first;

				SizeType f=lrs_.super().flavor2Index(f1,f2,ne1,ne2,j1,j2);
				flavorsOldInverse_[i]=flavorsOldInverse[f];
			}
		}

		void reorderMap()
		{
			if (flavorsOldInverse_.size()==0) return;
			typename PsimagLite::Vector<SizeType>::Type perm(flavorsOldInverse_.size());
			PsimagLite::Sort<typename PsimagLite::Vector<SizeType>::Type > sort;
			sort.sort(flavorsOldInverse_,perm);
			typename PsimagLite::Vector<PairType>::Type r(reducedEffective_.size());
			PsimagLite::Matrix<SizeType> reducedInverse(reducedInverse_.n_row(),reducedInverse_.n_col());

			for (SizeType i=0;i<reducedEffective_.size();i++) {
				r[i]=reducedEffective_[perm[i]];
				reducedInverse(r[i].first,r[i].second)=i;
			}
			reducedEffective_=r;
			reducedInverse_=reducedInverse;
		}

		SizeType m_;
		const LeftRightSuperType& lrs_;
		typename PsimagLite::Vector<PsimagLite::Matrix<SparseElementType> >::Type lfactor_;
		PsimagLite::Matrix<SparseElementType> lfactorHamiltonian_;
		SparseMatrixType hamiltonian2_,hamiltonian3_;
		typename PsimagLite::Vector<PairType>::Type reducedEffective_;
		PsimagLite::Matrix<SizeType> reducedInverse_;
		typename PsimagLite::Vector<SizeType>::Type flavorsOldInverse_;
		ClebschGordanType& cgObject_;

	}; // class

//	template<typename LeftRightSuperType,
//		typename ReflectionSymmetryType,
//		typename ConcurrencyType>
//	const typename LeftRightSuperType::BasisType* Su2Reduced<
//	LeftRightSuperType,ReflectionSymmetryType,ConcurrencyType>::basis1Ptr_=0;
} // namespace Dmrg
/*@}*/
#endif

