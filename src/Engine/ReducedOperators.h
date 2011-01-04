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

/*! \file ReducedOperators.h
 *
 *  FIXME
 *
 */
#ifndef REDUCEDOP_IMPL_H
#define REDUCEDOP_IMPL_H

#include "Utils.h"
#include "Su2SymmetryGlobals.h"
#include "Operator.h"

namespace Dmrg {
	template<typename OperatorType,typename DmrgBasisType>
	class ReducedOperators {
			typedef typename OperatorType::SparseMatrixType SparseMatrixType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef typename DmrgBasisType::RealType RealType;
			typedef typename OperatorType::PairType PairType;
			typedef ClebschGordanCached<RealType> ClebschGordanType;
			typedef psimag::Matrix<SparseElementType> DenseMatrixType;
			typedef Su2SymmetryGlobals<RealType> Su2SymmetryGlobalsType;
			
		public:
			ReducedOperators(const DmrgBasisType* thisBasis,size_t dof,size_t orbitals) 
				: thisBasis_(thisBasis),useSu2Symmetry_(DmrgBasisType::useSu2Symmetry()),
					     dof_(dof),nOrbitals_(orbitals),
					cgObject_(&(Su2SymmetryGlobalsType::clebschGordanObject))
			{
			}
			
			template<typename IoInputter>
			ReducedOperators(IoInputter& io,size_t level,const DmrgBasisType* thisBasis,size_t dof,size_t orbitals)
			 : thisBasis_(thisBasis),useSu2Symmetry_(DmrgBasisType::useSu2Symmetry()),
					     dof_(dof),nOrbitals_(orbitals),
					cgObject_(&(Su2SymmetryGlobalsType::clebschGordanObject))
			{
				if (!useSu2Symmetry_) return;
				io.read(reducedOperators_,"#OPERATORS");
			}

			const OperatorType& getReducedOperatorByIndex(int i) const 
			{
				return reducedOperators_[i];	
			}

			const OperatorType& getReducedOperatorByIndex(char modifier,int i) const 
			{
				size_t s = size_t(i/dof_);
				size_t tmp = i%dof_;
				size_t offset = reducedOperators_[s*dof_].su2Related.offset;
				size_t g = tmp%offset;
				size_t i0 = s*dof_+g;

				if (modifier=='N') return getReducedOperatorByIndex(i0);
				else return getReducedOperatorByIndex(i0+offset);
			}

			const SparseMatrixType& hamiltonian() const { return reducedHamiltonian_; }

			void setOperators(const std::vector<OperatorType>& ops)
			{
				if (!useSu2Symmetry_) return;
				findBasisInverse();		
				reducedOperators_.resize(ops.size());
				
				for (size_t i=0;i<ops.size();i++) {
					size_t angularMomentum =ops[i].jm.first;
					std::vector<const OperatorType*> opSrc(angularMomentum+1);
					if (ops[i].su2Related.source.size()==0) continue;
					OperatorType myOp[angularMomentum+1];
					size_t transposeCounter=0;
					
					for (size_t m=0;m<=angularMomentum;m++) {
						int tr = ops[i].su2Related.transpose[m];
						if (tr>=0) {
							transposeConjugate(myOp[transposeCounter].data,
									ops[ops[i].su2Related.source[m]].data );
							myOp[transposeCounter].fermionSign=ops[ops[i].su2Related.source[m]].fermionSign;
							myOp[transposeCounter].jm=typename OperatorType::PairType(angularMomentum,
								angularMomentum-ops[ops[i].su2Related.source[m]].jm.second);	
							myOp[transposeCounter].angularFactor =
									-ops[ops[i].su2Related.source[m]].angularFactor;
							opSrc[m]=&myOp[transposeCounter];
							transposeCounter++;
						} else {
							opSrc[m]=&ops[ops[i].su2Related.source[m]];
						}
					}
					createReducedOperator(reducedOperators_[i].data,opSrc);
					reducedOperators_[i].fermionSign=opSrc[0]->fermionSign;
					reducedOperators_[i].jm=opSrc[0]->jm;
					reducedOperators_[i].angularFactor=opSrc[0]->angularFactor;

					reducedOperators_[i].su2Related.offset = ops[i].su2Related.offset;
					size_t i1 = i +  ops[i].su2Related.offset;
					createReducedConj(angularMomentum,reducedOperators_[i1].data,reducedOperators_[i].data);
					reducedOperators_[i1].fermionSign=opSrc[0]->fermionSign;
					reducedOperators_[i1].jm=opSrc[0]->jm;
					reducedOperators_[i1].angularFactor=opSrc[0]->angularFactor;
				}
			}

			void setMomentumOfOperators(const std::vector<size_t>& momentum)
			{
				momentumOfOperators_=momentum;
				
			}

			void setHamiltonian(const SparseMatrixType& hamiltonian)
			{
				if (!useSu2Symmetry_) return;
				findBasisInverse();	
				size_t angularMomentum =0;
				std::vector<const OperatorType*> opSrc(angularMomentum+1);

				OperatorType myOp; 
				myOp.data = hamiltonian;
				myOp.fermionSign=1;
				myOp.jm=typename OperatorType::PairType(0,0);	
				myOp.angularFactor = 1.0;
				opSrc[0]=&myOp;
				createReducedOperator(reducedHamiltonian_,opSrc);
			}

			void setToProduct(const DmrgBasisType& basis2,const DmrgBasisType& basis3,size_t x,const DmrgBasisType* thisBasis)
			{
				if (!useSu2Symmetry_) return;
				thisBasis_ = thisBasis;
				reducedOperators_.resize(x);

				j1Max_=basis2.jMax();
				j2Max_=basis3.jMax();

				// build all lfactors
				lfactorLeft_.resize(momentumOfOperators_.size());
				for (size_t i=0;i<momentumOfOperators_.size();i++) {
					buildLfactor(lfactorLeft_[i],true,basis2,basis3,momentumOfOperators_[i]); // left
				}
				buildLfactor(lfactorHamLeft_,true,basis2,basis3,0);
				
				lfactorRight_.resize(momentumOfOperators_.size());
				for (size_t i=0;i<momentumOfOperators_.size();i++) {
					buildLfactor(lfactorRight_[i],false,basis2,basis3,momentumOfOperators_[i]); // right
				}
				buildLfactor(lfactorHamRight_,false,basis2,basis3,0);
				calcReducedMapping(basis2,basis3);
				cacheFlavorIndex(basis2,basis3);
				calcFastBasis(fastBasisLeft_,basis2,basis3,true,thisBasis_->reducedSize());
				calcFastBasis(fastBasisRight_,basis2,basis3,false,thisBasis_->reducedSize());
			}

			void prepareTransform(const DenseMatrixType& ftransform,const DmrgBasisType* thisBasis)
			{
				if (!useSu2Symmetry_) return;
				size_t nr=thisBasis->reducedSize();
				size_t nold = thisBasis_->reducedSize();
				ftransform_.resize(nold,nr);

				for (size_t i=0;i<nold;i++) {
					size_t ii =thisBasis_->reducedIndex(i); // old
					for (size_t j=0;j<nr;j++) {
						size_t jj = thisBasis->reducedIndex(j); //new 
						ftransform_(i,j)=ftransform(ii,jj);
					}
				}
				thisBasis_ = thisBasis;
			}

			void changeBasis(size_t k)
			{
				if (!useSu2Symmetry_) return;
				changeBasis(reducedOperators_[k].data);
			}

			void changeBasisHamiltonian()
			{
				if (!useSu2Symmetry_) return;
				changeBasis(reducedHamiltonian_);
			}

			void externalProduct(size_t ind,const DmrgBasisType& basis2,const DmrgBasisType& basis3,bool order,
					     const OperatorType& myOp)
			{
				if (!useSu2Symmetry_) return;
				size_t n = thisBasis_->reducedSize();

				const DmrgBasisType* basisA = &basis2;
				const DmrgBasisType* basisB = &basis3;
				size_t angularMomentum = myOp.jm.first;
				int fermionSign = myOp.fermionSign;
				const SparseMatrixType& A = myOp.data;
				std::vector<size_t> jvals;

				if (!order) {
					basisA = &basis3;
					basisB = &basis2;
				}

				int ki = utils::isInVector(momentumOfOperators_,angularMomentum);
				if (ki<0) throw std::runtime_error("Operator has unknown momentum\n");
				psimag::Matrix<SparseElementType> B(n,n);
				externalProd_(B,basisA,basisB,A,ki,order,fermionSign);
				fullMatrixToCrsMatrix(reducedOperators_[ind].data,B);
				reducedOperators_[ind].fermionSign = fermionSign;
				reducedOperators_[ind].jm = myOp.jm;
				reducedOperators_[ind].angularFactor = myOp.angularFactor;
				reducedOperators_[ind].su2Related = myOp.su2Related;
			}

			void outerProductHamiltonian(const DmrgBasisType& basis2,const DmrgBasisType& basis3,
					     const SparseMatrixType& h2, const SparseMatrixType& h3)
			{
				if (!useSu2Symmetry_) return;

				size_t n = thisBasis_->reducedSize();
				const DmrgBasisType* basisA = &basis2;
				const DmrgBasisType* basisB = &basis3;
				psimag::Matrix<SparseElementType> B2(n,n);
				externalProd_(B2,basisA,basisB,h2,-1,true,1);

				basisA = &basis3;
				basisB = &basis2;
				psimag::Matrix<SparseElementType> B3(n,n);
				externalProd_(B3,basisA,basisB,h3,-1,false,1);

				B2 += B3;

				fullMatrixToCrsMatrix(reducedHamiltonian_,B2);
			}

			void reorder(size_t k,const std::vector<size_t>& permutation)
			{
				if (!useSu2Symmetry_) return;
				for (size_t i=0;i<permutation.size();i++) 
					if (permutation[i]!=i) 
						throw std::runtime_error("reorderHamiltonian: permutation not the identity!\n");
			}

			void reorderHamiltonian(const std::vector<size_t>& permutation)
			{
				if (!useSu2Symmetry_) return;
				for (size_t i=0;i<permutation.size();i++) 
					if (permutation[i]!=i) 
						throw std::runtime_error("reorderHamiltonian: permutation not the identity!\n");
			}

			size_t size() const { return reducedOperators_.size(); }

			template<typename ConcurrencyType>
			void gather(ConcurrencyType& concurrency)
			{
				Dmrg::gather(reducedOperators_,concurrency);
			}

			template<typename ConcurrencyType>
			void broadcast( ConcurrencyType& concurrency)
			{
				Dmrg::broadcast(reducedOperators_,concurrency);
			}

			template<typename IoOutputter>
			void save(IoOutputter& io,const std::string& s) const
			{
				io.printVector(reducedOperators_,"#OPERATORS");
			}

		private:
			const DmrgBasisType* thisBasis_;
			bool useSu2Symmetry_;
			size_t dof_;
			size_t nOrbitals_;
			ClebschGordanType* cgObject_;
			std::vector<size_t> momentumOfOperators_;
			std::vector<size_t> basisrinverse_;
			std::vector<OperatorType> reducedOperators_;
			SparseMatrixType reducedHamiltonian_;
			size_t j1Max_,j2Max_;
			std::vector<std::vector<SparseElementType> > lfactorLeft_;
			std::vector<std::vector<SparseElementType> > lfactorRight_;
			std::vector<SparseElementType> lfactorHamLeft_,lfactorHamRight_;
			psimag::Matrix<int> reducedMapping_;
			std::vector<std::vector<size_t> > fastBasisLeft_,fastBasisRight_;
			psimag::Matrix<size_t> flavorIndexCached_;
			DenseMatrixType ftransform_;

			SparseElementType lfactor(int ki,bool order,size_t j1,size_t j2,size_t j1prime,size_t jProd,size_t jProdPrime) const
			{
				size_t ix=lfactorIndex(order,j1,j2,j1prime,jProd,jProdPrime);
				if (ki<0) {
					if (order) return lfactorHamLeft_[ix];
					return lfactorHamRight_[ix];
				}
				if (order) return lfactorLeft_[ki][ix];
				return lfactorRight_[ki][ix];
			}

			size_t lfactorIndex(bool order,size_t j1,size_t j2,size_t j1prime,size_t jProd,size_t jProdPrime) const
			{
				if (order) return lfactorIndexLeft(j1,j2,j1prime,jProd,jProdPrime);
				return lfactorIndexRight(j1,j2,j1prime,jProd,jProdPrime);
			}

			size_t lfactorIndexLeft(size_t j1,size_t j2,size_t j1prime,size_t jProd,size_t jProdPrime) const
			{
				size_t ix = j1+j2*j1Max_+j1prime*j1Max_*j2Max_+jProd*j1Max_*j1Max_*j2Max_+
						jProdPrime*thisBasis_->jMax()*j1Max_*j1Max_*j2Max_;
				return ix;	
			}

			size_t lfactorIndexRight(size_t j2,size_t j1,size_t j2prime,size_t jProd,size_t jProdPrime) const
			{
				size_t ix = j2+j1*j2Max_+j2prime*j2Max_*j1Max_+jProd*j2Max_*j2Max_*j1Max_+
						jProdPrime*thisBasis_->jMax()*j2Max_*j2Max_*j1Max_;
				return ix;	
			}

			void calcReducedMapping(const DmrgBasisType& basis2,const DmrgBasisType& basis3)
			{
				size_t fMax=0;
				const std::vector<size_t>& flavorsOld=thisBasis_->flavorsOld();
				
				for (size_t i=0;i<thisBasis_->reducedSize();i++) {
					size_t ii = thisBasis_->reducedIndex(i);
					if (ii>=flavorsOld.size()) throw std::runtime_error("Ouch!!\n");
					size_t f = flavorsOld[ii];
					if (f>fMax) fMax=f;
				}
				fMax++;
				reducedMapping_.resize(thisBasis_->jMax(),fMax);
				for (size_t i=0;i<reducedMapping_.n_row();i++) 
					for (size_t j=0;j<reducedMapping_.n_col();j++)
						reducedMapping_(i,j)= -1; 
				for (size_t i=0;i<thisBasis_->reducedSize();i++) {
					size_t ii = thisBasis_->reducedIndex(i);
					size_t j = thisBasis_->jmValue(ii).first;
					size_t f = flavorsOld[ii];
					reducedMapping_(j,f)=i;
				}
			}

			void findBasisInverse()
			{
				basisrinverse_.resize(thisBasis_->size());
				for (size_t i=0;i<thisBasis_->size();i++) {
					size_t f = thisBasis_->getFlavor(i);
					size_t j = thisBasis_->jmValue(i).first;
					basisrinverse_[i]=findJf(j,f);
				}
			}

			size_t findJf(size_t j,size_t f)
			{
				for (size_t i=0;i<thisBasis_->reducedSize();i++) {
					PairType jm=thisBasis_->jmValue(thisBasis_->reducedIndex(i));
					if (jm.first!=j) continue;
					if (thisBasis_->getFlavor(thisBasis_->reducedIndex(i))!=f) continue;
					return i;	
				}
				return 0; //bogus
			}

			void createReducedConj(size_t k1 ,SparseMatrixType& opDest,const SparseMatrixType& opSrc)
			{
				size_t n=opSrc.rank();
				transposeConjugate(opDest,opSrc);
				for (size_t i=0;i<n;i++) {
					PairType jm = thisBasis_->jmValue(thisBasis_->reducedIndex(i));
					for (int k=opDest.getRowPtr(i);k<opDest.getRowPtr(i+1);k++) {
						size_t j=opDest.getCol(k);
						PairType jmPrime = thisBasis_->jmValue(thisBasis_->reducedIndex(j));
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

			void createReducedOperator(SparseMatrixType& opDest,const std::vector<const OperatorType*>& opSrc)
			{
				size_t n = thisBasis_->reducedSize();
				DenseMatrixType opDest1(n,n);
				for (size_t i=0;i<opSrc.size();i++) 
					createReducedOperator(opDest1,*opSrc[i]); //,PairType(k,mu1),mysign);
				fullMatrixToCrsMatrix(opDest,opDest1);
			}

			void createReducedOperator(DenseMatrixType& opDest1,const OperatorType& opSrc)
			{
				for (size_t i=0;i<opSrc.data.rank();i++) {
					PairType jm = thisBasis_->jmValue(i);
					for (int l=opSrc.data.getRowPtr(i);l<opSrc.data.getRowPtr(i+1);l++) {
						size_t iprime = opSrc.data.getCol(l);
						PairType jmPrime = thisBasis_->jmValue(iprime);
						RealType divisor = opSrc.angularFactor*(jmPrime.first+1);
						opDest1(basisrinverse_[i],basisrinverse_[iprime]) += 
									opSrc.data(i,iprime)*cgObject_->operator()(jmPrime,jm,opSrc.jm)/divisor;
					}
				}
			}

			void changeBasis(SparseMatrixType &v)
			{	
				SparseMatrixType tmpMatrix=transformFullFast(v,ftransform_);
				v= tmpMatrix;
			}

			void buildLfactor(
					std::vector<SparseElementType>& lfactor,
					bool order,const DmrgBasisType& basis2,
     					const DmrgBasisType& basis3,
	  				size_t k)
			{
				size_t jMax=j1Max_;
				if (!order) jMax = j2Max_;

				const DmrgBasisType* basisA = &basis2;
				const DmrgBasisType* basisB = &basis3;
				if (!order) {
					basisA = &basis3;
					basisB = &basis2;
				}

				lfactor.resize(j1Max_*j2Max_*jMax*thisBasis_->jMax()*thisBasis_->jMax());				
				for (size_t i1=0;i1<basisA->jVals();i1++) {
					for (size_t i2=0;i2<basisB->jVals();i2++) {
						for (size_t i1prime=0;i1prime<basisA->jVals();i1prime++) {
							for (size_t i=0;i<thisBasis_->jVals();i++) {
								for (size_t iprime=0;iprime<thisBasis_->jVals();iprime++) {
									size_t ix = lfactorIndex(order,basisA->jVals(i1),
										basisB->jVals(i2),basisA->jVals(i1prime),
										thisBasis_->jVals(i),thisBasis_->jVals(iprime));
									
									lfactor[ix]=calcLfactor(order,basisA->jVals(i1),basisB->jVals(i2),
										basisA->jVals(i1prime),thisBasis_->jVals(i),
										thisBasis_->jVals(iprime),k);
								}
							}
						}
					}
				}
			}

			SparseElementType calcLfactor(bool order,size_t j1,size_t j2,size_t j1prime,size_t jProd,size_t jProdPrime,
					size_t k) const
			{
				if (order) return 	 calcLfactorLeft(j1,j2,j1prime,jProd,jProdPrime,k);
				return calcLfactorRight(j1,j2,j1prime,jProd,jProdPrime,k);
			}

			SparseElementType calcLfactorLeft(size_t j1,size_t j2,size_t j1prime,size_t jProd,size_t jProdPrime,size_t k) const
			{
				SparseElementType sum=0;
				for (size_t m1=0;m1<=j1;m1++) {
					PairType jm1(j1,m1);
					for (size_t m2=0;m2<=j2;m2++) {
						PairType jm2(j2,m2);
						for (size_t mCapital=0;mCapital<=k;mCapital++) {
							PairType jmCapital(k,mCapital);
							// get m
							int x = j1+j2-jProd;
							if (x%2!=0) continue;
							x/=2;
							if (m1+m2-x<0) continue;
							size_t m = m1 + m2 - x;
							if (m>jProd) continue;
							PairType jm(jProd,m);

							// get m'=m+M
							x = -jProdPrime + k +jProd;
							if (x%2!=0) continue;
							x/=2;
							if (m +mCapital-x<0) continue;
							size_t mPrime = m +mCapital-x;
							if (mPrime>jProdPrime) continue;
							PairType jmPrime(jProdPrime,mPrime);

							// get m1prime=mPrime - m2
							x = j1prime + j2 - jProdPrime;
							if (x%2!=0) continue;
							x/=2;
							if (mPrime + x - m2 < 0) continue;
							size_t m1prime = mPrime + x - m2;
							if (m1prime > j1prime) continue;
							PairType jm1prime(j1prime,m1prime);

							// m1prime = m1 + M  
							x = k - j1prime +j1;
							if (x%2!=0) continue;
							x/=2;
							if (m1 - x + mCapital < 0) continue;
							if (m1 - x + mCapital != m1prime) continue;

							sum +=  cgObject_->operator()(jmPrime,jmCapital,jm)*
								cgObject_->operator()(jm,jm1,jm2)*
								cgObject_->operator()(jm1prime,jmCapital,jm1)*
								cgObject_->operator()(jmPrime,jm1prime,jm2);
						}
					}
				}
				return sum/RealType(jProdPrime+1);
			}

			SparseElementType calcLfactorRight(size_t j2,size_t j1,size_t j2prime,size_t jProd,size_t jProdPrime,size_t k) const
			{
				SparseElementType sum=0;
				for (size_t m1=0;m1<=j1;m1++) {
					PairType jm1(j1,m1);
					for (size_t m2=0;m2<=j2;m2++) {
						PairType jm2(j2,m2);
						for (size_t mCapital=0;mCapital<=k;mCapital++) {
							PairType jmCapital(k,mCapital);
							// get m=m1+m2
							int x = j1+j2-jProd;
							if (x%2!=0) continue;
							x/=2;
							if (m1+m2-x<0) continue;
							size_t m = m1 + m2 - x;
							if (m>jProd) continue;
							PairType jm(jProd,m);

							// get m'=m+mCapital
							x = -jProdPrime + k +jProd;
							if (x%2!=0) continue;
							x/=2;
							if (m - x + mCapital<0) continue;
							size_t mPrime = m - x +mCapital;
							if (mPrime>jProdPrime) continue;
							PairType jmPrime(jProdPrime,mPrime);

							// get m2prime
							x = j2prime + j1 - jProdPrime;
							if (x%2!=0) continue;
							x/=2;
							if (mPrime + x - m1 < 0) continue;
							size_t m2prime = mPrime + x - m1;
							if (m2prime > j2prime) continue;
							PairType jm2prime(j2prime,m2prime);

							// check m2prime
							x = k - j2prime + j2;
							if (x%2!=0) continue;
							x/=2;
							if (m2 - x + mCapital < 0) continue;
							if (m2 - x + mCapital != m2prime) continue;

							sum +=  cgObject_->operator()(jmPrime,jmCapital,jm)*
								cgObject_->operator()(jm,jm1,jm2)*
								cgObject_->operator()(jm2prime,jmCapital,jm2)*
								cgObject_->operator()(jmPrime,jm1,jm2prime);
						}
					}
				}
				return sum/RealType(jProdPrime+1);
			}

			bool checkJvalues(size_t jProd,size_t j1,size_t j2) const
			{
				if (jProd>j1+j2) return false;
				if (j1>j2 && jProd<j1-j2) return false;
				if (j2>j1 && jProd<j2-j1) return false;	
				return true;
			}

			void externalProd_(psimag::Matrix<SparseElementType>& B,const DmrgBasisType* basisA,
					  const DmrgBasisType* basisB,const SparseMatrixType& A,int ki,bool order,int fermionSign)
			{
				size_t n = B.n_row();
				size_t angularMomentum = 0;
				if (ki>=0) angularMomentum = momentumOfOperators_[ki];
				const std::vector<std::vector<size_t> >* fastBasis = &fastBasisLeft_;
				if (!order) fastBasis = &fastBasisRight_;
				for (size_t i0=0;i0<fastBasis->size();i0++) {
					const std::vector<size_t>& twopairs = (*fastBasis)[i0];
					size_t i = twopairs[0];
					size_t iprime = twopairs[1];;
					size_t i1 = twopairs[2];
					size_t i2 = twopairs[3];

					size_t ii1 = basisA->reducedIndex(i1);
					size_t j1 = basisA->jmValue(ii1).first;

					size_t ii2 = basisB->reducedIndex(i2);
					size_t j2 = basisB->jmValue(ii2).first;
					size_t ne2 = basisB->electrons(ii2);

					size_t jProd = thisBasis_->jVals(i);
					size_t jProdPrime = thisBasis_->jVals(iprime);

					int ix = twopairs[4];
					RealType fsign = 1;
					if (!order && (ne2%2 !=0)) fsign = fermionSign;

					for (int k=A.getRowPtr(i1);k<A.getRowPtr(i1+1);k++) {
						size_t i1prime = A.getCol(k);
						size_t ii1prime = basisA->reducedIndex(i1prime);
						size_t j1prime  = basisA->jmValue(ii1prime).first;

						size_t fprime=0;
						if (order) fprime = flavorIndexCached_(i1prime,i2);
						else fprime = flavorIndexCached_(i2,i1prime);

						SparseElementType tmp = lfactor(ki,order,j1prime,j2,j1,jProdPrime,jProd);
						if (tmp==static_cast<SparseElementType>(0)) continue;
						int iy = reducedMapping_(jProdPrime,fprime);
						if (iy<0 || iy>=int(n)) continue;

						// magic correction::
						tmp *= sqrt(jProd+1.0)/sqrt(jProdPrime+1.0);
						tmp *= sqrt(j1prime+1.0)/sqrt(j1+1.0);

						tmp *= fsign * A.getValue(k);
						B(ix,iy)=tmp;
					}
				}
			}

			void calcFastBasis(std::vector<std::vector<size_t> >& fastBasis,const DmrgBasisType& basis2,
					  const DmrgBasisType& basis3,bool order,size_t n)
			{
				const DmrgBasisType* basisA = &basis2;
				const DmrgBasisType* basisB = &basis3;

				if (!order) {
					basisA = &basis3;
					basisB = &basis2;
				}
				fastBasis.clear();
				std::vector<size_t> twopairs(5);
				for (size_t i=0;i<thisBasis_->jVals();i++) {
					twopairs[0]=i;
					size_t jProd = thisBasis_->jVals(i);
					for (size_t iprime=0;iprime<thisBasis_->jVals();iprime++) {
						twopairs[1]=iprime;
						for (size_t i1=0;i1<basisA->reducedSize();i1++) {
							twopairs[2]=i1;
							size_t ii1 = basisA->reducedIndex(i1);
							size_t ne1 = basisA->electrons(ii1);
							size_t j1 = basisA->jmValue(ii1).first;
							size_t f1 = basisA->getFlavor(ii1);
							
							for (size_t i2=0;i2<basisB->reducedSize();i2++) {
								
								twopairs[3]=i2;
								size_t ii2 = basisB->reducedIndex(i2);
								size_t ne2 = basisB->electrons(ii2);
								size_t j2 = basisB->jmValue(ii2).first;
								size_t f2 = basisB->getFlavor(ii2);
								if (jProd>j1+j2) continue;
								if (j1>j2 && jProd<j1-j2) continue;
								if (j2>j1 && jProd<j2-j1) continue;

								size_t f = 0;
								if (order) f = thisBasis_->flavor2Index(f1,f2,ne1,ne2,j1,j2);
								else f = thisBasis_->flavor2Index(f2,f1,ne2,ne1,j2,j1);

								// jProd,f --> ix
								int ix = reducedMapping_(jProd,f);
								if (ix<0 || size_t(ix)>=n) continue;
								twopairs[4]=ix;

								fastBasis.push_back(twopairs);
							}
						}
					}
				}
			}

			void cacheFlavorIndex(const DmrgBasisType& basis2,
					  const DmrgBasisType& basis3)
			{
				flavorIndexCached_.resize(basis2.reducedSize(),basis3.reducedSize());
				for (size_t i1=0;i1<basis2.reducedSize();i1++) {
					size_t ii1 = basis2.reducedIndex(i1);
					size_t ne1 = basis2.electrons(ii1);
					size_t j1 = basis2.jmValue(ii1).first;
					size_t f1 = basis2.getFlavor(ii1);
					for (size_t i2=0;i2<basis3.reducedSize();i2++) {
						size_t ii2 = basis3.reducedIndex(i2);
						size_t j2  = basis3.jmValue(ii2).first;
						size_t f2 = basis3.getFlavor(ii2);
						size_t ne2 = basis3.electrons(ii2);
						flavorIndexCached_(i1,i2)=thisBasis_->flavor2Index(f1,f2,ne1,ne2,j1,j2);
					}
				}
			}
 	}; // ReducedOperators
} // namespace Dmrg

/*@}*/
#endif
