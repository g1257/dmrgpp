/*
Copyright (c) 2009-2015, 2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file ReducedOperators.h
 *
 *  FIXME
 *
 */
#ifndef REDUCEDOP_IMPL_H
#define REDUCEDOP_IMPL_H

#include "Su2SymmetryGlobals.h"
#include "Operator.h"
#include "ChangeOfBasis.h"
#include "BlockOffDiagMatrix.h"
#include "../KronUtil/MatrixDenseOrSparse.h"

namespace Dmrg {
template<typename BasisType>
class ReducedOperators {

	typedef typename BasisType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	// typedef BlockOffDiagMatrix<MatrixDenseOrSparse<SparseMatrixType> > BlockOffDiagMatrixType;
	typedef Operator<SparseMatrixType> OperatorType_;
	typedef typename BasisType::RealType RealType;
	typedef typename OperatorType_::PairType PairType;
	typedef ClebschGordanCached<RealType> ClebschGordanType;
	typedef PsimagLite::Matrix<SparseElementType> DenseMatrixType;
	typedef Su2SymmetryGlobals<RealType> Su2SymmetryGlobalsType;
	typedef typename PsimagLite::Vector<OperatorType_>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<const OperatorType_*>::Type
	VectorPointerOperatorType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef ChangeOfBasis<SparseMatrixType, DenseMatrixType> ChangeOfBasisType;

public:

	typedef typename ChangeOfBasisType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef OperatorType_ OperatorType;

	ReducedOperators(const BasisType* thisBasis)
	    : thisBasis_(thisBasis),
	      useSu2Symmetry_(BasisType::useSu2Symmetry()),
	      cgObject_(&(Su2SymmetryGlobalsType::clebschGordanObject))
	{
	}

	template<typename IoInputter>
	ReducedOperators(IoInputter& io,
	                 SizeType,
	                 const BasisType* thisBasis,
	                 typename PsimagLite::EnableIf<
	                 PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	    : thisBasis_(thisBasis),
	      useSu2Symmetry_(BasisType::useSu2Symmetry()),
	      cgObject_(&(Su2SymmetryGlobalsType::clebschGordanObject))
	{
		if (!useSu2Symmetry_) return;
		io.read(reducedOperators_,"Operators");
	}

	template<typename IoInputter>
	void read(IoInputter& io)
	{
		if (!useSu2Symmetry_) return;
		io.read(reducedOperators_,"Operators");
	}

	const OperatorType& getReducedOperatorByIndex(int i) const
	{
		return reducedOperators_[i];
	}

	const OperatorType& getReducedOperatorByIndex(char modifier,
	                                              const PairType& p) const
	{
		SizeType s = SizeType(p.first/p.second);
		SizeType tmp = p.first%p.second;
		SizeType offset = reducedOperators_[s*p.second].su2Related.offset;
		SizeType g = tmp%offset;
		SizeType i0 = s*p.second+g;

		if (modifier=='N') return getReducedOperatorByIndex(i0);
		else return getReducedOperatorByIndex(i0+offset);
	}

	const SparseMatrixType& hamiltonian() const { return reducedHamiltonian_; }

	void setOperators(const typename PsimagLite::Vector<OperatorType>::Type& ops)
	{
		if (!useSu2Symmetry_) return;
		findBasisInverse();
		reducedOperators_.resize(ops.size());

		for (SizeType i=0;i<ops.size();i++) {
			SizeType angularMomentum =ops[i].jm.first;
			VectorPointerOperatorType
			        opSrc(angularMomentum+1);
			if (ops[i].su2Related.source.size()==0) continue;
			VectorOperatorType myOp(angularMomentum+1);
			SizeType transposeCounter=0;

			for (SizeType m=0;m<=angularMomentum;m++) {
				int tr = ops[i].su2Related.transpose[m];
				if (tr>=0) {
					transposeConjugate(myOp[transposeCounter].data,
					                   ops[ops[i].su2Related.source[m]].data);
					myOp[transposeCounter].fermionOrBoson =
					        ops[ops[i].su2Related.source[m]].fermionOrBoson;
					myOp[transposeCounter].jm=typename OperatorType::PairType(
					            angularMomentum,
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
			reducedOperators_[i].fermionOrBoson=opSrc[0]->fermionOrBoson;
			reducedOperators_[i].jm=opSrc[0]->jm;
			reducedOperators_[i].angularFactor=opSrc[0]->angularFactor;

			reducedOperators_[i].su2Related.offset = ops[i].su2Related.offset;
			SizeType i1 = i +  ops[i].su2Related.offset;
			createReducedConj(angularMomentum,
			                  reducedOperators_[i1].data,
			                  reducedOperators_[i].data);
			reducedOperators_[i1].fermionOrBoson=opSrc[0]->fermionOrBoson;
			reducedOperators_[i1].jm=opSrc[0]->jm;
			reducedOperators_[i1].angularFactor=opSrc[0]->angularFactor;
		}
	}

	void setMomentumOfOperators(const PsimagLite::Vector<SizeType>::Type& momentum)
	{
		momentumOfOperators_=momentum;
	}

	void setHamiltonian(const SparseMatrixType& hamiltonian)
	{
		if (!useSu2Symmetry_) return;
		findBasisInverse();
		SizeType angularMomentum =0;
		VectorPointerOperatorType opSrc(angularMomentum+1);

		OperatorType myOp;
		myOp.data = hamiltonian;
		myOp.fermionOrBoson = ProgramGlobals::BOSON;
		myOp.jm=typename OperatorType::PairType(0,0);
		myOp.angularFactor = 1.0;
		opSrc[0]=&myOp;
		createReducedOperator(reducedHamiltonian_,opSrc);
	}

	void setToProduct(const BasisType& basis2,
	                  const BasisType& basis3,
	                  SizeType x,
	                  const BasisType* thisBasis)
	{
		if (!useSu2Symmetry_) return;
		thisBasis_ = thisBasis;
		reducedOperators_.resize(x);

		j1Max_=basis2.jMax();
		j2Max_=basis3.jMax();

		// build all lfactors
		lfactorLeft_.resize(momentumOfOperators_.size());
		for (SizeType i=0;i<momentumOfOperators_.size();i++) {
			buildLfactor(lfactorLeft_[i],true,basis2,basis3,momentumOfOperators_[i]);
		}

		buildLfactor(lfactorHamLeft_,true,basis2,basis3,0);

		lfactorRight_.resize(momentumOfOperators_.size());
		for (SizeType i=0;i<momentumOfOperators_.size();i++) {
			buildLfactor(lfactorRight_[i],false,basis2,basis3,momentumOfOperators_[i]);
		}

		buildLfactor(lfactorHamRight_,false,basis2,basis3,0);
		calcReducedMapping(basis2,basis3);
		cacheFlavorIndex(basis2,basis3);
		calcFastBasis(fastBasisLeft_,basis2,basis3,true,thisBasis_->reducedSize());
		calcFastBasis(fastBasisRight_,basis2,basis3,false,thisBasis_->reducedSize());
	}

	void prepareTransform(const BlockDiagonalMatrixType& ftransform1,
	                      const BasisType* thisBasis)
	{

		if (!useSu2Symmetry_) {
			changeOfBasis_.update(ftransform1);
			return;
		}

		SparseMatrixType ftransform;
		ftransform1.toSparse(ftransform);

		SizeType nr=thisBasis->reducedSize();
		SizeType nold = thisBasis_->reducedSize();

		PsimagLite::Matrix<SparseElementType> fm(nold,nr);

		PsimagLite::Vector<int>::Type inverseP(ftransform.cols(),-1);
		for (SizeType j=0;j<nr;j++) {
			SizeType jj = thisBasis->reducedIndex(j); //new
			assert(jj<inverseP.size());
			inverseP[jj] =  j;
		}
		for (SizeType i=0;i<nold;i++) {
			SizeType ii =thisBasis_->reducedIndex(i); // old
			for (int k = ftransform.getRowPtr(ii);k<ftransform.getRowPtr(ii+1);k++) {
				SizeType jj = ftransform.getCol(k);
				assert(jj<inverseP.size());
				int j = inverseP[jj];
				if (j<0) continue;
				fm(i,j)=ftransform.getValue(k);
			}
		}

		thisBasis_ = thisBasis;
		fullMatrixToCrsMatrix(su2Transform_, fm);
		transposeConjugate(su2TransformT_,su2Transform_);
	}

	void changeBasis(SizeType k)
	{
		if (!useSu2Symmetry_) return;
		changeBasis(reducedOperators_[k].data);
	}

	void changeBasisHamiltonian(SparseMatrixType& hamiltonian,
	                            const BlockDiagonalMatrixType& transform)
	{
		hamiltonian.checkValidity();
		ChangeOfBasisType::changeBasis(hamiltonian,transform);
		if (useSu2Symmetry_)
			changeBasis(reducedHamiltonian_);
	}

	void externalProduct(SizeType ind,
	                     const BasisType& basis2,
	                     const BasisType& basis3,
	                     bool order,
	                     const OperatorType& myOp)
	{
		if (!useSu2Symmetry_) return;
		SizeType n = thisBasis_->reducedSize();

		const BasisType* basisA = &basis2;
		const BasisType* basisB = &basis3;
		SizeType angularMomentum = myOp.jm.first;
		ProgramGlobals::FermionOrBosonEnum fermionOrBoson = myOp.fermionOrBoson;
		const SparseMatrixType& A = myOp.data;
		PsimagLite::Vector<SizeType>::Type jvals;

		if (!order) {
			basisA = &basis3;
			basisB = &basis2;
		}

		int ki = PsimagLite::indexOrMinusOne(momentumOfOperators_, angularMomentum);
		if (ki < 0)
			err("Operator has unknown momentum\n");

		PsimagLite::Matrix<SparseElementType> B(n,n);
		externalProd_(B, basisA, basisB, A, ki, order, fermionOrBoson);
		fullMatrixToCrsMatrix(reducedOperators_[ind].data,B);
		reducedOperators_[ind].fermionOrBoson = fermionOrBoson;
		reducedOperators_[ind].jm = myOp.jm;
		reducedOperators_[ind].angularFactor = myOp.angularFactor;
		reducedOperators_[ind].su2Related = myOp.su2Related;
	}

	void outerProductHamiltonian(const BasisType& basis2,
	                             const BasisType& basis3,
	                             const SparseMatrixType& h2,
	                             const SparseMatrixType& h3)
	{
		if (!useSu2Symmetry_) return;

		SizeType n = thisBasis_->reducedSize();
		const BasisType* basisA = &basis2;
		const BasisType* basisB = &basis3;
		PsimagLite::Matrix<SparseElementType> B2(n,n);
		externalProd_(B2,basisA,basisB,h2,-1,true,1);

		basisA = &basis3;
		basisB = &basis2;
		PsimagLite::Matrix<SparseElementType> B3(n,n);
		externalProd_(B3,basisA,basisB,h3,-1,false,1);

		B2 += B3;

		fullMatrixToCrsMatrix(reducedHamiltonian_,B2);
	}

	void reorder(SizeType,const PsimagLite::Vector<SizeType>::Type& permutation)
	{
		if (!useSu2Symmetry_) return;
		for (SizeType i=0;i<permutation.size();i++) {
			if (permutation[i]==i) continue;
			PsimagLite::String msg("ReducedOperators: reorder: permutation");
			throw PsimagLite::RuntimeError(msg + " not the identity!\n");
		}
	}

	void reorderHamiltonian(const PsimagLite::Vector<SizeType>::Type& permutation)
	{
		if (!useSu2Symmetry_) return;
		for (SizeType i=0;i<permutation.size();i++) {
			if (permutation[i]==i) continue;
			PsimagLite::String msg("ReducedOperators: reorderHamiltonian: ");
			throw PsimagLite::RuntimeError(msg + "permutation not the identity!\n");
		}
	}

	SizeType size() const { return reducedOperators_.size(); }

	void gather()
	{
		PsimagLite::MPI::pointByPointGather(reducedOperators_);
	}

	void bcast()
	{
		for (SizeType i = 0; i < reducedOperators_.size(); i++)
			Dmrg::bcast(reducedOperators_[i]);
	}

	void write(PsimagLite::IoNg::Out& io,
	           const PsimagLite::String&,
	           PsimagLite::IoNgSerializer::WriteMode mode) const
	{
		io.write(reducedOperators_, "Operators", mode);
	}

	void changeBasis(SparseMatrixType &v)
	{
		if (!useSu2Symmetry_)
			return changeOfBasis_(v);

		SparseMatrixType tmp;
		multiply(tmp,v,su2Transform_);
		multiply(v,su2TransformT_,tmp);
	}

	void clear()
	{
		momentumOfOperators_.clear();
		basisrinverse_.clear();
		reducedOperators_.clear();
		reducedHamiltonian_.clear();
		lfactorLeft_.clear();
		lfactorRight_.clear();
		lfactorHamLeft_.clear();
		lfactorHamRight_.clear();
		reducedMapping_.clear();
		fastBasisLeft_.clear();
		fastBasisRight_.clear();
		flavorIndexCached_.clear();
		changeOfBasis_.clear();
		su2Transform_.clear();
		su2TransformT_.clear();
	}

private:

	SparseElementType lfactor(int ki,
	                          bool order,
	                          SizeType j1,
	                          SizeType j2,
	                          SizeType j1prime,
	                          SizeType jProd,
	                          SizeType jProdPrime) const
	{
		SizeType ix=lfactorIndex(order,j1,j2,j1prime,jProd,jProdPrime);
		if (ki<0) {
			if (order) return lfactorHamLeft_[ix];
			return lfactorHamRight_[ix];
		}
		if (order) return lfactorLeft_[ki][ix];
		return lfactorRight_[ki][ix];
	}

	SizeType lfactorIndex(bool order,
	                      SizeType j1,
	                      SizeType j2,
	                      SizeType j1prime,
	                      SizeType jProd,
	                      SizeType jProdPrime) const
	{
		if (order) return lfactorIndexLeft(j1,j2,j1prime,jProd,jProdPrime);
		return lfactorIndexRight(j1,j2,j1prime,jProd,jProdPrime);
	}

	SizeType lfactorIndexLeft(SizeType j1,
	                          SizeType j2,
	                          SizeType j1prime,
	                          SizeType jProd,
	                          SizeType jProdPrime) const
	{
		SizeType ix = j1+j2*j1Max_+j1prime*j1Max_*j2Max_+jProd*j1Max_*j1Max_*j2Max_+
		        jProdPrime*thisBasis_->jMax()*j1Max_*j1Max_*j2Max_;
		return ix;
	}

	SizeType lfactorIndexRight(SizeType j2,
	                           SizeType j1,
	                           SizeType j2prime,
	                           SizeType jProd,
	                           SizeType jProdPrime) const
	{
		SizeType ix = j2+j1*j2Max_+j2prime*j2Max_*j1Max_+jProd*j2Max_*j2Max_*j1Max_+
		        jProdPrime*thisBasis_->jMax()*j2Max_*j2Max_*j1Max_;
		return ix;
	}

	void calcReducedMapping(const BasisType&,const BasisType&)
	{
		SizeType fMax=0;
		const PsimagLite::Vector<SizeType>::Type& flavorsOld=thisBasis_->flavorsOld();

		for (SizeType i=0;i<thisBasis_->reducedSize();i++) {
			SizeType ii = thisBasis_->reducedIndex(i);
			if (ii>=flavorsOld.size()) throw PsimagLite::RuntimeError("Ouch!!\n");
			SizeType f = flavorsOld[ii];
			if (f>fMax) fMax=f;
		}
		fMax++;
		reducedMapping_.resize(thisBasis_->jMax(), fMax, -1);
		for (SizeType i=0;i<thisBasis_->reducedSize();i++) {
			SizeType ii = thisBasis_->reducedIndex(i);
			SizeType j = thisBasis_->jmValue(ii).first;
			SizeType f = flavorsOld[ii];
			reducedMapping_(j,f)=i;
		}
	}

	void findBasisInverse()
	{
		basisrinverse_.resize(thisBasis_->size());
		for (SizeType i=0;i<thisBasis_->size();i++) {
			SizeType f = thisBasis_->getFlavor(i);
			SizeType j = thisBasis_->jmValue(i).first;
			basisrinverse_[i]=findJf(j,f);
		}
	}

	SizeType findJf(SizeType j,SizeType f)
	{
		for (SizeType i=0;i<thisBasis_->reducedSize();i++) {
			PairType jm=thisBasis_->jmValue(thisBasis_->reducedIndex(i));
			if (jm.first!=j) continue;
			if (thisBasis_->getFlavor(thisBasis_->reducedIndex(i))!=f)
				continue;
			return i;
		}
		return 0; //bogus
	}

	void createReducedConj(SizeType k1,
	                       SparseMatrixType& opDest,
	                       const SparseMatrixType& opSrc)
	{
		//			SizeType n=opSrc.rank();
		transposeConjugate(opDest,opSrc);
		for (SizeType i=0;i<opSrc.rows();i++) {
			PairType jm = thisBasis_->jmValue(thisBasis_->reducedIndex(i));
			for (int k=opDest.getRowPtr(i);k<opDest.getRowPtr(i+1);k++) {
				SizeType j=opDest.getCol(k);
				PairType jmPrime = thisBasis_->jmValue(thisBasis_->reducedIndex(j));
				SparseElementType factor=SparseElementType(jm.first+1)/
				        SparseElementType(jmPrime.first+1);
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
	                           const VectorPointerOperatorType& opSrc)
	{
		SizeType n = thisBasis_->reducedSize();
		DenseMatrixType opDest1(n,n);
		for (SizeType i=0;i<opSrc.size();i++)
			createReducedOperator(opDest1,*opSrc[i]);
		fullMatrixToCrsMatrix(opDest,opDest1);
	}

	void createReducedOperator(DenseMatrixType& opDest1,const OperatorType& opSrc)
	{
		DenseMatrixType opSrcDense;
		crsMatrixToFullMatrix(opSrcDense, opSrc.data);
		for (SizeType i=0;i<opSrc.data.rows();i++) {
			PairType jm = thisBasis_->jmValue(i);
			for (int l=opSrc.data.getRowPtr(i);l<opSrc.data.getRowPtr(i+1);l++) {
				SizeType iprime = opSrc.data.getCol(l);
				PairType jmPrime = thisBasis_->jmValue(iprime);
				RealType divisor = opSrc.angularFactor*(jmPrime.first+1);
				opDest1(basisrinverse_[i],basisrinverse_[iprime]) +=
				        opSrcDense(i,iprime)*
				        cgObject_->operator()(jmPrime,jm,opSrc.jm)/divisor;
			}
		}
	}

	void buildLfactor(VectorType& lfactor,
	                  bool order,const BasisType& basis2,
	                  const BasisType& basis3,
	                  SizeType k)
	{
		SizeType jMax=j1Max_;
		if (!order) jMax = j2Max_;

		const BasisType* basisA = &basis2;
		const BasisType* basisB = &basis3;
		if (!order) {
			basisA = &basis3;
			basisB = &basis2;
		}

		lfactor.resize(j1Max_*j2Max_*jMax*thisBasis_->jMax()*thisBasis_->jMax());
		for (SizeType i1=0;i1<basisA->jVals();i1++) {
			for (SizeType i2=0;i2<basisB->jVals();i2++) {
				for (SizeType i1prime=0;i1prime<basisA->jVals();i1prime++) {
					for (SizeType i=0;i<thisBasis_->jVals();i++) {
						for (SizeType iprime=0;iprime<thisBasis_->jVals();iprime++) {
							SizeType ix = lfactorIndex(order,basisA->jVals(i1),
							                           basisB->jVals(i2),
							                           basisA->jVals(i1prime),
							                           thisBasis_->jVals(i),
							                           thisBasis_->jVals(iprime));

							lfactor[ix]=calcLfactor(order,basisA->jVals(i1),
							                        basisB->jVals(i2),
							                        basisA->jVals(i1prime),
							                        thisBasis_->jVals(i),
							                        thisBasis_->jVals(iprime),k);
						}
					}
				}
			}
		}
	}

	SparseElementType calcLfactor(bool order,
	                              SizeType j1,
	                              SizeType j2,
	                              SizeType j1prime,
	                              SizeType jProd,
	                              SizeType jProdPrime,
	                              SizeType k) const
	{
		if (order) return 	 calcLfactorLeft(j1,j2,j1prime,jProd,jProdPrime,k);
		return calcLfactorRight(j1,j2,j1prime,jProd,jProdPrime,k);
	}

	SparseElementType calcLfactorLeft(SizeType j1,
	                                  SizeType j2,
	                                  SizeType j1prime,
	                                  SizeType jProd,
	                                  SizeType jProdPrime,
	                                  SizeType k) const
	{
		SparseElementType sum=0;
		for (SizeType m1=0;m1<=j1;m1++) {
			PairType jm1(j1,m1);
			for (SizeType m2=0;m2<=j2;m2++) {
				PairType jm2(j2,m2);
				for (SizeType mCapital=0;mCapital<=k;mCapital++) {
					PairType jmCapital(k,mCapital);
					// get m
					int x = j1+j2-jProd;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(m1 + m2 - x) < 0) continue;
					SizeType m = m1 + m2 - x;
					if (m>jProd) continue;
					PairType jm(jProd,m);

					// get m'=m+M
					x = -jProdPrime + k +jProd;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(m + mCapital - x) < 0) continue;
					SizeType mPrime = m +mCapital-x;
					if (mPrime>jProdPrime) continue;
					PairType jmPrime(jProdPrime,mPrime);

					// get m1prime=mPrime - m2
					x = j1prime + j2 - jProdPrime;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(mPrime + x - m2) < 0) continue;
					SizeType m1prime = mPrime + x - m2;
					if (m1prime > j1prime) continue;
					PairType jm1prime(j1prime,m1prime);

					// m1prime = m1 + M
					x = k - j1prime +j1;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(m1 - x + mCapital) < 0) continue;
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

	SparseElementType calcLfactorRight(SizeType j2,
	                                   SizeType j1,
	                                   SizeType j2prime,
	                                   SizeType jProd,
	                                   SizeType jProdPrime,
	                                   SizeType k) const
	{
		SparseElementType sum=0;
		for (SizeType m1=0;m1<=j1;m1++) {
			PairType jm1(j1,m1);
			for (SizeType m2=0;m2<=j2;m2++) {
				PairType jm2(j2,m2);
				for (SizeType mCapital=0;mCapital<=k;mCapital++) {
					PairType jmCapital(k,mCapital);
					// get m=m1+m2
					int x = j1+j2-jProd;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(m1 + m2 - x) < 0) continue;
					SizeType m = m1 + m2 - x;
					if (m>jProd) continue;
					PairType jm(jProd,m);

					// get m'=m+mCapital
					x = -jProdPrime + k +jProd;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(m - x + mCapital) < 0) continue;
					SizeType mPrime = m - x +mCapital;
					if (mPrime>jProdPrime) continue;
					PairType jmPrime(jProdPrime,mPrime);

					// get m2prime
					x = j2prime + j1 - jProdPrime;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(mPrime + x - m1) < 0) continue;
					SizeType m2prime = mPrime + x - m1;
					if (m2prime > j2prime) continue;
					PairType jm2prime(j2prime,m2prime);

					// check m2prime
					x = k - j2prime + j2;
					if (x%2!=0) continue;
					x/=2;
					if (static_cast<int>(m2 - x + mCapital) < 0) continue;
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

	bool checkJvalues(SizeType jProd,SizeType j1,SizeType j2) const
	{
		if (jProd>j1+j2) return false;
		if (j1>j2 && jProd<j1-j2) return false;
		if (j2>j1 && jProd<j2-j1) return false;
		return true;
	}

	void externalProd_(PsimagLite::Matrix<SparseElementType>& B,
	                   const BasisType* basisA,
	                   const BasisType* basisB,
	                   const SparseMatrixType& A,
	                   int ki,
	                   bool order,
	                   int fermionSign)
	{
		SizeType n = B.rows();
		//SizeType angularMomentum = 0;
		//if (ki>=0) angularMomentum = momentumOfOperators_[ki];
		const typename PsimagLite::Vector<VectorSizeType>::Type* fastBasis =
		        &fastBasisLeft_;
		if (!order) fastBasis = &fastBasisRight_;

		VectorSizeType basisBElectrons;
		basisB->su2ElectronsBridge(basisBElectrons);

		for (SizeType i0=0;i0<fastBasis->size();i0++) {
			const PsimagLite::Vector<SizeType>::Type& twopairs = (*fastBasis)[i0];
			SizeType i = twopairs[0];
			SizeType iprime = twopairs[1];;
			SizeType i1 = twopairs[2];
			SizeType i2 = twopairs[3];

			SizeType ii1 = basisA->reducedIndex(i1);
			SizeType j1 = basisA->jmValue(ii1).first;

			SizeType ii2 = basisB->reducedIndex(i2);
			SizeType j2 = basisB->jmValue(ii2).first;
			assert(ii2 < basisBElectrons.size());
			SizeType ne2 = basisBElectrons[ii2];

			SizeType jProd = thisBasis_->jVals(i);
			SizeType jProdPrime = thisBasis_->jVals(iprime);

			int ix = twopairs[4];
			RealType fsign = 1;
			if (!order && (ne2%2 !=0)) fsign = fermionSign;

			for (int k=A.getRowPtr(i1);k<A.getRowPtr(i1+1);k++) {
				SizeType i1prime = A.getCol(k);
				SizeType ii1prime = basisA->reducedIndex(i1prime);
				SizeType j1prime  = basisA->jmValue(ii1prime).first;

				SizeType fprime=0;
				if (order) fprime = flavorIndexCached_(i1prime,i2);
				else fprime = flavorIndexCached_(i2,i1prime);

				SparseElementType tmp = lfactor(ki,
				                                order,
				                                j1prime,
				                                j2,
				                                j1,
				                                jProdPrime,
				                                jProd);
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

	void calcFastBasis(PsimagLite::Vector<VectorSizeType>::Type& fastBasis,
	                   const BasisType& basis2,
	                   const BasisType& basis3,
	                   bool order,
	                   SizeType n)
	{
		const BasisType* basisA = &basis2;
		const BasisType* basisB = &basis3;

		if (!order) {
			basisA = &basis3;
			basisB = &basis2;
		}

		VectorSizeType basisAElectrons;
		basisA->su2ElectronsBridge(basisAElectrons);
		VectorSizeType basisBElectrons;
		basisB->su2ElectronsBridge(basisBElectrons);

		fastBasis.clear();
		PsimagLite::Vector<SizeType>::Type twopairs(5);
		for (SizeType i=0;i<thisBasis_->jVals();i++) {
			twopairs[0]=i;
			SizeType jProd = thisBasis_->jVals(i);
			for (SizeType iprime=0;iprime<thisBasis_->jVals();iprime++) {
				twopairs[1]=iprime;
				for (SizeType i1=0;i1<basisA->reducedSize();i1++) {
					twopairs[2]=i1;
					SizeType ii1 = basisA->reducedIndex(i1);
					assert(ii1 < basisAElectrons.size());
					SizeType ne1 = basisAElectrons[ii1];
					SizeType j1 = basisA->jmValue(ii1).first;
					SizeType f1 = basisA->getFlavor(ii1);

					for (SizeType i2=0;i2<basisB->reducedSize();i2++) {

						twopairs[3]=i2;
						SizeType ii2 = basisB->reducedIndex(i2);
						assert(ii2 < basisBElectrons.size());
						SizeType ne2 = basisBElectrons[ii2];
						SizeType j2 = basisB->jmValue(ii2).first;
						SizeType f2 = basisB->getFlavor(ii2);
						if (jProd>j1+j2) continue;
						if (j1>j2 && jProd<j1-j2) continue;
						if (j2>j1 && jProd<j2-j1) continue;

						SizeType f = 0;
						if (order) f = thisBasis_->flavor2Index(f1,f2,ne1,ne2,j1,j2);
						else f = thisBasis_->flavor2Index(f2,f1,ne2,ne1,j2,j1);

						// jProd,f --> ix
						int ix = reducedMapping_(jProd,f);
						if (ix<0 || SizeType(ix)>=n) continue;
						twopairs[4]=ix;

						fastBasis.push_back(twopairs);
					}
				}
			}
		}
	}

	void cacheFlavorIndex(const BasisType& basis2,const BasisType& basis3)
	{
		VectorSizeType basis2Electrons;
		basis2.su2ElectronsBridge(basis2Electrons);
		VectorSizeType basis3Electrons;
		basis2.su2ElectronsBridge(basis3Electrons);

		flavorIndexCached_.resize(basis2.reducedSize(), basis3.reducedSize());
		for (SizeType i1=0;i1<basis2.reducedSize();i1++) {
			SizeType ii1 = basis2.reducedIndex(i1);
			assert(ii1 < basis2Electrons.size());
			SizeType ne1 = basis2Electrons[ii1];
			SizeType j1 = basis2.jmValue(ii1).first;
			SizeType f1 = basis2.getFlavor(ii1);
			for (SizeType i2=0;i2<basis3.reducedSize();i2++) {
				SizeType ii2 = basis3.reducedIndex(i2);
				SizeType j2  = basis3.jmValue(ii2).first;
				SizeType f2 = basis3.getFlavor(ii2);
				assert(ii2 < basis3Electrons.size());
				SizeType ne2 = basis3Electrons[ii2];
				flavorIndexCached_(i1,i2)=thisBasis_->flavor2Index(f1,f2,ne1,ne2,j1,j2);
			}
		}
	}

	const BasisType* thisBasis_;
	bool useSu2Symmetry_;
	ClebschGordanType* cgObject_;
	PsimagLite::Vector<SizeType>::Type momentumOfOperators_;
	PsimagLite::Vector<SizeType>::Type basisrinverse_;
	typename PsimagLite::Vector<OperatorType>::Type reducedOperators_;
	SparseMatrixType reducedHamiltonian_;
	SizeType j1Max_;
	SizeType j2Max_;
	VectorVectorType lfactorLeft_;
	VectorVectorType lfactorRight_;
	VectorType lfactorHamLeft_,lfactorHamRight_;
	PsimagLite::Matrix<int> reducedMapping_;
	PsimagLite::Vector<PsimagLite::Vector<SizeType>::Type>::Type fastBasisLeft_;
	PsimagLite::Vector<PsimagLite::Vector<SizeType>::Type>::Type fastBasisRight_;
	PsimagLite::Matrix<SizeType> flavorIndexCached_;
	ChangeOfBasisType changeOfBasis_;
	SparseMatrixType su2Transform_;
	SparseMatrixType su2TransformT_;
}; // ReducedOperators
}// namespace Dmrg

/*@}*/
#endif

