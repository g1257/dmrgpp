/*
Copyright (c) 2009-2011, UT-Battelle, LLC
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

/*! \file Immm.h
 *
 *  DOC TBW FIXME
 *
 */
#ifndef IMMM_HEADER_H
#define IMMM_HEADER_H
#include "ModelBase.h"
#include "ParametersImmm.h"
#include "HilbertSpaceImmm.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "LinkProductImmm.h"
#include "ProgramGlobals.h"
#include <cassert>
#include "ModelCommon.h"

namespace Dmrg {
template<typename ModelBaseType>
class Immm : public ModelBaseType {

	typedef unsigned int long long WordType;

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef LinkProductImmm<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef  HilbertSpaceImmm<WordType> HilbertSpaceImmmType;
	typedef typename HilbertSpaceImmmType::HilbertState HilbertState;
	typedef typename PsimagLite::Vector<HilbertState>::Type HilbertBasisType;
	typedef	 typename ModelBaseType::MyBasis MyBasis;
	typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::BasisDataType BasisDataType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceImmmType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceImmmType::SPIN_DOWN;
	static const SizeType NUMBER_OF_SPINS=HilbertSpaceImmmType::NUMBER_OF_SPINS;

	enum AtomEnum {ATOM_COPPER, ATOM_OXYGEN};

	enum {ORBITALS_COPPER = 1, ORBITALS_OXYGEN = 2};

	Immm(const SolverParamsType& solverParams,
	     InputValidatorType& io,
	     GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry),
	      copperEach_(4),
	      hilbertSpace_(ORBITALS_OXYGEN)
	{
		statesCopper_ = (1<<ORBITALS_COPPER * NUMBER_OF_SPINS);
		statesOxygen_ = (1<<(ORBITALS_OXYGEN * NUMBER_OF_SPINS));

		switch (modelParameters_.minOxygenElectrons) {
		case 0:
			break;
		case 1:
			statesOxygen_--;
			break;
		case 2:
			statesOxygen_ -= 5;
			break;
		case 3:
			statesOxygen_ = 5;
			break;
		case 4:
			statesOxygen_ = 1;
			break;
		default:
			throw PsimagLite::RuntimeError("Immm: invalid param minOxygenElectrons\n");
		}
	}

	SizeType memResolv(PsimagLite::MemResolv& mres,
	                   SizeType x,
	                   PsimagLite::String msg = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType site) const
	{
		AtomEnum atom = atomAtSite(site);
		return (atom == ATOM_OXYGEN) ? statesOxygen_ : statesCopper_;
	}

	void print(std::ostream& os) const { operator<<(os,modelParameters_); }

	//! find creation operator matrices for (i,sigma) in the natural basis,
	//! find quantum numbers and number of electrons
	//! for each state in the basis
	void setNaturalBasis(typename PsimagLite::Vector<OperatorType>::Type& creationMatrix,
	                     SparseMatrixType& hamiltonian,
	                     BasisDataType& q,
	                     const BlockType& block,
	                     const RealType& time)  const
	{
		assert(block.size()==1);
		HilbertBasisType natBasis;
		typename PsimagLite::Vector<SizeType>::Type qvector;
		setNaturalBasis(natBasis,qvector,block);

		setOperatorMatrices(creationMatrix,block);

		//! Set symmetry related
		setSymmetryRelated(q,natBasis,block.size());

		//! set hamiltonian
		this->calcHamiltonian(hamiltonian,creationMatrix,block,time);
	}

	//! set creation matrices for sites in block
	void setOperatorMatrices(typename PsimagLite::Vector<OperatorType>::Type& creationMatrix,
	                         const BlockType& block) const
	{
		assert(block.size()==1);
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;
		typename PsimagLite::Vector<SizeType>::Type qvector;
		setNaturalBasis(natBasis,qvector,block);

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		creationMatrix.clear();

		SparseMatrixType nmatrix(natBasis.size(),natBasis.size());
		SparseMatrixType cTranspose;

		SizeType total = NUMBER_OF_SPINS * orbitalsAtSite(0);
		for (SizeType sigma=0;sigma<total;sigma++) {
			if (!isAllowedThisDofFull(1<<sigma,block[0])) continue;
			findOperatorMatrices(tmpMatrix,block[0],sigma,natBasis);
			typename OperatorType::Su2RelatedType su2related;

			OperatorType myOp(tmpMatrix,-1,typename OperatorType::PairType(0,0),1,su2related);
			creationMatrix.push_back(myOp);

			transposeConjugate(cTranspose,tmpMatrix);
			nmatrix += multiplyTc(cTranspose,cTranspose);
		}

		// add n_i
		typename OperatorType::Su2RelatedType su2related2;
		OperatorType nOp(nmatrix,1,typename OperatorType::PairType(0,0),1,su2related2);
		creationMatrix.push_back(nOp);
	}

	MatrixType naturalOperator(const PsimagLite::String& what,
	                           SizeType site,
	                           SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatrices(creationMatrix,block);
		SizeType orbitals = orbitalsAtSite(site);
		SizeType orbital = dof % orbitals;
		SizeType spin = dof / orbitals;
		if (what=="+" or what=="i") {
			return cDaggerCi(block,SPIN_UP,SPIN_DOWN);
		}

		if (what=="-") {
			return cDaggerCi(block,SPIN_DOWN,SPIN_UP);
		}

		if (what=="z") {
			return nUpOrDown(block,SPIN_UP)-nUpOrDown(block,SPIN_DOWN);
		}

		if (what=="n") {
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,creationMatrix[creationMatrix.size()-1].data);
			return tmp;
		}

		if (what=="c") {
			SparseMatrixType tmp2;
			transposeConjugate(tmp2,creationMatrix[orbital + spin*orbitals].data);
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,tmp2);
			return tmp;
		}

		if (what=="nup") {
			return nUpOrDown(block,SPIN_UP);
		}

		if (what=="ndown") {
			return nUpOrDown(block,SPIN_DOWN);
		}

		if (what == "o") {
			SparseMatrixType tmp2;
			transposeConjugate(tmp2,creationMatrix[orbital + spin*orbitals].data);
			SparseMatrixType tmp3 = creationMatrix[orbital + spin*orbitals].data * tmp2;
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,tmp3);
			return tmp;
		}

		if (what=="d") { // delta = c^\dagger * c^dagger
			throw PsimagLite::RuntimeError("delta unimplemented for this model\n");
		}
		std::cerr<<"what="<<what<<"\n";
		throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setNaturalBasis(HilbertBasisType& basis,
	                     typename PsimagLite::Vector<SizeType>::Type& q,
	                     const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		assert(block.size()==1);
		SizeType dof =  NUMBER_OF_SPINS * orbitalsAtSite(0);
		HilbertState total = (1<<dof);

		HilbertBasisType  basisTmp;
		for (HilbertState a=0;a<total;a++) {
			if (!isAllowedThisDof(a,block[0])) continue;
			basisTmp.push_back(a);
		}

		// reorder the natural basis (needed for MULTIPLE BANDS)
		findQuantumNumbers(q,basisTmp,block[0]);

		this->orderBasis(basis, q, basisTmp);
	}

	void findElectrons(typename PsimagLite::Vector<SizeType>::Type& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType site) const
	{
		electrons.resize(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			// nup
			SizeType nup = hilbertSpace_.electronsWithGivenSpin(basis[i],site,SPIN_UP);
			// ndown
			SizeType ndown = hilbertSpace_.electronsWithGivenSpin(basis[i],site,SPIN_DOWN);
			electrons[i] = nup + ndown;
		}
	}

	virtual SizeType maxElectronsOneSpin() const
	{
		return NUMBER_OF_SPINS * ORBITALS_OXYGEN * geometry_.numberOfSites() + 1;
	}

private:

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	RealType sign(HilbertState const &ket, SizeType site,SizeType sigma) const
	{
		int value=0;
		SizeType dof = NUMBER_OF_SPINS * orbitalsAtSite(0);
		for (SizeType alpha=0;alpha<dof;alpha++)
			value += hilbertSpace_.calcNofElectrons(ket,0,site,alpha);

		// add electron on site 0 if needed
		if (site>0) value += hilbertSpace_.electrons(ket);

		//order for sign is: up a (sigma==0), down a (sigma==2), up b (sigma==1), down b(sigma==3)
		unsigned int x = hilbertSpace_.get(ket,site);
		switch (sigma) {
		case 0:
			break;
		case 1:
			if (x & 1) value++;
			if (x & 4) value++;
			break;
		case 2:
			if (x&1) value++;
			break;
		case 3:
			if (x&1) value++;
			if (x&4) value++;
			if (x&2) value++;
			break;

		}
		if (value==0 || value%2==0) return 1.0;

		return FERMION_SIGN;
	}

	//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void findOperatorMatrices(SparseMatrixType& creationMatrix,
	                          SizeType site,
	                          SizeType sigma,
	                          const HilbertBasisType& natBasis) const
	{
		HilbertState bra,ket;
		SizeType n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<n;ii++) {
			bra=ket=natBasis[ii];

			if (hilbertSpace_.isNonZero(ket,0,sigma)) {

			} else {
				hilbertSpace_.create(bra,0,sigma);
				int jj = PsimagLite::isInVector(natBasis,bra);
				assert(jj>=0);
				assert(ii!=SizeType(jj));
				cm(ii,jj) =sign(ket,0,sigma);
			}
		}
		// here reinterpret for SU(2) if needed

		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,cm);
		transposeConjugate(creationMatrix,temp);
	}

	void findQuantumNumbers(typename PsimagLite::Vector<SizeType>::Type& q,
	                        const HilbertBasisType& basis,
	                        SizeType site) const
	{
		BasisDataType qq;
		setSymmetryRelated(qq,basis,site);
		MyBasis::findQuantumNumbers(q,qq);
	}

	void setSymmetryRelated(BasisDataType& q,
	                        const HilbertBasisType& basis,
	                        SizeType site) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		typename PsimagLite::Vector<PairType>::Type jmvalues;
		typename PsimagLite::Vector<SizeType>::Type flavors;
		PairType jmSaved = calcJmvalue<PairType>(basis[0]);
		jmSaved.first++;
		jmSaved.second++;

		typename PsimagLite::Vector<SizeType>::Type electronsUp(basis.size());
		typename PsimagLite::Vector<SizeType>::Type electronsDown(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair = calcJmvalue<PairType>(basis[i]);

			jmvalues.push_back(jmpair);

			SizeType flavor = 0; // na  + 3*nb;

			flavors.push_back(flavor);
			jmSaved = jmpair;

			// nup
			electronsUp[i] = hilbertSpace_.electronsWithGivenSpin(basis[i],site,SPIN_UP);
			// ndown
			electronsDown[i] = hilbertSpace_.electronsWithGivenSpin(basis[i],site,SPIN_DOWN);
		}
		q.jmValues=jmvalues;
		q.flavors = flavors;
		q.electrons = electronsUp + electronsDown;
		q.szPlusConst = electronsUp;
	}

	//! Not implemented, su(2) symmetry won't work
	template<typename PairType>
	PairType calcJmvalue(const HilbertState& ket) const
	{
		PairType jm(0,0);
		return jm;
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType& cm,
	                                const BlockType& block,
	                                RealType time,
	                                RealType factorForDiagonals=1.0) const
	{
		// on-site potential:
		SizeType site = block[0];
		SizeType linSize = geometry_.numberOfSites();

		SizeType siteCorrected  = 0;
		for (SizeType i=0;i<site;i++) {
			if (orbitalsAtSite(i) == ORBITALS_OXYGEN) siteCorrected++;
		}

		SizeType total = NUMBER_OF_SPINS * orbitalsAtSite(site);
		for (SizeType dof=0;dof<total;dof++) {
			SizeType norb = orbitalsAtSite(site);
			assert(norb==1 || norb==2);
			SizeType orb = dof % norb;
			SizeType siteCorrected2 = (orb==0) ? site : siteCorrected;
			SizeType index = siteCorrected2+orb*linSize;
			assert(index<modelParameters_.potentialV.size());
			SparseElementType value = modelParameters_.potentialV[index];
			SparseMatrixType tmpMatrix =value * n(cm[dof].data);
			hmatrix += factorForDiagonals * tmpMatrix;
		}

		// on-site U only for Cu sites, for now:
		if (total != NUMBER_OF_SPINS) return;
		SparseElementType tmp = factorForDiagonals * modelParameters_.hubbardU[site];
		hmatrix +=  tmp * nbar(cm[0].data) * nbar(cm[1].data);
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

		return tmpMatrix;
	}

	SparseMatrixType nbar(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,cdagger,c);

		return tmpMatrix;
	}

	void diagTest(const SparseMatrixType& fullm,const PsimagLite::String& str) const
	{
		if (fullm.rank()!=256) return;
		MatrixType fullm2;
		crsMatrixToFullMatrix(fullm2,fullm);
		typename PsimagLite::Vector<SparseElementType>::Type eigs(fullm2.n_row());
		PsimagLite::diag(fullm2,eigs,'V');
		std::cout<<str<<" diagTest size="<<fullm.rank()<<" eigs[0]="<<eigs[0]<<"\n";
		std::cout<<fullm;
	}

	MatrixType cDaggerCi(const typename PsimagLite::Vector<SizeType>::Type& block,
	                     SizeType spin1,
	                     SizeType spin2) const
	{
		assert(block.size()==1);
		SizeType site = block[0];
		SizeType norb = orbitalsAtSite(site);
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatrices(creationMatrix,block);
		assert(creationMatrix.size()>0);
		SizeType rank = creationMatrix[0].data.row();
		MatrixType tmp(rank,rank);
		assert(norb*2-1<creationMatrix.size());
		assert(spin1<2);
		assert(spin2<2);
		for (SizeType orb=0;orb<norb;orb++)
			tmp += multiplyTc(creationMatrix[orb+spin1*norb].data,
			        creationMatrix[orb+spin2*norb].data);
		return tmp;
	}

	MatrixType nUpOrDown(const typename PsimagLite::Vector<SizeType>::Type& block,
	                     SizeType spin) const
	{
		return cDaggerCi(block,spin,spin);
	}

	bool isAllowedThisDof(SizeType alpha,SizeType site) const
	{
		if (modelParameters_.minOxygenElectrons > 0)
			return isAllowedThisDofRestricted(alpha,site);
		else
			return isAllowedThisDofFull(alpha,site);
	}

	bool isAllowedThisDofFull(SizeType alpha,SizeType site) const
	{
		SizeType norb1 = orbitalsAtSite(site);
		if (norb1 == ORBITALS_OXYGEN) return true;
		return ((alpha & 10) == 0);
	}

	bool isAllowedThisDofRestricted(SizeType alpha,SizeType site) const
	{
		SizeType norb = orbitalsAtSite(site);
		if (norb == ORBITALS_COPPER)
			return ((alpha & 10) == 0);

		bool b1 = (alpha == 7 || alpha == 11 || alpha >= 13);
		bool b2 = (alpha == 3 || alpha == 5  || alpha == 9 ||
		           alpha == 6 || alpha == 10 || alpha == 12);
		switch (modelParameters_.minOxygenElectrons) {
		case 1:
			return (alpha > 0);
		case 2:
			return (b1 | b2);
		case 3:
			return b1;
		case 4:
			return (alpha == 15);
		}

		throw PsimagLite::RuntimeError("Immm: isAllowedThisDofRestricted\n");
	}

	AtomEnum atomAtSite(SizeType site) const
	{
		SizeType tmp = (site + 1) % copperEach_;
		return (tmp == 0) ? ATOM_COPPER : ATOM_OXYGEN;
	}

	SizeType orbitalsAtSite(SizeType site) const
	{
		return (atomAtSite(site) == ATOM_COPPER) ? ORBITALS_COPPER : ORBITALS_OXYGEN;
	}

	ParametersImmm<RealType>  modelParameters_;
	GeometryType const &geometry_;
	SizeType copperEach_;
	HilbertSpaceImmmType hilbertSpace_;
	SizeType statesCopper_;
	SizeType statesOxygen_;
};     //class Immm

} // namespace Dmrg
/*@}*/
#endif // IMMM_HEADER_H

