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

#ifndef IMMM_HEADER_H
#define IMMM_HEADER_H
#include "ModelBase.h"
#include "ParametersImmm.h"
#include "HilbertSpaceImmm.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include <cassert>

namespace Dmrg {
template<typename ModelBaseType>
class Immm : public ModelBaseType {

	typedef unsigned int long WordType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef  HilbertSpaceImmm<WordType> HilbertSpaceImmmType;
	typedef typename HilbertSpaceImmmType::HilbertState HilbertState;
	typedef typename PsimagLite::Vector<HilbertState>::Type HilbertBasisType;
	typedef	 typename ModelBaseType::MyBasis MyBasis;
	typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceImmmType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceImmmType::SPIN_DOWN;
	static const SizeType NUMBER_OF_SPINS=HilbertSpaceImmmType::NUMBER_OF_SPINS;

	enum AtomEnum {ATOM_COPPER, ATOM_OXYGEN};

	enum {ORBITALS_COPPER = 1, ORBITALS_OXYGEN = 2};

	Immm(const SolverParamsType& solverParams,
	     InputValidatorType& io,
	     SuperGeometryType const &geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      superGeometry_(geometry),
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
			err("Immm: invalid param minOxygenElectrons\n");
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		io.write(label + "/copperEach_", copperEach_);
		hilbertSpace_.write(label, io);
		io.write(label + "/statesCopper_", statesCopper_);
		io.write(label + "/statesOxygen_", statesOxygen_);
	}

	virtual SizeType maxElectronsOneSpin() const
	{
		return NUMBER_OF_SPINS * ORBITALS_OXYGEN * superGeometry_.numberOfSites() + 1;
	}

	SizeType differentTypesOfAtoms() const { return 2; }

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		for (SizeType site = 0; site < differentTypesOfAtoms(); ++site) {
			OpsLabelType& c = this->createOpsLabel("c", site);
			OpsLabelType& nop = this->createOpsLabel("nop", site);
			OpsLabelType& splus = this->createOpsLabel("splus", site);
			OpsLabelType& sminus = this->createOpsLabel("sminus", site);
			OpsLabelType& sz = this->createOpsLabel("sz", site);
			OpsLabelType& o = this->createOpsLabel("o", site);
			OpsLabelType& nup = this->createOpsLabel("nup", site);
			OpsLabelType& ndown = this->createOpsLabel("ndown", site);

			this->makeTrackable("c");
			this->makeTrackable("nop");

			BlockType block(1, site);
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			setOperatorMatricesInternal(creationMatrix, qns, block);
			//		SizeType orbitals = orbitalsAtSite(site);
			//		SizeType orbital = dof % orbitals;
			//		SizeType spin = dof / orbitals;

			SizeType total = creationMatrix.size();

			splus.push(cDaggerCi(block,SPIN_UP,SPIN_DOWN));

			sminus.push(cDaggerCi(block,SPIN_DOWN,SPIN_UP));

			{ // S^z
				MatrixType tmp1;
				crsMatrixToFullMatrix(tmp1,nUpOrDown(block,SPIN_UP).getCRS());
				MatrixType tmp2;
				crsMatrixToFullMatrix(tmp2,nUpOrDown(block,SPIN_DOWN).getCRS());
				tmp1 -= tmp2;
				SparseMatrixType tmp(tmp1);
				typename OperatorType::Su2RelatedType su2Related;
				typename OperatorType::PairType pairZero(0, 0);
				sz.push(OperatorType(tmp,
				                     ProgramGlobals::FermionOrBosonEnum::BOSON,
				                     pairZero,
				                     1.0,
				                     su2Related));
			}

			for (SizeType dof = 0; dof < total; ++dof)
				c.push(creationMatrix[dof]);

			assert(creationMatrix.size() > 0);
			nop.push(creationMatrix[creationMatrix.size()-1]);

			nup.push(nUpOrDown(block,SPIN_UP));

			ndown.push(nUpOrDown(block,SPIN_DOWN));

			for (SizeType dof = 0; dof < total; ++dof) {
				SparseMatrixType tmp2;
				transposeConjugate(tmp2,creationMatrix[dof].getCRS());
				SparseMatrixType tmp3 = creationMatrix[dof].getCRS() * tmp2;
				typename OperatorType::Su2RelatedType su2Related;
				o.push(OperatorType(tmp3,
				                    ProgramGlobals::FermionOrBosonEnum::BOSON,
				                    typename OperatorType::PairType(0,0),
				                    1.0,
				                    su2Related));
			}
		}
	}

	void fillModelLinks()
	{
		err("Immm is broken\n");
	}

private:

	//! set creation matrices for sites in block
	void setOperatorMatricesInternal(VectorOperatorType& creationMatrix,
	                                 VectorQnType& qns,
	                                 const BlockType& block) const
	{
		assert(block.size()==1);
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;
		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis, block[0]);

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		creationMatrix.clear();

		MatrixType nmatrix(natBasis.size(),natBasis.size());

		SizeType total = NUMBER_OF_SPINS * orbitalsAtSite(0);
		for (SizeType sigma=0;sigma<total;sigma++) {
			if (!isAllowedThisDofFull(1<<sigma,block[0])) continue;
			findOperatorMatrices(tmpMatrix,block[0],sigma,natBasis);
			typename OperatorType::Su2RelatedType su2related;

			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::FERMION,
			                  typename OperatorType::PairType(0, 0),
			                  1,
			                  su2related);
			creationMatrix.push_back(myOp);

			nmatrix += multiplyTc(tmpMatrix,tmpMatrix);
		}

		// add n_i
		typename OperatorType::Su2RelatedType su2related2;
		OperatorType nOp(SparseMatrixType(nmatrix),
		                 ProgramGlobals::FermionOrBosonEnum::BOSON,
		                 typename OperatorType::PairType(0,0),
		                 1
		                 ,su2related2);
		creationMatrix.push_back(nOp);
	}

	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		assert(block.size()==1);
		SizeType dof =  NUMBER_OF_SPINS * orbitalsAtSite(0);
		HilbertState total = (1<<dof);

		basis.clear();
		for (HilbertState a = 0; a < total; ++a) {
			if (!isAllowedThisDof(a,block[0])) continue;
			basis.push_back(a);
		}
	}

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
	                          SizeType,
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
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
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

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        SizeType site) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		VectorSizeType other(2, 0);
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = calcJmvalue<PairType>(basis[i]);
			SizeType flavor = 0; // na  + 3*nb;
			// nup
			other[1] = hilbertSpace_.electronsWithGivenSpin(basis[i],site,SPIN_UP);
			// ndown
			SizeType electronsDown = hilbertSpace_.electronsWithGivenSpin(basis[i],site,SPIN_DOWN);
			other[0] = electronsDown + other[1];
			bool sign = other[0] & 1;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	//! Not implemented, su(2) symmetry won't work
	template<typename PairType>
	PairType calcJmvalue(const HilbertState&) const
	{
		PairType jm(0,0);
		return jm;
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType time) const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		// on-site potential:
		SizeType site = block[0];
		SizeType linSize = superGeometry_.numberOfSites();

		SizeType siteCorrected  = 0;
		for (SizeType i=0;i<site;i++) {
			if (orbitalsAtSite(i) == ORBITALS_OXYGEN) siteCorrected++;
		}

		SizeType total = NUMBER_OF_SPINS * orbitalsAtSite(site);
		for (SizeType dof=0;dof<total;dof++) {
			const OperatorType& cmop = ModelBaseType::naturalOperator("c", site, dof);
			const SparseMatrixType& cm = cmop.getCRS();
			SizeType norb = orbitalsAtSite(site);
			assert(norb==1 || norb==2);
			SizeType orb = dof % norb;
			SizeType siteCorrected2 = (orb==0) ? site : siteCorrected;
			SizeType index = siteCorrected2+orb*linSize;
			assert(index<modelParameters_.potentialV.size());
			SparseElementType value = modelParameters_.potentialV[index];
			SparseMatrixType tmpMatrix =value * n(cm);
			hmatrix += tmpMatrix;
		}

		// on-site U only for Cu sites, for now:
		if (total != NUMBER_OF_SPINS) return;

		const OperatorType& cup = ModelBaseType::naturalOperator("c", site, 0);
		const OperatorType& cdown = ModelBaseType::naturalOperator("c", site, 1);
		hmatrix +=  modelParameters_.hubbardU[site] * nbar(cup.getCRS()) * nbar(cdown.getCRS());
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

	OperatorType cDaggerCi(const typename PsimagLite::Vector<SizeType>::Type& block,
	                       SizeType spin1,
	                       SizeType spin2) const
	{
		assert(block.size()==1);
		SizeType site = block[0];
		SizeType norb = orbitalsAtSite(site);
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		VectorQnType qns;
		setOperatorMatricesInternal(creationMatrix, qns, block);
		assert(creationMatrix.size()>0);
		SizeType rank = creationMatrix[0].getCRS().rows();
		MatrixType tmp(rank,rank);
		assert(norb*2-1<creationMatrix.size());
		assert(spin1<2);
		assert(spin2<2);
		for (SizeType orb=0;orb<norb;orb++)
			tmp += multiplyTc(creationMatrix[orb+spin1*norb].getCRS(),
			        creationMatrix[orb+spin2*norb].getCRS());

		typename OperatorType::Su2RelatedType su2Related;
		return OperatorType(SparseMatrixType(tmp),
		                    ProgramGlobals::FermionOrBosonEnum::BOSON,
		                    typename OperatorType::PairType(0,0),
		                    1.0,
		                    su2Related);
	}

	OperatorType nUpOrDown(const typename PsimagLite::Vector<SizeType>::Type& block,
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

	ParametersImmm<RealType, QnType> modelParameters_;
	const SuperGeometryType& superGeometry_;
	SizeType copperEach_;
	HilbertSpaceImmmType hilbertSpace_;
	SizeType statesCopper_;
	SizeType statesOxygen_;
}; //class Immm

} // namespace Dmrg
/*@}*/
#endif // IMMM_HEADER_H

