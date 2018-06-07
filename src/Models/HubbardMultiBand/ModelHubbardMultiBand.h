/*
Copyright (c) 2009-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

/*! \file ModelHubbardMultiBand.h
 *
 *  TBW
 *
 */
#ifndef MODEL_HUBBARD_MULTI_BAND_H
#define MODEL_HUBBARD_MULTI_BAND_H
#include "ModelBase.h"
#include "ParametersHubbardMultiBand.h"
#include "../Models/FeAsModel/HilbertSpaceFeAs.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "../Models/FeAsModel/LinkProductFeAs.h"
#include "ProgramGlobals.h"
#include "ModelCommon.h"

namespace Dmrg {
template<typename ModelBaseType>
class ModelHubbardMultiBand : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef TargetQuantumElectrons<RealType> TargetQuantumElectronsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef  HilbertSpaceFeAs<HilbertState> HilbertSpaceFeAsType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef LinkProductFeAs<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef	 typename ModelBaseType::MyBasis MyBasis;
	typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ParametersHubbardMultiBand<ComplexOrRealType> ParamsModelType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

	ModelHubbardMultiBand(const SolverParamsType& solverParams,
	                      InputValidatorType& io,
	                      GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry),
	      spinSquared_(spinSquaredHelper_,
	                   modelParameters_.orbitals,
	                   2*modelParameters_.orbitals)
	{
		ProgramGlobals::init(modelParameters_.orbitals*geometry_.numberOfSites() + 1);
		SizeType v1 = 2*modelParameters_.orbitals*geometry.numberOfSites();
		SizeType v2 = v1*modelParameters_.orbitals;
		if (modelParameters_.potentialV.size() != v1 &&
		        modelParameters_.potentialV.size() != v2) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "potentialV length must be 2*orbitals times the number of sites or";
			str += " 2*orbitals*orbitals times the number of sites\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		LinkProductType::setOrbitals(modelParameters_.orbitals);
		HilbertSpaceFeAsType::setOrbitals(modelParameters_.orbitals);

		VectorSizeType block(1,0);
		int sitesTimesDof=2*modelParameters_.orbitals;
		HilbertState total = (1<<sitesTimesDof);
		basis_.resize(total);
		for (HilbertState a = 0; a < total; ++a)
			basis_[a] = a;

		HilbertBasisType basisTmp = basis_;
		setSymmetryRelatedInternal(qq_, basis_, 1);
		ModelBaseType::orderBasis(basis_, basisTmp, qq_);

		setOperatorMatricesInternal(creationMatrix_, block);

		cacheInteractionOp();
	}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType) const
	{
		return basis_.size();
	}

	void print(std::ostream& os) const { operator<<(os,modelParameters_); }

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         const BlockType&) const
	{
		creationMatrix = creationMatrix_;
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		SizeType orbital = dof % modelParameters_.orbitals;
		SizeType spin = dof/modelParameters_.orbitals;
		assert(creationMatrix_.size()>0);
		SizeType nrow = creationMatrix_[0].data.rows();
		PsimagLite::String what2 = what;

		if (what2 == "splus") {
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < modelParameters_.orbitals; ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp(nrow,nrow);
			tmp += multiplyTc(creationMatrix_[dof].data,
			                  creationMatrix_[dof+modelParameters_.orbitals].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp2,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what2 == "sminus") {
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < modelParameters_.orbitals; ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp(nrow,nrow);
			tmp += multiplyTc(creationMatrix_[dof+modelParameters_.orbitals].data,
			        creationMatrix_[dof].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp2,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what2 == "z" || what2 == "sz") { // S^z
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < modelParameters_.orbitals; ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp(nrow,nrow);
			MatrixType tmp2(nrow,nrow);

			tmp += multiplyTc(creationMatrix_[dof].data,creationMatrix_[dof].data);
			tmp2 += multiplyTc(creationMatrix_[dof+modelParameters_.orbitals].data,
			        creationMatrix_[dof+modelParameters_.orbitals].data);

			tmp = 0.5*(tmp-tmp2);
			SparseMatrixType tmp3(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp3,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}
		if (what2=="n") {
			VectorSizeType allowed(2*modelParameters_.orbitals,0);
			for (SizeType x = 0; x < allowed.size(); ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp =
			        multiplyTc(creationMatrix_[dof].data,creationMatrix_[dof].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp2,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what2=="c") {
			VectorSizeType allowed(2*modelParameters_.orbitals,0);
			for (SizeType x = 0; x < allowed.size(); ++x) allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			VectorOperatorType cm = creationMatrix_;
			cm[orbital + spin*modelParameters_.orbitals].conjugate();
			return cm[orbital + spin*modelParameters_.orbitals];
		}

		if (what2=="d") { // delta = c^\dagger * c^dagger
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < allowed.size(); ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			SparseMatrixType atmp;
			multiply(atmp,
			         creationMatrix_[orbital+orbital+modelParameters_.orbitals].data,
			        creationMatrix_[orbital].data);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(atmp,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		PsimagLite::String str("ModelHubbardMultiBand: naturalOperator: no label ");
		str += what + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setBasis(HilbertBasisType& basis,
	              SymmetryElectronsSzType& qq,
	              const VectorSizeType&) const
	{
		basis = basis_;
		qq = qq_;
	}

	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType) const
	{
		electrons.resize(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			// nup
			SizeType nup = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                            SPIN_UP);
			// ndown
			SizeType ndown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                              SPIN_DOWN);
			electrons[i] = nup + ndown;
		}
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType& cm,
	                                const BlockType& block,
	                                RealType time,
	                                RealType factorForDiagonals=1.0) const
	{
		SizeType n=block.size();

		for (SizeType i = 0; i < n; ++i) {

			addHoppingOnSite(hmatrix,cm,i,factorForDiagonals,block[i]);

			addInteraction(hmatrix,cm,i,factorForDiagonals,block[i]);

			addPotentialV(hmatrix,
			              cm,
			              i,
			              block[i],
			              factorForDiagonals,
			              modelParameters_.potentialV);

		}
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

private:

	void setOperatorMatricesInternal(VectorOperatorType& creationMatrix,
	                                 const BlockType& block) const
	{
		HilbertBasisType natBasis = basis_;
		SparseMatrixType tmpMatrix;

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		creationMatrix.clear();
		SizeType dofs = 2*modelParameters_.orbitals;
		for (SizeType i=0;i<block.size();i++) {
			for (SizeType sigma=0;sigma<dofs;sigma++) {
				findOperatorMatrices(tmpMatrix,i,sigma,natBasis);

				SizeType m=0;
				int asign=1;
				if (sigma>modelParameters_.orbitals-1) {
					m=1;
					asign= -1;
				}
				typename OperatorType::Su2RelatedType su2related;
				if (sigma <modelParameters_.orbitals) {
					su2related.source.push_back(i*dofs+sigma);
					su2related.source.push_back(i*dofs+sigma + modelParameters_.orbitals);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = modelParameters_.orbitals;
				}

				OperatorType myOp(tmpMatrix,
				                  -1,
				                  typename OperatorType::PairType(1,m),
				                  asign,
				                  su2related);
				creationMatrix.push_back(myOp);
			}
		}
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to
	//! basis state ket
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	RealType sign(HilbertState const &ket, int i,SizeType sigma) const
	{
		int value=0;
		SizeType dofs=2*modelParameters_.orbitals;
		for (SizeType alpha=0;alpha<dofs;alpha++)
			value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);
		// add electron on site 0 if needed
		if (i>0) value += HilbertSpaceFeAsType::electrons(ket);

		//order for sign is: a up, b up, a down, b down, etc
		unsigned int x = HilbertSpaceFeAsType::get(ket,i);
		int spin = sigma/modelParameters_.orbitals;
		SizeType orb = sigma % modelParameters_.orbitals;

		for (SizeType j=0;j<orb;j++) {
			for (SizeType k=0;k<2;k++) {
				SizeType ind = j + k * modelParameters_.orbitals;
				int mask = (1<<ind);
				if (x & mask) value++;
			}
		}
		if (spin==SPIN_DOWN) {
			int mask = (1<<orb);
			if (x & mask) value++;
		}
		if (value==0 || value%2==0) return 1.0;

		return FERMION_SIGN;
	}

	//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void findOperatorMatrices(SparseMatrixType& creationMatrix,
	                          int i,
	                          int sigma,
	                          const HilbertBasisType& natBasis) const
	{
		int n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertState ket = natBasis[ii];
			HilbertState bra = ket;

			if (HilbertSpaceFeAsType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceFeAsType::create(bra,i,sigma);
				int jj = PsimagLite::isInVector(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices: internal error\n");
				if (ii==SizeType(jj)) {
					std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra<<" sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("Creation operator cannot be diagonal\n");
				}

				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,cm);
		transposeConjugate(creationMatrix,temp);
	}

	void setSymmetryRelatedInternal(SymmetryElectronsSzType& q,
	                                const HilbertBasisType& basis,
	                                int n) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		typename PsimagLite::Vector<PairType>::Type jmvalues;
		VectorSizeType flavors;

		VectorSizeType electronsUp(basis.size());
		VectorSizeType electrons(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair(0,0);

			jmvalues.push_back(jmpair);

			flavors.push_back(0.0);

			// nup
			electronsUp[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                              SPIN_UP);
			// ndown
			SizeType electronsDown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                                      SPIN_DOWN);
			electrons[i] = electronsDown + electronsUp[i];
		}

		q.set(jmvalues,flavors,electrons,electronsUp);
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// does not work for 6 or 9
	template<typename PairType>
	PairType calcJmvalue(const HilbertState& ket) const
	{
		SizeType site0=0;
		SizeType site1=0;

		spinSquared_.doOnePairOfSitesA(ket,site0,site1);
		spinSquared_.doOnePairOfSitesB(ket,site0,site1);
		spinSquared_.doDiagonal(ket,site0,site1);

		RealType sz = spinSquared_.spinZ(ket,site0);
		PairType jm= spinSquaredHelper_.getJmPair(sz);
		return jm;
	}

	void addInteraction(SparseMatrixType &hmatrix,
	                    const VectorOperatorType&,
	                    SizeType,
	                    RealType factorForDiagonals,
	                    SizeType actualSite) const
	{
		assert(actualSite < modelParameters_.hubbardU.size());
		hmatrix += factorForDiagonals*modelParameters_.hubbardU[actualSite]*qx_;
	}

	void addHoppingOnSite(SparseMatrixType& hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      RealType factorForDiagonals,
	                      SizeType actualSite) const
	{
		SizeType orbitals = modelParameters_.orbitals;
		SizeType dof = 2*orbitals;
		SparseMatrixType tmpMatrix;
		SparseMatrixType tmp;
		for (SizeType orb1 = 0; orb1 < orbitals; ++orb1) {
			const SparseMatrixType& m1up = cm[orb1 + SPIN_UP*orbitals + i*dof].data;
			const SparseMatrixType& m1down = cm[orb1 + SPIN_DOWN*orbitals + i*dof].data;

			for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
				ComplexOrRealType val = getOnSiteHopping(actualSite, orb1, orb2);
				if (val == 0.0) continue;

				const SparseMatrixType& m2up = cm[orb2 + SPIN_UP*orbitals + i*dof].data;
				SparseMatrixType m2t;
				transposeConjugate(m2t, m2up);
				multiply(tmpMatrix, m1up, m2t);
				tmp += val*tmpMatrix;

				const SparseMatrixType& m2down = cm[orb2 + SPIN_DOWN*orbitals + i*dof].data;
				transposeConjugate(m2t, m2down);
				multiply(tmpMatrix, m1down, m2t);

				tmp += val*tmpMatrix;
			}
		}

		hmatrix += factorForDiagonals*tmp;
	}

	ComplexOrRealType getOnSiteHopping(SizeType actualSite, SizeType orb1, SizeType orb2) const
	{
		const typename PsimagLite::Vector<PsimagLite::Matrix<ComplexOrRealType> >::Type& v =
		        modelParameters_.hopOnSite;
		if (v.size() == 0)
			err("getOnSiteHopping");
		if (v.size() > 1 && actualSite >= v.size())
			err("getOnSiteHopping too small\n");
		return (v.size() == 1) ? v[0](orb1, orb2) : v[actualSite](orb1, orb2);
		}

		void addPotentialV(SparseMatrixType &hmatrix,
		const VectorOperatorType& cm,
		SizeType i,
		SizeType actualIndexOfSite,
		RealType factorForDiagonals,
		const typename PsimagLite::Vector<RealType>::Type& V) const
		{
		SizeType v1 = 2*modelParameters_.orbitals*geometry_.numberOfSites();
		SizeType v2 = v1*modelParameters_.orbitals;
		if (V.size() != v1 && V.size() != v2) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "potentialV[T] length must be 2*orbitals times the number of sites or";
			str += " 2*orbitals*orbitals times the number of sites\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		if (V.size() == v1) {
			for (SizeType orb=0;orb<modelParameters_.orbitals;orb++)
				addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb,factorForDiagonals,V);
		}

		if (V.size() == v2) {
			for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {
				for (SizeType orb2=0;orb2<modelParameters_.orbitals;orb2++) {
					addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb,orb2,factorForDiagonals,V);
				}
			}

			return;
		}
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorOperatorType& cm,
	                   SizeType i,
	                   SizeType actualIndexOfSite,
	                   SizeType orbital,
	                   RealType factorForDiagonals,
	                   const typename PsimagLite::Vector<RealType>::Type& V) const
	{
		int dof = 2*modelParameters_.orbitals;
		SparseMatrixType nup = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data);
		SparseMatrixType ndown = n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data);

		SizeType linSize = geometry_.numberOfSites();

		SizeType iUp = actualIndexOfSite + (orbital + 0*modelParameters_.orbitals)*linSize;
		hmatrix += factorForDiagonals * V[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orbital + 1*modelParameters_.orbitals)*linSize;
		hmatrix += factorForDiagonals * V[iDown] * ndown;
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorOperatorType& cm,
	                   SizeType i,
	                   SizeType actualIndexOfSite,
	                   SizeType orb,
	                   SizeType orb2,
	                   RealType factorForDiagonals,
	                   const typename PsimagLite::Vector<RealType>::Type& V) const
	{
		int dof=2*modelParameters_.orbitals;
		SizeType orbitalsSquared = modelParameters_.orbitals*modelParameters_.orbitals;

		SparseMatrixType nup = nEx(cm[orb+SPIN_UP*modelParameters_.orbitals+i*dof].data,
		        cm[orb2+SPIN_UP*modelParameters_.orbitals+i*dof].data);
		SparseMatrixType ndown = nEx(cm[orb+SPIN_DOWN*modelParameters_.orbitals+i*dof].data,
		        cm[orb2+SPIN_DOWN*modelParameters_.orbitals+i*dof].data);

		SizeType linSize = geometry_.numberOfSites();

		SizeType iUp = actualIndexOfSite + (orb + orb2*modelParameters_.orbitals +
		                                    0*orbitalsSquared)*linSize;
		hmatrix += factorForDiagonals * V[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orb + orb2*modelParameters_.orbitals +
		                                      1*orbitalsSquared)*linSize;
		hmatrix += factorForDiagonals * V[iDown] * ndown;
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

		return tmpMatrix;
	}

	SparseMatrixType nEx(const SparseMatrixType& c1,
	                     const SparseMatrixType& c2) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c2);
		multiply(tmpMatrix,c1,cdagger);

		return tmpMatrix;
	}

	void cacheInteractionOp()
	{
		VectorSizeType block(1, 0);
		VectorOperatorType cm;
		setOperatorMatricesInternal(cm, block);

		SparseMatrixType m0m1;
		SparseMatrixType m2m3;
		SparseMatrixType m2m3t;
		SparseMatrixType term;

		SizeType orbitals = modelParameters_.orbitals;
		for (SizeType k0 = 0; k0 < orbitals; ++k0) {
			const SparseMatrixType& m0 = cm[k0 + SPIN_DOWN*orbitals].data;
			for (SizeType k1 = 0; k1 < orbitals; ++k1) {
				const SparseMatrixType& m1 = cm[k1 + SPIN_UP*orbitals].data;
				multiply(m0m1, m0, m1);
				for (SizeType k2 = 0; k2 < orbitals; ++k2) {
					const SparseMatrixType& m2 = cm[k2 + SPIN_UP*orbitals].data;
					SizeType k3 = (k0 + k1 + k2) % orbitals;

					const SparseMatrixType& m3 = cm[k3 + SPIN_DOWN*orbitals].data;
					multiply(m2m3, m3, m2);
					transposeConjugate(m2m3t, m2m3);
					multiply(term, m2m3t, m0m1);
					qx_ += 0.5*term;
				}
			}
		}
	}

	ParamsModelType  modelParameters_;
	const GeometryType& geometry_;
	SpinSquaredHelper<RealType,HilbertState> spinSquaredHelper_;
	SpinSquared<SpinSquaredHelper<RealType,HilbertState> > spinSquared_;
	HilbertBasisType basis_;
	SymmetryElectronsSzType qq_;
	VectorSizeType q_;
	VectorOperatorType creationMatrix_;
	SparseMatrixType qx_;
}; //class ModelHubbardMultiBand
} // namespace Dmrg
/*@}*/
#endif

