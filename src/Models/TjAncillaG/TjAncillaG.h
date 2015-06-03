/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file TjAncillaG.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef DMRG_TJ_ANCILLAG_H
#define DMRG_TJ_ANCILLAG_H
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "../Models/Tj1Orb/Tj1Orb.h"
#include "../Models/TjAncillaG/LinkProductTjAncillaG.h"
#include "../Models/Tj1Orb/ParametersModelTj1Orb.h"
#include "ModelCommon.h"

namespace Dmrg {
//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
template<typename ModelBaseType>
class TjAncillaG : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef ModelHubbard<ModelBaseType> ModelHubbardType;
	typedef Tj1Orb<ModelBaseType> Tj1OrbType;
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
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef LinkProductTjAncillaG<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ModelHubbardType::HilbertState HilbertStateType;
	typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelHubbardType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename OperatorType::PairType PairType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename PsimagLite::Vector<HilbertStateType>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	static const int DEGREES_OF_FREEDOM=2;
	static const int NUMBER_OF_ORBITALS = 1;
	static const int FERMION_SIGN = -1;

	enum {SPIN_UP, SPIN_DOWN};

	TjAncillaG(const SolverParamsType& solverParams,
	       InputValidatorType& io,
	       GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry),
	      tj1orb_(solverParams,io,geometry),
	      offset_(DEGREES_OF_FREEDOM+3), // c^\dagger_up, c^\dagger_down, S+, Sz, n
	      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM)
	{}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType site) const
	{
		return tj1orb_.hilbertSize(site);
	}

	//! find creation operator matrices for (i,sigma) in the natural basis,
	//! find quantum numbers and number of electrons
	//! for each state in the basis
	virtual void setNaturalBasis(VectorOperatorType& creationMatrix,
	                             SparseMatrixType &hamiltonian,
	                             SymmetryElectronsSzType& q,
	                             const BlockType& block,
	                             const RealType& time) const
	{

		HilbertBasisType natBasis;
		VectorSizeType quantumNumbs;
		tj1orb_.setNaturalBasis(natBasis,quantumNumbs,block);

		setOperatorMatrices(creationMatrix,block);

		//! Set symmetry related
		setSymmetryRelated(q,natBasis,block.size());

		//! set hamiltonian
		this->calcHamiltonian(hamiltonian,creationMatrix,block,time);
	}

	virtual void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                                 const BlockType& block) const
	{
		return tj1orb_.setOperatorMatrices(creationMatrix,block);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType&,
	                                const BlockType& block,
	                                RealType,
	                                RealType factorForDiagonals=1.0) const
	{
		SizeType n=block.size();

		SizeType linSize = geometry_.numberOfSites();
		for (SizeType i=0;i<n;i++) {
			// potentialV
			SparseMatrixType nup(naturalOperator("nup",i,0));
			SparseMatrixType ndown(naturalOperator("ndown",i,0));
			SparseMatrixType m = nup;
			assert(block[i]+linSize<modelParameters_.potentialV.size());
			m *= modelParameters_.potentialV[block[i]];
			m += modelParameters_.potentialV[block[i]+linSize]*ndown;
			hmatrix += factorForDiagonals * m;
		}
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setNaturalBasis(HilbertBasisType& basis,
	                     VectorSizeType& q,
	                     const VectorSizeType& block) const
	{
		assert(block.size()==1);
		HilbertStateType a=0;
		int sitesTimesDof=DEGREES_OF_FREEDOM;
		HilbertStateType total = (1<<sitesTimesDof);
		total--;

		HilbertBasisType  basisTmp;
		for (a=0;a<total;a++) basisTmp.push_back(a);

		// reorder the natural basis (needed for MULTIPLE BANDS)
		findQuantumNumbers(q,basisTmp,1);
		VectorSizeType iperm(q.size());

		PsimagLite::Sort<VectorSizeType > sort;
		sort.sort(q,iperm);
		basis.clear();
		for (a=0;a<total;a++) basis.push_back(basisTmp[iperm[a]]);
	}

	/** \cppFunction{!PTEX_THISFUNCTION} returns the operator in the
	 *unmangled (natural) basis of one-site */
	MatrixType naturalOperator(const PsimagLite::String& what,
	                                                      SizeType site,
	                                                      SizeType dof) const
	{
		return tj1orb_.naturalOperator(what,site,dof);
	}

	//! find total number of electrons for each state in the basis
	void findElectrons(VectorSizeType& electrons,
	                   const VectorHilbertStateType& basis,
	                   SizeType other) const
	{
		tj1orb_.findElectrons(electrons,basis,other);
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

	void print(std::ostream& os) const
	{
		os<<modelParameters_;
	}

private:

	void findQuantumNumbers(VectorSizeType& q,
	                        const HilbertBasisType& basis,int n) const
	{
		SymmetryElectronsSzType qq;
		setSymmetryRelated(qq,basis,n);
		qq.findQuantumNumbers(q, MyBasis::useSu2Symmetry());
	}

	void setSymmetryRelated(SymmetryElectronsSzType& q,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		assert(n==1);

		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		typename PsimagLite::Vector<PairType>::Type jmvalues;
		VectorSizeType flavors;
		PairType jmSaved = calcJmvalue<PairType>(basis[0]);
		jmSaved.first++;
		jmSaved.second++;

		VectorSizeType szPlusConst(basis.size());
		VectorSizeType bogus(basis.size(),0);
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair = calcJmvalue<PairType>(basis[i]);

			jmvalues.push_back(jmpair);
			// nup
			SizeType ups = HilbertSpaceHubbardType::getNofDigits(basis[i],SPIN_UP);
			// ndown
			SizeType downs = HilbertSpaceHubbardType::getNofDigits(basis[i],SPIN_DOWN);

			flavors.push_back(0);
			jmSaved = jmpair;
			szPlusConst[i] = (ups + geometry_.numberOfSites()) - downs;
		}

		q.set(jmvalues,flavors,szPlusConst,bogus);
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// Reinterprets 6 and 9
	template<typename PairType>
	PairType calcJmvalue(const HilbertStateType& ket) const
	{
		return calcJmValueAux<PairType>(ket);
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// does not work for 6 or 9
	template<typename PairType>
	PairType calcJmValueAux(const HilbertStateType& ket) const
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

	//serializr start class TjAncillaG
	//serializr vptr
	//serializr normal modelParameters_
	ParametersModelTj1Orb<RealType>  modelParameters_;
	//serializr ref geometry_ end
	const GeometryType &geometry_;
	Tj1OrbType tj1orb_;
	//serializr normal offset_
	SizeType offset_;
	//serializr normal spinSquaredHelper_
	SpinSquaredHelper<RealType,HilbertStateType> spinSquaredHelper_;
	//serializr normal spinSquared_
	SpinSquared<SpinSquaredHelper<RealType,HilbertStateType> > spinSquared_;

};	//class TjAncillaG

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_ANCILLAG_H

