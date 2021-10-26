/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

/*! \file FermionSpinless.h
 *
 *  TBW
 *
 */
#ifndef DMRG_FERMION_SPINLESS
#define DMRG_FERMION_SPINLESS
#include <cassert>
#include "Sort.h" // in PsimagLite
#include "ParametersFermionSpinless.h"
#include "HilbertSpaceFermionSpinless.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"

namespace Dmrg {
//! Model Hubbard for DMRG solver, inherits from ModelBase and implements its interface:
template<typename ModelBaseType>
class FermionSpinless : public ModelBaseType {

	static const int FERMION_SIGN = -1;
	static const int DEGREES_OF_FREEDOM=1;
	static const int NUMBER_OF_ORBITALS=1;

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef unsigned int long WordType;
	typedef  HilbertSpaceFermionSpinless<WordType> HilbertSpaceType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename HilbertSpaceType::HilbertState HilbertState;
	typedef typename PsimagLite::Vector<HilbertState>::Type VectorHilbertStateType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename PsimagLite::Vector<HilbertState>::Type HilbertBasisType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	FermionSpinless(const SolverParamsType& solverParams,
	                InputValidatorType& io,
	                const SuperGeometryType& geometry,
	                PsimagLite::String extra)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      offset_(DEGREES_OF_FREEDOM),
	      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM),
	      hasDelta_(extra == "WithDelta")
	{
		if (extra != "" && extra != "WithDelta")
			err("FermionSpinLess can only be followed by WithDelta and not " + extra + "\n");

		assert(geometry.numberOfSites() == modelParameters_.potentialV.size());
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		io.write(label + "/offset_", offset_);
		spinSquaredHelper_.write(label, io);
		spinSquared_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType time)  const
	{
		SizeType n = block.size();
		if (n != 1)
			err("addDiagonalsInNaturalBasis: block.size() != 1\n");

		const SizeType site = block[0];
		const OperatorType& niupop = ModelBaseType::naturalOperator("n", site, 0);
		const SparseMatrixType& niup = niupop.getCRS();

		// V_iup term
		RealType tmp = modelParameters_.potentialV[site];
		hmatrix += tmp*niup;
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;
		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis);

		//! Set the operators c^\daggger_{i\sigma} in the natural basis
		for (SizeType i=0;i<block.size();i++) {
			for (int sigma=0;sigma<DEGREES_OF_FREEDOM;sigma++) {
				tmpMatrix = findOperatorMatrices(i,sigma,natBasis);
				int asign= 1;
				if (sigma>0) asign= 1;
				typename OperatorType::Su2RelatedType su2related;
				if (sigma==0) {
					su2related.source.push_back(i*offset_);
					su2related.source.push_back(i*offset_+1);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = NUMBER_OF_ORBITALS;
				}

				OperatorType myOp(tmpMatrix,
				                  ProgramGlobals::FermionOrBosonEnum::FERMION,
				                  typename OperatorType::PairType(1,1-sigma),
				                  asign,
				                  su2related);

				this->createOpsLabel("c").push(myOp);
			}

			tmpMatrix = findOperatorMatrices(i,natBasis);
			RealType angularFactor= 1;
			typename OperatorType::Su2RelatedType su2related;
			su2related.offset = 1; //check FIXME
			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  typename OperatorType::PairType(0,0),
			                  angularFactor,
			                  su2related);

			this->createOpsLabel("n").push(myOp);
		}

		this->makeTrackable("c");
		this->makeTrackable("n");
	}

	void fillModelLinks()
	{
		ModelTermType& hop = ModelBaseType::createTerm("hopping");
		ModelTermType& ninj = ModelBaseType::createTerm("ninj");

		OpForLinkType cup("c");
		hop.push(cup, 'N', cup, 'C', typename ModelTermType::Su2Properties(1, 1, 0));

		OpForLinkType n("n");
		ninj.push(n, 'N', n, 'N');

		if (hasDelta_) {
			ModelTermType& cicj = ModelBaseType::createTerm("delta");
			cicj.push(cup, 'N', cup, 'N');
		}
	}

	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		int sitesTimesDof=DEGREES_OF_FREEDOM*block.size();
		HilbertState total = (1<<sitesTimesDof);

		basis.resize(total);
		for (HilbertState a = 0; a < total; ++a)
			basis[a] = a;
	}

	// Calculate fermionic sign when applying operator
	// c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceType::HilbertState const &ket,
	              int i,
	              int sigma) const
	{
		int value=0;
		value += HilbertSpaceType::calcNofElectrons(ket,0,i,0);
		int tmp1 = HilbertSpaceType::get(ket,0) &1;
		if (i>0 && tmp1>0) value++;

		return (value%2==0) ? 1.0 : FERMION_SIGN;
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(int i,
	                                      int sigma,
	                                      const HilbertBasisType& natBasis) const
	{
		typename HilbertSpaceType::HilbertState bra,ket;
		int n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			if (HilbertSpaceType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceType::create(bra,i,sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				assert(jj >= 0);
				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis) const
	{
		const SizeType localSymms = ModelBaseType::targetQuantum().sizeOfOther();
		if (localSymms == 0) {
			if (hasDelta_) {
				PsimagLite::String msg(__FILE__);
				msg += ": You should be using one local symmetry, not zero\n";
				std::cerr<<msg;
				std::cerr<<msg;
			}
		} else if (localSymms == 1) {
			if (hasDelta_) {
				PsimagLite::String msg(__FILE__);
				err(msg + ": You should be using zero local symmetry, not one\n");
			}
		} else {
			PsimagLite::String msg(__FILE__);
			err(msg + ": Two many local symmetries in input file\n");
		}

		const bool isCanonical = (localSymms == 1);

		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		qns.resize(basis.size(), QnType::zero());
		VectorSizeType other;
		if (isCanonical) other.resize(1);
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = calcJmValue<PairType>(basis[i]);
			SizeType electrons = HilbertSpaceType::getNofDigits(basis[i],0);
			SizeType flavor = electrons;

			bool sign = electrons & 1;
			if (other.size() == 1)
				other[0] = electrons;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(int i,
	                                      const VectorHilbertStateType& natBasis) const
	{

		SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertState ket=natBasis[ii];
			cm(ii,ii) = 0.0;
			for (int sigma=0;sigma<DEGREES_OF_FREEDOM;sigma++)
				if (HilbertSpaceType::isNonZero(ket,i,sigma))
					cm(ii,ii) += 1.0;
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// does not work for 6 or 9
	template<typename PairType>
	PairType calcJmValue(const HilbertState& ket) const
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

	ParametersFermionSpinless<RealType, QnType>  modelParameters_;
	SizeType offset_;
	SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
	SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
	const bool hasDelta_;
};	//class FermionSpinless

} // namespace Dmrg
/*@}*/
#endif

