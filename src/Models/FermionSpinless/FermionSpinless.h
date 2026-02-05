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
#include "CrsMatrix.h"
#include "HilbertSpaceFermionSpinless.h"
#include "ParametersFermionSpinless.h"
#include "ProgramGlobals.h"
#include "Sort.h" // in PsimagLite
#include "SpinSquared.h"
#include "SpinSquaredHelper.h"
#include "VerySparseMatrix.h"
#include <cassert>

namespace Dmrg {
//! Model Hubbard for DMRG solver, inherits from ModelBase and implements its interface:
template <typename ModelBaseType> class FermionSpinless : public ModelBaseType {

	static const int FERMION_SIGN       = -1;
	static const int DEGREES_OF_FREEDOM = 1;
	static const int NUMBER_OF_ORBITALS = 1;

public:

	using ModelHelperType        = typename ModelBaseType::ModelHelperType;
	using SuperGeometryType      = typename ModelBaseType::SuperGeometryType;
	using LeftRightSuperType     = typename ModelBaseType::LeftRightSuperType;
	using LinkType               = typename ModelBaseType::LinkType;
	using OperatorsType          = typename ModelHelperType::OperatorsType;
	using OperatorType           = typename OperatorsType::OperatorType;
	using RealType               = typename ModelHelperType::RealType;
	using QnType                 = typename ModelBaseType::QnType;
	using VectorQnType           = typename QnType::VectorQnType;
	using SparseMatrixType       = typename ModelBaseType::SparseMatrixType;
	using ComplexOrRealType      = typename SparseMatrixType::value_type;
	using WordType               = unsigned int long;
	using HilbertSpaceType       = HilbertSpaceFermionSpinless<WordType>;
	using VectorOperatorType     = typename ModelBaseType::VectorOperatorType;
	using VectorSizeType         = typename ModelBaseType::VectorSizeType;
	using OpsLabelType           = typename ModelBaseType::OpsLabelType;
	using BlockType              = typename ModelBaseType::BlockType;
	using SolverParamsType       = typename ModelBaseType::SolverParamsType;
	using VectorType             = typename ModelBaseType::VectorType;
	using HilbertState           = typename HilbertSpaceType::HilbertState;
	using VectorHilbertStateType = typename PsimagLite::Vector<HilbertState>::Type;
	using InputValidatorType     = typename ModelBaseType::InputValidatorType;
	using BasisType              = typename ModelBaseType::MyBasis;
	using MyBasisWithOperators   = typename ModelBaseType::BasisWithOperatorsType;
	using HilbertBasisType       = typename PsimagLite::Vector<HilbertState>::Type;
	using OpForLinkType          = typename ModelBaseType::OpForLinkType;
	using ModelTermType          = typename ModelBaseType::ModelTermType;
	using VectorPairRealBoolType = typename PsimagLite::Vector<std::pair<RealType, bool>>::Type;

	FermionSpinless(const SolverParamsType&  solverParams,
	                InputValidatorType&      io,
	                const SuperGeometryType& geometry,
	                PsimagLite::String       extra)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , offset_(DEGREES_OF_FREEDOM)
	    , spinSquared_(spinSquaredHelper_, NUMBER_OF_ORBITALS, DEGREES_OF_FREEDOM)
	    , hasDelta_(extra == "WithDelta")
	    , hasCalcMu_(false)
	    , tau_(0)
	    , mu_(0)
	    , previousTimeStep_(0)
	    , vectorMu_(geometry.numberOfSites())
	{
		if (extra != "" && extra != "WithDelta")
			err("FermionSpinLess can only be followed by WithDelta and not " + extra
			    + "\n");

		const SizeType n = geometry.numberOfSites();

		if (geometry.numberOfSites() != modelParameters_.potentialV.size())
			err("potentialV must have exactly " + ttos(n) + " entries.\n");

		bool hasTau = false;
		try {
			io.readline(tau_, "TSPTau=");
			hasTau = true;
		} catch (std::exception&) { }

		try {
			io.readline(mu_, "TSPMu=");
			hasCalcMu_ = true;
		} catch (std::exception&) { }

		bool b1 = (hasTau && !hasCalcMu_);
		bool b2 = (!hasTau && hasCalcMu_);
		if (b1 || b2)
			err("FermionSpinless: Both or none of TSPTau= and TSPMu= must appear\n");
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

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType&  block,
	                                RealType          time) const
	{
		static bool firstCall = true;
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		SizeType n = block.size();
		if (n != 1)
			err("addDiagonalsInNaturalBasis: block.size() != 1\n");

		const SizeType          site   = block[0];
		const OperatorType&     niupop = ModelBaseType::naturalOperator("n", site, 0);
		const SparseMatrixType& niup   = niupop.getCRS();

		// V_iup term
		RealType tmp = modelParameters_.potentialV[site];
		hmatrix += tmp * niup;

		if (hasCalcMu_) {
			const SizeType nsites = ModelBaseType::superGeometry().numberOfSites();
			const SizeType currentTimeStep = time / tau_;
			const RealType effectiveMu
			    = calcMu(site, currentTimeStep, nsites, tau_, mu_);
			hmatrix += effectiveMu * niup;
			if (previousTimeStep_ != currentTimeStep || firstCall) {
				printMuOfR(time);
				previousTimeStep_ = currentTimeStep;
				firstCall         = false;
				for (SizeType i = 0; i < vectorMu_.size(); ++i)
					vectorMu_[site].second = false;

				vectorMu_[site].first  = effectiveMu;
				vectorMu_[site].second = true;
			} else {
				vectorMu_[site].first  = effectiveMu;
				vectorMu_[site].second = true;
			}
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType         site = 0;
		BlockType        block(1, site);
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;
		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis);

		//! Set the operators c^\daggger_{i\sigma} in the natural basis
		for (SizeType i = 0; i < block.size(); i++) {
			for (int sigma = 0; sigma < DEGREES_OF_FREEDOM; sigma++) {
				tmpMatrix = findOperatorMatrices(i, sigma, natBasis);
				int asign = 1;
				if (sigma > 0)
					asign = 1;
				typename OperatorType::Su2RelatedType su2related;
				if (sigma == 0) {
					su2related.source.push_back(i * offset_);
					su2related.source.push_back(i * offset_ + 1);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = NUMBER_OF_ORBITALS;
				}

				OperatorType myOp(tmpMatrix,
				                  ProgramGlobals::FermionOrBosonEnum::FERMION,
				                  typename OperatorType::PairType(1, 1 - sigma),
				                  asign,
				                  su2related);

				this->createOpsLabel("c").push(myOp);
			}

			tmpMatrix = findOperatorMatrices(i, natBasis);
			RealType                              angularFactor = 1;
			typename OperatorType::Su2RelatedType su2related;
			su2related.offset = 1; // check FIXME
			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  typename OperatorType::PairType(0, 0),
			                  angularFactor,
			                  su2related);

			this->createOpsLabel("n").push(myOp);
		}

		this->makeTrackable("c");
		this->makeTrackable("n");
	}

	void fillModelLinks()
	{
		ModelTermType& hop  = ModelBaseType::createTerm("hopping");
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

	void setBasis(HilbertBasisType& basis, const VectorSizeType& block) const
	{
		int          sitesTimesDof = DEGREES_OF_FREEDOM * block.size();
		HilbertState total         = (1 << sitesTimesDof);

		basis.resize(total);
		for (HilbertState a = 0; a < total; ++a)
			basis[a] = a;
	}

	// Calculate fermionic sign when applying operator
	// c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceType::HilbertState const& ket, int i, int sigma) const
	{
		int value = 0;
		value += HilbertSpaceType::calcNofElectrons(ket, 0, i, 0);
		int tmp1 = HilbertSpaceType::get(ket, 0) & 1;
		if (i > 0 && tmp1 > 0)
			value++;

		return (value % 2 == 0) ? 1.0 : FERMION_SIGN;
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType
	findOperatorMatrices(int i, int sigma, const HilbertBasisType& natBasis) const
	{
		typename HilbertSpaceType::HilbertState                   bra, ket;
		int                                                       n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n, n);

		for (SizeType ii = 0; ii < natBasis.size(); ii++) {
			bra = ket = natBasis[ii];
			if (HilbertSpaceType::isNonZero(ket, i, sigma)) {

			} else {
				HilbertSpaceType::create(bra, i, sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis, bra);
				assert(jj >= 0);
				cm(ii, jj) = sign(ket, i, sigma);
			}
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	void setSymmetryRelated(VectorQnType& qns, const HilbertBasisType& basis) const
	{
		const SizeType localSymms = ModelBaseType::targetQuantum().sizeOfOther();
		if (localSymms == 0) {
			if (!hasDelta_) {
				PsimagLite::String msg(__FILE__);
				msg += ": You should be using one local symmetry, not zero\n";
				std::cerr << msg;
				std::cerr << msg;
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
		using PairType = std::pair<SizeType, SizeType>;
		qns.resize(basis.size(), QnType::zero());
		VectorSizeType other;
		if (isCanonical)
			other.resize(1);
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair    = calcJmValue<PairType>(basis[i]);
			SizeType electrons = HilbertSpaceType::getNofDigits(basis[i], 0);
			SizeType flavor    = electrons;

			bool sign = electrons & 1;
			if (other.size() == 1)
				other[0] = electrons;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(int i, const VectorHilbertStateType& natBasis) const
	{

		SizeType                                                  n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n, n);

		for (SizeType ii = 0; ii < natBasis.size(); ii++) {
			HilbertState ket = natBasis[ii];
			cm(ii, ii)       = 0.0;
			for (int sigma = 0; sigma < DEGREES_OF_FREEDOM; sigma++)
				if (HilbertSpaceType::isNonZero(ket, i, sigma))
					cm(ii, ii) += 1.0;
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// does not work for 6 or 9
	template <typename PairType> PairType calcJmValue(const HilbertState& ket) const
	{
		SizeType site0 = 0;
		SizeType site1 = 0;

		spinSquared_.doOnePairOfSitesA(ket, site0, site1);
		spinSquared_.doOnePairOfSitesB(ket, site0, site1);
		spinSquared_.doDiagonal(ket, site0, site1);

		RealType sz = spinSquared_.spinZ(ket, site0);
		PairType jm = spinSquaredHelper_.getJmPair(sz);

		return jm;
	}

	static RealType calcMu(SizeType site, SizeType n, SizeType N1, RealType tau, RealType mu)
	{
		//  (nm is the number of steps to increase the onsite
		//   chemical potential at site i)
		static SizeType nm = 0;
		static SizeType rm = 1; //( steps changes when n=1000)

		if (site + 1 == N1 || site == 0) {
			++nm;

			if (n % 1000 == 0 && n > 0) {
				++rm;
				nm = 0;

				std::stringstream msg;
				msg << "FermionSpinless::calcMu(): nm=0 "
				    << "site=" << site << " n=" << n << " tau=" << tau
				    << " mu=" << mu << " rm=" << rm << "\n";
				std::cout << msg.str();
			}
		}

		if (site + rm < N1) { //(chemical potential for  i <= end site - rm)
			return -mu;
		}

		if (site + rm > N1) { //(chemical potential for  i > end site -rm+1)
			return -64.0;
		}

		assert(site + rm == N1); //(chemical potential for  i = end site -rm+1)
		return -64.0 * tau * nm;
	}

	void printMuOfR(RealType time) const
	{
		std::cout << "timeStep=" << previousTimeStep_ << " time=" << time << " ";
		std::cout << "potentialV=[";
		for (SizeType i = 0; i < vectorMu_.size(); ++i) {
			RealType value = (vectorMu_[i].second) ? vectorMu_[i].first : -100;
			std::cout << value;
			if (i + 1 < vectorMu_.size())
				std::cout << ", ";
		}

		std::cout << "];\n";
	}
	ParametersFermionSpinless<RealType, QnType>        modelParameters_;
	SizeType                                           offset_;
	SpinSquaredHelper<RealType, WordType>              spinSquaredHelper_;
	SpinSquared<SpinSquaredHelper<RealType, WordType>> spinSquared_;
	const bool                                         hasDelta_;
	bool                                               hasCalcMu_;
	RealType                                           tau_;
	RealType                                           mu_;
	mutable SizeType                                   previousTimeStep_;
	mutable VectorPairRealBoolType                     vectorMu_;

}; // class FermionSpinless

} // namespace Dmrg
/*@}*/
#endif
