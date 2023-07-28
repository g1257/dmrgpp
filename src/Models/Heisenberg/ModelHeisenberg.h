/*
Copyright (c) 2009, 2017-2019, UT-Battelle, LLC
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

/*! \file ModelHeisenberg.h
 *
 *  An implementation of the Quantum Heisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_MODEL_HEISENBERG_HEADER_H
#define DMRG_MODEL_HEISENBERG_HEADER_H

#include "../../Engine/ProgramGlobals.h"
#include "../../Engine/SpinSquared.h"
#include "../../Engine/SpinSquaredHelper.h"
#include "../../Engine/Utils.h"
#include "../../Engine/VerySparseMatrix.h"
#include "Aklt.h"
#include "CrsMatrix.h"
#include "ParametersModelHeisenberg.h"
#include <algorithm>

namespace Dmrg
{

template <typename ModelBaseType>
class ModelHeisenberg : public ModelBaseType
{

	static const int NUMBER_OF_ORBITALS = 1;
	static const int DEGREES_OF_FREEDOM = 2; // spin up and down

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename ModelBaseType::VectorQnType VectorQnType;
	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef unsigned int long WordType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef typename PsimagLite::Vector<SizeType>::Type HilbertBasisType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ModelBaseType::MyBasis MyBasis;
	typedef typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef Aklt<ModelBaseType> AkltType;
	using ModelParametersType = ParametersModelHeisenberg<RealType, QnType>;

	ModelHeisenberg(const SolverParamsType& solverParams,
	    InputValidatorType& io,
	    const SuperGeometryType& geometry,
	    PsimagLite::String additional)
	    : ModelBaseType(solverParams,
		geometry,
		io)
	    , modelParameters_(io)
	    , superGeometry_(geometry)
	    , additional_(additional)
	    , spinSquared_(spinSquaredHelper_, NUMBER_OF_ORBITALS, DEGREES_OF_FREEDOM)
	    , aklt_(*this, additional)
	{
		SizeType n = superGeometry_.numberOfSites();
		SizeType md = modelParameters_.anisotropyD.size();
		SizeType me = modelParameters_.anisotropyE.size();

		ModelParametersType::checkMagneticField(modelParameters_.magneticFieldX.size(), 'X', n);

		ModelParametersType::checkMagneticField(modelParameters_.magneticFieldZ.size(), 'Z', n);

		if (md > 0 && md != n) {
			PsimagLite::String msg("ModelHeisenberg: If provided, ");
			msg += " AnisotropyD must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (me > 0 && me != n) {
			PsimagLite::String msg("ModelHeisenberg: If provided, ");
			msg += " AnisotropyE must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		spinSquaredHelper_.write(label, io);
		spinSquared_.write(label, io);
	}

	/* PSIDOC Heisenberg::addDiagonalsInNaturalBasis
	 We describe only the addition of a Zeeman term to the Heisenberg model here; note
	 that this function is more complicated.
	 Please look carefully at the following C++ lines:
	 PSIDOCCOPY $FirstFunctionBelow::MagneticField
	 */
	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	    const BlockType& block,
	    RealType time) const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		SizeType linSize = superGeometry_.numberOfSites();
		SizeType n = block.size();

		for (SizeType i = 0; i < n; ++i) {

			// PSIDOCMARK_BEGIN MagneticField
			SizeType site = block[i];

			const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);
			const OperatorType& splus = ModelBaseType::naturalOperator("splus", site, 0);

			addMagneticField(hmatrix, 'X', site);
			addMagneticField(hmatrix, 'Z', site);

			// PSIDOCMARK_END

			// anisotropyD
			if (modelParameters_.anisotropyD.size() == linSize) {
				SparseMatrixType Szsquare;
				RealType tmp = modelParameters_.anisotropyD[site];
				multiply(Szsquare, sz.getCRS(), sz.getCRS());
				hmatrix += tmp * Szsquare;
			}

			// anisotropyE
			if (modelParameters_.anisotropyE.size() == linSize) {

				SparseMatrixType splusSquared;

				RealType tmp = 0.5 * modelParameters_.anisotropyE[site];
				multiply(splusSquared, splus.getCRS(), splus.getCRS());
				hmatrix += tmp * splusSquared;

				SparseMatrixType sminusSquared;
				transposeConjugate(sminusSquared, splusSquared);
				hmatrix += tmp * sminusSquared;
			}
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		SizeType total = utils::powUint(modelParameters_.twiceTheSpin + 1, block.size());
		HilbertBasisType natBasis(total);
		for (SizeType i = 0; i < total; ++i)
			natBasis[i] = i;

		setSymmetryRelated(qns, natBasis, block.size());

		for (SizeType i = 0; i < block.size(); i++) {
			// Set the operators S^+_i in the natural basis
			SparseMatrixType tmpMatrix = findSplusMatrices(i, natBasis);

			typename OperatorType::Su2RelatedType su2related;
			su2related.source.push_back(i * DEGREES_OF_FREEDOM);
			su2related.source.push_back(i * DEGREES_OF_FREEDOM + NUMBER_OF_ORBITALS);
			su2related.source.push_back(i * DEGREES_OF_FREEDOM);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(1);
			su2related.offset = NUMBER_OF_ORBITALS;

			OperatorType myOp(tmpMatrix,
			    ProgramGlobals::FermionOrBosonEnum::BOSON,
			    PairType(2, 2),
			    -1,
			    su2related);
			this->createOpsLabel("splus").push(myOp);
			if (additional_ != "2")
				this->makeTrackable("splus");

			myOp.dagger();
			this->createOpsLabel("sminus").push(myOp);

			// Set the operators S^z_i in the natural basis
			tmpMatrix = findSzMatrices(i, natBasis);
			typename OperatorType::Su2RelatedType su2related2;
			OperatorType myOp2(tmpMatrix,
			    ProgramGlobals::FermionOrBosonEnum::BOSON,
			    PairType(2, 1),
			    1.0 / sqrt(2.0),
			    su2related2);
			this->createOpsLabel("sz").push(myOp2);
			this->makeTrackable("sz");

			// Set the operators S^x_i in the natural basis
			tmpMatrix = findSxOrSyBarMatrices(i, natBasis, "sx");
			typename OperatorType::Su2RelatedType su2related3;
			OperatorType myOp3(tmpMatrix,
			    ProgramGlobals::FermionOrBosonEnum::BOSON,
			    PairType(2, 1),
			    1.0 / sqrt(2.0),
			    su2related3);
			this->createOpsLabel("sx").push(myOp3);

			if (additional_ == "Anisotropic" || additional_ == "2")
				this->makeTrackable("sx");

			// Set the operators S^ybar_i in the natural basis
			tmpMatrix = findSxOrSyBarMatrices(i, natBasis, "sybar");
			OperatorType syBarOp(tmpMatrix,
			    ProgramGlobals::FermionOrBosonEnum::BOSON,
			    PairType(0, 0),
			    1,
			    su2related3);
			this->createOpsLabel("sybar").push(syBarOp);

			if (additional_ == "2")
				this->makeTrackable("sybar");

			aklt_.fillLabeledOperators(i, myOp.getCRS(), myOp2.getCRS());

			tmpMatrix = findMaximal(i, natBasis);
			OperatorType myOp4(tmpMatrix,
			    ProgramGlobals::FermionOrBosonEnum::BOSON,
			    PairType(2, 1),
			    1.0 / sqrt(2.0),
			    su2related3);
			this->createOpsLabel("maximal").push(myOp4);
		}
	}

	void fillModelLinks()
	{
		if (additional_ == "2") {
			fillModelLinks2();
			return;
		}

		if (BasisType::useSu2Symmetry())
			err("SU(2) no longer supported\n");

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");

		OpForLinkType splus("splus");

		auto valueModiferTerm0 = [](ComplexOrRealType& value) { value *= 0.5; };

		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0);

		connectionSzSz();

		if (additional_ == "Anisotropic") {

			ModelTermType& sxsx = ModelBaseType::createTerm("sxsx");

			OpForLinkType sx("sx");

			sxsx.push(sx, 'N', sx, 'N', typename ModelTermType::Su2Properties(2, 1, 0));
		}

		aklt_.fillModelLinks();
	}

private:

	void addMagneticField(SparseMatrixType& hmatrix,
	    char c,
	    SizeType site) const
	{
		assert(c == 'X' || c == 'Z');

		const SizeType linSize = superGeometry_.numberOfSites();
		const VectorRealType& v = (c == 'X') ? modelParameters_.magneticFieldX
						     : modelParameters_.magneticFieldZ;

		if (v.size() != linSize)
			return;

		assert(site < v.size());
		RealType tmp = v[site];
		const OperatorType& sminus = ModelBaseType::naturalOperator("sminus", site, 0);

		const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);

		if (c == 'Z') {
			hmatrix += tmp * sz.getCRS();
			return;
		}

		assert(c == 'X');
		const OperatorType& splus = ModelBaseType::naturalOperator("splus", site, 0);
		constexpr RealType zeroPointFive = 0.5;
		hmatrix += zeroPointFive * tmp * splus.getCRS();
		hmatrix += zeroPointFive * tmp * sminus.getCRS();
	}

	void fillModelLinks2()
	{
		ModelTermType& sxsx = ModelBaseType::createTerm("sxsx");
		OpForLinkType sx("sx");
		sxsx.push(sx, 'N', sx, 'N');

		auto sybarsybarModifier = [](ComplexOrRealType& value) { value *= -1.0; };
		ModelTermType& sybarsybar = ModelBaseType::createTerm("sybarsybar");
		OpForLinkType sybar("sybar");
		sybarsybar.push(sybar, 'N', sybar, 'N', sybarsybarModifier);

		connectionSzSz();
	}

	void connectionSzSz()
	{
		ModelTermType& szsz = ModelBaseType::createTerm("szsz");

		OpForLinkType sz("sz");
		szsz.push(sz, 'N', sz, 'N', typename ModelTermType::Su2Properties(2, 0.5));
	}

	//! Find S^+_site in the natural basis natBasis
	SparseMatrixType findSplusMatrices(SizeType site,
	    const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total, total);
		RealType j = 0.5 * modelParameters_.twiceTheSpin;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType bits = 1 + ProgramGlobals::logBase2(modelParameters_.twiceTheSpin);
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);

			assert(ketsite == ket);
			SizeType brasite = ketsite + 1;
			if (brasite >= modelParameters_.twiceTheSpin + 1)
				continue;

			SizeType bra = ket & (~mask);
			assert(bra == 0);
			brasite <<= (site * bitsForOneSite);
			bra |= brasite;
			assert(bra == brasite);

			RealType m = ketsite - j;
			RealType x = j * (j + 1) - m * (m + 1);
			assert(x >= 0);

			cm(ket, bra) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(SizeType site,
	    const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total, total);
		RealType j = 0.5 * modelParameters_.twiceTheSpin;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType bits = ProgramGlobals::logBase2(modelParameters_.twiceTheSpin) + 1;
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);
			assert(ketsite == ket);
			RealType m = ketsite - j;
			cm(ket, ket) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find Maximal_i in the natural basis natBasis
	SparseMatrixType findMaximal(SizeType site,
	    const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total, total);
		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType bits = ProgramGlobals::logBase2(modelParameters_.twiceTheSpin) + 1;
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);
			assert(ketsite == ket);
			if (ketsite != modelParameters_.twiceTheSpin)
				continue;
			cm(ket, ket) = 1;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findSxOrSyBarMatrices(SizeType site,
	    const HilbertBasisType& natBasis,
	    PsimagLite::String what) const
	{
		if (what != "sx" && what != "sybar")
			err("findSxOrSyBarMatrices: don't know how to calculate " + what + "\n");

		const RealType sign = (what == "sx") ? 1.0 : -1.0;
		SparseMatrixType Splus_temp = findSplusMatrices(site, natBasis);
		SparseMatrixType Sminus_temp, Sx;
		transposeConjugate(Sminus_temp, Splus_temp);
		const RealType tmp = 0.5;

		Sx = tmp * Splus_temp;
		Sx += sign * tmp * Sminus_temp;

		return Sx;
	}

	void setSymmetryRelated(VectorQnType& qns,
	    const HilbertBasisType& basis,
	    int n) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType, SizeType> PairType;

		bool isCanonical = (ModelBaseType::targetQuantum().sizeOfOther() == 1);
		if (isCanonical && additional_ == "Anisotropic")
			err(PsimagLite::String(__FILE__) + ": Anisotropic sub-model CANNOT be canonical. Please " + "delete the TargetSzPlusConst= from the input file\n");

		if (isCanonical && !modelParameters_.magneticFieldX.empty())
			err(PsimagLite::String(__FILE__) + ": MagneticFieldX CANNOT be canonical. Please " + "delete the TargetSzPlusConst= from the input file\n");

		VectorSizeType other;
		if (isCanonical)
			other.resize(1, 0);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(modelParameters_.twiceTheSpin, basis[i]);
			if (isCanonical)
				other[0] = getSzPlusConst(basis[i], n);
			SizeType flavor = 1;
			qns[i] = QnType(false, other, jmpair, flavor);
		}
	}

	SizeType getSzPlusConst(SizeType ket, SizeType n) const
	{
		if (n == 1)
			return ket;

		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);

		SizeType sum = 0;
		for (SizeType site = 0; site < n; ++site) {
			SizeType mask = modelParameters_.twiceTheSpin;
			mask <<= (site * bitsForOneSite);
			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);
			sum += ketsite;
		}

		return sum;
	}

	ModelParametersType modelParameters_;
	const SuperGeometryType& superGeometry_;
	SpinSquaredHelper<RealType, WordType> spinSquaredHelper_;
	PsimagLite::String additional_;
	SpinSquared<SpinSquaredHelper<RealType, WordType>> spinSquared_;
	AkltType aklt_;
}; // class ModelHeisenberg

} // namespace Dmrg
/*@}*/
#endif // DMRG_MODEL_HEISENBERG_HEADER_H
