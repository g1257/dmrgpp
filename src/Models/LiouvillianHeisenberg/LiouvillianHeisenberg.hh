/*
Copyright (c) 2009, 2017-2026, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 6+]
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

/*! \file LiouvillianHeisenberg.hh
 *
 *  An implementation of the LiouvilleHeisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_LIOUVILLIANHEISENBERG_H
#define DMRG_LIOUVILLIANHEISENBERG_H

#include "CrsMatrix.h"
#include "ImaginaryUnitOrFail.hh"
#include "ParamsLiouvillianHeisenberg.hh"
#include "ProgramGlobals.h"
#include "Utils.h"
#include "VerySparseMatrix.h"
#include <algorithm>
#include <utility>

namespace Dmrg {

template <typename ModelBaseType> class LiouvillianHeisenberg : public ModelBaseType {
public:

	using ModelHelperType      = typename ModelBaseType::ModelHelperType;
	using BasisType            = typename ModelHelperType::BasisType;
	using SuperGeometryType    = typename ModelBaseType::SuperGeometryType;
	using LeftRightSuperType   = typename ModelBaseType::LeftRightSuperType;
	using LinkType             = typename ModelBaseType::LinkType;
	using OperatorsType        = typename ModelHelperType::OperatorsType;
	using RealType             = typename ModelHelperType::RealType;
	using VectorType           = typename ModelBaseType::VectorType;
	using QnType               = typename ModelBaseType::QnType;
	using VectorQnType         = typename ModelBaseType::VectorQnType;
	using BlockType            = typename ModelBaseType::BlockType;
	using SolverParamsType     = typename ModelBaseType::SolverParamsType;
	using SparseMatrixType     = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType    = typename SparseMatrixType::value_type;
	using WordType             = unsigned int long;
	using InputValidatorType   = typename ModelBaseType::InputValidatorType;
	using MatrixType           = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorSizeType       = typename PsimagLite::Vector<SizeType>::Type;
	using VectorRealType       = typename ModelBaseType::VectorRealType;
	using ModelTermType        = typename ModelBaseType::ModelTermType;
	using HilbertBasisType     = typename PsimagLite::Vector<SizeType>::Type;
	using OperatorType         = typename OperatorsType::OperatorType;
	using PairType             = typename OperatorType::PairType;
	using VectorOperatorType   = typename PsimagLite::Vector<OperatorType>::Type;
	using MyBasis              = typename ModelBaseType::MyBasis;
	using MyBasisWithOperators = typename ModelBaseType::BasisWithOperatorsType;
	using OpsLabelType         = typename ModelBaseType::OpsLabelType;
	using OpForLinkType        = typename ModelBaseType::OpForLinkType;
	using ModelParametersType  = ParamsLiouvillianHeisenberg<RealType, QnType>;

	LiouvillianHeisenberg(const SolverParamsType&  solverParams,
	                      InputValidatorType&      io,
	                      const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , superGeometry_(geometry)
	    , is_separate_(modelParameters_.ancillas == ModelParametersType::Ancillas::SEPARATE)

	{
		SizeType n = superGeometry_.numberOfSites();

		ModelParametersType::checkMagneticField(
		    modelParameters_.magneticFieldZ.size(), 'Z', n);
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType&  block,
	                                RealType          time) const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		assert(block.size() == 1);
		SizeType site = block[0];
		addMagneticFieldZ(hmatrix, modelParameters_.magneticFieldZ, site);

		// physical ancilla jump connection
		const auto iter = jump_mapping_.find(site);
		if (iter == jump_mapping_.end()) {
			// ATTENTION: Early exit here
			return;
		}

		SizeType index = iter->second;
		assert(index < jump_operators_.size());

		hmatrix += jump_operators_[index];
	}

	bool isHermitian() const final { return false; }

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType hilbert_small = modelParameters_.twiceTheSpin + 1; // = 2
		SizeType hilbert_big
		    = hilbert_small * hilbert_small; // physical and ancilla together
		SizeType total = (is_separate_) ? hilbert_small : hilbert_big;

		std::vector<SizeType> perm(total);
		{
			// Do not use natBasis outside this block
			HilbertBasisType natBasis(total);
			for (SizeType i = 0; i < total; ++i) {
				natBasis[i] = i;
				perm[i]     = i;
			}

			setSymmetryRelated(qns, natBasis);
		}

		std::vector<RealType> signs(hilbert_small, 1);

		// Small size operators
		SparseMatrixType splus_small = findSplusMatrices(); // 2-state space only
		SparseMatrixType sz_small    = findSzMatrices(); // 2-state space only

		// Set the operators S^+ in the natural basis
		SparseMatrixType tmpMatrix;
		std::string      label;

		if (is_separate_) {
			label     = "splus";
			tmpMatrix = splus_small;
		} else {
			// Physical S+
			label = "splus_p";
			externalProduct(tmpMatrix, splus_small, hilbert_small, signs, false, perm);
		}

		tmpMatrix.checkValidity();
		OperatorType myOp1(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
		this->createOpsLabel(label).push(myOp1);
		this->makeTrackable(label);

		// Physical S-
		myOp1.dagger();
		label = (is_separate_) ? "sminus" : "sminus_p";
		this->createOpsLabel(label).push(myOp1);

		if (!is_separate_) {
			// Ancilla S+ and S-
			externalProduct(tmpMatrix, splus_small, hilbert_small, signs, true, perm);
			tmpMatrix.checkValidity();
			OperatorType myOp2(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("splus_a").push(myOp2);
			this->makeTrackable("splus_a");

			myOp2.dagger();
			this->createOpsLabel("sminus_a").push(myOp2);
		}

		// Set the operators S^z in the natural basis
		if (is_separate_) {
			label     = "sz";
			tmpMatrix = sz_small;
		} else {
			// Physical Sz
			label = "sz_p";
			externalProduct(tmpMatrix, sz_small, hilbert_small, signs, false, perm);
		}

		tmpMatrix.checkValidity();
		OperatorType myOp3(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
		this->createOpsLabel(label).push(myOp3);
		this->makeTrackable(label);

		if (!is_separate_) {
			// Ancilla Sz
			externalProduct(tmpMatrix, sz_small, hilbert_small, signs, true, perm);
			tmpMatrix.checkValidity();
			OperatorType myOp4(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("sz_a").push(myOp4);
			this->makeTrackable("sz_a");
		}

		cacheJumpOperators();
	}

	void fillModelLinks()
	{
		// physical H
		// Multiply by -i
		ComplexOrRealType sqrt_minus_one = ImginaryUnitOrFail<ComplexOrRealType>::value();
		std::string       addition       = (is_separate_) ? "" : "_p";
		connectionsSpSm("spsm" + addition, "splus" + addition, -sqrt_minus_one);
		connectionSzSz("szsz" + addition, "sz" + addition, -sqrt_minus_one);

		if (!is_separate_) {
			// ancilla H^T
			// Multiply by +i
			connectionsSpSm("spsm_a", "splus_a", sqrt_minus_one);
			connectionSzSz("szsz_a", "sz_a", sqrt_minus_one);
		}

		// if (is_separate) the minus sign will be chosen from the input file
	}

private:

	void addMagneticFieldZ(SparseMatrixType&            hmatrix,
	                       const std::vector<RealType>& v,
	                       SizeType                     site) const
	{
		const SizeType linSize = superGeometry_.numberOfSites();
		if (v.size() != linSize) {
			// no magnetic field supplied
			if (v.size() == 0)
				return;

			// magnetic field supplied with wrong size
			// should have been caught earlier
			assert(false);
		}

		assert(site < v.size());
		RealType            tmp   = v[site];
		std::string         label = (is_separate_) ? "sz" : "sz_p";
		const OperatorType& sz_p  = ModelBaseType::naturalOperator(label, site, 0);
		hmatrix += tmp * sz_p.getCRS();

		if (!is_separate_) {
			const OperatorType& sz_a = ModelBaseType::naturalOperator("sz_a", site, 0);
			hmatrix += tmp * sz_a.getCRS();
		}
	}

	void connectionsSpSm(const std::string&       connection_name,
	                     const std::string&       op_name,
	                     const ComplexOrRealType& factor)
	{
		ModelTermType& spsm = ModelBaseType::createTerm(connection_name);
		OpForLinkType  splus(op_name);

		auto valueModiferTerm0 = [&factor = std::as_const(factor)](ComplexOrRealType& value)
		{ value *= (0.5 * factor); };
		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0);
	}

	void connectionSzSz(const std::string&       connection_name,
	                    const std::string&       op_name,
	                    const ComplexOrRealType& factor)
	{
		ModelTermType& szsz = ModelBaseType::createTerm(connection_name);

		OpForLinkType sz(op_name);
		auto valueModiferTerm0 = [&factor = std::as_const(factor)](ComplexOrRealType& value)
		{ value *= factor; };
		szsz.push(sz, 'N', sz, 'N', valueModiferTerm0);
	}

	//! Find S^+_site in the physical "small" basis only
	SparseMatrixType findSplusMatrices() const
	{
		// We assume there's only one kind of sites
		// FIXME TODO for SDHS
		// SDHS = site-dependent Hilbert spaces
		SizeType site = 0;

		SizeType total = modelParameters_.twiceTheSpin + 1;
		assert(total == 2);
		MatrixType cm(total, total);
		RealType   j              = 0.5 * modelParameters_.twiceTheSpin;
		SizeType   bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType   bits = 1 + ProgramGlobals::logBase2(modelParameters_.twiceTheSpin);
		SizeType   mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = ii;

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

	//! Find S^z_i in the physical small basis only
	SparseMatrixType findSzMatrices() const
	{
		// We assume there's only one kind of sites
		// FIXME TODO for SDHS
		// SDHS = site-dependent Hilbert spaces
		SizeType site = 0;

		SizeType total = modelParameters_.twiceTheSpin + 1;
		assert(total == 2);
		MatrixType cm(total, total);
		RealType   j              = 0.5 * modelParameters_.twiceTheSpin;
		SizeType   bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType   bits = ProgramGlobals::logBase2(modelParameters_.twiceTheSpin) + 1;
		SizeType   mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = ii;

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);
			assert(ketsite == ket);
			RealType m   = ketsite - j;
			cm(ket, ket) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	// There are no symmetries due to the s^- x s^- term and similar for s^+ s^+
	void setSymmetryRelated(VectorQnType& qns, const HilbertBasisType& basis) const
	{
		using PairType = std::pair<SizeType, SizeType>;

		VectorSizeType other(1);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			SizeType flavor = 1;
			qns[i]          = QnType(false, other, PairType(0, 0), flavor);
		}
	}

	void cacheJumpOperators()
	{
		SizeType n = superGeometry_.numberOfSites();
		if (modelParameters_.bath_gamma.size() != n) {
			err("bath_gamma.size need to be equal to number of sites\n");
		}

		if (modelParameters_.bath_f.size() != n) {
			err("bath_f.size need to be equal to number of sites\n");
		}

		for (SizeType site = 0; site < n; ++site) {
			if (modelParameters_.bath_gamma[site] == 0)
				continue;

			ComplexOrRealType factor_a
			    = modelParameters_.bath_gamma[site] * modelParameters_.bath_f[site];
			ComplexOrRealType factor_b = modelParameters_.bath_gamma[site] - factor_a;

			MatrixType dense = (is_separate_)
			    ? jumpOperatorWhenSeparate(site, factor_a, factor_b)
			    : jumpOperatorWhenTogether(site, factor_a, factor_b);

			SparseMatrixType tmpMatrix(dense);
			jump_mapping_[site] = jump_operators_.size();
			jump_operators_.push_back(tmpMatrix);
		}
	}

	// When separate, in the input file, if a physical site has f_i its ancilla ought to
	// have 1 - f_i
	MatrixType jumpOperatorWhenSeparate(SizeType                 site,
	                                    const ComplexOrRealType& factor_a,
	                                    const ComplexOrRealType& factor_b) const
	{
		SparseMatrixType sm = ModelBaseType::naturalOperator("sminus", site, 0).getCRS();
		SparseMatrixType sp = ModelBaseType::naturalOperator("splus", site, 0).getCRS();

		SparseMatrixType p0 = sm * sp;
		SparseMatrixType p1 = sp * sm;

		MatrixType p0_dense = p0.toDense();
		MatrixType p1_dense = p1.toDense();

		SizeType rank = p0_dense.rows();
		assert(rank == p0_dense.cols());
		assert(rank == 2);
		MatrixType tmp_dense(rank, rank);
		for (SizeType i = 0; i < rank; ++i) {
			for (SizeType j = 0; j < rank; ++j) {
				tmp_dense(i, j) = -0.5 * factor_a * p1_dense(i, j)
				    - 0.5 * factor_b * p0_dense(i, j);
			}
		}

		return tmp_dense;
	}

	MatrixType jumpOperatorWhenTogether(SizeType                 site,
	                                    const ComplexOrRealType& factor_a,
	                                    const ComplexOrRealType& factor_b) const
	{
		SparseMatrixType sm_p
		    = ModelBaseType::naturalOperator("sminus_p", site, 0).getCRS();
		SparseMatrixType sm_a
		    = ModelBaseType::naturalOperator("sminus_a", site, 0).getCRS();
		SparseMatrixType sp_p = ModelBaseType::naturalOperator("splus_p", site, 0).getCRS();
		SparseMatrixType sp_a = ModelBaseType::naturalOperator("splus_a", site, 0).getCRS();

		SparseMatrixType p0_p = sm_p * sp_p;
		SparseMatrixType p0_a = sm_a * sp_a;
		SparseMatrixType p1_p = sp_p * sm_p;
		SparseMatrixType p1_a = sp_a * sm_a;

		SizeType hilbert = sm_p.rows();
		assert(hilbert == 4); // because of the ancillas
		SizeType physical_hilbert = sqrt(hilbert);
		assert(physical_hilbert == 2); // spin 1/2
		std::vector<RealType> signs(physical_hilbert, 1);
		std::vector<SizeType> perm(hilbert);
		for (SizeType i = 0; i < hilbert; ++i) {
			perm[i] = i;
		}

		SparseMatrixType sp = findSplusMatrices();
		SparseMatrixType sm;
		transposeConjugate(sm, sp);

		SparseMatrixType sm_cross_sm;
		SparseMatrixType sp_cross_sp;
		externalProduct(sm_cross_sm, sm, sm, signs, true, perm);
		externalProduct(sp_cross_sp, sp, sp, signs, true, perm);

		MatrixType       sm_cross_sm_dense = sm_cross_sm.toDense();
		MatrixType       sp_cross_sp_dense = sp_cross_sp.toDense();
		SparseMatrixType p0                = p0_p;
		p0 += p0_a;
		SparseMatrixType p1 = p1_p;
		p1 += p1_a;
		MatrixType p0_dense = p0.toDense();
		MatrixType p1_dense = p1.toDense();

		SizeType rank = p0_dense.rows();
		assert(rank == p0_dense.cols());
		assert(rank == 4);
		MatrixType tmp_dense(rank, rank);
		for (SizeType i = 0; i < rank; ++i) {
			for (SizeType j = 0; j < rank; ++j) {
				tmp_dense(i, j)
				    = factor_a * (sm_cross_sm_dense(i, j) - 0.5 * p1_dense(i, j))
				    + factor_b * (sp_cross_sp_dense(i, j) - 0.5 * p0_dense(i, j));
			}
		}

		return tmp_dense;
	}

	ModelParametersType           modelParameters_;
	const SuperGeometryType&      superGeometry_;
	bool                          is_separate_;
	std::vector<SparseMatrixType> jump_operators_;
	std::map<SizeType, SizeType>  jump_mapping_;
}; // class LiouvillianHeisenberg

} // namespace Dmrg
/*@}*/
#endif // DMRG_LIOUVILLIANHEISENBERG_H
