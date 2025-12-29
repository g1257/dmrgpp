/*
Copyright (c) 2009, 2017-2025, UT-Battelle, LLC
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

/*! \file LiouvillianHeisenberg.hh
 *
 *  An implementation of the LiouvilleHeisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_LIOUVILLIANHEISENBERG_H
#define DMRG_LIOUVILLIANHEISENBERG_H

#include "../../Engine/ProgramGlobals.h"
#include "../../Engine/Utils.h"
#include "../../Engine/VerySparseMatrix.h"
#include "CrsMatrix.h"
#include "ParamsLiouvillianHeisenberg.hh"
#include <algorithm>

namespace Dmrg
{

template <typename ModelBaseType>
class LiouvillianHeisenberg : public ModelBaseType
{
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
	using ModelParametersType = ParamsLiouvillianHeisenberg<RealType, QnType>;

	LiouvillianHeisenberg(const SolverParamsType& solverParams,
	    InputValidatorType& io,
	    const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams,
		  geometry,
		  io)
	    , modelParameters_(io)
	    , superGeometry_(geometry)
	{
		SizeType n = superGeometry_.numberOfSites();

		ModelParametersType::checkMagneticField(modelParameters_.magneticFieldZ.size(), 'Z', n);

		cacheJumpOperators();
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
	    const BlockType& block,
	    RealType time) const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		constexpr SizeType oneSitePhysical = 2;
		SizeType total = oneSitePhysical * oneSitePhysical;

		std::vector<SizeType> perm(total);
		for (SizeType i = 0; i < total; ++i) {
			perm[i] = i;
		}

		std::vector<RealType> signs(oneSitePhysical, 1);

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

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		SizeType oneSitePhysical = modelParameters_.twiceTheSpin + 1;
		// Contains ancilla
		SizeType total = utils::powUint(oneSitePhysical * oneSitePhysical, block.size());
		HilbertBasisType natBasis(total);
		std::vector<SizeType> perm(total);
		for (SizeType i = 0; i < total; ++i) {
			natBasis[i] = i;
			perm[i] = i;
		}

		setSymmetryRelated(qns, natBasis, block.size());

		std::vector<RealType> signs(oneSitePhysical, 1);

		for (SizeType i = 0; i < block.size(); i++) {
			// Set the operators S^+_i in the natural basis
			SparseMatrixType tmpMatrix1 = findSplusMatrices(i, natBasis);
			SparseMatrixType tmpMatrix;

			// Physical S+ and S-
			externalProduct(tmpMatrix, tmpMatrix1, oneSitePhysical, signs, false, perm);
			OperatorType myOp1(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("splus_p").push(myOp1);
			this->makeTrackable("splus_p");

			myOp1.dagger();
			this->createOpsLabel("sminus_p").push(myOp1);

			// Ancilla S+ and S-
			externalProduct(tmpMatrix, tmpMatrix1, oneSitePhysical, signs, true, perm);
			OperatorType myOp2(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("splus_a").push(myOp2);
			this->makeTrackable("splus_a");

			myOp2.dagger();
			this->createOpsLabel("sminus_a").push(myOp2);

			// Physical Sz
			tmpMatrix1 = findSzMatrices(i, natBasis);
			externalProduct(tmpMatrix, tmpMatrix1, oneSitePhysical, signs, false, perm);
			OperatorType myOp3(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("sz_p").push(myOp3);
			this->makeTrackable("sz_p");

			// Ancilla Sz
			externalProduct(tmpMatrix, tmpMatrix1, oneSitePhysical, signs, true, perm);
			OperatorType myOp4(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("sz_a").push(myOp4);
			this->makeTrackable("sz_a");
		}
	}

	void fillModelLinks()
	{
		if (BasisType::useSu2Symmetry())
			err("SU(2) no longer supported\n");

		// physical H
		// FIXME: Multiply by -i
		connectionsSpSm("spsm_p", "splus_p");
		connectionSzSz("szsz_p", "sz_p");

		// ancilla H^T
		// FIXME: Multiply by +i
		connectionsSpSm("spsm_a", "splus_a");
		connectionSzSz("szsz_a", "sz_a");
	}

private:

	void addMagneticFieldZ(SparseMatrixType& hmatrix, const std::vector<RealType>& v, SizeType site) const
	{
		const SizeType linSize = superGeometry_.numberOfSites();
		if (v.size() != linSize)
			return;

		assert(site < v.size());
		RealType tmp = v[site];
		const OperatorType& sz_p = ModelBaseType::naturalOperator("sz_p", site, 0);
		hmatrix += tmp * sz_p.getCRS();

		const OperatorType& sz_a = ModelBaseType::naturalOperator("sz_a", site, 0);
		hmatrix += tmp * sz_a.getCRS();
	}

	void connectionsSpSm(const std::string& connection_name, const std::string& op_name)
	{
		ModelTermType& spsm = ModelBaseType::createTerm(connection_name);
		OpForLinkType splus(op_name);
		auto valueModiferTerm0 = [](ComplexOrRealType& value) { value *= 0.5; };
		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0);
	}

	void connectionSzSz(const std::string& connection_name, const std::string& op_name)
	{
		ModelTermType& szsz = ModelBaseType::createTerm(connection_name);

		OpForLinkType sz(op_name);
		szsz.push(sz, 'N', sz, 'N');
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

	// There are no symmetries due to the s^- x s^- term and similar for s^+ s^+
	void setSymmetryRelated(VectorQnType& qns,
	    const HilbertBasisType& basis,
	    int) const
	{
		typedef std::pair<SizeType, SizeType> PairType;

		VectorSizeType other;
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			SizeType flavor = 1;
			qns[i] = QnType(false, other, PairType(0, 0), flavor);
		}
	}

	bool mustApplyJumpOps(SizeType site) const
	{
		return (site == 0 || site == superGeometry_.numberOfSites());
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
			cacheJumpOperator(site);
		}
	}

	void cacheJumpOperator(SizeType site)
	{
		assert(modelParameters_.bath_gamma.size() > site);
		assert(modelParameters_.bath_f.size() > site);
		ComplexOrRealType factor_a = modelParameters_.bath_gamma[site] * modelParameters_.bath_f[site];
		ComplexOrRealType factor_b = modelParameters_.bath_gamma[site] - factor_a;

		const SparseMatrixType& sm_p = ModelBaseType::naturalOperator("sm_p", site, 0).getCRS();
		const SparseMatrixType& sm_a = ModelBaseType::naturalOperator("sm_a", site, 0).getCRS();
		const SparseMatrixType& sp_p = ModelBaseType::naturalOperator("sp_p", site, 0).getCRS();
		const SparseMatrixType& sp_a = ModelBaseType::naturalOperator("sp_a", site, 0).getCRS();

		SparseMatrixType p0_p = sm_p * sp_p;
		SparseMatrixType p0_a = sm_a * sp_a;
		SparseMatrixType p1_p = sp_p * sm_p;
		SparseMatrixType p1_a = sp_a * sm_a;

		SparseMatrixType sm_a_cross_sm_b;
		SparseMatrixType sp_a_cross_sp_b;
		SizeType hilbert = sm_p.rows();
		assert(hilbert == 4); // because of the ancillas
		SizeType physical_hilbert = sqrt(hilbert);
		assert(physical_hilbert == 2); // spin 1/2
		std::vector<RealType> signs(physical_hilbert, 1);
		std::vector<SizeType> perm(hilbert);
		for (SizeType i = 0; i < hilbert; ++i) {
			perm[i] = i;
		}

		externalProduct(sm_a_cross_sm_b, sm_a, sm_p, signs, true, perm);
		externalProduct(sp_a_cross_sp_b, sp_a, sp_p, signs, true, perm);

		MatrixType sm_a_cross_sm_b_dense = sm_a_cross_sm_b.toDense();
		MatrixType sp_a_cross_sp_b_dense = sp_a_cross_sp_b.toDense();
		SparseMatrixType p0 = p0_p;
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
				tmp_dense(i, j) = factor_a * (sm_a_cross_sm_b_dense(i, j) - 0.5 * p1_dense(i, j)) + factor_b * (sp_a_cross_sp_b_dense(i, j) - 0.5 * p0_dense(i, j));
			}
		}

		SparseMatrixType tmpMatrix(tmp_dense);
		jump_mapping_[site] = jump_operators_.size();
		jump_operators_.push_back(tmpMatrix);
	}

	ModelParametersType modelParameters_;
	const SuperGeometryType& superGeometry_;
	std::vector<SparseMatrixType> jump_operators_;
	std::map<SizeType, SizeType> jump_mapping_;
}; // class LiouvillianHeisenberg

} // namespace Dmrg
/*@}*/
#endif // DMRG_LIOUVILLIANHEISENBERG_H
