#ifndef DMRGPP_MODEL_KONDO_H
#define DMRGPP_MODEL_KONDO_H
#include "ParametersKondo.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ModelBaseType>
class Kondo : public ModelBaseType {

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename ModelBaseType::VectorQnType VectorQnType;
	typedef typename ModelBaseType::BlockType VectorSizeType;
	typedef typename ModelBaseType::RealType RealType;
	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef ParametersKondo<RealType, QnType> ParametersKondoType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef std::pair<SizeType, SizeType> PairSizeType;

public:

	Kondo(const SolverParamsType& solverParams,
	      InputValidatorType& io,
	      const SuperGeometryType& geometry,
	      PsimagLite::String option)
	    : ModelBaseType(solverParams, geometry, io),
	      solverParams_(solverParams),
	      modelParams_(io, option == "Ex")
	{
		SizeType sitesTimesDof = 2 + modelParams_.twiceTheSpin;
		SizeType total = (1<<sitesTimesDof);
		basis_.resize(total);
		for (SizeType a = 0; a < total; ++a)
			basis_[a] = a;

		setSymmetryRelatedInternal(qn_, basis_);

		ModelBaseType::orderByQuantum(basis_, qn_);

		setOperatorMatricesInternal();
	}

	// For information purposes only. Write model parameters
	// String contains the group
	// Serializer object is second argument
	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParams_.write(label, io);
		io.write(label + "/basis_", basis_);
		io.write(label + "/ops_", ops_);
	}

	// Fill SparseMatrixType with the on-site Hamiltonian terms in the on-site basis
	// Give SparseMatrixType in the order you chose to give the
	// operators in setOperatorMatrices
	// The RealType contain the physical time in case your onsite terms
	// depend on it
	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const VectorSizeType& block,
	                                RealType)  const
	{
		SizeType n = block.size();
		assert(n == 1);
		SizeType linSize = ModelBaseType::superGeometry().numberOfSites();
		SparseMatrixType tmpMatrix;
		SparseMatrixType niup;
		SparseMatrixType nidown;

		for (SizeType i = 0; i < n; ++i) {
			SizeType ind = block[i];
			const OperatorType& cup = ModelBaseType::naturalOperator("c", ind, 0);
			const OperatorType& cdown = ModelBaseType::naturalOperator("c", ind, 1);
			assert(ops_.size() > 1 + i*4);

			// onsite U hubbard
			//n_i up
			transposeConjugate(tmpMatrix, cup.getCRS());
			multiply(niup, tmpMatrix, cup.getCRS());

			//n_i down
			transposeConjugate(tmpMatrix, cdown.getCRS());
			multiply(nidown, tmpMatrix, cdown.getCRS());

			multiply(tmpMatrix,niup,nidown);
			assert(ind < modelParams_.hubbardU.size());
			RealType tmp = modelParams_.hubbardU[ind];
			hmatrix += tmp*tmpMatrix;

			// V_iup term
			assert(ind + 1*linSize < modelParams_.potentialV.size());
			tmp = modelParams_.potentialV[ind + 0*linSize];
			hmatrix += tmp*niup;

			// V_idown term
			tmp = modelParams_.potentialV[ind + 1*linSize];
			hmatrix += tmp*nidown;

			// Kondo term
			assert(ind < modelParams_.kondoJ.size());
			hmatrix += modelParams_.kondoJ[ind] * kondoOnSite(ind, niup, nidown);

			if (!modelParams_.extended) continue; // EARLY CONTINUE HERE

			const OperatorType& Splus = ModelBaseType::naturalOperator("Splus", 0, 0);
			hmatrix += modelParams_.kondoHx*Splus.getCRS();
			transposeConjugate(tmpMatrix, Splus.getCRS());
			hmatrix += modelParams_.kondoHx*tmpMatrix;

			transposeConjugate(tmpMatrix, cup.getCRS());
			multiply(niup, tmpMatrix, cdown.getCRS());
			hmatrix += modelParams_.electronHx*niup;

			transposeConjugate(tmpMatrix, cdown.getCRS());
			multiply(nidown, tmpMatrix, cup.getCRS());
			hmatrix += modelParams_.electronHx*nidown;

			multiply(niup, cup.getCRS(), cdown.getCRS());
			hmatrix += modelParams_.pairingField*niup;

			transposeConjugate(nidown, niup);
			hmatrix += modelParams_.pairingField*nidown;
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		qns = qn_;
		assert(ops_.size() >= 5);
		OpsLabelType& c = this->createOpsLabel("c");
		for (SizeType sigma = 0; sigma < 2; ++sigma)
			c.push(ops_[sigma]);

		this->createOpsLabel("Splus").push(ops_[2]);
		this->createOpsLabel("Sz").push(ops_[3]);
		this->createOpsLabel("n").push(ops_[4]);

		this->makeTrackable("c");
		this->makeTrackable("Splus");
		this->makeTrackable("Sz");
		this->makeTrackable("n");

		{
			SparseMatrixType mup = ops_[0].getCRS();
			SparseMatrixType mupT;
			transposeConjugate(mupT, mup);
			SparseMatrixType mdown = ops_[1].getCRS();
			SparseMatrixType mdownT;
			transposeConjugate(mdownT, mdown);
			SparseMatrixType szMatrix = mdownT*mdown;
			szMatrix *= (-1.0);
			szMatrix += mupT*mup;
			szMatrix *= 0.5;

			PairSizeType zeroPair(0, 0);
			typename OperatorType::Su2RelatedType su2Related;
			this->createOpsLabel("sz").push(OperatorType(szMatrix,
			                                             ProgramGlobals::FermionOrBosonEnum::BOSON,
			                                             zeroPair,
			                                             1,
			                                             su2Related));
		}
	}

	void fillModelLinks()
	{
		const bool isSu2 = BasisType::useSu2Symmetry();
		ModelTermType& hop = ModelBaseType::createTerm("hopping");
		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");
		ModelTermType& szsz = ModelBaseType::createTerm("SzSz");
		ModelTermType& ninj = ModelBaseType::createTerm("ninj");

		for (SizeType spin = 0; spin < 2; ++spin) {
			OpForLinkType c("c", spin);
			typename ModelTermType::Su2Properties su2properties(1, (spin == 1) ? -1 : 1, spin);
			hop.push(c, 'N', c, 'C', su2properties);
		}

		auto valueModiferTerm0 = [isSu2](ComplexOrRealType& value)
		{ value *= (isSu2) ? -0.5 : 0.5;};
		auto valueModifierTermOther = [isSu2](ComplexOrRealType& value)
		{ if (isSu2) value = -value;};

		OpForLinkType splus("Splus");
		OpForLinkType sz("Sz");
		OpForLinkType n("n");

		typename ModelTermType::Su2Properties su2properties(2, -1, 2);
		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0, su2properties);

		if (!isSu2) {
			typename ModelTermType::Su2Properties su2properties(2, 0.5);
			szsz.push(sz, 'N', sz, 'N', su2properties);
		} else {
			typename ModelTermType::Su2Properties su2properties(2, -1, 2);
			spsm.push(splus, 'N', splus, 'C', valueModifierTermOther, su2properties);
		}

		ninj.push(n, 'N', n, 'N');
	}

private:

	void setSymmetryRelatedInternal(VectorQnType& qns,
	                                const VectorSizeType& basis) const
	{
		qns.resize(basis.size(), QnType::zero());
		SizeType nsym = (modelParams_.extended) ? 0 : 2;
		VectorSizeType other(nsym);

		// bit 0 <--- up electron
		// bit 1 <--- down electron
		// bit 2 <--- localized spin down
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairSizeType jmpair(0,0);

			// nup
			SizeType electronsUp = (basis[i] & 1) ? 1 : 0;
			// ndown
			SizeType electronsDown = (basis[i] & 2) ? 1 : 0;

			SizeType electrons = electronsDown + electronsUp;

			SizeType bosonicSz = basis[i];
			bosonicSz >>= 2; // delete electronic part

			if (other.size() >= 2) {
				other[0] = electrons;
				other[1] = electronsUp + bosonicSz;
			}

			bool sign = electrons & 1;
			qns[i] = QnType(sign, other, jmpair, 0);
		}
	}

	void setOperatorMatricesInternal()
	{
		typename OperatorType::Su2RelatedType su2related;

		for (SizeType sigma = 0; sigma < 2; ++sigma) {
			SparseMatrixType tmpMatrix = findCmatrix(sigma, basis_);

			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::FERMION,
			                  typename OperatorType::PairType(0, 0),
			                  1,
			                  su2related);

			ops_.push_back(myOp);
		}

		// now the S+ and Sz for local spins
		SparseMatrixType m = findSplusMatrix(basis_);
		OperatorType sp(m,
		                ProgramGlobals::FermionOrBosonEnum::BOSON,
		                typename OperatorType::PairType(0, 0),
		                1,
		                su2related);
		ops_.push_back(sp);

		m = findSzMatrix(basis_);
		OperatorType sz(m,
		                ProgramGlobals::FermionOrBosonEnum::BOSON,
		                typename OperatorType::PairType(0, 0),
		                1,
		                su2related);
		ops_.push_back(sz);

		m = findNmatrix(basis_);
		OperatorType nm(m,
		                ProgramGlobals::FermionOrBosonEnum::BOSON,
		                typename OperatorType::PairType(0, 0),
		                1,
		                su2related);
		ops_.push_back(nm);

	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findCmatrix(SizeType sigma,
	                             const VectorSizeType& basis) const
	{
		SizeType n = basis.size();
		const ComplexOrRealType zero = 0.0;
		MatrixType cm(n, n, zero);
		SizeType mask = (1 << sigma);

		for (SizeType i = 0; i < n; ++i) {
			SizeType ket = basis[i];
			SizeType bra = ket;
			if (ket & mask) {
			} else {
				bra ^= mask;
				int j = PsimagLite::indexOrMinusOne(basis, bra);
				if (j < 0)
					err("findCmatrix\n");
				cm(i, j) = sign(ket, sigma);
			}
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}


	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	RealType sign(const SizeType& ket, SizeType sigma) const
	{
		return (sigma == 1 && (ket & 1)) ? -1 : 1;
	}

	//! Find S^+ in the natural basis natBasis
	SparseMatrixType findSplusMatrix(const VectorSizeType& basis) const
	{
		SizeType site = 0;
		SizeType total = basis.size();
		MatrixType cm(total, total, 0);
		RealType j = 0.5*modelParams_.twiceTheSpin;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParams_.twiceTheSpin);
		SizeType bits = 1 + ProgramGlobals::logBase2(modelParams_.twiceTheSpin);
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;

		for (SizeType i = 0; i < total; ++i) {
			SizeType ket = basis[i];
			// save and discard electronic part
			SizeType electronic = ket & 3;
			ket >>= 2;

			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);

			assert(ketsite == ket);
			SizeType brasite = ketsite + 1;
			if (brasite >= modelParams_.twiceTheSpin + 1) continue;

			SizeType bra = ket & (~mask);
			assert(bra == 0);
			brasite <<= (site*bitsForOneSite);
			bra |= brasite;
			assert(bra == brasite);

			RealType m = ketsite - j;
			RealType x = j*(j+1)-m*(m+1);
			assert(x>=0);

			// restore electronic part to bra
			bra <<= 2;
			bra |= electronic;
			int j = PsimagLite::indexOrMinusOne(basis, bra);
			if (j < 0)
				err("findCmatrix\n");
			cm(i, j) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrix(const VectorSizeType& basis) const
	{
		SizeType site = 0;
		SizeType total = basis.size();
		MatrixType cm(total, total, 0);
		RealType j = 0.5*modelParams_.twiceTheSpin;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(modelParams_.twiceTheSpin);
		SizeType bits = ProgramGlobals::logBase2(modelParams_.twiceTheSpin) + 1;
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site*bitsForOneSite);

		for (SizeType i = 0; i < total; ++i) {
			SizeType ket = basis[i];
			// discard electronic part
			ket >>= 2;

			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);
			assert(ketsite == ket);
			RealType m = ketsite - j;

			cm(i, i) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find n in the natural basis natBasis
	SparseMatrixType findNmatrix(const VectorSizeType& basis) const
	{
		const ComplexOrRealType zero = 0.0;
		SizeType n = basis.size();
		MatrixType cm(n, n, zero);
		VectorSizeType mask(2, 0);
		mask[0] = 1;
		mask[1] = 2;
		for (SizeType i = 0; i < n; ++i) {
			SizeType ket = basis[i];
			for (SizeType sigma = 0; sigma < 2; ++sigma) {
				if (ket & mask[sigma]) cm(i, i) += 1.0;
			}
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	SparseMatrixType kondoOnSite(SizeType,
	                             const SparseMatrixType& niup,
	                             const SparseMatrixType& nidown) const
	{
		// cdu[d] is actually cu[d] not cu[d] dagger.
		const SparseMatrixType& cdu = ops_[0].getCRS();
		const SparseMatrixType& cdd = ops_[1].getCRS();
		const SparseMatrixType& Sp = ops_[2].getCRS();
		const SparseMatrixType& Sz = ops_[3].getCRS();

		SparseMatrixType Sm;
		transposeConjugate(Sm, Sp);

		SparseMatrixType sz = niup;
		const RealType minusOne = -1.0;
		const RealType zeroPointFive = 0.5;
		sz += minusOne*nidown;
		sz *= zeroPointFive;

		SparseMatrixType cu;
		transposeConjugate(cu, cdu);

		SparseMatrixType cd;
		transposeConjugate(cd, cdd);

		SparseMatrixType sp;
		multiply(sp, cdu, cd);

		SparseMatrixType sm;
		transposeConjugate(sm, sp);
#ifndef NDEBUG
		SparseMatrixType smtest;
		multiply(smtest, cdd, cu);
		assert(smtest == sm);
#endif

		SparseMatrixType m = sp*Sm;
		m += sm*Sp;
		m *= zeroPointFive;
		m += sz*Sz;
		return m;
	}

	const SolverParamsType& solverParams_;
	ParametersKondoType modelParams_;
	VectorSizeType basis_;
	VectorQnType qn_;
	VectorOperatorType ops_;
};
}

#endif
