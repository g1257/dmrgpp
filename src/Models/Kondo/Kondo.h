#ifndef DMRGPP_MODEL_KONDO_H
#define DMRGPP_MODEL_KONDO_H
#include "ParametersKondo.h"
#include "LinkProductKondo.h"

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
	typedef ParametersKondo<RealType, QnType> ParametersKondoType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef LinkProductKondo<ModelHelperType, GeometryType> LinkProductType;
	typedef typename LinkProductType::PairSizeType PairSizeType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;

public:

	Kondo(const SolverParamsType& solverParams,
	      InputValidatorType& io,
	      GeometryType const &geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    new LinkProductType(io),
	                    io),
	      solverParams_(solverParams),
	      geometry_(geometry),
	      modelParams_(io)
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

	// START Functions that each model MUST implement

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
		SizeType linSize = geometry_.numberOfSites();
		SparseMatrixType tmpMatrix;
		SparseMatrixType niup;
		SparseMatrixType nidown;

		for (SizeType i = 0; i < n; ++i) {
			SizeType ind = block[i];
			const VectorOperatorType& ops = ModelBaseType::trackableOps(ind);
			assert(ops_.size() > 1 + i*4);
			// onsite U hubbard
			//n_i up
			SizeType sigma =0; // up sector
			transposeConjugate(tmpMatrix, ops[sigma + i*4].data);
			multiply(niup, tmpMatrix, ops[sigma + i*4].data);
			//n_i down
			sigma = 1; // down sector
			transposeConjugate(tmpMatrix, ops[sigma + i*4].data);
			multiply(nidown,tmpMatrix, ops[sigma + i*4].data);

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
			hmatrix += modelParams_.kondoJ[ind] * kondoOnSite(ops, niup, nidown);
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		qns = qn_;
		assert(ops_.size() >= 6);
		OpsLabelType& c = this->createOpsLabel("c");
		for (SizeType sigma = 0; sigma < 2; ++sigma)
			c.push(ops_[sigma]);

		this->createOpsLabel("Splus").push(ops_[2]);
		this->createOpsLabel("Sz").push(ops_[3]);
		this->createOpsLabel("n").push(ops_[4]);

		this->makeTrackableOrderMatters("c");
		this->makeTrackableOrderMatters("Splus");
		this->makeTrackableOrderMatters("Sz");
		this->makeTrackableOrderMatters("n");

		{
			SparseMatrixType mup = ops_[0].data;
			SparseMatrixType mupT;
			transposeConjugate(mupT, mup);
			SparseMatrixType mdown = ops_[1].data;
			SparseMatrixType mdownT;
			transposeConjugate(mdownT, mdown);
			SparseMatrixType szMatrix = mdownT*mdown;
			szMatrix *= (-1.0);
			szMatrix += mupT*mup;
			szMatrix *= 0.5;

			PairSizeType zeroPair(0, 0);
			typename OperatorType::Su2RelatedType su2Related;
			this->createOpsLabel("sz").push(OperatorType(szMatrix,
			                                             1,
			                                             zeroPair,
			                                             1,
			                                             su2Related));
		}
	}

	// END ^^^^^^^^^^^Functions that each model needs to implement

private:

	void setSymmetryRelatedInternal(VectorQnType& qns,
	                                const VectorSizeType& basis) const
	{
		qns.resize(basis.size(), QnType::zero());
		VectorSizeType other(2, 0);

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

			other[0] = electrons;
			other[1] = electronsUp + bosonicSz;

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
			                  -1,
			                  typename OperatorType::PairType(0, 0),
			                  1,
			                  su2related);

			ops_.push_back(myOp);
		}

		// now the S+ and Sz for local spins
		SparseMatrixType m = findSplusMatrix(basis_);
		OperatorType sp(m,
		                1,
		                typename OperatorType::PairType(0, 0),
		                1,
		                su2related);
		ops_.push_back(sp);

		m = findSzMatrix(basis_);
		OperatorType sz(m,
		                1,
		                typename OperatorType::PairType(0, 0),
		                1,
		                su2related);
		ops_.push_back(sz);

		m = findNmatrix(basis_);
		OperatorType nm(m,
		                1,
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

	SparseMatrixType kondoOnSite(const VectorOperatorType& ops,
	                             const SparseMatrixType& niup,
	                             const SparseMatrixType& nidown) const
	{
		const SparseMatrixType& cdu = ops[0].data;
		const SparseMatrixType& cdd = ops[1].data;
		const SparseMatrixType& Sp = ops[2].data;
		const SparseMatrixType& Sz = ops[3].data;

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
	const GeometryType& geometry_;
	ParametersKondoType modelParams_;
	VectorSizeType basis_;
	VectorQnType qn_;
	VectorOperatorType ops_;
};
}

#endif
