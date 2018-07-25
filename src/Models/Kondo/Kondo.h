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

public:

	Kondo(const SolverParamsType& solverParams,
	      InputValidatorType& io,
	      GeometryType const &geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    new LinkProductType,
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

		setOperatorMatricesInternal();
	}

	// START Functions that each model MUST implement

	// For information purposes only. Write model parameters
	// String contains the group
	// Serializer object is second argument
	void write(PsimagLite::String,
	           PsimagLite::IoNg::Out::Serializer&) const
	{
		err("write unimplemented\n");
	}

	// Given an operator what with degree of freedom dof
	// create the one-site matrix for this operator
	// and create a OperatorType object from it and return it
	// site MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		err("Unimplemented naturalOperator\n");
		SizeType h = qn_.size();
		SparseMatrixType m(h, h);
		m.makeDiagonal(h, 1.0);
		typename OperatorType::Su2RelatedType su2Related;
		return OperatorType(m,
		                    1.0,
		                    typename OperatorType::PairType(0,0),
		                    1.0,
		                    su2Related);
	}

	// Return the size of the one-site Hilbert space basis for this model
	// site MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	SizeType hilbertSize(SizeType site) const
	{
		return qn_.size();
	}

	// Fill the VectorOperatorType with operators that need to be kept
	// track by the DMRG++ Engine.
	// Fill VectorQnType with the qns of the one site basis in the order
	// you chose to give the operators
	// You can check that block.size() == 1 or throw otherwise
	// The contents of block MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	void setOperatorMatrices(VectorOperatorType& ops,
	                         VectorQnType& qn,
	                         const VectorSizeType& block) const
	{
		assert(block.size() == 1);
		ops = ops_;
		qn = qn_;
	}

	// Fill SparseMatrixType with the on-site Hamiltonian terms in the on-site basis
	// Give SparseMatrixType in the order you chose to give the
	// operators in setOperatorMatrices
	// The RealType contain the physical time in case your onsite terms
	// depend on it
	void addDiagonalsInNaturalBasis(SparseMatrixType&,
	                                const VectorOperatorType&,
	                                const VectorSizeType& block,
	                                RealType)  const
	{
		err("addDiagonalsInNaturalBasis unimplemented\n");
	}

	// END ^^^^^^^^^^^Functions that each model needs to implement

private:

	void setSymmetryRelatedInternal(VectorQnType& qns,
	                                const VectorSizeType& basis) const
	{
		qns.resize(basis.size(), ModelBaseType::QN_ZERO);
		SizeType mode = ModelBaseType::targetQuantum().qn.other.size();
		assert(mode == 1); // Sz only as other
		VectorSizeType other(1, 0);

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

			SizeType bosonicSz = (basis[i] & 4) ? 1 : 0;

			other[0] = electronsUp + bosonicSz;

			qns[i] = QnType(electrons, other, jmpair, 0);
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

	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findCmatrix(SizeType sigma,
	                             const VectorSizeType& basis) const
	{
		SizeType n = basis.size();
		MatrixType cm(n, n, static_cast<ComplexOrRealType>(0.0));
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

	const SolverParamsType& solverParams_;
	const GeometryType& geometry_;
	ParametersKondoType modelParams_;
	VectorSizeType basis_;
	VectorQnType qn_;
	VectorOperatorType ops_;
};
}


#endif
