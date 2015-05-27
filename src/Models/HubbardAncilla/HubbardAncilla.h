/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

/*! \file HubbardAncilla.h
 *
 *
 */
#ifndef DMRG_HUBBARD_ANCILLA_H
#define DMRG_HUBBARD_ANCILLA_H
#include "ModelCommon.h"
#include "LinkProductHubbardAncilla.h"
#include "ParametersHubbardAncilla.h"

namespace Dmrg {

template<typename ModelBaseType>
class HubbardAncilla : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef ModelFeBasedSc<ModelBaseType> ModelFeAsType;
	typedef typename ModelFeAsType::MatrixType MatrixType;
	typedef typename ModelFeAsType::HilbertState HilbertState;
	typedef typename ModelFeAsType::HilbertBasisType HilbertBasisType;
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
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef LinkProductHubbardAncilla<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef	 typename ModelBaseType::MyBasis MyBasis;
	typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename MyBasis::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;

	static const SizeType SPIN_UP = ModelFeAsType::SPIN_UP;
	static const SizeType SPIN_DOWN = ModelFeAsType::SPIN_DOWN;

	HubbardAncilla(const SolverParamsType& solverParams,
	               InputValidatorType& io,
	               GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry),
	      modelFeAs_(solverParams,io,geometry),
	      orbitals_(modelParameters_.orbitals)
	{
		if (orbitals_ != 2)
			throw PsimagLite::RuntimeError("HubbardAncilla: only for 2 orbitals\n");
	}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType site) const { return modelFeAs_.hilbertSize(site); }

	void print(std::ostream& os) const { modelFeAs_.print(os); }

	//! find creation operator matrices for (i,sigma) in the natural basis,
	//! find quantum numbers and number of electrons
	//! for each state in the basis
	void setNaturalBasis(VectorOperatorType& creationMatrix,
	                     SparseMatrixType &hamiltonian,
	                     SymmetryElectronsSzType &q,
	                     BlockType const &block,
	                     const RealType& time)  const
	{
		blockIsSize1OrThrow(block);

		modelFeAs_.setNaturalBasis(creationMatrix,hamiltonian,q,block,time);

		// add \Gamma^+_i to creationMatrix
		setGammaMatrix(creationMatrix,block);

		// add ONSITE J S^+_i S^-_i + S^-_i S^+_i to Hamiltonian
		addSplusSminus(hamiltonian,creationMatrix,block);
	}

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         BlockType const &block) const
	{
		blockIsSize1OrThrow(block);

		modelFeAs_.setOperatorMatrices(creationMatrix,block);

		// add \Gamma^+_i to creationMatrix
		setGammaMatrix(creationMatrix,block);
	}

	PsimagLite::Matrix<SparseElementType> naturalOperator(const PsimagLite::String& what,
	                                                      SizeType site,
	                                                      SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatrices(creationMatrix,block);

		if (what=="g") {
			PsimagLite::Matrix<SparseElementType> tmp;
			SizeType x = 2*orbitals_;
			crsMatrixToFullMatrix(tmp,creationMatrix[x].data);
			return tmp;
		}

		return modelFeAs_.naturalOperator(what,site,dof);
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setNaturalBasis(typename PsimagLite::Vector<HilbertState>  ::Type&basis,
	                     typename PsimagLite::Vector<SizeType>::Type& q,
	                     const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		modelFeAs_.setNaturalBasis(basis,q,block);
	}

	void findElectrons(typename PsimagLite::Vector<SizeType>::Type& electrons,
	                   const typename PsimagLite::Vector<HilbertState>::Type& basis,
	                   SizeType site) const
	{
		modelFeAs_.findElectrons(electrons,basis,site);
	}

	virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                        const VectorOperatorType& cm,
	                                        const BlockType& block,
	                                        RealType time,
	                                        RealType factorForDiagonals=1.0)  const
	{
		modelFeAs_.addDiagonalsInNaturalBasis(hmatrix,cm,block,time,factorForDiagonals);
		// add ONSITE J S^+_i S^-_i + S^-_i S^+_i to Hamiltonian
		SparseMatrixType tmpMatrix;
		addSplusSminus(tmpMatrix,cm,block);
		hmatrix += factorForDiagonals*tmpMatrix;
	}

	virtual SizeType maxElectronsOneSpin() const
	{
		return 2 * modelParameters_.orbitals * geometry_.numberOfSites() + 1;
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelFeAs_.targetQuantum();
	}

private:

	void setGammaMatrix(VectorOperatorType& creationMatrix,
	                    const BlockType& block)  const
	{
		blockIsSize1OrThrow(block);
		assert(creationMatrix.size() == 2*modelParameters_.orbitals);

		MatrixType m1 = creationMatrix[0].data.toDense();
		MatrixType m2 = creationMatrix[1].data.toDense();
		MatrixType m3 = creationMatrix[2].data.toDense();
		MatrixType m4 = creationMatrix[3].data.toDense();
		MatrixType m = m1*transposeConjugate(m2) + m3*transposeConjugate(m4);

		Su2RelatedType su2related;
		SparseMatrixType tmpMatrix(m);
		OperatorType myOp(tmpMatrix,
		                  1,
		                  typename OperatorType::PairType(0,0),
		                  1,
		                  su2related);

		creationMatrix.push_back(myOp);
	}

	void addSplusSminus(SparseMatrixType &hmatrix,
	                    const VectorOperatorType& cm,
	                    const BlockType& block) const
	{
		blockIsSize1OrThrow(block);
		MatrixType m1 = cm[0].data.toDense();
		MatrixType m2 = cm[1].data.toDense();
		MatrixType m3 = cm[2].data.toDense();
		MatrixType m4 = cm[3].data.toDense();
		MatrixType m = m1*transposeConjugate(m3)*m4*transposeConjugate(m2);
		MatrixType m5 = modelParameters_.ancillaJ*(m + transposeConjugate(m));
		fullMatrixToCrsMatrix(hmatrix,m5);
	}

	void blockIsSize1OrThrow(const BlockType& block) const
	{
		if (block.size()==1) return;
		throw PsimagLite::RuntimeError("FeAsBasedExtended:: blocks must be of size 1\n");
	}

	ParametersHubbardAncilla<RealType>  modelParameters_;
	const GeometryType& geometry_;
	ModelFeAsType modelFeAs_;
	SizeType orbitals_;
}; //class HubbardAncilla

} // namespace Dmrg
/*@}*/
#endif // DMRG_HUBBARD_ANCILLA_H

