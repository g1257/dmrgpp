/*
Copyright (c) 2009, UT-Battelle, LLC
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

/*! \file FeBasedScExtedned.h
 *
 *  An implementation of a Hubbard model for Fe-based superconductors
 *  to use with the DmrgSolver
 *  This extends the FeAsBasedSc model to include JNN and JNNN couplings
 *  FIXME: Merge into FeAsBasedSc
 *
 */
#ifndef FEAS_BASED_SC_EX
#define FEAS_BASED_SC_EX
#include "../Models/FeAsModel/ModelFeBasedSc.h"

namespace Dmrg {

template<typename ModelBaseType>
class FeAsBasedScExtended : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef ModelFeBasedSc<ModelBaseType> ModelFeAsType;
	typedef typename ModelFeAsType::HilbertState HilbertState;
	typedef typename ModelFeAsType::HilbertBasisType HilbertBasisType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename OperatorType::Su2RelatedType Su2RelatedType;
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename BasisType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const SizeType SPIN_UP = ModelFeAsType::SPIN_UP;
	static const SizeType SPIN_DOWN = ModelFeAsType::SPIN_DOWN;

	FeAsBasedScExtended(const SolverParamsType& solverParams,
	                    InputValidatorType& io,
	                    const GeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      geometry_(geometry),
	      modelFeAs_(solverParams,io,geometry),
	      orbitals_(modelParameters_.orbitals)
	{}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		modelFeAs_.write(label, io);
		io.write(label + "/orbitals_", orbitals_);
	}

	virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                        const BlockType& block,
	                                        RealType time)  const
	{
		modelFeAs_.addDiagonalsInNaturalBasis(hmatrix, block, time);
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		modelFeAs_.fillLabeledOperators(qns);
		SizeType site = 0;
		BlockType block(1, site);
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatricesInternal(creationMatrix, qns, block);

		OpsLabelType& splus = this->createOpsLabel("naturalSplus");
		SizeType x = 2*orbitals_;
		assert(x < creationMatrix.size());
		splus.push(creationMatrix[x]);

		OpsLabelType& sminus = this->createOpsLabel("naturalSminus");
		x = 2*orbitals_;
		assert(x < creationMatrix.size());
		creationMatrix[x].dagger();
		sminus.push(creationMatrix[x]);
		creationMatrix[x].dagger();

		OpsLabelType& sz = this->createOpsLabel("naturalSz");
		x = 2*orbitals_ + 1;
		assert(x < creationMatrix.size());
		sz.push(creationMatrix[x]);

		this->makeTrackable("naturalSplus");
		this->makeTrackable("naturalSz");
	}

	void fillModelLinks()
	{
		modelFeAs_.fillModelLinks();

		bool isSu2 = BasisType::useSu2Symmetry();

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");

		OpForLinkType splus("naturalSplus");

		auto valueModiferTerm0 = [isSu2](ComplexOrRealType& value)
		{ value *= (isSu2) ? -0.5 : 0.5;};

		spsm.push(splus, 'N', splus, 'C', 2, -1, 2, valueModiferTerm0);

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");

		if (!isSu2) {
			OpForLinkType sz("naturalSz");
			szsz.push(sz, 'N', sz, 'N', 2, 0.5);
		} else {
			auto valueModifierTermOther = [isSu2](ComplexOrRealType& value)
			{ if (isSu2) value = -value;};
			spsm.push(splus, 'N', splus, 'C', 2, -1, 2, valueModifierTermOther);
		}
	}

private:

	//! set creation matrices for sites in block
	void setOperatorMatricesInternal(VectorOperatorType& creationMatrix,
	                                VectorQnType& qns,
	                                const BlockType& block) const
	{
		blockIsSize1OrThrow(block);

		modelFeAs_.setQns(qns);
		modelFeAs_.setOperatorMatricesInternal(creationMatrix, block);

		// add S^+_i to creationMatrix
		setSplus(creationMatrix,block);

		// add S^z_i to creationMatrix
		setSz(creationMatrix,block);
	}

	// add S^+_i to creationMatrix
	void setSplus(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
	              const BlockType& block) const
	{
		SparseMatrixType m;
		cDaggerC(m,creationMatrix,block,1.0,SPIN_UP,SPIN_DOWN);
		Su2RelatedType su2related;
		SizeType offset = 2*orbitals_;
		su2related.source.push_back(offset);
		su2related.source.push_back(offset+1);
		su2related.source.push_back(offset);
		su2related.transpose.push_back(-1);
		su2related.transpose.push_back(-1);
		su2related.transpose.push_back(1);
		su2related.offset = 1;

		OperatorType sPlus(m,
		                   ProgramGlobals::BOSON,
		                   typename OperatorType::PairType(2, 2),
		                   -1,
		                   su2related);
		creationMatrix.push_back(sPlus);
	}

	// add S^z_i to creationMatrix
	void setSz(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
	           const BlockType& block) const
	{
		SparseMatrixType m1,m2;
		cDaggerC(m1,creationMatrix,block,0.5,SPIN_UP,SPIN_UP);
		cDaggerC(m2,creationMatrix,block,-0.5,SPIN_DOWN,SPIN_DOWN);
		Su2RelatedType su2related2;
		SparseMatrixType m = m1;
		m += m2;
		OperatorType sz(m,
		                ProgramGlobals::BOSON,
		                typename OperatorType::PairType(2 ,1),
		                1.0/sqrt(2.0),
		                su2related2);
		creationMatrix.push_back(sz);
	}

	// add S^+_i to creationMatrix
	void cDaggerC(SparseMatrixType& sum,
	              const typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
	              const BlockType&,
	              RealType value,
	              SizeType spin1,
	              SizeType spin2) const
	{
		SparseMatrixType tmpMatrix,tmpMatrix2;
		for (SizeType orbital=0;orbital<orbitals_;orbital++) {
			transposeConjugate(tmpMatrix2,
			                   creationMatrix[orbital+spin2*orbitals_].data);
			multiply(tmpMatrix,
			         creationMatrix[orbital+spin1*orbitals_].data,
			        tmpMatrix2);

			if (orbital == 0) sum = value*tmpMatrix;
			else sum += value*tmpMatrix;
		}
	}

	// add J_{ij} S^+_i S^-_j + S^-_i S^+_j to Hamiltonia
	void addSplusSminus(SparseMatrixType &,
	                    const typename PsimagLite::Vector<OperatorType> ::Type&,
	                    const BlockType&) const
	{
		// nothing if block.size == 1
	}

	// add J_{ij} S^z_i S^z_j to Hamiltonian
	void addSzSz(SparseMatrixType&,
	             const typename PsimagLite::Vector<OperatorType> ::Type&,
	             const BlockType&) const
	{
		// nothing if block.size == 1
	}

	void blockIsSize1OrThrow(const BlockType& block) const
	{
		if (block.size()==1) return;
		throw PsimagLite::RuntimeError("FeAsBasedExtended:: blocks must be of size 1\n");
	}

	ParametersModelFeAs<RealType, QnType>  modelParameters_;
	const GeometryType& geometry_;
	ModelFeAsType modelFeAs_;
	SizeType orbitals_;
}; //class FeAsBasedScExtended

} // namespace Dmrg
/*@}*/
#endif // FEAS_BASED_SC_EX

