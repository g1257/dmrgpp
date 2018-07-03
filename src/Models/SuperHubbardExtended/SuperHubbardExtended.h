/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file SuperHubbardExtended.h
 *
 *  Hubbard + V_{ij} n_i n_j
 *
 */
#ifndef SUPER_HUBBARD_EXTENDED_H
#define SUPER_HUBBARD_EXTENDED_H
#include "../Models/ExtendedHubbard1Orb/ExtendedHubbard1Orb.h"
#include "LinkProdSuperHubbardExtended.h"
#include "ModelCommon.h"

namespace Dmrg {
//! Extended Hubbard for DMRG solver, uses ExtendedHubbard1Orb by containment
template<typename ModelBaseType>
class SuperHubbardExtended : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef ExtendedHubbard1Orb<ModelBaseType> ModelHubbardType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef TargetQuantumElectrons<RealType> TargetQuantumElectronsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef LinkProdSuperHubbardExtended<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelHubbardType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
	typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename PsimagLite::Vector<HilbertState>::Type VectorHilbertStateType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeTypeType;

	SuperHubbardExtended(const SolverParamsType& solverParams,
	                     InputValidatorType& io,
	                     GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry),
	      extendedHubbard_(solverParams,io,geometry)
	{}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType site) const
	{
		return extendedHubbard_.hilbertSize(site);
	}

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         const BlockType& block) const
	{
		extendedHubbard_.setOperatorMatrices(creationMatrix,block);
		// add S+ and Sz to creationMatrix
		setSpinMatrices(creationMatrix,block);
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		SizeType matrixIndex = findMatrixIndex(what);

		if (matrixIndex < 3)
			return extendedHubbard_.naturalOperator(what,site,dof);

		BlockType block;
		block.resize(1);
		block[0]=site;
		VectorOperatorType creationMatrix;
		setSpinMatrices(creationMatrix,block);
		assert(creationMatrix.size() > 2);
		return creationMatrix[matrixIndex-3];
	}

	//! find total number of electrons for each state in the basis
	void findElectrons(VectorSizeTypeType& electrons,
	                   const VectorHilbertStateType& basis,
	                   SizeType site) const
	{
		extendedHubbard_.findElectrons(electrons,basis,site);
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
		        io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		extendedHubbard_.write(label, io);
	}

	virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                        const VectorOperatorType& cm,
	                                        const BlockType& block,
	                                        RealType time,
	                                        RealType factorForDiagonals=1.0)  const
	{
		extendedHubbard_.addDiagonalsInNaturalBasis(hmatrix,
		                                            cm,
		                                            block,
		                                            time,
		                                            factorForDiagonals);
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

	void setBasis(HilbertBasisType& basis,
	              SymmetryElectronsSzType& qq,
	              const VectorSizeType& block) const
	{
		extendedHubbard_.setBasis(basis, qq, block);
	}

private:

	SizeType findMatrixIndex(const PsimagLite::String& what) const
	{
		if (what == "+") return 3;
		if (what == "z") return 4;
		return 0;
	}

	//! Full hamiltonian from creation matrices cm
	void addSiSj(SparseMatrixType &,
	             const VectorOperatorType&,
	             const BlockType& block) const
	{
		// Assume block.size()==1 and then problem solved!!
		// there are no connection if there's only one site ;-)
		assert(block.size()==1);
	}

	void setSpinMatrices(VectorOperatorType& creationMatrix,
	                     const BlockType& block) const
	{
		assert(block.size()==1);

		SparseMatrixType sPlus = extendedHubbard_.naturalOperator("+",0,0).data;
		RealType angularFactor= 1;
		typename OperatorType::Su2RelatedType su2related;
		su2related.offset = 1; //check FIXME
		OperatorType sPlusOp(sPlus,1,typename OperatorType::PairType(0,0),angularFactor,su2related);
		creationMatrix.push_back(sPlusOp);

		SparseMatrixType sz = extendedHubbard_.naturalOperator("z",0,0).data;
		sz *= 0.5;
		OperatorType szOp(sz,1,typename OperatorType::PairType(0,0),angularFactor,su2related);
		creationMatrix.push_back(szOp);
	}

	ParametersModelHubbard<RealType>  modelParameters_;
	const GeometryType &geometry_;
	ModelHubbardType extendedHubbard_;
};	//class SuperHubbardExtended

} // namespace Dmrg
/*@}*/
#endif // SUPER_HUBBARD_EXTENDED_H
