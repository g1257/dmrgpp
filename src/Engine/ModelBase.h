/*
Copyright (c) 2009-2018, UT-Battelle, LLC
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

/*! \file ModelBase.h
 *
 *
 */

#ifndef MODEL_BASE_H
#define MODEL_BASE_H

#include "ReflectionOperatorEmpty.h"
#include "Vector.h"
#include "Sort.h"
#include "MemResolv.h"
#include "TargetQuantumElectrons.h"
#include "Io/IoSerializerStub.h"
#include "ModelCommon.h"
#include "NotReallySort.h"
#include "LabeledOperators.h"

namespace Dmrg {

template<typename ModelHelperType_,
         typename ParametersType_,
         typename InputValidatorType_,
         typename GeometryType_>
class ModelBase  {

public:

	typedef ParametersType_ ParametersType;
	typedef InputValidatorType_ InputValidatorType;
	typedef ModelHelperType_ ModelHelperType;
	typedef GeometryType_ GeometryType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::BasisType MyBasis;
	typedef typename ModelHelperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef ReflectionOperatorEmpty<LeftRightSuperType> ReflectionSymmetryType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename MyBasis::QnType QnType;
	typedef TargetQuantumElectrons<RealType, QnType> TargetQuantumElectronsType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::SparseElementType ComplexOrRealType;
	typedef ModelCommon<ParametersType, GeometryType, ModelHelperType> ModelCommonType;
	typedef typename ModelCommonType::HamiltonianConnectionType HamiltonianConnectionType;
	typedef typename ModelCommonType::VectorLinkType VectorLinkType;
	typedef typename ModelCommonType::LinkProductBaseType LinkProductBaseType;
	typedef typename ModelCommonType::VectorType VectorType;
	typedef ParametersType SolverParamsType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef VectorSizeType HilbertBasisType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename QnType::PairSizeType PairSizeType;
	typedef typename BasisWithOperatorsType::VectorBoolType VectorBoolType;

	ModelBase(const ParametersType& params,
	          const GeometryType_& geometry,
	          const LinkProductBaseType* lpb,
	          InputValidatorType& io)
	    : modelCommon_(params, geometry, lpb),
	      targetQuantum_(io)
	{}

	void postCtor()
	{
		VectorOperatorType cm;
		VectorQnType qns;
		BlockType block(1, 0);
		setOperatorMatrices(cm, qns, block);
		modelCommon_.postCtor(cm);
	}

	virtual ~ModelBase() {}

	// START Functions that each model MUST implement

	// For information purposes only. Write model parameters
	// String contains the group
	// Serializer object is second argument
	virtual void write(PsimagLite::String,
	                   PsimagLite::IoNg::Out::Serializer&) const = 0;

	// Given an operator what with degree of freedom dof
	// create the one-site matrix for this operator
	// and create a OperatorType object from it and return it
	// site MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	virtual OperatorType naturalOperator(const PsimagLite::String& what,
	                                     SizeType site,
	                                     SizeType dof) const = 0;

	// Return the size of the one-site Hilbert space basis for this model
	// site MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	virtual SizeType hilbertSize(SizeType site) const = 0;

	// Fill the VectorOperatorType with operators that need to be kept
	// track by the DMRG++ Engine.
	// Fill VectorQnType with the qns of the one site basis in the order
	// you chose to give the operators
	// You can check that block.size() == 1 or throw otherwise
	// The contents of block MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	virtual void setOperatorMatrices(VectorOperatorType&,
	                                 VectorQnType&,
	                                 const BlockType& block) const = 0;

	// Fill SparseMatrixType with the on-site Hamiltonian terms in the on-site basis
	// Give SparseMatrixType in the order you chose to give the
	// operators in setOperatorMatrices
	// The RealType contain the physical time in case your onsite terms
	// depend on it
	virtual void addDiagonalsInNaturalBasis(SparseMatrixType&,
	                                        const VectorOperatorType&,
	                                        const BlockType& block,
	                                        RealType)  const = 0;

	// END ^^^^^^^^^^^Functions that each model needs to implement

	virtual const LinkProductBaseType& linkProduct() const
	{
		return modelCommon_.linkProduct();
	}

	virtual void findOddElectronsOfOneSite(VectorBoolType& oddElectrons,
	                                       SizeType site) const
	{
		typename PsimagLite::Vector<SizeType>::Type block(1, site);
		typename PsimagLite::Vector<OperatorType>::Type cm;
		VectorQnType qq;
		setOperatorMatrices(cm, qq, block);
		SizeType n = qq.size();
		oddElectrons.resize(n);
		for (SizeType i = 0; i < n; ++i)
			oddElectrons[i] = qq[i].oddElectrons;
	}

	virtual void matrixVectorProduct(VectorType& x,
	                                 const VectorType& y,
	                                 const HamiltonianConnectionType& hc) const
	{
		return modelCommon_.matrixVectorProduct(x, y, hc);
	}

	virtual void addHamiltonianConnection(SparseMatrixType &matrix,
	                                      const LeftRightSuperType& lrs,
	                                      RealType currentTime) const
	{
		return modelCommon_.addHamiltonianConnection(matrix,lrs,currentTime);
	}

	virtual void fullHamiltonian(SparseMatrixType& matrix,
	                             const HamiltonianConnectionType& hc) const
	{
		return modelCommon_.fullHamiltonian(matrix, hc);
	}

	//! Full hamiltonian from creation matrices cm
	virtual void calcHamiltonian(SparseMatrixType &hmatrix,
	                             const VectorOperatorType& cm,
	                             const BlockType& block,
	                             RealType time)  const
	{
		hmatrix.makeDiagonal(cm[0].data.rows());

		modelCommon_.addConnectionsInNaturalBasis(hmatrix,cm,block,time);

		addDiagonalsInNaturalBasis(hmatrix,cm,block,time);
	}

	virtual SizeType maxElectronsOneSpin() const
	{
		SizeType tmp = hilbertSize(0);
		tmp = static_cast<SizeType>(log(tmp)/log(2.0));
		SizeType maxElectrons = static_cast<SizeType>(tmp/2);
		if (tmp & 1) maxElectrons++;

		return maxElectrons*modelCommon_.geometry().numberOfSites() + 1;
	}

	virtual PsimagLite::String symmName() const { return "undefined"; }

	virtual bool instrospect() const
	{
		std::cerr<<"Instrospection not available for this model\n";
		return false;
	}

	void printBasis(SizeType site) const
	{
		BlockType block(1, site);
		typename PsimagLite::Vector<OperatorType>::Type cm;
		VectorQnType qq;
		setOperatorMatrices(cm, qq, block);
		std::cout<<"block="<<block;
		std::cout<<"qq="<<qq;

		SizeType n = cm.size();
		for (SizeType i = 0; i < n; ++i) {
			std::cout<<"Matrix "<<i<<"\n";
			std::cout<<cm[i];
		}
	}

	const GeometryType& geometry() const { return modelCommon_.geometry(); }

	const ParametersType& params() const { return modelCommon_.params(); }

	static void checkNaturalOperatorDof(SizeType dof,
	                                    PsimagLite::String label,
	                                    const VectorSizeType& allowed)
	{
		if (std::find(allowed.begin(),allowed.end(),dof) != allowed.end()) return;
		PsimagLite::String str("For this model and label=");
		str += label + " dof=" + ttos(dof) + " is not allowed\n";
		str += "Allowed dofs are ";
		for (SizeType i = 0; i < allowed.size(); ++i)
			str += ttos(i) + " ";
		str += "\n";
		throw PsimagLite::RuntimeError(str);
	}

	const TargetQuantumElectronsType& targetQuantum() const
	{
		return targetQuantum_;
	}

	static void orderByQuantum(VectorSizeType& basis, VectorQnType& qn)
	{
		VectorSizeType newBasis;
		VectorSizeType partition;
		VectorQnType qns;
		NotReallySort notReallySort;
		notReallySort(newBasis,
		              qns,
		              partition,
		              basis,
		              qn,
		              false,
		              10,
		              ProgramGlobals::VERBOSE_NO);

		SizeType n = partition.size();
		assert(n > 0);
		--n;
		assert(n == qns.size());
		for (SizeType i = 0; i < n; ++i)
			for (SizeType j = partition[i]; j < partition[i + 1]; ++j)
				qn[j] = qns[i];

		basis = newBasis;
	}

private:

	ModelCommonType modelCommon_;
	TargetQuantumElectronsType targetQuantum_;
};     //class ModelBase
} // namespace Dmrg
/*@}*/
#endif
