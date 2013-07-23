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

/*! \file ModelBase.h
 *
 *
 */

#ifndef MODEL_BASE_H
#define MODEL_BASE_H

#include "ReflectionOperatorEmpty.h"
#include "LinkProductStruct.h"

namespace Dmrg {

template<typename ModelHelperType_,
         typename ParametersType,
         typename InputValidatorType_,
         typename GeometryType_>
class ModelBase  {

public:

	typedef InputValidatorType_ InputValidatorType;
	typedef typename PsimagLite::Vector<unsigned int long long>::Type HilbertBasisType;
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
	typedef typename MyBasis::BasisDataType BasisDataType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::SparseElementType ComplexOrRealType;
	typedef LinkProductStruct<ComplexOrRealType> LinkProductStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef ParametersType SolverParamsType;
	typedef typename ModelHelperType::LinkType LinkType;

	ModelBase(const ParametersType& params,
	          InputValidatorType& io,
	          const GeometryType& geometry)
	    : params_(params),
	      geometry_(geometry)
	{}

	virtual ~ModelBase() {}

	const GeometryType& geometry() const { return geometry_; }

	const ParametersType& params() const { return params_; }

	virtual void setNaturalBasis(VectorOperatorType& creationMatrix,
	                             SparseMatrixType &hamiltonian,
	                             BasisDataType& q,
	                             const BlockType& block,
	                             const RealType& time) const = 0;

	virtual PsimagLite::Matrix<ComplexOrRealType>
	naturalOperator(const PsimagLite::String& what,
	                SizeType site,
	                SizeType dof) const = 0;

	virtual void findElectrons(typename PsimagLite::Vector<SizeType> ::Type& electrons,
	                           const HilbertBasisType& basis,
	                           SizeType site) const = 0;

	virtual void print(std::ostream& os) const = 0;

	virtual void matrixVectorProduct(VectorType& x,
	                                 const VectorType& y,
	                                 ModelHelperType const &modelHelper) const = 0;

	virtual void addHamiltonianConnection(SparseMatrixType &matrix,
	                                      const LeftRightSuperType& lrs) const = 0;

	virtual void hamiltonianConnectionProduct(VectorType& x,
	                                          const VectorType& y,
	                                          ModelHelperType const &modelHelper) const = 0;

	virtual void fullHamiltonian(SparseMatrixType& matrix,
	                             const ModelHelperType& modelHelper) const = 0;

	virtual SizeType hilbertSize(SizeType site) const = 0;

	virtual void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                                 const BlockType& block) const = 0;

	virtual void setNaturalBasis(HilbertBasisType& basis,
	                             typename PsimagLite::Vector<SizeType>::Type& q,
	                             const BlockType& block) const = 0;

	virtual void calcHamiltonian(SparseMatrixType &hmatrix,
	                             const VectorOperatorType& cm,
	                             const BlockType& block,
	                             RealType time,
	                             RealType factorForDiagonals) const = 0;

	virtual SizeType getLinkProductStruct(LinkProductStructType** lps,
	                                      const ModelHelperType& modelHelper) const = 0;

	virtual LinkType getConnection(const SparseMatrixType** A,
	                               const SparseMatrixType** B,
			                       SizeType ix,
			                       const LinkProductStructType& lps,
			                       const ModelHelperType& modelHelper) const = 0;

	virtual void findElectronsOfOneSite(BlockType& electrons,SizeType site) const
	{
		typename PsimagLite::Vector<SizeType>::Type block(1,site);
		HilbertBasisType basis;
		typename PsimagLite::Vector<SizeType>::Type quantumNumbs;
		setNaturalBasis(basis,quantumNumbs,block);
		findElectrons(electrons,basis,site);
	}

	virtual void hamiltonianOnLink(SparseMatrixType& hmatrix,
	                               const BlockType& block,
	                               const RealType& time,
	                               RealType factorForDiagonals) const
	{
		assert(block.size()==2);

		typename PsimagLite::Vector<OperatorType>::Type cm;
		setOperatorMatrices(cm,block);
		calcHamiltonian(hmatrix,cm,block,time,factorForDiagonals);
	}

private:

	const ParametersType& params_;
	const GeometryType& geometry_;

};     //class ModelBase

template<typename ModelHelperType,
         typename ParametersType,
         typename InputValidatorType,
         typename GeometryType>
std::ostream& operator<<(std::ostream& os,
                         const ModelBase<ModelHelperType,
                                         ParametersType,
                                         InputValidatorType,
                                         GeometryType>& model)
{
	model.print(os);
	return os;
}
} // namespace Dmrg
/*@}*/
#endif
