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
#include "ModelCommonBase.h"
#include "Vector.h"
#include "Sort.h"
#include "MemResolv.h"
#include "TargetQuantumElectrons.h"
#include "Io/IoSerializerStub.h"

namespace Dmrg {

template<typename ModelHelperType_,
         typename ParametersType,
         typename InputValidatorType_,
         typename GeometryType_>
class ModelBase  {

public:

	typedef InputValidatorType_ InputValidatorType;
	typedef typename PsimagLite::Vector<unsigned int long>::Type HilbertBasisType;
	typedef ModelHelperType_ ModelHelperType;
	typedef GeometryType_ GeometryType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelHelperType::RealType RealType;
	typedef TargetQuantumElectrons<RealType> TargetQuantumElectronsType;
	typedef typename ModelHelperType::BasisType MyBasis;
	typedef typename ModelHelperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef ReflectionOperatorEmpty<LeftRightSuperType> ReflectionSymmetryType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::SparseElementType ComplexOrRealType;
	typedef ModelCommonBase<ModelHelperType,ParametersType,GeometryType> ModelCommonBaseType;
	typedef typename ModelCommonBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelCommonBaseType::VectorType VectorType;
	typedef ParametersType SolverParamsType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	ModelBase(InputValidatorType&,
	          ModelCommonBaseType* modelCommon)
	    : modelCommon_(modelCommon)
	{}

	virtual ~ModelBase()
	{
		delete modelCommon_;
	}

	virtual OperatorType naturalOperator(const PsimagLite::String& what,
	                                     SizeType site,
	                                     SizeType dof) const = 0;

	virtual void findElectrons(VectorSizeType& electrons,
	                           const HilbertBasisType& basis,
	                           SizeType site) const = 0;

	virtual void print(std::ostream& os) const = 0;

	virtual SizeType hilbertSize(SizeType site) const = 0;

	virtual void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                                 const BlockType& block) const = 0;

	virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                        const VectorOperatorType& cm,
	                                        const BlockType& block,
	                                        RealType time,
	                                        RealType factorForDiagonals=1.0)  const = 0;

	virtual void findElectronsOfOneSite(BlockType& electrons,SizeType site) const
	{
		typename PsimagLite::Vector<SizeType>::Type block(1, site);
		HilbertBasisType basis;
		SymmetryElectronsSzType qq;
		setBasis(basis, qq, block);
		findElectrons(electrons, basis, site);
	}

	virtual void hamiltonianOnLink(SparseMatrixType& hmatrix,
	                               const BlockType& block,
	                               RealType time,
	                               RealType factorForDiagonals) const
	{
		typename PsimagLite::Vector<OperatorType>::Type cm;
		setOperatorMatrices(cm,block);
		calcHamiltonian(hmatrix,cm,block,time,factorForDiagonals,true);
	}

	virtual void matrixVectorProduct(VectorType& x,
	                                 const VectorType& y,
	                                 ModelHelperType const &modelHelper) const
	{
		return modelCommon_->matrixVectorProduct(x,y,modelHelper);
	}

	virtual void addHamiltonianConnection(SparseMatrixType &matrix,
	                                      const LeftRightSuperType& lrs,
	                                      RealType currentTime) const
	{
		return modelCommon_->addHamiltonianConnection(matrix,lrs,currentTime);
	}

	virtual void hamiltonianConnectionProduct(VectorType& x,
	                                          const VectorType& y,
	                                          ModelHelperType const &modelHelper) const
	{
		return modelCommon_->hamiltonianConnectionProduct(x,y,modelHelper);
	}

	virtual void fullHamiltonian(SparseMatrixType& matrix,
	                             const ModelHelperType& modelHelper) const
	{
		return modelCommon_->fullHamiltonian(matrix,modelHelper);
	}

	virtual SizeType getLinkProductStruct(const ModelHelperType& modelHelper) const
	{
		return modelCommon_->getLinkProductStruct(modelHelper);
	}

	virtual LinkType getConnection(const SparseMatrixType** A,
	                               const SparseMatrixType** B,
	                               SizeType ix,
	                               const ModelHelperType& modelHelper) const
	{
		return modelCommon_->getConnection(A,B,ix,modelHelper);
	}

	//! Full hamiltonian from creation matrices cm
	virtual void calcHamiltonian(SparseMatrixType &hmatrix,
	                             const VectorOperatorType& cm,
	                             const BlockType& block,
	                             RealType time,
	                             RealType factorForDiagonals=1.0,
	                             bool sysEnvOnly=false)  const
	{
		hmatrix.makeDiagonal(cm[0].data.rows());

		modelCommon_->addConnectionsInNaturalBasis(hmatrix,cm,block,sysEnvOnly,time);

		addDiagonalsInNaturalBasis(hmatrix,cm,block,time,factorForDiagonals);
	}

	virtual SizeType maxElectronsOneSpin() const
	{
		SizeType tmp = hilbertSize(0);
		tmp = static_cast<SizeType>(log(tmp)/log(2.0));
		SizeType maxElectrons = static_cast<SizeType>(tmp/2);
		if (tmp & 1) maxElectrons++;

		return maxElectrons*modelCommon_->geometry().numberOfSites() + 1;
	}

	virtual void setBasis(HilbertBasisType& basis,
	                      SymmetryElectronsSzType& qq,
	                      const VectorSizeType& block) const = 0;

	void printBasis(SizeType site) const
	{
		BlockType block(1, site);
		HilbertBasisType natBasis;
		SymmetryElectronsSzType qq;
		setBasis(natBasis, qq, block);
		VectorSizeType q;
		qq.findQuantumNumbers(q, MyBasis::useSu2Symmetry());
		std::cout<<"block="<<block;
		std::cout<<"natBasis="<<natBasis;
		std::cout<<"quantumNumbs="<<q;
	}

	const GeometryType& geometry() const { return modelCommon_->geometry(); }

	const ParametersType& params() const { return modelCommon_->params(); }

	virtual const TargetQuantumElectronsType& targetQuantum() const = 0;

	virtual SizeType memResolv(PsimagLite::MemResolv& mres,
	                           SizeType x,
	                           PsimagLite::String msg = "") const = 0;

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

	virtual PsimagLite::String symmName() const { return "undefined"; }

protected:

	void orderBasis(HilbertBasisType& basis,
	                const HilbertBasisType& basisUnordered,
	                const SymmetryElectronsSzType& qq) const
	{
		VectorSizeType q;
		qq.findQuantumNumbers(q, MyBasis::useSu2Symmetry());
		orderBasis(basis, q, basisUnordered);
	}

private:

	void orderBasis(HilbertBasisType& basis,
	                VectorSizeType& q,
	                const HilbertBasisType& basisTmp) const
	{
		// reorder the natural basis
		VectorSizeType iperm(q.size());
		PsimagLite::Sort<VectorSizeType> sort;
		sort.sort(q,iperm);

		SizeType total = basisTmp.size();
		VectorSizeType basis2(total);
		for (SizeType a=0;a<total;a++)
			basis2[a] = basisTmp[iperm[a]];

		// Ensure deterministic order for the natural basis
		SizeType offset = 0;
		VectorSizeType symmetryBlock;

		basis.resize(total);
		for (SizeType a=0;a<total;a++) {
			if (a>0 && q[a] != q[a-1]) {
				iperm.resize(symmetryBlock.size());
				sort.sort(symmetryBlock,iperm);

				for (SizeType k = 0; k < symmetryBlock.size(); ++k)
					basis[k + offset] = symmetryBlock[k];

				offset += symmetryBlock.size();
				symmetryBlock.clear();
			}

			symmetryBlock.push_back(basis2[a]);
		}

		if (symmetryBlock.size() == 0) return;

		iperm.resize(symmetryBlock.size());
		sort.sort(symmetryBlock,iperm);

		for (SizeType k = 0; k < symmetryBlock.size(); ++k)
			basis[k + offset] = symmetryBlock[k];

		offset += symmetryBlock.size();
		symmetryBlock.clear();
	}

	ModelCommonBaseType* modelCommon_;

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
