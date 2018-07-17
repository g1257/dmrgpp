/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file ExtendedHubbard1Orb.h
 *  FIXME: Merge into Hubbard
 *  Hubbard + V_{ij} n_i n_j
 *
 */
#ifndef EXTENDED_HUBBARD_1ORB_H
#define EXTENDED_HUBBARD_1ORB_H
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "LinkProdExtendedHubbard1Orb.h"

namespace Dmrg {
//! Extended Hubbard for DMRG solver, uses ModelHubbard by containment
template<typename ModelBaseType>
class ExtendedHubbard1Orb : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef ModelHubbard<ModelBaseType> ModelHubbardType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef LinkProdExtendedHubbard1Orb<ModelHelperType, GeometryType> LinkProductType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelHubbardType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
	typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename PsimagLite::Vector<HilbertState>::Type VectorHilbertStateType;

	ExtendedHubbard1Orb(const SolverParamsType& solverParams,
	                    InputValidatorType& io,
	                    GeometryType const &geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    new LinkProductType,
	                    io),
	      modelParameters_(io),
	      geometry_(geometry),
	      modelHubbard_(solverParams,io,geometry, 1)
	{}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const { return 0; }

	SizeType hilbertSize(SizeType site) const
	{
		return modelHubbard_.hilbertSize(site);
	}

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType&creationMatrix,
	                         VectorQnType& qns,
	                         const BlockType& block) const
	{
		modelHubbard_.setOperatorMatrices(creationMatrix, qns, block);
		// add ni to creationMatrix
		setNi(creationMatrix, block);
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		VectorQnType qns;
		setOperatorMatrices(creationMatrix,qns, block);

		if (what=="n") {
			VectorSizeType allowed(1,0);
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			PsimagLite::Matrix<ComplexOrRealType> tmp;
			return creationMatrix[2];
		} else {
			return modelHubbard_.naturalOperator(what,site,dof);
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
		        io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		modelHubbard_.write(label, io);
	}

	virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                        const VectorOperatorType& cm,
	                                        const BlockType& block,
	                                        RealType time,
	                                        RealType factorForDiagonals=1.0)  const
	{
		// remove the n matrices before sending to Hubbard
		VectorOperatorType cmCorrected;
		SizeType k = 0;
		for (SizeType i = 0; i < block.size(); ++i) {
			for (SizeType j = 0; j < 2; ++j) {
				cmCorrected.push_back(cm[k++]);
			}
			k++;
		}

		modelHubbard_.addDiagonalsInNaturalBasis(hmatrix,
		                                         cmCorrected,
		                                         block,
		                                         time,
		                                         factorForDiagonals);
	}

private:

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(int i,
	                                      const VectorHilbertStateType& natBasis) const
	{

		SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertState ket=natBasis[ii];
			cm(ii,ii) = 0.0;
			for (SizeType sigma=0;sigma<2;sigma++)
				if (HilbertSpaceHubbardType::isNonZero(ket,i,sigma))
					cm(ii,ii) += 1.0;
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	void setNi(VectorOperatorType& creationMatrix,
	           const BlockType& block) const
	{
		VectorOperatorType creationMatrix2 = creationMatrix;
		creationMatrix.clear();
		VectorHilbertStateType natBasis;
		modelHubbard_.setBasis(natBasis, block);
		SizeType operatorsPerSite = utils::exactDivision(creationMatrix2.size(),
		                                                 block.size());
		SizeType k = 0;

		for (SizeType i = 0; i < block.size(); ++i) {
			SparseMatrixType tmpMatrix = findOperatorMatrices(i,natBasis);
			RealType angularFactor= 1;
			typename OperatorType::Su2RelatedType su2related;
			su2related.offset = 1; //check FIXME
			OperatorType myOp(tmpMatrix,
			                  1,
			                  typename OperatorType::PairType(0,0),
			                  angularFactor,
			                  su2related);

			for (SizeType j = 0; j < operatorsPerSite; ++j)
				creationMatrix.push_back(creationMatrix2[k++]);

			creationMatrix.push_back(myOp);
		}
	}

	ParametersModelHubbard<RealType, QnType>  modelParameters_;
	const GeometryType &geometry_;
	ModelHubbardType modelHubbard_;
};	//class ExtendedHubbard1Orb

} // namespace Dmrg
/*@}*/
#endif // EXTENDED_HUBBARD_1ORB_H

