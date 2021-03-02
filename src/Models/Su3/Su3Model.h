/*
Copyright (c) 2009-2021, UT-Battelle, LLC
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

/*! \file Su3Model.h
 *
 *  An implementation of the Quantum Heisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_SU3_MODEL_H
#define DMRG_SU3_MODEL_H

#include <algorithm>
#include "ParametersSu3.h"
#include "CrsMatrix.h"
#include "../../Engine/VerySparseMatrix.h"
#include "../../Engine/ProgramGlobals.h"
#include "../../Engine/Utils.h"
#include "Su3RepresentationP1.h"

namespace Dmrg {

template<typename ModelBaseType>
class Su3Model : public ModelBaseType {

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::RealType RealType;
	typedef	typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename ModelBaseType::VectorQnType VectorQnType;
	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef unsigned int long WordType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef typename PsimagLite::Vector<SizeType>::Type HilbertBasisType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef ParametersSu3<RealType, QnType> ParametersSu3Type;
	typedef Su3RepresentationBase<ComplexOrRealType> Su3RepresentationBaseType;
	typedef Su3RepresentationP1<ComplexOrRealType> Su3RepresentationP1Type;

	Su3Model(const SolverParamsType& solverParams,
	                InputValidatorType& io,
	                const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      superGeometry_(geometry),
	      modelParameters_(io),
	      su3Rep_(nullptr)
	{
		if (modelParameters_.p == 1) {
			su3Rep_ = new Su3RepresentationP1Type();
		} else {
			err("Implementation for p = " + ttos(modelParameters_.p) +
			    " has not been implemented\n");
		}
	}

	~Su3Model()
	{
		delete su3Rep_;
		su3Rep_ = nullptr;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

	// m*T3(i)*T3(i) + m*T8(i)*T8(i)
	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType& block,
	                                RealType) const
	{
		assert(block.size() == 1);

		MatrixType m3;
		su3Rep_->getMatrix(m3, 2);

		MatrixType m8;
		su3Rep_->getMatrix(m8, 7);

		MatrixType m = m3*m3;
		m *= modelParameters_.mass;

		SparseMatrixType mSparse(m);

		hmatrix += mSparse;
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		SizeType total = su3Rep_->size();
		HilbertBasisType natBasis(total);
		for (SizeType i = 0; i < total; ++i) natBasis[i] = i;

		setSymmetryRelated(qns, natBasis, block.size());

		for (SizeType a = 0; a < 8; ++a) {

			MatrixType m;

			su3Rep_->getMatrix(m, a);

			SparseMatrixType sparseMatrix(m);

			typename OperatorType::Su2RelatedType su2related;

			OperatorType myOp(sparseMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("T" + ttos(a + 1)).push(myOp);
			this->makeTrackable("T" + ttos(a + 1));
		}
	}

	void fillModelLinks()
	{
		ModelTermType& jOne = ModelBaseType::createTerm("jOne");
		for (SizeType a = 0; a < 8; ++a) {

			OpForLinkType aOpForLink("T" + ttos(a + 1));

			jOne.push(aOpForLink, 'N', aOpForLink, 'C');
		}

		ModelTermType& jTwo = ModelBaseType::createTerm("jTwo");
		for (SizeType a = 0; a < 8; ++a) {

			OpForLinkType aOpForLink("T" + ttos(a + 1));

			jTwo.push(aOpForLink, 'N', aOpForLink, 'C');
		}
	}

private:

	// \sum_i T3(i) and \sum_i T4(i) are conserved separately
	// We delegate to the representation the actual values and computation
	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		typedef std::pair<SizeType,SizeType> PairType;

		VectorSizeType other(2);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0, 0);
			other[0] = su3Rep_->t3OfState(i);
			other[1] = su3Rep_->t8OfState(i);
			SizeType flavor = 1;
			qns[i] = QnType(false, other, jmpair, flavor);
		}
	}

	const SuperGeometryType& superGeometry_;
	ParametersSu3Type modelParameters_;
	Su3RepresentationBaseType* su3Rep_;
}; // class Su3Model

} // namespace Dmrg
/*@}*/
#endif //DMRG_SU3_MODEL_H

