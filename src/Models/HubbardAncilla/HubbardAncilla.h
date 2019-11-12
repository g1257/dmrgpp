/*
Copyright (c) 2009-2015-2019, UT-Battelle, LLC
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

/*! \file HubbardAncilla.h
 *
 *  An implementation of a FeBasedSc model for Fe-based superconductors to use with
 *  the DmrgSolver
 *
 */
#ifndef DMRG_HUBBARD_ANCILLA_H
#define DMRG_HUBBARD_ANCILLA_H
#include "ModelBase.h"
#include "ParametersHubbardAncilla.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include "Geometry/GeometryDca.h"
#include <cstdlib>
#include "HelperHubbardAncilla.h"

namespace Dmrg {
template<typename ModelBaseType>
class HubbardAncilla : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
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
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ParametersHubbardAncilla<RealType, QnType> ParametersHubbardAncillaType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef HelperHubbardAncilla<ModelBaseType, ParametersHubbardAncillaType>
	HelperHubbardAncillaType;
	typedef typename HelperHubbardAncillaType::HilbertSpaceFeAsType HilbertSpaceFeAsType;

	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;
	static SizeType const ORBITALS = 2;

	HubbardAncilla(const SolverParamsType& solverParams,
	               InputValidatorType& io,
	               GeometryType const &geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      geometry_(geometry),
	      helperHubbardAncilla_(geometry_, modelParameters_)
	{}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		helperHubbardAncilla_.write(label1, io, this->params().model);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType t) const
	{
		helperHubbardAncilla_.addDiagonalsInNaturalBasis(hmatrix, block, t);
	}

	virtual PsimagLite::String oracle() const
	{
		const RealType ne = ModelBaseType::targetQuantum().qn(0).other[0];
		const RealType nup = ModelBaseType::targetQuantum().qn(0).other[1];
		const RealType ndown = ne - nup;
		const RealType n = ModelBaseType::geometry().numberOfSites();
		RealType energy = -nup*(n - nup) - ndown*(n - ndown);
		return ModelBaseType::oracle(energy, " -Nup*(L-Nup) -Ndown*(L-Ndown)");
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0; // FIXME for Immm SDHS
		BlockType block(1, site);
		HilbertBasisType natBasis;
		HelperHubbardAncillaType::setBasis(natBasis, block);
		HelperHubbardAncillaType::setSymmetryRelated(qns, natBasis);

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& ll = this->createOpsLabel("l");

		this->makeTrackable("c");
		this->makeTrackable("l");

		const bool hot = helperHubbardAncilla_.isHot();

		SizeType dofs = 2*ORBITALS;
		for (SizeType i=0;i<block.size();i++) {
			VectorSparseMatrixType vm;
			HelperHubbardAncillaType::findAllMatrices(vm,i,natBasis);
			for (SizeType sigma = 0; sigma < dofs; ++sigma) {
				if (!hot && (sigma & 1)) continue;
				MatrixType tmp;
				HelperHubbardAncillaType::findOperatorMatrices(tmp,i,sigma,natBasis);
				SparseMatrixType tmpMatrix(tmp);
				SizeType m=0;
				int asign=1;
				if (sigma>ORBITALS-1) {
					m=1;
					asign= -1;
				}

				typename OperatorType::Su2RelatedType su2related;
				if (sigma <ORBITALS) {
					su2related.source.push_back(i*dofs+sigma);
					su2related.source.push_back(i*dofs+sigma + ORBITALS);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = ORBITALS;
				}

				OperatorType myOp(tmpMatrix,
				                  ProgramGlobals::FermionOrBosonEnum::FERMION,
				                  typename OperatorType::PairType(1,m),
				                  asign,
				                  su2related);
				c.push(myOp);
			}

			HelperHubbardAncillaType::setLambdaMatrices(ll, vm);
		}
	}

	void fillModelLinks()
	{
		ModelTermType& hop = ModelBaseType::createTerm("hopping");
		ModelTermType& ll = ModelBaseType::createTerm("LambdaLambda");

		const SizeType orbitals = (helperHubbardAncilla_.isHot()) ? 2 : 1;

		for (SizeType spin = 0; spin < 2; ++spin) {
			for (SizeType orb = 0; orb < orbitals; ++orb) {
				OpForLinkType c("c", orb + spin*orbitals, orb);

				hop.push(c, 'N', c, 'C', 1, (spin == 1) ? -1 : 1, spin);
			}

			OpForLinkType l("l", spin);

			ll.push(l, 'N', l, 'C');
		}
	}

private:

	//serializr normal modelParameters_
	ParametersHubbardAncillaType  modelParameters_;
	//serializr ref geometry_ start
	const GeometryType& geometry_;
	HelperHubbardAncillaType helperHubbardAncilla_;
}; //class HubbardAncilla
} // namespace Dmrg
/*@}*/
#endif

