/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

/*! \file HubbardAncillaExtended.h
 *
 *  An implementation of a FeBasedSc model for Fe-based superconductors to use with
 *  the DmrgSolver
 *
 */
#ifndef DMRG_HUBBARD_ANCILLA_EXTENDED_H
#define DMRG_HUBBARD_ANCILLA_EXTENDED_H
#include "ModelBase.h"
#include "ParametersHubbardAncillaExtended.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include "Geometry/GeometryDca.h"
#include <cstdlib>
#include "../HubbardAncilla/HelperHubbardAncilla.h"

namespace Dmrg {
template<typename ModelBaseType>
class HubbardAncillaExtended : public ModelBaseType {

public:

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
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef ParametersHubbardAncillaExtended<RealType, QnType> ParametersHubbardAncillaType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef HelperHubbardAncilla<ModelBaseType, ParametersHubbardAncillaType>
	HelperHubbardAncillaType;
	typedef typename HelperHubbardAncillaType::HilbertSpaceFeAsType HilbertSpaceFeAsType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;
	static SizeType const ORBITALS = 2;

	HubbardAncillaExtended(const SolverParamsType& solverParams,
	                       InputValidatorType& io,
	                       GeometryType const &geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
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

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		HilbertBasisType natBasis;
		HelperHubbardAncillaType::setBasis(natBasis, block);
		helperHubbardAncilla_.setSymmetryRelated(qns, natBasis);

		//! Set the operators c^\dagger_{i\gamma\sigma} in the natural basis
		SizeType dofs = 2*ORBITALS;
		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& d = this->createOpsLabel("d");
		OpsLabelType& splus = this->createOpsLabel("splus");
		OpsLabelType& sminus = this->createOpsLabel("sminus");
		OpsLabelType& sz = this->createOpsLabel("sz");
		OpsLabelType& p = this->createOpsLabel("p");
		OpsLabelType& nop = this->createOpsLabel("n");

		this->makeTrackable("c");
		this->makeTrackable("d");
		this->makeTrackable("splus");
		this->makeTrackable("sz");
		this->makeTrackable("p");
		this->makeTrackable("n");

		const bool hot = helperHubbardAncilla_.isHot();
		for (SizeType i=0;i<block.size();i++) {
			VectorSparseMatrixType vm;
			HelperHubbardAncillaType::findAllMatrices(vm,i,natBasis);
			for (SizeType sigma=0;sigma<dofs;++sigma) {
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


			HelperHubbardAncillaType::setLambdaMatrices(d, vm);
			setSplus(splus, sminus, i, 0, vm);
			if (hot) setSplus(splus, sminus, i, 1, vm);
			setSz(sz, i, 0, vm);
			if (hot) setSz(sz, i, 1, vm);
			setPair(p, i, 0, vm);
			if (hot) setPair(p, i, 1, vm);
			setN(nop, i ,0 ,vm);
			if (hot) setN(nop, i, 1, vm);
		}

		//		if (hot_)
		//			assert(creationMatrix.size() == 14);
		//		else
		//			assert(creationMatrix.size() == 8);
	}

	void fillModelLinks()
	{
		const SizeType orbitals = (helperHubbardAncilla_.isHot()) ? 2 : 1;
		const bool isSu2 = BasisType::useSu2Symmetry();
		ModelTermType& hop = ModelBaseType::createTerm("hopping"); // diagonal in orbital
		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");
		ModelTermType& szsz = ModelBaseType::createTerm("szsz");
		ModelTermType& ninj = ModelBaseType::createTerm("ninj");
		ModelTermType& ll = ModelBaseType::createTerm("LambdaLambda");
		ModelTermType& pp = ModelBaseType::createTerm("PairPair");

		for (SizeType spin = 0; spin < 2; ++spin) {
			for (SizeType orb = 0; orb < orbitals; ++orb) {
				OpForLinkType c("c", orb + spin*orbitals, orb);

				hop.push(c, 'N', c, 'C', 1, (spin == 1) ? -1 : 1, spin);
			}

			OpForLinkType d("d", spin);

			ll.push(d, 'N', d, 'C');
		}

		auto valueModiferTerm0 = [isSu2](ComplexOrRealType& value)
		{ value *= (isSu2) ? -0.5 : 0.5;};
		auto valueModifierTermOther = [isSu2](ComplexOrRealType& value)
		{ if (isSu2) value = -value;};

		for (SizeType orb = 0; orb < orbitals; ++orb) {
			OpForLinkType splus("splus", orb, orb);
			OpForLinkType sz("sz", orb, orb);
			OpForLinkType n("n", orb, orb);
			OpForLinkType pair("p", orb, orb);

			spsm.push(splus, 'N', splus, 'C', 2, -1, 2, valueModiferTerm0);

			if (!isSu2)
				szsz.push(sz, 'N', sz, 'N', 2, 0.5);
			else
				spsm.push(splus, 'N', splus, 'C', 2, -1, 2, valueModifierTermOther);

			ninj.push(n, 'N', n, 'N');

			pp.push(pair, 'N', pair, 'C', 1, 1, 0,
			        [](ComplexOrRealType& value) {value *= (-1.0);});
		}
	}

private:

	void setSplus(OpsLabelType& splusop,
	              OpsLabelType& sminus,
	              SizeType,
	              SizeType orbital,
	              const VectorSparseMatrixType& vm) const
	{
		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType m;
		transposeConjugate(m,vm[2+orbital]);
		SparseMatrixType splus;
		multiply(splus,vm[0+orbital],m);

		OperatorType myOp(splus,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  typename OperatorType::PairType(0,0),
		                  -1,
		                  su2related);
		splusop.push(myOp);

		myOp.dagger();
		sminus.push(myOp);
	}

	void setSz(OpsLabelType& szop,
	           SizeType,
	           SizeType orbital,
	           const VectorSparseMatrixType& vm) const
	{
		typename OperatorType::Su2RelatedType su2related;

		SparseMatrixType cm1(vm[0+orbital]);
		SparseMatrixType cm2(vm[2+orbital]);
		SparseMatrixType n1 = HelperHubbardAncillaType::n(cm1);
		SparseMatrixType n2 = HelperHubbardAncillaType::n(cm2);
		MatrixType dn1;
		MatrixType dn2;
		crsMatrixToFullMatrix(dn1,n1);
		crsMatrixToFullMatrix(dn2,n2);

		SizeType n = dn1.rows();
		MatrixType szmatrix(n,n);

		for (SizeType i = 0; i < n; ++i)
			szmatrix(i,i) = static_cast<RealType>(0.5)*(dn1(i,i) - dn2(i,i));

		OperatorType sz(SparseMatrixType(szmatrix),
		                ProgramGlobals::FermionOrBosonEnum::BOSON,
		                typename OperatorType::PairType(0,0),
		                1,
		                su2related);

		szop.push(sz);
	}

	void setPair(OpsLabelType& p,
	             SizeType,
	             SizeType orbital,
	             const VectorSparseMatrixType& vm) const
	{
		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType pair;
		multiply(pair,vm[0+orbital],vm[2+orbital]);

		OperatorType myOp(pair,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  typename OperatorType::PairType(0,0),
		                  1,
		                  su2related);
		p.push(myOp);
	}

	void setN(OpsLabelType& nopop,
	          SizeType,
	          SizeType orbital,
	          const VectorSparseMatrixType& vm) const
	{
		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType cm1(vm[0+orbital]);
		SparseMatrixType cm2(vm[2+orbital]);
		SparseMatrixType n1 = HelperHubbardAncillaType::n(cm1);
		SparseMatrixType n2 = HelperHubbardAncillaType::n(cm2);
		MatrixType dn1;
		MatrixType dn2;
		crsMatrixToFullMatrix(dn1,n1);
		crsMatrixToFullMatrix(dn2,n2);
		SizeType n = dn1.rows();
		MatrixType nmatrix(n,n);

		for (SizeType i = 0; i < n; ++i) {
			nmatrix(i,i) = 1.0*(dn1(i,i) + dn2(i,i));
		}

		OperatorType nop(SparseMatrixType(nmatrix),
		                 ProgramGlobals::FermionOrBosonEnum::BOSON,
		                 typename OperatorType::PairType(0,0),
		                 1,
		                 su2related);

		nopop.push(nop);
	}

	ParametersHubbardAncillaType modelParameters_;
	const GeometryType& geometry_;
	HelperHubbardAncillaType helperHubbardAncilla_;
}; //class HubbardAncilla
} // namespace Dmrg
/*@}*/
#endif

