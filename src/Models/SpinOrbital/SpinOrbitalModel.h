/*
Copyright (c) 2009, 2017-2021, UT-Battelle, LLC
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

/*! \file SpinOrbitalModel.h
 *
 *  TBW
 *
 */

#ifndef DMRG_SPIN_ORBITAL_MODEL_H
#define DMRG_SPIN_ORBITAL_MODEL_H

#include <algorithm>
#include "CrsMatrix.h"
#include "../../Engine/VerySparseMatrix.h"
#include "../../Engine/ProgramGlobals.h"
#include "../../Engine/Utils.h"
#include "ParametersSpinOrbital.h"

namespace Dmrg {

template<typename ModelBaseType>
class SpinOrbitalModel : public ModelBaseType {

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
	typedef typename ModelBaseType::ModelLinksType ModelLinksType;
	typedef ParametersSpinOrbital<RealType, QnType> ParametersSpinOrbitalType;

	//typedef Aklt<ModelBaseType> AkltType;

	SpinOrbitalModel(const SolverParamsType& solverParams,
	              InputValidatorType& io,
	              const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParams_(io),
	      superGeometry_(geometry)
	{
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParams_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType& block,
	                                RealType) const
	{
		SizeType linSize = superGeometry_.numberOfSites();
		SizeType n = block.size();

		for (SizeType i = 0; i < n; ++i) {

			SizeType site = block[i];

			const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);
			const OperatorType& splus = ModelBaseType::naturalOperator("splus", site, 0);

			if (modelParams_.magneticFieldV.size() == linSize) {

				RealType tmp = modelParams_.magneticFieldV[site];
				const OperatorType& sminus = ModelBaseType::naturalOperator("sminus", site, 0);

				if (modelParams_.magneticFieldDirection == "z") {
					hmatrix += tmp*sz.getCRS();
				} else if (modelParams_.magneticFieldDirection == "x") {
					static const RealType zeroPointFive = 0.5;
					hmatrix += zeroPointFive*tmp*splus.getCRS();
					hmatrix += zeroPointFive*tmp*sminus.getCRS();
				}
			}
		}
	}

protected:

	// Trackable Operators are
	// 0   S+
	// 1   Sz
	// 2   L+
	// 3   Lz
	// 4   L+2
	// 5   L+Lz
	// 6   Lz2
	// 7   S+L+
	// 8   S+Lz
	// 9   SzL+
	// 10  SzLz
	// more will be needed for J3L term
	//
	//
	// Basis is (S, L)
	// 0                  -s     -l
	// 1                  -s     -l + 1
	// ...
	// (2s + 1)(2l + 1)    s      l
	//
	void fillLabeledOperators(VectorQnType& qns)
	{
		HilbertBasisType natBasis;
		setBasis(natBasis, site);

		setSymmetryRelated(qns, natBasis);

		// Set the operators S^+_i in the natural basis
		SparseMatrixType tmpMatrix = findSplusMatrices(site, natBasis);

		typename OperatorType::Su2RelatedType su2related;

		PsimagLite::String border = (typeOfSite == SiteType::SITE_BORDER) ? "B" : "";

		OperatorType myOp(tmpMatrix,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  PairType(2, 2),
		                  -1,
		                  su2related);
		this->createOpsLabel("splus" + border, indOfKindOfSite).push(myOp);
		this->makeTrackable("splus" + border);

		myOp.dagger();
		this->createOpsLabel("sminus" + border, indOfKindOfSite).push(myOp);

		// Set the operators S^z_i in the natural basis
		tmpMatrix = findSzMatrices(site, natBasis);
		typename OperatorType::Su2RelatedType su2related2;
		OperatorType myOp2(tmpMatrix,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(2, 1),
		                   1.0/sqrt(2.0),
		                   su2related2);
		this->createOpsLabel("sz" + border, indOfKindOfSite).push(myOp2);
		this->makeTrackable("sz" + border);

		// Set the operators S^x_i in the natural basis
		tmpMatrix = findSxMatrices(site, natBasis);
		typename OperatorType::Su2RelatedType su2related3;
		OperatorType myOp3(tmpMatrix,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(2, 1),
		                   1.0/sqrt(2.0),
		                   su2related3);
		this->createOpsLabel("sx" + border, indOfKindOfSite).push(myOp3);

		tmpMatrix = findMaximal(site, natBasis);
		OperatorType myOp4(tmpMatrix,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(2, 1),
		                   1.0/sqrt(2.0),
		                   su2related3);
		this->createOpsLabel("maximal" + border, indOfKindOfSite).push(myOp4);

		return natBasis.size();
	}

	// Connectors and Factors (if factor is missing, 1 is assumed)
	// HC means Hermitian conjugate
	// Class is a letter starting with a
	//
	// a        0 0      .5
	// a        1 1
	//
	// b        2 2      .5
	// b        3 3
	//
	// c        4 4      .5
	// c        5 5
	// c        2 5      .5
	// c (HC)   5 2      .5
	// c        6 6
	//
	// d        7 7      .25
	// d        8 8      .5
	// d        9 9      .5
	// d       10 10
	//
	// more will be needed for J3L term
	void fillModelLinks()
	{
		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");

		OpForLinkType splus("splus");

		auto valueModiferTerm0 = [](ComplexOrRealType& value) { value *=  0.5;};

		typename ModelTermType::Su2Properties su2properties(2, -1, 2);
		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0, su2properties);

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");

		OpForLinkType sz("sz");
		szsz.push(sz, 'N', sz, 'N', typename ModelTermType::Su2Properties(2, 0.5));

		OpForLinkType splusB("splusB");

		const bool wantsHermit = true;
		ModelTermType& spsmB1 = ModelBaseType::createTerm("SplusSminusB1",
		                                                  wantsHermit,
		                                                  "SplusSminus");
		spsmB1.push(splusB, 'N', splus, 'C', valueModiferTerm0, su2properties);

		ModelTermType& spsmB2 = ModelBaseType::createTerm("SplusSminusB2",
		                                                  wantsHermit,
		                                                  "SplusSminus");
		spsmB2.push(splus, 'N', splusB, 'C', valueModiferTerm0, su2properties);

		OpForLinkType szB("szB");

		ModelTermType& szszB1 = ModelBaseType::createTerm("szszB1",
		                                                  wantsHermit,
		                                                  "szsz");
		szszB1.push(szB, 'N', sz, 'N', typename ModelTermType::Su2Properties(2, 0.5));

		ModelTermType& szszB2 = ModelBaseType::createTerm("szszB2",
		                                                  wantsHermit,
		                                                  "szsz");
		szszB2.push(sz, 'N', szB, 'N', typename ModelTermType::Su2Properties(2, 0.5));
	}

private:

	void setBasis(HilbertBasisType& natBasis, SizeType site) const
	{
		const SizeType total1 = modelParams_.twiceS + 1;
		const SizeType total2 = modelParams_.twiceL + 1;
		const SizeType total = total1*total2;
		natBasis.resize(total);
		for (SizeType i = 0; i < total; ++i) natBasis[i] = i;
	}

	//! Find S^+_site in the natural basis natBasis
	MatrixType findSplusMatrices(const HilbertBasisType& natBasis, SizeType orbital) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total, total);
		const SizeType twiceTheSpin = (orbital == 0) ? modelParams_.twiceS : modelParams_.twiceL;
		RealType j = 0.5*twiceTheSpin;

		for (SizeType ii = 0; ii < total; ++ii) {
			SizeType ket = natBasis[ii];

			int bra = plusToKet(ket, orbital);

			if (bra < 0) continue;

			RealType mPlusJ = mPlusJ(ket, orbital);
			RealType m = mPlusJ - j;
			RealType x = j*(j + 1) - m*(m + 1);
			assert(x >= 0);

			cm(ket, bra) = sqrt(x);
		}

		return cm;
	}

	//! Find S^z_i in the natural basis natBasis
	MatrixType findSzMatrices(const HilbertBasisType& natBasis, SizeType orbital) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total, total);
		const SizeType twiceTheSpin = (orbital == 0) ? modelParams_.twiceS : modelParams_.twiceL;
		RealType j = 0.5*twiceTheSpin;

		for (SizeType ii = 0; ii < total; ++ii) {
			SizeType ket = natBasis[ii];

			RealType mPlusJ = mPlusJ(ket, orbital);
			RealType m = mPlusJ - j;
			cm(ket ,ket) = m;
		}

		return cm;
	}

	// Here the conserved quantity is (sz + s) + (lz + l)
	// Note that sz + s is always an integer, even if s is half integer
	// Note that lz + l is always an integer, even if l is half integer
	void setSymmetryRelated(VectorQnType& qns, const HilbertBasisType& basis) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;

		bool isCanonical = (ModelBaseType::targetQuantum().sizeOfOther() == 1);

		if (!isCanonical)
			err(PsimagLite::String(__FILE__) + ": must have exactly one symmetry\n");

		VectorSizeType other;
		other.resize(1, 0);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0, 0);
			SizeType mOfSpinPlusJ = mPlusJ(ket, 0);
			SizeType mOfOrbitalPlusJ = mPlusJ(ket, 1);
			other[0] = mOfSpinPlusJ + mOfOrbitalPlusJ;
			SizeType flavor = 1;
			qns[i + offset] = QnType(false, other, jmpair, flavor);
		}
	}

	ParametersSpinOrbitalType modelParams_;
	const SuperGeometryType& superGeometry_;
}; // class SpinOrbitalModel

} // namespace Dmrg
/*@}*/
#endif //DMRG_SPIN_ORBITAL_MODEL_H

