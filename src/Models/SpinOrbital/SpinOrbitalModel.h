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
#include "Sort.h"

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
		HilbertBasisType natBasis;
		setBasis(natBasis);

		MatrixType splus = findSplusMatrices(natBasis, 0);
		MatrixType lplus = findSplusMatrices(natBasis, 1);
		MatrixType lminus;
		transposeConjugate(lminus, lplus);
		MatrixType tmp = splus*lminus;
		MatrixType tmp2;
		transposeConjugate(tmp2, tmp);
		MatrixType tmp3 = tmp;
		tmp3 += tmp2;
		tmp3 *= 0.5;
		MatrixType sz = findSzMatrices(natBasis, 0);
		MatrixType lz = findSzMatrices(natBasis, 1);
		tmp = sz*lz;
		sDotL_ = tmp3 + tmp;
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
		assert(block.size() == 1);
//		SizeType site = block[0]; // lambda1 and lambda2 have no site depedence
		MatrixType tmp = sDotL_;
		tmp *= modelParams_.lambda1;

		MatrixType sdotlSquared = modelParams_.lambda2*sDotL_*sDotL_;

		tmp += sdotlSquared;

		SparseMatrixType tmp3(tmp);
		hmatrix += tmp3;
	}

protected:

	// Trackable Operators are
	// 0   S+
	// 1   Sz
	// 2   L+
	// 3   Lz
	// 4   L+2
	// 5   L+L-
	// 6   L-L+
	// 7   L+Lz
	// 8   Lz2
	// 9   S+L+
	// 10  S+L-
	// 11  S+Lz
	// 12  L+Sz
	// 13  SzLz
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
		setBasis(natBasis);

		setSymmetryRelated(qns, natBasis);

		// this creates trackables 0 to 3
		for (SizeType orbital = 0; orbital < 2; ++orbital) {
			PsimagLite::String sOrL = (orbital == 0) ? "s" : "l";
			// Set the operators S^+_i in the natural basis
			MatrixType tmpMatrix = findSplusMatrices(natBasis, orbital);
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(tmpMatrix),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel(sOrL + "plus").push(myOp);
			this->makeTrackable(sOrL + "plus");

			myOp.dagger();
			this->createOpsLabel(sOrL + "minus").push(myOp);

			// Set the operators S^z_i in the natural basis
			tmpMatrix = findSzMatrices(natBasis, orbital);
			typename OperatorType::Su2RelatedType su2related2;
			OperatorType myOp2(SparseMatrixType(tmpMatrix),
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1,
			                   su2related2);
			this->createOpsLabel(sOrL + "z").push(myOp2);
			this->makeTrackable(sOrL + "z");
		}

		MatrixType lplus = findSplusMatrices(natBasis, 1); //lplus
		MatrixType lminus;
		transposeConjugate(lminus, lplus);
		MatrixType lpluslminus = lplus*lminus;
		MatrixType lminuslplus = lminus*lplus;
		MatrixType lplusSquared = lplus*lplus;
		MatrixType lz = findSzMatrices(natBasis, 1); // lz
		MatrixType lplusLz = lplus*lz;
		MatrixType lzSquared = lz*lz;

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(lplusSquared),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("lplusSquared").push(myOp);
			this->makeTrackable("lplusSquared");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(lpluslminus),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("lpluslminus").push(myOp);
			this->makeTrackable("lpluslminus");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(lminuslplus),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("lminuslplus").push(myOp);
			this->makeTrackable("lminuslplus");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(lplusLz),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("lplusLz").push(myOp);
			this->makeTrackable("lplusLz");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(lzSquared),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("lzSquared").push(myOp);
			this->makeTrackable("lzSquared");
		}

		MatrixType splus = findSplusMatrices(natBasis, 0); //splus
		MatrixType sz = findSzMatrices(natBasis, 0); // sz
		MatrixType spluslplus = splus*lplus;
		MatrixType spluslminus = splus*lminus;
		MatrixType splusLz = splus*lz;
		MatrixType lplusSz = lplus*sz;
		MatrixType szlz = sz*lz;

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(spluslplus),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("spluslplus").push(myOp);
			this->makeTrackable("spluslplus");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(spluslminus),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("spluslminus").push(myOp);
			this->makeTrackable("spluslminus");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(splusLz),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("splusLz").push(myOp);
			this->makeTrackable("splusLz");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(lplusSz),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("lplusSz").push(myOp);
			this->makeTrackable("lplusSz");
		}

		{
			typename OperatorType::Su2RelatedType su2related;
			OperatorType myOp(SparseMatrixType(szlz),
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1,
			                  su2related);
			this->createOpsLabel("szlz").push(myOp);
			this->makeTrackable("szlz");
		}
	}

	// Connectors and Factors (if factor is missing, 1 is assumed)
	// HC means Hermitian conjugate
	// Class is a letter starting with a
	//
	// a        S+ S+      N C    .5
	// a        Sz Sz      N N
	//
	// b        L+ L+      N C    .5
	// b        Lz Lz      N N
	//
	// c        L+2 L+2    N C    .25
	// c        L+L- L+L-  N N    .25
	// c        L-L+ L-L+  N N    .25
	// c        L+Lz L+Lz  N C
	// c        L+Lz L+    N C    0.5
	// c (HC)   of the above      0.5
	// c        L+   L+Lz  N C    0.5
	// c (HC)   of the above      0.5
	// c        Lz2 Lz2    N N
	//
	// d        S+L+ S+L+  N C      .25
	// d  (HC) of the above
	// d        S+L- S+L-  N C      .25
	// d  (HC) of the above
	// d        S+Lz S+Lz  N C      .5
	// d  (HC) of the above
	// d        L+Sz L+Sz  N C      .5
	// d  (HC) of the above
	// d        SzLz SzLz  N N
	//
	// more will be needed for J3L term
	void fillModelLinks()
	{
		auto valueModiferTerm0 = [](ComplexOrRealType& value) { value *=  0.5;};
		auto valueModiferTerm1 = [](ComplexOrRealType& value) { value *=  0.25;};

		// this creates connections a and b
		for (SizeType orbital = 0; orbital < 2; ++orbital) {
			PsimagLite::String sOrL = (orbital == 0) ? "s" : "l";
			ModelTermType& sdotS = ModelBaseType::createTerm(sOrL + "Dot" + sOrL);

			OpForLinkType splus(sOrL + "plus");

			sdotS.push(splus, 'N', splus, 'C', valueModiferTerm0);

			OpForLinkType sz(sOrL + "z");
			sdotS.push(sz, 'N', sz, 'N');
		}

		// this creates connections in c listed above
		ModelTermType& ldotLSquared = ModelBaseType::createTerm("ldotLSquared", false);

		OpForLinkType lplusSquared("lplusSquared");
		ldotLSquared.push(lplusSquared, 'N', lplusSquared, 'C', valueModiferTerm1);
		ldotLSquared.push(lplusSquared, 'C', lplusSquared, 'N', valueModiferTerm1);

		OpForLinkType lpluslminus("lpluslminus");
		OpForLinkType lminuslplus("lminuslplus");
		ldotLSquared.push(lpluslminus, 'N', lminuslplus, 'N', valueModiferTerm1);
		ldotLSquared.push(lminuslplus, 'N', lpluslminus, 'N', valueModiferTerm1);

		OpForLinkType lplusLz("lplusLz");

		ldotLSquared.push(lplusLz, 'N', lplusLz, 'C');
		ldotLSquared.push(lplusLz, 'C', lplusLz, 'N');

		OpForLinkType lplus("lplus");
		ldotLSquared.push(lplusLz, 'N', lplus, 'C', valueModiferTerm0);
		ldotLSquared.push(lplus, 'N', lplusLz, 'C', valueModiferTerm0);

		ldotLSquared.push(lplusLz, 'C', lplus, 'N', valueModiferTerm0);
		ldotLSquared.push(lplus, 'C', lplusLz, 'N', valueModiferTerm0);

		OpForLinkType lzSquared("lzSquared");
		ldotLSquared.push(lzSquared, 'N', lzSquared, 'N');

		// this creates connections in d listed above
		ModelTermType& ldotLsDotS = ModelBaseType::createTerm("ldotLsDotS");

		OpForLinkType spluslplus("spluslplus");
		ldotLsDotS.push(spluslplus, 'N', spluslplus, 'C', valueModiferTerm1);

		OpForLinkType spluslminus("spluslminus");
		ldotLsDotS.push(spluslminus, 'N', spluslminus, 'C', valueModiferTerm1);

		OpForLinkType splusLz("splusLz");
		ldotLsDotS.push(splusLz, 'N', splusLz, 'C', valueModiferTerm0);

		OpForLinkType lplusSz("lplusSz");
		ldotLsDotS.push(lplusSz, 'N', lplusSz, 'C', valueModiferTerm0);

		OpForLinkType szlz("szlz");
		ldotLsDotS.push(szlz, 'N', szlz, 'N');
	}

private:

	void setBasis(HilbertBasisType& natBasis) const
	{
		const SizeType total1 = modelParams_.twiceS + 1;
		const SizeType total2 = modelParams_.twiceL + 1;
		const SizeType total = total1*total2;
		natBasis.resize(total);
		for (SizeType i = 0; i < total; ++i) natBasis[i] = i;

		orderOneSiteBasis(natBasis);
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

			if (mPlusJ(ket, orbital) >= twiceTheSpin)
				continue;

			SizeType mPlusj0 = mPlusJ(ket, 0);
			SizeType mPlusj1 = mPlusJ(ket, 1);
			SizeType bra = (orbital == 0) ? packM(mPlusj0 + 1, mPlusj1)
			                              : packM(mPlusj0, mPlusj1 + 1);

			RealType mPlusj = mPlusJ(ket, orbital);
			RealType m = mPlusj - j;
			RealType x = j*(j + 1) - m*(m + 1);
			assert(x >= 0);

			cm(bra, ket) = sqrt(x);
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

			RealType mPlusj = mPlusJ(ket, orbital);
			RealType m = mPlusj - j;
			cm(ket ,ket) = m;
		}

		return cm;
	}

	// ket = sz' + lz'*(2s + 1)
	SizeType mPlusJ(SizeType ket, SizeType orbital) const
	{
		div_t q = div(ket, modelParams_.twiceS + 1);
		assert(static_cast<SizeType>(q.rem) <= modelParams_.twiceS);
		assert(static_cast<SizeType>(q.quot) <= modelParams_.twiceL);
		return (orbital == 0) ? q.rem : q.quot;
	}

	SizeType packM(SizeType szp, SizeType lzp) const
	{
		assert(szp <= modelParams_.twiceS);
		assert(lzp <= modelParams_.twiceL);
		return szp + lzp*(modelParams_.twiceS + 1);
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

		bool isCanonical = (ModelBaseType::targetQuantum().sizeOfOther() == 2);

		if (!isCanonical)
			err(PsimagLite::String(__FILE__) + ": must have exactly one symmetry\n");

		VectorSizeType other(2);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0, 0);
			SizeType mOfSpinPlusJ = mPlusJ(basis[i], 0);
			SizeType mOfOrbitalPlusJ = mPlusJ(basis[i], 1);
			other[0] = mOfSpinPlusJ;
		        other[1] = mOfOrbitalPlusJ;
			SizeType flavor = 1;
			qns[i] = QnType(false, other, jmpair, flavor);
		}
	}

	void orderOneSiteBasis(HilbertBasisType& basis) const
	{
		const SizeType n = basis.size();
		VectorSizeType symm(n);
		for (SizeType i = 0; i < n; ++i) {
			SizeType mOfSpinPlusJ = mPlusJ(basis[i], 0);
			SizeType mOfOrbitalPlusJ = mPlusJ(basis[i], 1);
			symm[i] = mOfSpinPlusJ + mOfOrbitalPlusJ*(modelParams_.twiceS + 1);
		}

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(n);
		sort.sort(symm, iperm);

		HilbertBasisType basisSorted(n);
		for (SizeType i = 0; i < n; ++i) {
			basisSorted[i] = basis[iperm[i]];
		}

		basisSorted.swap(basis);
	}

	ParametersSpinOrbitalType modelParams_;
	const SuperGeometryType& superGeometry_;
	MatrixType sDotL_;
}; // class SpinOrbitalModel

} // namespace Dmrg
/*@}*/
#endif //DMRG_SPIN_ORBITAL_MODEL_H

