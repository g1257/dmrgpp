/*
Copyright (c) 2009, 2017-2019, UT-Battelle, LLC
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

/*! \file HeisenbergMix.h
 *
 *  An implementation of the Quantum Heisenberg Model
 *  with spin 1 and with spin 1/2 on the border sites
 *  (to be used with DmrgSolver)
 *
 */

#ifndef DMRG_HEISENBERG_MIX_H
#define DMRG_HEISENBERG_MIX_H

#include <algorithm>
#include "ParametersModelHeisenberg.h"
#include "CrsMatrix.h"
#include "../../Engine/VerySparseMatrix.h"
#include "../../Engine/ProgramGlobals.h"
#include "../../Engine/Utils.h"
#include "Aklt.h"

namespace Dmrg {

template<typename ModelBaseType>
class HeisenbergMix : public ModelBaseType {

	static const int NUMBER_OF_ORBITALS=1;
	static const int DEGREES_OF_FREEDOM=2; // spin up and down

public:

	enum class SiteType {SITE_MIDDLE, SITE_BORDER};

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
	typedef Aklt<ModelBaseType> AkltType;
	typedef typename ModelLinksType::AtomKindBase AtomKindBaseType;
	using ModelParametersType = ParametersModelHeisenberg<RealType, QnType>;

	class AtomKind : public AtomKindBaseType {

	public:

		AtomKind(SizeType n) : n_(n) {}

		virtual SizeType siteToAtomKind(SizeType site) const
		{
			return (site == 0 || site == n_ - 1) ? 1 : 0;
		}

		virtual SizeType kindsOfAtoms() const { return 2; }

	private:

		SizeType n_;
	};

	HeisenbergMix(const SolverParamsType& solverParams,
	              InputValidatorType& io,
	              const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      superGeometry_(geometry),
	      atomKind_(nullptr)
	{
		SizeType n = superGeometry_.numberOfSites();
		SizeType md = modelParameters_.anisotropyD.size();
		SizeType me = modelParameters_.anisotropyE.size();

		ModelParametersType::checkMagneticField(modelParameters_.magneticFieldX.size(), 'X', n);

		ModelParametersType::checkMagneticField(modelParameters_.magneticFieldZ.size(), 'Z', n);

		if (md > 0 && md != n) {
			PsimagLite::String msg("HeisenbergMix: If provided, ");
			msg += " AnisotropyD must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (me > 0 && me != n) {
			PsimagLite::String msg("HeisenbergMix: If provided, ");
			msg += " AnisotropyE must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (BasisType::useSu2Symmetry()) {
			PsimagLite::String msg("HeisenbergMix: SU(2) symmetry ");
			msg += "is not implemented yet.\n";
			throw PsimagLite::RuntimeError(msg);
		}
	}

	~HeisenbergMix()
	{
		delete atomKind_;
		atomKind_ = nullptr;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

	/* PSIDOC HeisenbergMix::addDiagonalsInNaturalBasis
	 We describe only the addition of a Zeeman term to the Heisenberg model here; note
	 that this function is more complicated.
	 Please look carefully at the following C++ lines:
	 PSIDOCCOPY $FirstFunctionBelow::MagneticField
	 */
	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType time) const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		SizeType linSize = superGeometry_.numberOfSites();
		SizeType n = block.size();

		for (SizeType i = 0; i < n; ++i) {

			// PSIDOCMARK_BEGIN MagneticField
			SizeType site = block[i];

			const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);
			const OperatorType& splus = ModelBaseType::naturalOperator("splus", site, 0);

			addMagneticField(hmatrix, 'X', site);
			addMagneticField(hmatrix, 'Z', site);

			// PSIDOCMARK_END

			// anisotropyD
			if (modelParameters_.anisotropyD.size() == linSize) {
				SparseMatrixType Szsquare;
				RealType tmp = modelParameters_.anisotropyD[site];
				multiply(Szsquare, sz.getCRS(), sz.getCRS());
				hmatrix += tmp*Szsquare;

			}

			// anisotropyE
			if (modelParameters_.anisotropyE.size() == linSize) {

				SparseMatrixType splusSquared;

				RealType tmp = 0.5*modelParameters_.anisotropyE[site];
				multiply(splusSquared, splus.getCRS(), splus.getCRS());
				hmatrix += tmp*splusSquared;

				SparseMatrixType sminusSquared;
				transposeConjugate(sminusSquared, splusSquared);
				hmatrix += tmp*sminusSquared;
			}
		}
	}

protected:

	virtual const AtomKindBaseType& getAtomKind()
	{
		if (!atomKind_)
			atomKind_ = new AtomKind(superGeometry_.numberOfSites());
		return *atomKind_;
	}

	void fillLabeledOperators(VectorQnType& qns)
	{
		qns.clear();
		SizeType offset = fillLabeledOperators(qns, 0, SiteType::SITE_MIDDLE);
		fillLabeledOperators(qns, offset, SiteType::SITE_BORDER);
	}

	// this function can be private
	SizeType fillLabeledOperators(VectorQnType& qns, SizeType offset, SiteType typeOfSite)
	{
		SizeType site = (typeOfSite == SiteType::SITE_MIDDLE) ? 1 : 0;
		SizeType indOfKindOfSite = (typeOfSite == SiteType::SITE_MIDDLE) ? 0 : 1;
		HilbertBasisType natBasis;
		setBasis(natBasis, site);

		setSymmetryRelated(qns, natBasis, offset, typeOfSite);

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

	void fillModelLinks()
	{
		if (BasisType::useSu2Symmetry())
			err("SU(2) not supported\n");

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");

		OpForLinkType splus("splus");

		auto valueModiferTerm0 = [](ComplexOrRealType& value) { value *=  0.5;};

		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0);

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");

		OpForLinkType sz("sz");
		szsz.push(sz, 'N', sz, 'N');

		OpForLinkType splusB("splusB");

		const bool wantsHermit = true;
		ModelTermType& spsmB1 = ModelBaseType::createTerm("SplusSminusB1",
		                                                  wantsHermit,
		                                                  "SplusSminus");
		spsmB1.push(splusB, 'N', splus, 'C', valueModiferTerm0);

		ModelTermType& spsmB2 = ModelBaseType::createTerm("SplusSminusB2",
		                                                  wantsHermit,
		                                                  "SplusSminus");
		spsmB2.push(splus, 'N', splusB, 'C', valueModiferTerm0);

		OpForLinkType szB("szB");

		ModelTermType& szszB1 = ModelBaseType::createTerm("szszB1",
		                                                  wantsHermit,
		                                                  "szsz");
		szszB1.push(szB, 'N', sz, 'N');

		ModelTermType& szszB2 = ModelBaseType::createTerm("szszB2",
		                                                  wantsHermit,
		                                                  "szsz");
		szszB2.push(sz, 'N', szB, 'N');
	}

private:

	void addMagneticField(SparseMatrixType& hmatrix,
	                      char c,
	                      SizeType site) const
	{
		assert(c == 'X' || c == 'Z');

		const SizeType linSize = superGeometry_.numberOfSites();
		const VectorRealType& v = (c == 'X') ? modelParameters_.magneticFieldX
		                                     : modelParameters_.magneticFieldZ;

		if (v.size() != linSize) return;

		assert(site < v.size());
		RealType tmp = v[site];
		const OperatorType& sminus = ModelBaseType::naturalOperator("sminus", site, 0);

		const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);

		if (c == 'Z') {
			hmatrix += tmp*sz.getCRS();
			return;
		}

		assert(c == 'X');
		const OperatorType& splus = ModelBaseType::naturalOperator("splus", site, 0);
		constexpr RealType zeroPointFive = 0.5;
		hmatrix += zeroPointFive*tmp*splus.getCRS();
		hmatrix += zeroPointFive*tmp*sminus.getCRS();
	}

	// atomKind_ is not setup here yet, so do not call it
	SizeType getSpin(SizeType site) const
	{
		const SizeType n = superGeometry_.numberOfSites();
		return (site == 0 || site == n - 1) ? modelParameters_.twiceTheSpinBorder
		                                    : modelParameters_.twiceTheSpin;
	}

	void setBasis(HilbertBasisType& natBasis, SizeType site) const
	{
		const SizeType twiceTheSpin = getSpin(site);
		const SizeType total = twiceTheSpin + 1;
		natBasis.resize(total);
		for (SizeType i = 0; i < total; ++i) natBasis[i] = i;
	}

	//! Find S^+_site in the natural basis natBasis
	SparseMatrixType findSplusMatrices(SizeType site,
	                                   const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		const SizeType twiceTheSpin = getSpin(site);
		RealType j = 0.5*twiceTheSpin;
		SizeType bits = 1 + ProgramGlobals::logBase2(twiceTheSpin);
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;

			assert(ketsite == ket);
			SizeType brasite = ketsite + 1;
			if (brasite >= twiceTheSpin + 1) continue;

			SizeType bra = ket & (~mask);
			assert(bra == 0);
			bra |= brasite;
			assert(bra == brasite);

			RealType m = ketsite - j;
			RealType x = j*(j+1)-m*(m+1);
			assert(x>=0);

			cm(ket,bra) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(SizeType site,
	                                const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		const SizeType twiceTheSpin = getSpin(site);
		RealType j = 0.5*twiceTheSpin;
		SizeType bits = ProgramGlobals::logBase2(twiceTheSpin) + 1;
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			assert(ketsite == ket);
			RealType m = ketsite - j;
			cm(ket,ket) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find Maximal_i in the natural basis natBasis
	SparseMatrixType findMaximal(SizeType site,
	                             const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		const SizeType twiceTheSpin = getSpin(site);
		SizeType bits = ProgramGlobals::logBase2(twiceTheSpin) + 1;
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			assert(ketsite == ket);
			if (ketsite != twiceTheSpin)
				continue;
			cm(ket, ket) = 1;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}


	SparseMatrixType findSxMatrices(SizeType site,
	                                const HilbertBasisType& natBasis) const
	{
		SparseMatrixType Splus_temp = findSplusMatrices(site, natBasis);
		SparseMatrixType Sminus_temp,Sx;
		transposeConjugate(Sminus_temp,Splus_temp);
		RealType tmp=0.5;

		Sx = tmp*Splus_temp;
		Sx += tmp*Sminus_temp;

		return Sx;
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        SizeType offset,
	                        SiteType) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;

		bool isCanonical = (ModelBaseType::targetQuantum().sizeOfOther() == 1);

		if (isCanonical && !modelParameters_.magneticFieldX.empty())
			err(PsimagLite::String(__FILE__) +
			    ": MagneticFieldX CANNOT be canonical. Please " +
			    "delete the TargetSzPlusConst= from the input file\n");

		VectorSizeType other;
		if (isCanonical) other.resize(1, 0);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size() + offset, QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0, 0);
			if (isCanonical)
				other[0] = basis[i];
			SizeType flavor = 1;
			qns[i + offset] = QnType(false, other, jmpair, flavor);
		}
	}

	ModelParametersType modelParameters_;
	const SuperGeometryType& superGeometry_;
	const AtomKind* atomKind_;
}; // class HeisenbergMix

} // namespace Dmrg
/*@}*/
#endif //DMRG_HEISENBERG_MIX_H

