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

/*! \file HolsteinSpinlessThin.h
 *
 *  An implementation of a Hubbard Holstein model to use with the DmrgSolver
 *
 */
#ifndef DMRG_HOLSTEIN_SPINLESS_THIN_H
#define DMRG_HOLSTEIN_SPINLESS_THIN_H
#include "../HubbardHolsteinSpinless/ParametersHubbardHolsteinSpinless.h"
#include "../HubbardOneBand/HilbertSpaceHubbard.h"
#include "CrsMatrix.h"
#include "Geometry/GeometryDca.h"
#include "ModelBase.h"
#include "ProgramGlobals.h"
#include "SpinSquared.h"
#include "SpinSquaredHelper.h"
#include "VerySparseMatrix.h"
#include <cstdlib>
#include <numeric>

namespace Dmrg
{
template <typename ModelBaseType>
class HolsteinSpinlessThin : public ModelBaseType
{

public:

	enum class SiteType { SITE_BOSON,
		SITE_FERMION };

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
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
	typedef unsigned int long WordType;
	typedef HilbertSpaceHubbard<WordType> HilbertSpaceHubbardType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::MyBasis BasisType;
	typedef typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ParametersHubbardHolsteinSpinless<RealType, QnType> ParametersHolsteinSpinlessThinType;
	typedef std::pair<SizeType, SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef typename ModelBaseType::ModelLinksType ModelLinksType;
	typedef typename ModelLinksType::AtomKindBase AtomKindBaseType;

	class AtomKind : public AtomKindBaseType
	{

	public:

		virtual SizeType siteToAtomKind(SizeType site) const
		{
			return (site & 1);
		}

		virtual SizeType kindsOfAtoms() const { return 2; }
	};

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP = HilbertSpaceHubbardType::SPIN_UP;

	HolsteinSpinlessThin(const SolverParamsType& solverParams,
	    InputValidatorType& io,
	    const SuperGeometryType& geometry,
	    PsimagLite::String additional)
	    : ModelBaseType(solverParams,
		geometry,
		io)
	    , modelParameters_(io)
	    , isSsh_(additional == "SSH")
	    , atomKind_(0)
	{
		if (isSsh_)
			err("SSH not supported in thin version yet!\n");
	}

	~HolsteinSpinlessThin()
	{
		delete atomKind_;
		atomKind_ = nullptr;
	}

	void print(std::ostream& os) const { operator<<(os, modelParameters_); }

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	    const BlockType& block,
	    RealType time) const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		assert(block.size() == 1);
		const SizeType actualSite = block[0];
		SiteType kindOfSite = determineKindOfSiteFromSite(actualSite);

		if (kindOfSite == SiteType::SITE_FERMION) {
			HilbertBasisType fermionicBasis;
			setBasis(fermionicBasis, SiteType::SITE_FERMION);
			VectorSparseMatrixType cm;
			findAllFermionicMatrices(cm, fermionicBasis);
			addPotentialFV(hmatrix, cm, actualSite);
		} else {
			assert(kindOfSite == SiteType::SITE_BOSON);
			HilbertBasisType bosonicBasis;
			setBasis(bosonicBasis, SiteType::SITE_BOSON);
			SparseMatrixType tmpMatrix = findPhononadaggerMatrix(bosonicBasis);
			addPotentialPhononV(hmatrix, tmpMatrix, actualSite);
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

protected:

	virtual const AtomKindBaseType& getAtomKind()
	{
		if (!atomKind_)
			atomKind_ = new AtomKind();
		return *atomKind_;
	}

	void fillLabeledOperators(VectorQnType& qns)
	{
		OpsLabelType& c = this->createOpsLabel("c", 1); // 1 == fermionic site
		OpsLabelType& a = this->createOpsLabel("a", 0); // 0 == bosonic site
		this->makeTrackable("c");
		if (modelParameters_.numberphonons > 0)
			this->makeTrackable("a");

		// qns for bosons are all the same
		HilbertBasisType bosonicBasis;
		setBasis(bosonicBasis, SiteType::SITE_BOSON);
		qns.clear();
		setSymmetryRelated(qns, bosonicBasis, 0, SiteType::SITE_BOSON);

		HilbertBasisType fermionicBasis;
		setBasis(fermionicBasis, SiteType::SITE_FERMION);
		setSymmetryRelated(qns, fermionicBasis, bosonicBasis.size(), SiteType::SITE_FERMION);

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		SparseMatrixType tmpMatrix;
		SparseMatrixType nmatrix;
		for (SizeType sigma = 0; sigma < 1; ++sigma) {
			tmpMatrix = findOperatorMatrices(sigma, fermionicBasis);
			int asign = 1;
			if (sigma > 0)
				asign = 1;
			typename OperatorType::Su2RelatedType su2related;
			if (sigma == 0) {
				su2related.source.push_back(0);
				su2related.source.push_back(1);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(-1);
				su2related.offset = 1;
			}

			OperatorType myOp(tmpMatrix,
			    ProgramGlobals::FermionOrBosonEnum::FERMION,
			    typename OperatorType::PairType(1, 1 - sigma),
			    asign,
			    su2related);

			c.push(myOp);

			if (sigma == 0)
				nmatrix = n(tmpMatrix);
		}

		OpsLabelType& n = this->createOpsLabel("n", 1); // 1 == fermionic site
		typename OperatorType::Su2RelatedType su2related;
		OperatorType myOp(nmatrix,
		    ProgramGlobals::FermionOrBosonEnum::BOSON,
		    typename OperatorType::PairType(0, 0),
		    1,
		    su2related);

		n.push(myOp);

		if (modelParameters_.numberphonons == 0)
			return;

		OpsLabelType& disp = this->createOpsLabel("x", 0);
		this->makeTrackable("n");
		this->makeTrackable("x");

		tmpMatrix = findPhononadaggerMatrix(bosonicBasis);

		typename OperatorType::Su2RelatedType su2related2;
		su2related2.source.push_back(0 * 2);
		su2related2.source.push_back(0 * 2 + 1);
		su2related2.source.push_back(0 * 2);
		su2related2.transpose.push_back(-1);
		su2related2.transpose.push_back(-1);
		su2related2.transpose.push_back(1);
		su2related2.offset = 1;
		OperatorType myOp2(tmpMatrix,
		    ProgramGlobals::FermionOrBosonEnum::BOSON,
		    PairType(2, 2),
		    -1,
		    su2related2);
		a.push(myOp2);

		SparseMatrixType tmp2;
		transposeConjugate(tmp2, tmpMatrix);
		tmp2 += tmpMatrix;
		typename OperatorType::Su2RelatedType su2Related3;
		disp.push(OperatorType(tmp2,
		    ProgramGlobals::FermionOrBosonEnum::BOSON,
		    typename OperatorType::PairType(0, 0),
		    1.0,
		    su2Related3));
	}

	void fillModelLinks()
	{
		ModelTermType& hopf = ModelBaseType::createTerm("HoppingFermionic");

		OpForLinkType cup("c", 0);
		typename ModelTermType::Su2Properties su2properties(1, 1, 0);
		hopf.push(cup, 'C', cup, 'N', su2properties);

		if (modelParameters_.numberphonons > 0) {
			ModelTermType& hopb = ModelBaseType::createTerm("HoppingBosonic");

			OpForLinkType a("a");
			hopb.push(a, 'C', a, 'N');

			const bool wantsHerm = false;
			ModelTermType& phononFermion = ModelBaseType::createTerm("PhononFermion", wantsHerm);

			OpForLinkType n("n");
			OpForLinkType x("x");
			phononFermion.push(x, 'N', n, 'N');
		}
	}

private:

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setBasis(HilbertBasisType& basis, SiteType kindOfSite) const
	{
		const HilbertState total = (kindOfSite == SiteType::SITE_BOSON) ? (modelParameters_.numberphonons + 1) : 2;

		basis.resize(total);

		for (SizeType i = 0; i < total; ++i)
			basis[i] = i;
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceHubbardType::HilbertState const& ket,
	    int i,
	    int) const
	{
		int value = 0;
		value += HilbertSpaceHubbardType::calcNofElectrons(ket, 0, i, 0);
		int tmp1 = HilbertSpaceHubbardType::get(ket, 0) & 1;
		if (i > 0 && tmp1 > 0)
			value++;
		return (value % 2 == 0) ? 1.0 : FERMION_SIGN;
	}

	void findAllFermionicMatrices(VectorSparseMatrixType& cm,
	    const HilbertBasisType& natBasis) const
	{
		cm.resize(1);
		cm[0] = findOperatorMatrices(SPIN_UP, natBasis);
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(SizeType sigma,
	    const HilbertBasisType& natBasis) const
	{
		typename HilbertSpaceHubbardType::HilbertState bra, ket;
		int n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n, n);

		for (SizeType ii = 0; ii < natBasis.size(); ii++) {
			bra = ket = natBasis[ii];
			if (HilbertSpaceHubbardType::isNonZero(ket, 0, sigma)) {

			} else {
				HilbertSpaceHubbardType::create(bra, 0, sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis, bra);
				if (jj < 0)
					throw PsimagLite::RuntimeError("findOperatorMatrices\n");
				cm(ii, jj) = sign(ket, 0, sigma);
			}
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	// only fermions, bosons have no symmetry
	void setSymmetryRelated(VectorQnType& qns,
	    const HilbertBasisType& basis,
	    SizeType offset,
	    SiteType typeOfSite) const
	{
		typedef std::pair<SizeType, SizeType> PairType;

		qns.resize(basis.size() + offset, QnType::zero());
		VectorSizeType other(1);
		for (SizeType i = 0; i < basis.size(); ++i) {

			PairType jmpair(0, 0);

			bool sign = false;

			if (typeOfSite == SiteType::SITE_FERMION) {
				// nup
				SizeType electronsUp = HilbertSpaceHubbardType::getNofDigits(basis[i], 0);
				other[0] = electronsUp;

				sign = other[0] & 1;
			}

			qns[i + offset] = QnType(sign, other, jmpair, other[0]);
		}
	}

	void addPotentialFV(SparseMatrixType& hmatrix,
	    const VectorSparseMatrixType& cm,
	    SizeType actualIndexOfSite) const
	{
		SparseMatrixType nup = n(cm[SPIN_UP]);
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialFV.size());
		hmatrix += modelParameters_.potentialFV[iUp] * nup;
	}

	void addPotentialPhononV(SparseMatrixType& hmatrix,
	    const SparseMatrixType& amatrix,
	    SizeType actualIndexOfSite) const
	{
		if (modelParameters_.numberphonons == 0)
			return;
		SparseMatrixType nphon = n(amatrix);
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialPV.size());
		hmatrix += modelParameters_.potentialPV[iUp] * nphon;
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger, c);
		multiply(tmpMatrix, cdagger, c);

		return tmpMatrix;
	}

	//! Find a^+_site in the natural basis natBasis
	SparseMatrixType findPhononadaggerMatrix(const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total, total);
		for (SizeType ii = 0; ii < total; ++ii) {
			const WordType ket = natBasis[ii];
			const SizeType nphon = ket;
			if (nphon >= modelParameters_.numberphonons)
				continue; //  too many phonons, cannot create
			const WordType bra = ket + 1;
			const RealType x = bra;
			cm(ii, bra) = sqrt(x);
		}

		// print cm matrix
		//		std::cerr<<cm;

		SparseMatrixType operatorMatrix(cm);
		//		SparseMatrixType temp;
		//		fullMatrixToCrsMatrix(temp,cm);
		//		transposeConjugate(operatorMatrix, temp);
		return operatorMatrix;
	}

	static SiteType determineKindOfSiteFromSite(SizeType site)
	{
		return (site & 1) ? SiteType::SITE_FERMION : SiteType::SITE_BOSON;
	}

	ParametersHolsteinSpinlessThinType modelParameters_;
	bool isSsh_;
	const AtomKind* atomKind_;
}; // class HolsteinThinSpinless
} // namespace Dmrg
/*@}*/
#endif
