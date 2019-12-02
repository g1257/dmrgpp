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

/*! \file HolsteinThin.h
 *
 *  An implementation of a Hubbard Holstein model to use with the DmrgSolver
 *
 */
#ifndef DMRG_HOLSTEIN_THIN_H
#define DMRG_HOLSTEIN_THIN_H
#include "ModelBase.h"
#include "../HubbardHolstein/ParametersHubbardHolstein.h"
#include "../HubbardOneBand/HilbertSpaceHubbard.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include "Geometry/GeometryDca.h"
#include <cstdlib>
#include <numeric>

namespace Dmrg {
template<typename ModelBaseType>
class HolsteinThin : public ModelBaseType {

public:

	enum class SiteType {SITE_BOSON, SITE_FERMION};

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
	typedef unsigned int long WordType;
	typedef HilbertSpaceHubbard<WordType> HilbertSpaceHubbardType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ParametersHubbardHolstein<RealType, QnType> ParametersHolsteinThinType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef typename ModelBaseType::ModelLinksType ModelLinksType;
	typedef typename ModelLinksType::AtomKindBase AtomKindBaseType;

	class AtomKind : public AtomKindBaseType {

	public:

		virtual SizeType siteToAtomKind(SizeType site) const
		{
			return (site & 1);
		}

		virtual SizeType kindsOfAtoms() const { return 2; }
	};

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP = HilbertSpaceHubbardType::SPIN_UP;
	static const int SPIN_DOWN = HilbertSpaceHubbardType::SPIN_DOWN;

	HolsteinThin(const SolverParamsType& solverParams,
	             InputValidatorType& io,
	             const GeometryType& geometry,
	             PsimagLite::String additional)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      geometry_(geometry),
	      isSsh_(additional == "SSH"),
	      atomKind_(0)
	{
		if (isSsh_)
			err("SSH not supported in thin version yet!\n");
	}

	~HolsteinThin()
	{
		delete atomKind_;
		atomKind_ = nullptr;
	}

	void print(std::ostream& os) const { operator<<(os,modelParameters_); }

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType) const
	{
		assert(block.size() == 1);
		const SizeType actualSite = block[0];
		SiteType kindOfSite = determineKindOfSiteFromSite(actualSite);

		if (kindOfSite == SiteType::SITE_FERMION) {
			HilbertBasisType fermionicBasis;
			setBasis(fermionicBasis, SiteType::SITE_FERMION);
			VectorSparseMatrixType cm;
			findAllFermionicMatrices(cm, fermionicBasis);
			addInteractionFU(hmatrix, cm, actualSite);
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
		MatrixType nmatrix;
		for (SizeType sigma = 0; sigma < 2; ++sigma) {
			tmpMatrix = findOperatorMatrices(sigma, fermionicBasis);
			int asign = 1;
			if (sigma > 0) asign= 1;
			typename OperatorType::Su2RelatedType su2related;
			if (sigma==0) {
				su2related.source.push_back(0);
				su2related.source.push_back(1);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(-1);
				su2related.offset = 1;
			}

			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::FERMION,
			                  typename OperatorType::PairType(1 ,1 - sigma),
			                  asign,
			                  su2related);

			c.push(myOp);

			if (sigma == 0)
				nmatrix = multiplyTc(tmpMatrix, tmpMatrix);
			else
				nmatrix += multiplyTc(tmpMatrix, tmpMatrix);
		}


		OpsLabelType& n = this->createOpsLabel("n", 1); // 1 == fermionic site
		SparseMatrixType nmatrixCrs;
		fullMatrixToCrsMatrix(nmatrixCrs, nmatrix);
		typename OperatorType::Su2RelatedType su2related;
		OperatorType myOp(nmatrixCrs,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  typename OperatorType::PairType(0, 0),
		                  1,
		                  su2related);

		n.push(myOp);

		if (modelParameters_.numberphonons == 0) return;

		this->makeTrackable("n");

		tmpMatrix = findPhononadaggerMatrix(bosonicBasis);

		typename OperatorType::Su2RelatedType su2related2;
		su2related2.source.push_back(0*2);
		su2related2.source.push_back(0*2+1);
		su2related2.source.push_back(0*2);
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
	}

	void fillModelLinks()
	{
		ModelTermType& hopf = ModelBaseType::createTerm("HoppingFermionic");

		OpForLinkType cup("c", 0);
		hopf.push(cup, 'C', cup, 'N', 1, 1, 0);

		OpForLinkType cdown("c", 1);
		hopf.push(cdown, 'C', cdown, 'N', 1, -1, 1);

		if (modelParameters_.numberphonons > 0) {
			ModelTermType& hopb = ModelBaseType::createTerm("HoppingBosonic");

			OpForLinkType a("a");
			hopb.push(a, 'C', a, 'N', 1, 1, 0);

			ModelTermType& phononFermion = ModelBaseType::createTerm("PhononFermion");

			OpForLinkType n("n");
			phononFermion.push(n, 'N', a, 'N', 1, 1, 0);

			ModelTermType& fermionPhonon = ModelBaseType::createTerm("FermionPhonon");

			fermionPhonon.push(a, 'N', n, 'N', 1, 1, 0);
		}
	}

private:

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setBasis(HilbertBasisType& basis, SiteType kindOfSite) const
	{
		const HilbertState total = (kindOfSite == SiteType::SITE_BOSON) ?
		            (modelParameters_.numberphonons + 1) : 4;

		basis.resize(total);

		for (SizeType i = 0; i < total; ++i)
			basis[i] = i;
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceHubbardType::HilbertState const &ket,
	              int i,
	              int sigma) const
	{
		int value=0;
		value += HilbertSpaceHubbardType::calcNofElectrons(ket,0,i,0);
		value += HilbertSpaceHubbardType::calcNofElectrons(ket,0,i,1);
		int tmp1 = HilbertSpaceHubbardType::get(ket,0) &1;
		int tmp2 = HilbertSpaceHubbardType::get(ket,0) &2;
		if (i>0 && tmp1>0) value++;
		if (i>0 && tmp2>0) value++;

		if (sigma==1) { // spin down
			if ((HilbertSpaceHubbardType::get(ket,i) &1)) value++;

		}
		if (value%2==0) return 1.0;

		return FERMION_SIGN;
	}

	void findAllFermionicMatrices(VectorSparseMatrixType& cm,
	                              const HilbertBasisType& natBasis) const
	{
		cm.resize(2);
		cm[0] = findOperatorMatrices(SPIN_UP, natBasis);
		cm[1] = findOperatorMatrices(SPIN_DOWN, natBasis);
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(SizeType sigma,
	                                      const HilbertBasisType& natBasis) const
	{
		typename HilbertSpaceHubbardType::HilbertState bra,ket;
		int n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			if (HilbertSpaceHubbardType::isNonZero(ket, 0, sigma)) {

			} else {
				HilbertSpaceHubbardType::create(bra, 0, sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices\n");
				cm(ii,jj) = sign(ket, 0, sigma);
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
		typedef std::pair<SizeType,SizeType> PairType;

		qns.resize(basis.size() + offset, QnType::zero());
		static const bool isCanonical = true;
		VectorSizeType other((isCanonical) ? 2 : 1, 0);
		for (SizeType i = 0; i < basis.size(); ++i) {

			PairType jmpair(0,0);

			bool sign = false;

			if (typeOfSite == SiteType::SITE_FERMION) {
				// nup
				SizeType electronsUp = HilbertSpaceHubbardType::getNofDigits(basis[i], 0);
				// ndown
				SizeType electronsDown = HilbertSpaceHubbardType::getNofDigits(basis[i], 1);

				other[0] = electronsUp + electronsDown;

				if (isCanonical) other[1] = electronsUp;
				sign = other[0] & 1;
			}

			qns[i + offset] = QnType(sign, other, jmpair, other[0]);
		}
	}

	//! Term is U \sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteractionFU(SparseMatrixType& hmatrix,
	                      const VectorSparseMatrixType& cm,
	                      SizeType actualSite) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType m1=cm[SPIN_UP];
		SparseMatrixType m2=cm[SPIN_DOWN];

		multiply(tmpMatrix,n(m1),n(m2));
		assert(actualSite < modelParameters_.hubbardFU.size());
		hmatrix += modelParameters_.hubbardFU[actualSite]*tmpMatrix;
	}

	void addPotentialFV(SparseMatrixType &hmatrix,
	                    const VectorSparseMatrixType& cm,
	                    SizeType actualIndexOfSite) const
	{
		SparseMatrixType nup = n(cm[SPIN_UP]);
		SparseMatrixType ndown = n(cm[SPIN_DOWN]);

		SizeType linSize = geometry_.numberOfSites();
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialFV.size());
		hmatrix += modelParameters_.potentialFV[iUp] * nup;
		SizeType iDown = actualIndexOfSite + linSize;
		assert(iDown < modelParameters_.potentialFV.size());
		hmatrix += modelParameters_.potentialFV[iDown] * ndown;
	}

	void addPotentialPhononV(SparseMatrixType& hmatrix,
	                         const SparseMatrixType& amatrix,
	                         SizeType actualIndexOfSite) const
	{
		if (modelParameters_.numberphonons == 0) return;
		SparseMatrixType nphon = n(amatrix);
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialPV.size());
		hmatrix += modelParameters_.potentialPV[iUp] * nphon;
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

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
		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,cm);
		transposeConjugate(operatorMatrix, temp);
		return operatorMatrix;
	}

	static SiteType determineKindOfSiteFromSite(SizeType site)
	{
		return (site & 1) ? SiteType::SITE_FERMION : SiteType::SITE_BOSON;
	}

	ParametersHolsteinThinType modelParameters_;
	const GeometryType& geometry_;
	bool isSsh_;
	const AtomKind* atomKind_;
}; //class HolsteinThin
} // namespace Dmrg
/*@}*/
#endif

