/*
Copyright (c) 2009-2018-2019, UT-Battelle, LLC
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

/*! \file Graphene.h
 *
 *  An implementation of a FeBasedSc model for Fe-based superconductors to use with
 *  the DmrgSolver
 *
 */
#ifndef DMRGPP_GRAPHENE_H
#define DMRGPP_GRAPHENE_H
#include "ModelBase.h"
#include "ParametersGraphene.h"
#include "../FeAsModel/HilbertSpaceFeAs.h"
//#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include <numeric>
#include "../../../PsimagLite/src/Vector.h"

namespace Dmrg {
template<typename ModelBaseType>
class Graphene : public ModelBaseType {

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
	typedef  HilbertSpaceFeAs<HilbertState> HilbertSpaceFeAsType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef	 typename ModelBaseType::MyBasis MyBasis;
	typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ParametersGraphene<ComplexOrRealType, QnType> ParametersGrapheneType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

	Graphene(const SolverParamsType& solverParams,
	         InputValidatorType& io,
	         GeometryType const &geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      geometry_(geometry)
	{
		ProgramGlobals::init(modelParameters_.orbitals*geometry_.numberOfSites() + 1);

		HilbertSpaceFeAsType::setOrbitals(modelParameters_.orbitals);
		statesPerSite_ = (1 << (modelParameters_.orbitals*2));

		VectorSizeType block(1,0);
		int sitesTimesDof = 2*modelParameters_.orbitals;
		HilbertState total = (1<<sitesTimesDof);
		basis_.resize(total);
		for (HilbertState a = 0; a < total; ++a) basis_[a] = a;

		setOperatorMatricesInternal(creationMatrix_, block);

		setSymmetryRelatedInternal(qq_,basis_,1);
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		io.write(label + "/statesPerSite_", statesPerSite_);
		io.write(label + "/basis_", basis_);
		io.write(label + "/qq_", qq_);
		io.write(label + "/creationMatrix_", creationMatrix_);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType& block,
	                                RealType time) const
	{
		SizeType n = block.size();
		assert(n == 1);
		for (SizeType i = 0; i < n; ++i) {
			addNiSquared(hmatrix);
			addPairHopping(hmatrix);
		}
	}

	void fillLabeledOperators(VectorQnType& qns)
	{
		qns = qq_;
		assert(creationMatrix_.size()>0);
		SizeType nrow = creationMatrix_[0].getCRS().rows();

		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& cc = this->createOpsLabel("C");
		for (SizeType dof = 0; dof < 2*modelParameters_.orbitals; ++dof) {
			VectorOperatorType cm = creationMatrix_;
			cc.push(creationMatrix_[dof]);
			cm[dof].dagger();
			c.push(cm[dof]);
		}

		this->makeTrackable("C");

		OpsLabelType& yop = this->createOpsLabel("yop");
		for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
			for (SizeType spin2 = 0; spin2 < 2; ++spin2) {
				MatrixType tmp =
				        multiplyTc(creationMatrix_[0+spin1*modelParameters_.orbitals].getCRS(),
				        creationMatrix_[1+spin2*modelParameters_.orbitals].getCRS());
				SparseMatrixType tmp2(tmp);
				typename OperatorType::Su2RelatedType su2Related;
				yop.push(OperatorType(tmp2,
				                      ProgramGlobals::FermionOrBosonEnum::BOSON,
				                      typename OperatorType::PairType(0,0),
				                      1.0,
				                      su2Related));
			}
		}

		this->makeTrackable("yop");

		OpsLabelType& nop = this->createOpsLabel("n");
		OpsLabelType& nn = this->createOpsLabel("nn");
		OpsLabelType& ntotal = this->createOpsLabel("ntotal");
		SparseMatrixType ntotalcrs;
		for (SizeType dof = 0; dof < 2*modelParameters_.orbitals; ++dof) {
			MatrixType tmp =
			        multiplyTc(creationMatrix_[dof].getCRS(),creationMatrix_[dof].getCRS());
			SparseMatrixType tmp2(tmp);
			ntotalcrs += tmp2;
			typename OperatorType::Su2RelatedType su2Related;
			nop.push(OperatorType(tmp2,
			                      ProgramGlobals::FermionOrBosonEnum::BOSON,
			                      typename OperatorType::PairType(0,0),
			                      1.0,
			                      su2Related));
			if (dof == 0 || dof == 2) // naup nadown
				nn.push(OperatorType(tmp2,
				                     ProgramGlobals::FermionOrBosonEnum::BOSON,
				                     typename OperatorType::PairType(0,0),
				                     1.0,
				                     su2Related));
		}

		this->makeTrackable("nn");


		typename OperatorType::Su2RelatedType su2Related;
		ntotal.push(OperatorType(ntotalcrs,
		                         ProgramGlobals::FermionOrBosonEnum::BOSON,
		                         typename OperatorType::PairType(0,0),
		                         1.0,
		                         su2Related));
		this->makeTrackable("ntotal");

		OpsLabelType& splus = this->createOpsLabel("splus");
		OpsLabelType& sminus = this->createOpsLabel("sminus");
		OpsLabelType& splus2 = this->createOpsLabel("splus2");
		for (SizeType dof = 0; dof < modelParameters_.orbitals; ++dof) {
			MatrixType tmp(nrow,nrow);
			tmp += multiplyTc(creationMatrix_[dof].getCRS(),
			                  creationMatrix_[dof + modelParameters_.orbitals].getCRS());
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			splus.push(OperatorType(tmp2,
			                        ProgramGlobals::FermionOrBosonEnum::BOSON,
			                        typename OperatorType::PairType(0,0),
			                        1.0,
			                        su2Related));
			SparseMatrixType tmp3;
			transposeConjugate(tmp3, tmp2);
			sminus.push(OperatorType(tmp3,
			                         ProgramGlobals::FermionOrBosonEnum::BOSON,
			                         typename OperatorType::PairType(0,0),
			                         1.0,
			                         su2Related));

			// s+ a
			if (dof == 0) splus2.push(OperatorType(tmp2,
			                                       ProgramGlobals::FermionOrBosonEnum::BOSON,
			                                       typename OperatorType::PairType(0,0),
			                                       1.0,
			                                       su2Related));
		}

		this->makeTrackable("splus2");


		OpsLabelType& sz = this->createOpsLabel("sz");
		for (SizeType dof = 0; dof < modelParameters_.orbitals; ++dof) {
			MatrixType tmp(nrow,nrow);
			MatrixType tmp2(nrow,nrow);

			tmp += multiplyTc(creationMatrix_[dof].getCRS(),creationMatrix_[dof].getCRS());
			tmp2 += multiplyTc(creationMatrix_[dof+modelParameters_.orbitals].getCRS(),
			        creationMatrix_[dof+modelParameters_.orbitals].getCRS());

			tmp = 0.5*(tmp-tmp2);
			SparseMatrixType tmp3(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			sz.push(OperatorType(tmp3,
			                     ProgramGlobals::FermionOrBosonEnum::BOSON,
			                     typename OperatorType::PairType(0,0),
			                     1.0,
			                     su2Related));
		}
	}

	/* PSIDOC FeAs::fillModelLinks
	  Write DOC here TBW FIXME TODO
	 */
	void fillModelLinks()
	{
		const SizeType orbitals = modelParameters_.orbitals;
		ModelTermType& hopA = ModelBaseType::createTerm("hoppingA");

		for (SizeType spin = 0; spin < 2; ++spin) {
			OpForLinkType c1("C", 0 + spin*orbitals);
			OpForLinkType c2("C", 0 + spin*orbitals);
			hopA.push(c1, 'N', c2, 'C', 1, (spin == 1) ? -1 : 1, spin);
		}

		ModelTermType& hopB = ModelBaseType::createTerm("hoppingB");

		for (SizeType spin = 0; spin < 2; ++spin) {
			OpForLinkType c1("C", 1 + spin*orbitals);
			OpForLinkType c2("C", 1 + spin*orbitals);
			hopB.push(c1, 'N', c2, 'C', 1, (spin == 1) ? -1 : 1, spin);
		}

		ModelTermType& ninj = ModelBaseType::createTerm("ninj");
		OpForLinkType ni("ntotal");
		ninj.push(ni, 'N', ni, 'N');

		ModelTermType& nana = ModelBaseType::createTerm("nana");
		for (SizeType spin = 0; spin < 2; ++spin) {
			OpForLinkType ni("nn", spin);
			nana.push(ni, 'N', ni, 'N');
		}

		ModelTermType& spasma = ModelBaseType::createTerm("spasma");
		OpForLinkType sp("splus2");
		spasma.push(sp, 'N', sp, 'C');

		ModelTermType& yoyo = ModelBaseType::createTerm("yoyo");
		for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
			for (SizeType spin2 = 0; spin2 < 2; ++spin2) {
				OpForLinkType yi("yop", spin1 + spin2*2);
				OpForLinkType yj("yop", spin2 + spin1*2);
				yoyo.push(yi, 'N', yj, 'C');
			}
		}
	}

	void setQns(VectorQnType& qns) const
	{
		qns = qq_;
	}

	void setOperatorMatricesInternal(VectorOperatorType& creationMatrix,
	                                 const BlockType& block) const
	{
		const HilbertBasisType& natBasis = basis_;
		SparseMatrixType tmpMatrix;

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		creationMatrix.clear();
		SizeType dofs = 2*modelParameters_.orbitals;
		for (SizeType i=0;i<block.size();i++) {
			for (SizeType sigma=0;sigma<dofs;sigma++) {
				findOperatorMatrices(tmpMatrix,i,sigma,natBasis);

				SizeType m=0;
				int asign=1;
				if (sigma>modelParameters_.orbitals-1) {
					m=1;
					asign= -1;
				}
				typename OperatorType::Su2RelatedType su2related;
				if (sigma <modelParameters_.orbitals) {
					su2related.source.push_back(i*dofs+sigma);
					su2related.source.push_back(i*dofs+sigma + modelParameters_.orbitals);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = modelParameters_.orbitals;
				}

				OperatorType myOp(tmpMatrix,
				                  ProgramGlobals::FermionOrBosonEnum::FERMION,
				                  typename OperatorType::PairType(1,m),
				                  asign,
				                  su2related);
				creationMatrix.push_back(myOp);
			}
		}
	}

private:

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to
	//! basis state ket
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	RealType sign(HilbertState const &ket, int i,SizeType sigma) const
	{
		int value=0;
		SizeType dofs=2*modelParameters_.orbitals;
		for (SizeType alpha=0;alpha<dofs;alpha++)
			value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);
		// add electron on site 0 if needed
		if (i>0) value += HilbertSpaceFeAsType::electrons(ket);

		//order for sign is: a up, b up, a down, b down, etc
		unsigned int x = HilbertSpaceFeAsType::get(ket,i);
		int spin = sigma/modelParameters_.orbitals;
		SizeType orb = sigma % modelParameters_.orbitals;

		for (SizeType j=0;j<orb;j++) {
			for (SizeType k=0;k<2;k++) {
				SizeType ind = j + k * modelParameters_.orbitals;
				int mask = (1<<ind);
				if (x & mask) value++;
			}
		}
		if (spin==SPIN_DOWN) {
			int mask = (1<<orb);
			if (x & mask) value++;
		}
		if (value==0 || value%2==0) return 1.0;

		return FERMION_SIGN;
	}

	//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void findOperatorMatrices(SparseMatrixType& creationMatrix,
	                          int i,
	                          int sigma,
	                          const HilbertBasisType& natBasis) const
	{
		HilbertState bra,ket;
		int n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];

			if (HilbertSpaceFeAsType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceFeAsType::create(bra,i,sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices: internal error\n");
				if (ii==SizeType(jj)) {
					std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra<<" sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("Creation operator cannot be diagonal\n");
				}
				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,cm);
		transposeConjugate(creationMatrix,temp);
	}

	void setSymmetryRelatedInternal(VectorQnType& qns,
	                                const HilbertBasisType& basis,
	                                int n) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;

		qns.resize(basis.size(), QnType::zero());
		VectorSizeType other(2, 0);
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0,0);

			SizeType na = HilbertSpaceFeAsType::calcNofElectrons(basis[i],0) +
			        HilbertSpaceFeAsType::calcNofElectrons(basis[i],0+2);
			SizeType nb = HilbertSpaceFeAsType::calcNofElectrons(basis[i],1) +
			        HilbertSpaceFeAsType::calcNofElectrons(basis[i],1+2);

			SizeType flavor = na  + 3*nb;

			// nup
			SizeType electronsUp = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                                    SPIN_UP);
			// ndown
			SizeType electronsDown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                                      SPIN_DOWN);
			SizeType electrons = electronsDown + electronsUp;

			other[0] = electrons;
			other[1] = electronsUp;
			bool sign = electrons & 1;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	//! Term is U*\sum_{\alpha}n_i^2
	void addNiSquared(SparseMatrixType& hmatrix) const
	{
		SizeType dofs = 2*modelParameters_.orbitals;
		SparseMatrixType nmatrix;
		const VectorOperatorType& cm = creationMatrix_;
		for (SizeType alpha = 0; alpha < dofs; ++alpha)
			nmatrix += n(cm[alpha].getCRS());

		SparseMatrixType nmatrixSquared = nmatrix*nmatrix;
		hmatrix += modelParameters_.hubbardU*nmatrixSquared;
	}

	// Term is P*\sum_i(0.25*ni_a ni_b +
	// 0.5*si_a^\dagger si_b + 0.5*si_b^\dagger si_a +
	// si_a^z si_b^z
	void addPairHopping(SparseMatrixType& hmatrix) const
	{
		SizeType orbitals = modelParameters_.orbitals;
		SparseMatrixType nmatrixA;
		SparseMatrixType nmatrixB;
		const VectorOperatorType& cm = creationMatrix_;
		for (SizeType spin = 0; spin < 2; ++spin) {
			nmatrixA += n(cm[0 + spin*orbitals].getCRS()); // orbital a
			nmatrixB += n(cm[1 + spin*orbitals].getCRS()); // orbital b
		}

		SparseMatrixType nmatrixAB = nmatrixA*nmatrixB;
		hmatrix += static_cast<ComplexOrRealType>(0.25)*modelParameters_.pairHopping*nmatrixAB;

		const SparseMatrixType& splusA =
		        ModelBaseType::naturalOperator("splus", 0, 0).getCRS(); // splus_a
		const SparseMatrixType& splusB =
		        ModelBaseType::naturalOperator("splus", 0, 1).getCRS(); // splus_b
		const SparseMatrixType& sminusA =
		        ModelBaseType::naturalOperator("sminus", 0, 0).getCRS(); // sminus_a
		const SparseMatrixType& sminusB =
		        ModelBaseType::naturalOperator("sminus", 0, 1).getCRS(); // sminus_b
		SparseMatrixType tmp = splusA*sminusB;
		tmp += splusB*sminusA;
		hmatrix += static_cast<ComplexOrRealType>(0.5)*modelParameters_.pairHopping*tmp;

		const SparseMatrixType& szA =
		        ModelBaseType::naturalOperator("sz", 0, 0).getCRS(); // sz_a
		const SparseMatrixType& szB =
		        ModelBaseType::naturalOperator("sz", 0, 1).getCRS(); // sz_b
		tmp = szA*szB;
		hmatrix += modelParameters_.pairHopping*tmp;
	}

	static SparseMatrixType n(const SparseMatrixType& c)
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger, c);
		multiply(tmpMatrix, c, cdagger);

		return tmpMatrix;
	}

	ParametersGrapheneType  modelParameters_;
	const GeometryType& geometry_;
	SizeType statesPerSite_;
	HilbertBasisType basis_;
	VectorQnType qq_;
	VectorOperatorType creationMatrix_;
}; //class Graphene
} // namespace Dmrg
/*@}*/
#endif

