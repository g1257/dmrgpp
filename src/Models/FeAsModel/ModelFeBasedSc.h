/*
Copyright (c) 2009-2018, UT-Battelle, LLC
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

/*! \file ModelFeBasedSc.h
 *
 *  An implementation of a FeBasedSc model for Fe-based superconductors to use with
 *  the DmrgSolver
 *
 */
#ifndef MODEL_FEAS_DMRG
#define MODEL_FEAS_DMRG
#include "ModelBase.h"
#include "ParametersModelFeAs.h"
#include "HilbertSpaceFeAs.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include "Geometry/GeometryDca.h"
#include "FeAsJzSymmetry.h"
#include <numeric>

namespace Dmrg {
template<typename ModelBaseType>
class ModelFeBasedSc : public ModelBaseType {

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
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ParametersModelFeAs<ComplexOrRealType, QnType> ParamsModelFeAsType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef FeAsJzSymmetry<HilbertBasisType,
	VectorOperatorType,
	PsimagLite::IsComplexNumber<ComplexOrRealType>::True> FeAsJzSymmetryType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

	ModelFeBasedSc(const SolverParamsType& solverParams,
	               InputValidatorType& io,
	               GeometryType const &geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      reinterpretX_(6),
	      reinterpretY_(9),
	      modelParameters_(io),
	      geometry_(geometry),
	      geometryDca_(geometry, modelParameters_.orbitals),
	      spinSquared_(spinSquaredHelper_,
	                   modelParameters_.orbitals,
	                   2*modelParameters_.orbitals),
	      reinterpret_(!modelParameters_.jzSymmetry),
	      feAsJzSymmetry_(modelParameters_.jzSymmetry)
	{
		ProgramGlobals::init(modelParameters_.orbitals*geometry_.numberOfSites() + 1);

		PsimagLite::String tspAlgo = "";
		try {
			io.readline(tspAlgo,"TSPAlgorithm=");
		} catch (std::exception&) {}

		if (tspAlgo == "SuzukiTrotter") reinterpret_ = false;

		SizeType v1 = 2*modelParameters_.orbitals*geometry.numberOfSites();
		SizeType v2 = v1*modelParameters_.orbitals;
		if (modelParameters_.potentialV.size() != v1 &&
		        modelParameters_.potentialV.size() != v2) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "potentialV length must be 2*orbitals times the number of sites or";
			str += " 2*orbitals*orbitals times the number of sites\n";
			err(str.c_str());
		}

		HilbertSpaceFeAsType::setOrbitals(modelParameters_.orbitals);
		statesPerSite_ = (1 << (modelParameters_.orbitals*2));

		VectorSizeType block(1,0);
		int sitesTimesDof = 2*modelParameters_.orbitals;
		HilbertState total = (1<<sitesTimesDof);
		basis_.resize(total);
		for (HilbertState a = 0; a < total; ++a) basis_[a] = a;

		if (basis_.size() == 16) {
			basis_[0] = 0;
			basis_[1] = 4;
			basis_[2] = 8;
			basis_[3] = 1;
			basis_[4] = 2;
			basis_[5] = 12;
			basis_[6] = 5;
			basis_[7] = 6;
			basis_[8] = 9;
			basis_[9] = 10;
			basis_[10] = 3;
			basis_[11] = 13;
			basis_[12] = 14;
			basis_[13] = 7;
			basis_[14] = 11;
			basis_[15] = 15;
		}

		if (basis_.size() == 64) {
			SizeType counter = 0;
			basis_[counter++] = 0;
			basis_[counter++] = 8;
			basis_[counter++] = 16;
			basis_[counter++] = 32;
			basis_[counter++] = 1;
			basis_[counter++] = 2;
			basis_[counter++] = 4;
			basis_[counter++] = 24;
			basis_[counter++] = 40;
			basis_[counter++] = 48;
			basis_[counter++] = 9;
			basis_[counter++] = 10;
			basis_[counter++] = 12;
			basis_[counter++] = 17;
			basis_[counter++] = 18;
			basis_[counter++] = 20; // 15
			basis_[counter++] = 33;
			basis_[counter++] = 34;
			basis_[counter++] = 36;
			basis_[counter++] = 3;
			basis_[counter++] = 5;
			basis_[counter++] = 6;
			basis_[counter++] = 56;
			basis_[counter++] = 25;
			basis_[counter++] = 26;
			basis_[counter++] = 28;
			basis_[counter++] = 41;
			basis_[counter++] = 42;
			basis_[counter++] = 44;
			basis_[counter++] = 49;
			basis_[counter++] = 50;
			basis_[counter++] = 52; // 31
			basis_[counter++] = 11;
			basis_[counter++] = 13;
			basis_[counter++] = 14;
			basis_[counter++] = 19;
			basis_[counter++] = 21;
			basis_[counter++] = 22;
			basis_[counter++] = 35;
			basis_[counter++] = 37;
			basis_[counter++] = 38;
			basis_[counter++] = 7;
			basis_[counter++] = 57;
			basis_[counter++] = 58;
			basis_[counter++] = 60;
			basis_[counter++] = 27;
			basis_[counter++] = 29;
			basis_[counter++] = 30; // 47
			basis_[counter++] = 43;
			basis_[counter++] = 45;
			basis_[counter++] = 46;
			basis_[counter++] = 51;
			basis_[counter++] = 53;
			basis_[counter++] = 54;
			basis_[counter++] = 15;
			basis_[counter++] = 23;
			basis_[counter++] = 39;
			basis_[counter++] = 59;
			basis_[counter++] = 61;
			basis_[counter++] = 62;
			basis_[counter++] = 31;
			basis_[counter++] = 47;
			basis_[counter++] = 55;
			basis_[counter++] = 63;

			assert(counter == 64);
		}

		SizeType sum = std::accumulate(basis_.begin(), basis_.end(), 0);
		SizeType n = basis_.size();
		if (sum != n*(n-1)/2)
			err("ModelFeBasedSc: basis set up wrong\n");

		setOperatorMatricesInternal(creationMatrix_, block);
		if (feAsJzSymmetry_.isEnabled() && !feAsJzSymmetry_.isSet())
			feAsJzSymmetry_.init(basis_, creationMatrix_);

		setSymmetryRelatedInternal(qq_,basis_,1);

		if (feAsJzSymmetry_.isEnabled()) {
			basis_[9] = 9;
			basis_[10] = 10;
			basis_[11] = 48;

			basis_[17] = 3;
			basis_[18] = 34;
			basis_[19] = 36;

			basis_[41] = 27;
			basis_[42] = 29;
			basis_[43] = 60;

			basis_[52] = 15;
			basis_[53] = 53;
			basis_[54] = 54;
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		io.write(label + "/reinterpretX_", reinterpretX_);
		io.write(label + "/reinterpretY_", reinterpretY_);
		modelParameters_.write(label, io);
		geometryDca_.write(label, io);
		spinSquaredHelper_.write(label, io);
		spinSquared_.write(label, io);
		io.write(label + "/reinterpret_", reinterpret_);
		io.write(label + "/statesPerSite_", statesPerSite_);
		io.write(label + "/basis_", basis_);
		io.write(label + "/qq_", qq_);
		io.write(label + "/creationMatrix_", creationMatrix_);
		feAsJzSymmetry_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType time) const
	{
		SizeType n=block.size();

		for (SizeType i = 0; i < n; ++i) {
			SizeType site = block[i];
			const VectorOperatorType& cm = ModelBaseType::trackableOps(site);
			addInteraction(hmatrix,cm,i,block[i]);
			addMagneticField(hmatrix,cm,i,block[i]);
			addSpinOrbit(hmatrix,cm,i);

			if (modelParameters_.potentialT.size()==0 || time==0) {
				addPotentialV(hmatrix,
				              cm,
				              i,
				              block[i],
				              modelParameters_.potentialV);
			} else {
				addPotentialV(hmatrix,
				              cm,
				              i,
				              block[i],
				              modelParameters_.potentialT);
			}
		}
	}

	void fillLabeledOperators(VectorQnType& qns)
	{
		qns = qq_;
		assert(creationMatrix_.size()>0);
		SizeType nrow = creationMatrix_[0].data.rows();

		OpsLabelType& splus = this->createOpsLabel("splus");
		OpsLabelType& sminus = this->createOpsLabel("sminus");
		for (SizeType dof = 0; dof < modelParameters_.orbitals; ++dof) {
			MatrixType tmp(nrow,nrow);
			tmp += multiplyTc(creationMatrix_[dof].data,
			                  creationMatrix_[dof + modelParameters_.orbitals].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			splus.push(OperatorType(tmp2,
			                        1.0,
			                        typename OperatorType::PairType(0,0),
			                        1.0,
			                        su2Related));
			SparseMatrixType tmp3;
			transposeConjugate(tmp3, tmp2);
			sminus.push(OperatorType(tmp3,
			                         1.0,
			                         typename OperatorType::PairType(0,0),
			                         1.0,
			                         su2Related));
		}

		OpsLabelType& sz = this->createOpsLabel("sz");
		for (SizeType dof = 0; dof < modelParameters_.orbitals; ++dof) {
			MatrixType tmp(nrow,nrow);
			MatrixType tmp2(nrow,nrow);

			tmp += multiplyTc(creationMatrix_[dof].data,creationMatrix_[dof].data);
			tmp2 += multiplyTc(creationMatrix_[dof+modelParameters_.orbitals].data,
			        creationMatrix_[dof+modelParameters_.orbitals].data);

			tmp = 0.5*(tmp-tmp2);
			SparseMatrixType tmp3(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			sz.push(OperatorType(tmp3,
			                     1.0,
			                     typename OperatorType::PairType(0,0),
			                     1.0,
			                     su2Related));
		}

		OpsLabelType& nop = this->createOpsLabel("n");
		for (SizeType dof = 0; dof < 2*modelParameters_.orbitals; ++dof) {
			MatrixType tmp =
			        multiplyTc(creationMatrix_[dof].data,creationMatrix_[dof].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			nop.push(OperatorType(tmp2,
			                      1.0,
			                      typename OperatorType::PairType(0,0),
			                      1.0,
			                      su2Related));
		}

		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& cc = this->createOpsLabel("C");
		for (SizeType dof = 0; dof < 2*modelParameters_.orbitals; ++dof) {
			VectorOperatorType cm = creationMatrix_;
			cc.push(creationMatrix_[dof]);
			cm[dof].dagger();
			c.push(cm[dof]);
		}

		OpsLabelType& d = this->createOpsLabel("d");
		for (SizeType dof = 0; dof < modelParameters_.orbitals; ++dof) {
			SizeType orbital = dof % modelParameters_.orbitals;
			SparseMatrixType atmp;
			multiply(atmp,
			         creationMatrix_[orbital + modelParameters_.orbitals].data,
			        creationMatrix_[orbital].data);
			typename OperatorType::Su2RelatedType su2Related;
			d.push(OperatorType(atmp,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related));
		}

		this->makeTrackable("C");
	}

	/* PSIDOC FeAs::fillModelLinks
	Here we have as many connections as  2*orbitals*orbitals:
	a up a up, a up b up, b up a up, b up, b up, etc
	and similarly for spin down.
	We'll do all of this in one single term, named ``hopping'' and put it in a C++ variable
	called \verb!hop!. So we'll do 2*orbitals*orbitals pushes into that \verb!hop!.
	Each push does the pair c1 with normal 'N', and c2 with transpose conjugate 'C'.
	(We shall not discuss the SU(2) related numbers.)
	We observe that c1 is a ``C'' trackable operator with degree of freedom
	orb1 + spin*orbitals; likewise c2 is a ``C'' trackable operator with orb2 + spin*orbitals.
	The orbitals are different but the spin is the same, as expected.
	In addition, lines (B) and (C) have a 3rd number for each operator.
	This 3rd number indicates the connector dependence, which is not on spin but only on orbital.
	That is why the 3rd number in (B) is orb1 and not orb1 + spin*orbitals, because
	the hoppings in the input file do not depend on spin.
	This is very important to note and could be a cause of confusion.
	By default this 3rd number (which can be omitted) is 0, and indicates no
	dependence of the connector on things other than site. The site dependence
	is handled by the geometry and must not be specified here.
	PSIDOCCOPY $FirstFunctionBelow

	In formulas, we can explain the distinction between the 2nd number and the 3rd number
	by writing
	\[
	t(orb1, orb2)\,\, c^\dagger_{orb1, spin1} (site1)\,\, c_{orb2, spin2} (site2)
	\]
	The first operator has label $orb1, spin1$
	that gets packed into $orb1 + spin1*orbitals$ and must be given as the 2nd argument
	for the first operator, and similarly for the second operator.
	In contrast, the indices orb1 and orb2 are in the connector dependence
	$t(orb1, orb2)$, and must be as
	3rd number of each operator, respectively.
	 */
	void fillModelLinks()
	{
		const SizeType orbitals = modelParameters_.orbitals;
		ModelTermType& hop = ModelBaseType::createTerm("hopping");//(A)
		for (SizeType spin = 0; spin < 2; ++spin) {
			for (SizeType orb1 = 0; orb1 < orbitals; ++orb1) {
				OpForLinkType c1("C", orb1 + spin*orbitals, orb1); // (B)
				for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
					OpForLinkType c2("C", orb2 + spin*orbitals, orb2); // (C)

					hop.push(c1, 'N', c2, 'C', 1, (spin == 1) ? -1 : 1, spin);
				}
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
				                  -1,
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

		reinterpret(cm,natBasis);

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
			if (n == 1) jmpair = calcJmvalue<PairType>(basis[i]);

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

			if (modelParameters_.spinOrbit.rows() > 0 && !modelParameters_.jzSymmetry)
				electronsUp = 0;

			feAsJzSymmetry_.setElectronsAndJz(electrons, electronsUp, i);

			other[0] = electrons;
			other[1] = electronsUp;
			bool sign = electrons & 1;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// Reinterprets 6 and 9
	template<typename PairType>
	PairType calcJmvalue(const HilbertState& ket) const
	{
		PairType jm(0,0);
		if (modelParameters_.orbitals!=2) return jm;
		SizeType x=reinterpretX_,y=reinterpretY_; // these states need reinterpretation

		if (ket==x) {
			jm=std::pair<SizeType,SizeType>(2,1);
		} else if (ket==y) {
			jm=std::pair<SizeType,SizeType>(0,0);
		} else jm=calcJmValueAux<PairType>(ket);

		return jm;
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// does not work for 6 or 9
	template<typename PairType>
	PairType calcJmValueAux(const HilbertState& ket) const
	{
		SizeType site0=0;
		SizeType site1=0;

		spinSquared_.doOnePairOfSitesA(ket,site0,site1);
		spinSquared_.doOnePairOfSitesB(ket,site0,site1);
		spinSquared_.doDiagonal(ket,site0,site1);

		RealType sz = spinSquared_.spinZ(ket,site0);
		PairType jm= spinSquaredHelper_.getJmPair(sz);
		return jm;
	}

	//! SU(2) symmetry related block
	//! Let |9> = |up a down b> and
	//! Let |6> = |up b down a>  then
	void reinterpret(MatrixType& cm,
	                 const HilbertBasisType& basis) const
	{
		if (!reinterpret_ || modelParameters_.orbitals!=2) return;

		int n  = cm.rows();
		if (n!=16)
			throw PsimagLite::RuntimeError("blocks.size must be 1, and basis.size 16\n");

		MatrixType cmCopy(n,n);
		int i,j;
		int x=PsimagLite::indexOrMinusOne(basis,reinterpretX_);
		int y=PsimagLite::indexOrMinusOne(basis,reinterpretY_);

		RealType factor = 0.7071067811865475244;
		for (i=0;i<n;i++) {
			if (i==x || i==y) continue;
			for (j=0;j<n;j++) {
				if (j==x || j==y) continue;
				cmCopy(i,j)=cm(i,j);
			}
		}
		for (j=0;j<n;j++) {
			if (j==x || j==y) continue;
			cmCopy(x,j)=factor*(cm(x,j)+cm(y,j));
			cmCopy(y,j)=factor*(cm(x,j)-cm(y,j));
		}
		for (i=0;i<n;i++) {
			if (i==x || i==y) continue;
			cmCopy(i,x)=factor*(cm(i,x)+cm(i,y));
			cmCopy(i,y)=factor*(cm(i,x)-cm(i,y));
		}

		const RealType sf= -1;
		const RealType zeroPointFive = 0.5;

		cmCopy(x,x) = zeroPointFive*(cm(x,x)+cm(x,y)+cm(y,x)+cm(y,y));
		cmCopy(x,y) = zeroPointFive*(cm(x,x)-sf*cm(x,y)+sf*cm(y,x)-cm(y,y));
		cmCopy(y,x) = zeroPointFive*(cm(x,x)+sf*cm(x,y)-sf*cm(y,x)-cm(y,y));
		cmCopy(y,y) = zeroPointFive*(cm(x,x)-cm(x,y)-cm(y,x)+cm(y,y));

		cm = cmCopy;
	}

	void addInteraction(SparseMatrixType &hmatrix,
	                    const VectorOperatorType& cm,
	                    SizeType i,
	                    SizeType actualSite) const
	{
		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_ORBITAL0) {
			return addInteractionAncilla(hmatrix,cm,i,actualSite);
		}

		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_KSPACE) {
			return addInteractionKspace(hmatrix,cm,i,actualSite);
		}

		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_IMPURITY) {
			return addInteractionImpurity(hmatrix,cm,i,actualSite);
		}

		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_CODE2) {
			return addInteractionUmatrix(hmatrix,cm,i);
		}

		addInteractionU1(hmatrix,cm,i);
		addInteractionU2(hmatrix,cm,i);
		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_PAPER33) {
			addInteractionJ1(hmatrix,cm,i);
			addInteractionJ2(hmatrix,cm,i);
		} else {
			addInteractionV(hmatrix,cm,i);
		}
	}

	RealType findHubbardU(SizeType index, SizeType orb1, SizeType orb2) const
	{
		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_PAPER33 ||
		        modelParameters_.feAsMode == ParamsModelFeAsType::INT_IMPURITY) {
			assert(index < modelParameters_.hubbardU.size());
			return modelParameters_.hubbardU[index];
		}

		assert(orb1 + orb2*modelParameters_.orbitals < modelParameters_.hubbardU.size());
		return modelParameters_.hubbardU[orb1 + orb2*modelParameters_.orbitals];
	}

	//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteractionU1(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType tmpMatrix;

		for (SizeType alpha=0;alpha<SizeType(modelParameters_.orbitals);alpha++) {
			SparseMatrixType m1=cm[alpha+SPIN_UP*modelParameters_.orbitals+i*dof].data;
			SparseMatrixType m2=cm[alpha+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;

			multiply(tmpMatrix,n(m1),n(m2));
			hmatrix += findHubbardU(0,alpha,alpha)*tmpMatrix;
		}
	}

	//! Term is U[1] n_{i BAND0 } n_{i BAND1}
	void addInteractionU2(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i) const
	{

		for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
			for (SizeType orb2=orb1+1;orb2<modelParameters_.orbitals;orb2++) {
				SparseMatrixType tmpMatrix;
				multiply(tmpMatrix, nSummedOverSpin(cm,i,orb1),nSummedOverSpin(cm,i,orb2));
				hmatrix += findHubbardU(1,orb1,orb2)*tmpMatrix;
			}
		}
	}

	//! Term is U[2] S_{i BAND0 } S_{i BAND1}
	void addInteractionJ1(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i) const
	{
		RealType val=0;
		RealType val2=2.0;
		RealType val3=4.0;

		for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
			for (SizeType orb2=orb1+1;orb2<modelParameters_.orbitals;orb2++) {
				SparseMatrixType tmpMatrix;

				multiply(tmpMatrix,
				         spinOperator(cm,i,orb1,0),
				         spinOperator(cm,i,orb2,1));
				val = modelParameters_.hubbardU[2]/val2;
				// this is -2*J
				hmatrix += val*tmpMatrix;

				multiply(tmpMatrix,
				         spinOperator(cm,i,orb1,1),
				         spinOperator(cm,i,orb2,0));
				val = modelParameters_.hubbardU[2]/val2;
				// this is -2*J
				hmatrix += val*tmpMatrix;

				multiply(tmpMatrix,
				         spinOperator(cm,i,orb1,2),
				         spinOperator(cm,i,orb2,2));
				val = modelParameters_.hubbardU[4]/val3;
				// this is -2*J
				hmatrix += val*tmpMatrix;
			}
		}
	}

	//! Term is U[3] \sum_{\alpha}\bar{n}_{i\alpha UP} \bar{n}_{i\alpha DOWN}
	//! where \bar{n}_{i\alpha \spin} = c^\dagger_{i\alpha\spin} c_{i\bar{\alpha}\bar{spin}}
	void addInteractionJ2(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i) const
	{
		SparseMatrixType tmpMatrix;

		for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
			for (SizeType orb2=0;orb2<modelParameters_.orbitals;orb2++) {
				if (orb1==orb2) continue;
				multiply(tmpMatrix,
				         nBar(cm,i,orb1,orb2,SPIN_UP),
				         nBar(cm,i,orb1,orb2,SPIN_DOWN));
				// -J
				hmatrix += modelParameters_.hubbardU[3]*tmpMatrix;
			}
		}
	}

	void addInteractionV(SparseMatrixType& hmatrix,
	                     const VectorOperatorType& cm,
	                     SizeType site) const
	{
		assert(modelParameters_.orbitals == 3);

		SizeType dofs = 2 * modelParameters_.orbitals;

		SparseMatrixType tmpMatrix1;
		SparseMatrixType tmpMatrix2;
		SparseMatrixType tmpMatrix;

		RealType value = modelParameters_.coulombV;

		for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
			for (SizeType spin2 = 0; spin2 < 2; ++spin2) {
				if (spin1 == spin2) continue;

				multiply(tmpMatrix1,
				         cm[2+spin1*modelParameters_.orbitals+site*dofs].data,
				        cm[0+spin2*modelParameters_.orbitals+site*dofs].data);

				SparseMatrixType cmTranspose1;
				transposeConjugate(cmTranspose1,
				                   cm[1+spin2*modelParameters_.orbitals+site*dofs].data);

				SparseMatrixType cmTranspose2;
				transposeConjugate(cmTranspose2,
				                   cm[1+spin1*modelParameters_.orbitals+site*dofs].data);

				multiply(tmpMatrix2,cmTranspose1,cmTranspose2);

				multiply(tmpMatrix,tmpMatrix1,tmpMatrix2);
				tmpMatrix1 = value*tmpMatrix;
				hmatrix += tmpMatrix1;

				transposeConjugate(tmpMatrix2,tmpMatrix1);
				hmatrix += tmpMatrix2;
			}
		}
	}

	void addSpinOrbit(SparseMatrixType &hmatrix,
	                  const VectorOperatorType& cm,
	                  SizeType i) const
	{
		if (modelParameters_.spinOrbit.rows()<4) return;
		SizeType orbitals = modelParameters_.orbitals;
		int dof=2*orbitals;
		SizeType nrow = cm[0].data.rows();
		MatrixType tmp(nrow,nrow);
		SparseMatrixType tmpMatrix;

		for (SizeType spin1=0;spin1<2;spin1++) {
			for (SizeType spin2=0;spin2<2;spin2++) {
				for (SizeType orb1=0;orb1<orbitals;orb1++) {
					for (SizeType orb2=0;orb2<orbitals;orb2++) {


						SparseMatrixType c1 = cm[orb1+spin1*orbitals+i*dof].data;

						SparseMatrixType c2 = cm[orb2+spin2*orbitals+i*dof].data;

						tmp += static_cast<RealType>(0.5)*
						        modelParameters_.spinOrbit(spin1*2+spin2,
						                                   orb1*orbitals + orb2)*
						        multiplyTc(c1,c2);
					}
				}
			}
		}

		fullMatrixToCrsMatrix(tmpMatrix,tmp);
		hmatrix += tmpMatrix;
	}

	void addMagneticField(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      SizeType actualIndexOfSite) const
	{
		if (modelParameters_.magneticField.rows()<3) return;
		assert(actualIndexOfSite<modelParameters_.magneticField.n_col());
		for (SizeType orb=0;orb<modelParameters_.orbitals;orb++)
			addMagneticField(hmatrix,cm,i,actualIndexOfSite,orb);
	}

	void addMagneticField(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      SizeType actualIndexOfSite,
	                      SizeType orbital) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType cup = cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data;
		SparseMatrixType cupTranspose;
		transposeConjugate(cupTranspose,cup);
		SparseMatrixType cdown = cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;
		SparseMatrixType A = cupTranspose * cdown;
		SparseMatrixType Atranspose;
		transposeConjugate(Atranspose,A);

		hmatrix += modelParameters_.magneticField(0,actualIndexOfSite) * A;

		hmatrix += modelParameters_.magneticField(1,actualIndexOfSite) * Atranspose;

		SparseMatrixType nup = n(cup);
		SparseMatrixType ndown = n(cdown);

		SparseMatrixType tmp = nup;
		const RealType f1 = (-1.0);
		tmp += f1*ndown;

		hmatrix += modelParameters_.magneticField(2,actualIndexOfSite) * tmp;

	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorOperatorType& cm,
	                   SizeType i,
	                   SizeType actualIndexOfSite,
	                   const typename PsimagLite::Vector<RealType>::Type& V) const
	{
		SizeType v1 = 2*modelParameters_.orbitals*geometry_.numberOfSites();
		SizeType v2 = v1*modelParameters_.orbitals;
		if (V.size() != v1 && V.size() != v2) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "potentialV[T] length must be 2*orbitals times the number of sites or";
			str += " 2*orbitals*orbitals times the number of sites\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		if (V.size() == v1) {
			for (SizeType orb=0;orb<modelParameters_.orbitals;orb++)
				addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb,V);
		}

		if (V.size() == v2) {
			for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {
				for (SizeType orb2=0;orb2<modelParameters_.orbitals;orb2++) {
					addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb,orb2,V);
				}
			}

			return;
		}
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorOperatorType& cm,
	                   SizeType i,
	                   SizeType actualIndexOfSite,
	                   SizeType orbital,
	                   const typename PsimagLite::Vector<RealType>::Type& V) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType nup = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data);
		SparseMatrixType ndown = n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data);

		SizeType linSize = geometry_.numberOfSites();

		SizeType iUp = actualIndexOfSite + (orbital + 0*modelParameters_.orbitals)*linSize;
		hmatrix += V[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orbital + 1*modelParameters_.orbitals)*linSize;
		hmatrix += V[iDown] * ndown;
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorOperatorType& cm,
	                   SizeType i,
	                   SizeType actualIndexOfSite,
	                   SizeType orb,
	                   SizeType orb2,
	                   const typename PsimagLite::Vector<RealType>::Type& V) const
	{
		int dof=2*modelParameters_.orbitals;
		SizeType orbitalsSquared = modelParameters_.orbitals*modelParameters_.orbitals;

		SparseMatrixType nup = nEx(cm[orb+SPIN_UP*modelParameters_.orbitals+i*dof].data,
		        cm[orb2+SPIN_UP*modelParameters_.orbitals+i*dof].data);
		SparseMatrixType ndown = nEx(cm[orb+SPIN_DOWN*modelParameters_.orbitals+i*dof].data,
		        cm[orb2+SPIN_DOWN*modelParameters_.orbitals+i*dof].data);

		SizeType linSize = geometry_.numberOfSites();

		SizeType iUp = actualIndexOfSite + (orb + orb2*modelParameters_.orbitals +
		                                    0*orbitalsSquared)*linSize;
		hmatrix += V[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orb + orb2*modelParameters_.orbitals +
		                                      1*orbitalsSquared)*linSize;
		hmatrix += V[iDown] * ndown;
	}

	SparseMatrixType nEx(const SparseMatrixType& c1, const SparseMatrixType& c2) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c2);
		multiply(tmpMatrix,c1,cdagger);

		return tmpMatrix;
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

		return tmpMatrix;
	}

	SparseMatrixType nBar(const VectorOperatorType& cm,
	                      SizeType i,
	                      SizeType orb1,
	                      SizeType orb2,
	                      SizeType spin) const
	{
		SizeType dofs = 2 * modelParameters_.orbitals;
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger=cm[orb1+spin*modelParameters_.orbitals+i*dofs].data;
		SparseMatrixType cbar;
		transposeConjugate(cbar,cm[orb2+(1-spin)*modelParameters_.orbitals+i*dofs].data);
		multiply(tmpMatrix,cdagger,cbar);
		return tmpMatrix;
	}

	SparseMatrixType nSummedOverSpin(const VectorOperatorType& cm,
	                                 SizeType i,
	                                 SizeType orbital) const
	{
		SizeType dofs = 2 * modelParameters_.orbitals;
		SparseMatrixType tmpMatrix = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dofs].data);
		tmpMatrix += n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dofs].data);
		return tmpMatrix;
	}

	SparseMatrixType spinOperator(const VectorOperatorType& cm,
	                              SizeType i,
	                              SizeType orbital,
	                              SizeType component) const
	{
		switch (component) {
		case 0: // S^+
			return spinOperatorAux(cm,i,orbital,SPIN_UP,SPIN_DOWN);
		case 1: // S^-
			return spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_UP);
		}

		SparseMatrixType tmpMatrix=spinOperatorAux(cm,i,orbital,SPIN_UP,SPIN_UP);
		SparseMatrixType tmpMatrix2=spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_DOWN);
		const RealType f1 = (-1.0);
		tmpMatrix += f1*tmpMatrix2;
		return tmpMatrix;
	}

	SparseMatrixType spinOperatorAux(const VectorOperatorType& cm,
	                                 SizeType i,
	                                 SizeType orbital,
	                                 SizeType spin1,
	                                 SizeType spin2) const
	{
		SizeType dofs = 2 * modelParameters_.orbitals;
		SparseMatrixType result,temp;
		transposeConjugate(temp,cm[orbital+spin2*modelParameters_.orbitals+i*dofs].data);
		multiply(result, // =
		         cm[orbital+spin1*modelParameters_.orbitals+i*dofs].data, // times
		        temp);

		return result;
	}

	//! only for feAsMode == 2
	void addInteractionUmatrix(SparseMatrixType &hmatrix,
	                           const VectorOperatorType& cm,
	                           SizeType site) const
	{
		const typename PsimagLite::Vector<RealType>::Type& U = modelParameters_.hubbardU;
		SizeType orbitals = modelParameters_.orbitals;
		SizeType dofs = orbitals * 2;
		for (SizeType interaction = 0; interaction < 2; ++interaction) {
			for (SizeType orb1=0;orb1<orbitals;orb1++) {
				for (SizeType orb2=0;orb2<orbitals;orb2++) {
					for (SizeType spin = 0; spin < 2; ++spin) {
						SizeType spin2 = (interaction == 0) ? spin : 1 - spin;
						const SparseMatrixType& cm1 = cm[orb1+spin*orbitals+site*dofs].data;
						const SparseMatrixType& cm2 = cm[orb1+spin2*orbitals+site*dofs].data;
						SizeType offset = orb1 + orb2*orbitals;
						if (interaction == 1) offset += orbitals * orbitals;
						SparseMatrixType tmpMatrix;

						multiply(tmpMatrix, n(cm1),n(cm2));
						assert(offset < U.size());

						hmatrix += U[offset]*tmpMatrix;
					}
				}
			}
		}
	}

	//! only for feAsMode == 4
	void addInteractionKspace(SparseMatrixType &hmatrix,
	                          const VectorOperatorType& cm,
	                          SizeType i,
	                          SizeType actualSite) const
	{
		if (actualSite > 0) return;

		SizeType orbs = modelParameters_.orbitals;
		SizeType dofs = orbs * 2;

		for (SizeType orb1=0;orb1<orbs;orb1++) {
			const SparseMatrixType& cm1 = cm[orb1+SPIN_UP*orbs+i*dofs].data;
			for (SizeType orb2=0;orb2<orbs;orb2++) {
				const SparseMatrixType& cm2 = cm[orb2+SPIN_UP*orbs+i*dofs].data;
				SparseMatrixType tmpMatrix(multiplyTc(cm1,cm2));
				for (SizeType orb3=0;orb3<orbs;orb3++) {
					const SparseMatrixType& cm3 = cm[orb3+SPIN_DOWN*orbs+i*dofs].data;

					SizeType orb4 = getMomentum(orb1, orb2, orb3);
					const SparseMatrixType& cm4 = cm[orb4+SPIN_DOWN*orbs+i*dofs].data;

					SparseMatrixType tmpMatrix2(multiplyTc(cm3,cm4));

					SparseMatrixType tmpMatrix3;
					multiply(tmpMatrix3, tmpMatrix2, tmpMatrix);

					hmatrix += modelParameters_.hubbardU[0]*tmpMatrix3;
				}
			}
		}
	}

	SizeType getMomentum(SizeType orb1, SizeType orb2, SizeType orb3) const
	{
		assert(orb1 < modelParameters_.orbitals);
		assert(orb2 < modelParameters_.orbitals);
		assert(orb3 < modelParameters_.orbitals);

		SizeType tmp = geometryDca_.kSum(orb3,orb1);
		SizeType orb4 = geometryDca_.kSustract(tmp,orb2);

		assert(orb4 < modelParameters_.orbitals);
		return orb4;
	}

	//! only for feAsMode == 3
	void addInteractionImpurity(SparseMatrixType &hmatrix,
	                            const VectorOperatorType& cm,
	                            SizeType i,
	                            SizeType actualSite) const
	{
		if (actualSite > 0) return;

		addInteractionU1(hmatrix,cm,i);
		addInteractionImp2(hmatrix,cm,i);
		addInteractionImp3(hmatrix,cm,i);
		addInteractionImp4(hmatrix,cm,i);
	}

	void addInteractionImp2(SparseMatrixType &hmatrix,
	                        const VectorOperatorType& cm,
	                        SizeType site) const
	{
		const typename PsimagLite::Vector<RealType>::Type& U = modelParameters_.hubbardU;
		SizeType orbitals = modelParameters_.orbitals;
		SizeType dofs = orbitals * 2;

		for (SizeType orb1=0;orb1<orbitals;orb1++) {
			for (SizeType orb2=orb1+1;orb2<orbitals;orb2++) {
				for (SizeType spin = 0; spin < 2; ++spin) {
					const SparseMatrixType& cm1 = cm[orb1+spin*orbitals+site*dofs].data;
					const SparseMatrixType& cm2 = cm[orb2+spin*orbitals+site*dofs].data;

					SparseMatrixType tmpMatrix;

					multiply(tmpMatrix, n(cm1),n(cm2));
					hmatrix += U[1]*tmpMatrix;
				}
			}
		}
	}

	void addInteractionImp3(SparseMatrixType &hmatrix,
	                        const VectorOperatorType& cm,
	                        SizeType site) const
	{
		const typename PsimagLite::Vector<RealType>::Type& U = modelParameters_.hubbardU;
		SizeType orbitals = modelParameters_.orbitals;
		SizeType dofs = orbitals * 2;

		for (SizeType orb1=0;orb1<orbitals;orb1++) {
			for (SizeType orb2=0;orb2<orbitals;orb2++) {

				if (orb1 == orb2) continue;

				const SparseMatrixType& cm1 = cm[orb1+0*orbitals+site*dofs].data;
				const SparseMatrixType& cm2 = cm[orb2+1*orbitals+site*dofs].data;

				SparseMatrixType tmpMatrix;

				multiply(tmpMatrix, n(cm1),n(cm2));
				hmatrix += U[2]*tmpMatrix;
			}
		}
	}

	void addInteractionImp4(SparseMatrixType &hmatrix,
	                        const VectorOperatorType& cm,
	                        SizeType site) const
	{
		const typename PsimagLite::Vector<RealType>::Type& U = modelParameters_.hubbardU;
		SizeType orbitals = modelParameters_.orbitals;
		SizeType dofs = orbitals * 2;

		for (SizeType type =0; type < 2; type++) {
			for (SizeType orb1=0;orb1<orbitals;orb1++) {
				for (SizeType orb2=0;orb2<orbitals;orb2++) {

					if (orb1 == orb2) continue;

					SizeType orb3 = (type == 0) ? orb2 : orb1;
					SizeType orb4 = (type == 0) ? orb1 : orb2;
					const SparseMatrixType& cm1 = cm[orb1+0*orbitals+site*dofs].data;
					const SparseMatrixType& cm2 = cm[orb2+0*orbitals+site*dofs].data;
					const SparseMatrixType& cm3 = cm[orb3+1*orbitals+site*dofs].data;
					const SparseMatrixType& cm4 = cm[orb4+1*orbitals+site*dofs].data;

					SparseMatrixType tmpMatrix(multiplyTc(cm1,cm2));
					SparseMatrixType tmpMatrix2(multiplyTc(cm3,cm4));
					hmatrix += U[3]*tmpMatrix*tmpMatrix2;
				}
			}
		}
	}

	//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteractionAncilla(SparseMatrixType &hmatrix,
	                           const VectorOperatorType& cm,
	                           SizeType i,
	                           SizeType actualSite) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType tmpMatrix;

		SizeType alpha = 0; // real sites, no ancilla
		SparseMatrixType m1=cm[alpha+SPIN_UP*modelParameters_.orbitals+i*dof].data;
		SparseMatrixType m2=cm[alpha+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;

		multiply(tmpMatrix,n(m1),n(m2));
		assert(actualSite < modelParameters_.hubbardU.size());
		hmatrix += modelParameters_.hubbardU[actualSite]*tmpMatrix;
	}

	void diagTest(const SparseMatrixType& fullm,const PsimagLite::String& str) const
	{
		if (fullm.rank()!=256) return;
		MatrixType fullm2;
		crsMatrixToFullMatrix(fullm2,fullm);
		typename PsimagLite::Vector<RealType>::Type eigs(fullm2.rows());
		PsimagLite::diag(fullm2,eigs,'V');
		std::cout<<str<<" diagTest size="<<fullm.rank()<<" eigs[0]="<<eigs[0]<<"\n";
		std::cout<<fullm;
	}

	//serializr start class ModelFeBasedSc
	//serializr vptr
	//serializr normal reinterpretX_
	HilbertState reinterpretX_;
	//serializr normal reinterpretY_
	HilbertState reinterpretY_;
	//serializr normal modelParameters_
	ParamsModelFeAsType  modelParameters_;
	//serializr ref geometry_ start
	const GeometryType& geometry_;
	GeometryDcaType geometryDca_;
	//serializr normal spinSquaredHelper_
	SpinSquaredHelper<RealType,HilbertState> spinSquaredHelper_;
	//serializr normal spinSquared_
	SpinSquared<SpinSquaredHelper<RealType,HilbertState> > spinSquared_;
	bool reinterpret_;
	//serializr normal statesPerSite_
	SizeType statesPerSite_;
	HilbertBasisType basis_;
	VectorQnType qq_;
	VectorOperatorType creationMatrix_;
	FeAsJzSymmetryType feAsJzSymmetry_;
}; //class ModelFeBasedSc
} // namespace Dmrg
/*@}*/
#endif

