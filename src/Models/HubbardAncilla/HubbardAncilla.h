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
#include "../FeAsModel/HilbertSpaceFeAs.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "LinkProductHubbardAncilla.h"
#include "ProgramGlobals.h"
#include "Geometry/GeometryDca.h"
#include <cstdlib>

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
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef HilbertSpaceFeAs<HilbertState> HilbertSpaceFeAsType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef LinkProductHubbardAncilla<ModelHelperType, GeometryType> LinkProductType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef ParametersHubbardAncilla<RealType, QnType> ParametersHubbardAncillaType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;
	static SizeType const ORBITALS  = 2;

	HubbardAncilla(const SolverParamsType& solverParams,
	               InputValidatorType& io,
	               GeometryType const &geometry)
	    : ModelBaseType(solverParams, geometry, new LinkProductType(io), io),
	      modelParameters_(io),
	      geometry_(geometry)
	{}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType&,
	                                const BlockType& block,
	                                RealType) const
	{
		SizeType n=block.size();
		HilbertBasisType natBasis;
		setBasis(natBasis, block);

		for (SizeType i=0;i<n;i++) {
			VectorSparseMatrixType cm;
			findAllMatrices(cm,i,natBasis);
			addInteraction(hmatrix, cm, block[i]);

			addPotentialV(hmatrix, cm, block[i]);
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0; // FIXME for Immm SDHS
		BlockType block(1, site);
		HilbertBasisType natBasis;
		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis);

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& ll = this->createOpsLabel("l");

		this->makeTrackableOrderMatters("c");
		this->makeTrackableOrderMatters("l");

		SizeType dofs = 2*ORBITALS;
		for (SizeType i=0;i<block.size();i++) {
			for (SizeType sigma=0;sigma<dofs;sigma+=2) {
				MatrixType tmp;
				findOperatorMatrices(tmp,i,sigma,natBasis);
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
				                  -1,
				                  typename OperatorType::PairType(1,m),
				                  asign,
				                  su2related);
				c.push(myOp);
			}

			setLambdaMatrices(ll, i, natBasis);
		}
	}

private:

	//! find all states in the natural basis for a block of n sites
	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		SizeType n = block.size();
		HilbertState total = (1<<(2*ORBITALS));
		total = pow(total,n);
		basis.resize(total);
		for (HilbertState a = 0; a < total; ++a) basis[a] = a;
	}

	//! set creation matrices for sites in block
	void setLambdaMatrices(OpsLabelType& ll,
	                       SizeType i,
	                       const HilbertBasisType& natBasis) const
	{
		VectorSparseMatrixType vm;
		findAllMatrices(vm,i,natBasis);

		typename OperatorType::Su2RelatedType su2related;
		for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
			SizeType spin2 = 1 - spin1;
			SparseMatrixType lambda;
			assert(1+spin2*ORBITALS < vm.size());
			multiply(lambda,vm[0+spin1*ORBITALS],vm[1+spin2*ORBITALS]);
			MatrixType dlambda;
			crsMatrixToFullMatrix(dlambda,lambda);
			correctLambda(dlambda,spin1,vm);

			OperatorType myOp(SparseMatrixType(dlambda),
			                  1,
			                  typename OperatorType::PairType(0,0),
			                  1,
			                  su2related);
			ll.push(myOp);
		}
	}

	void correctLambda(MatrixType& dlambda,
	                   SizeType spin1,
	                   VectorSparseMatrixType& vm) const
	{
		SizeType n = dlambda.rows();
		MatrixType corrector(n,n);
		computeCorrector(corrector,spin1,vm);

		MatrixType dlambda2 = dlambda;

		dlambda = dlambda2 * corrector;
	}

	void computeCorrector(MatrixType& corrector,
	                      SizeType spin1,
	                      VectorSparseMatrixType& vm) const
	{
		SizeType spin2 = 1 - spin1;
		SparseMatrixType cm1(vm[0+spin2*ORBITALS]);
		SparseMatrixType cm2(vm[1+spin1*ORBITALS]);
		SparseMatrixType n1 = n(cm1);
		SparseMatrixType n2 = n(cm2);
		MatrixType dn1;
		MatrixType dn2;
		crsMatrixToFullMatrix(dn1,n1);
		crsMatrixToFullMatrix(dn2,n2);

		SizeType n = corrector.rows();
		SparseElementType f1 = (-1.0);
		for (SizeType i = 0; i < n; ++i)
			corrector(i,i) = std::abs(dn1(i,i) + dn2(i,i) + f1);
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to
	//! basis state ket
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	RealType sign(HilbertState const &ket, int i,SizeType sigma) const
	{
		int value=0;
		SizeType dofs = 2*ORBITALS;
		for (SizeType alpha=0;alpha<dofs;alpha++)
			value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);
		// add electron on site 0 if needed
		if (i>0) value += HilbertSpaceFeAsType::electrons(ket);

		//order for sign is: a up, a down, b up, b down, etc
		unsigned int x = HilbertSpaceFeAsType::get(ket,i);
		int spin = sigma/ORBITALS;
		SizeType orb = sigma % ORBITALS;

		for (SizeType j=0;j<orb;j++) {
			for (SizeType k=0;k<2;k++) {
				SizeType ind = j + k * ORBITALS;
				int mask = (1<<ind);
				if (x & mask) value++;
			}
		}

		if (spin==SPIN_DOWN) {
			int mask = (1<<orb);
			if (x & mask) value++;
		}

		return (value==0 || value%2==0) ? 1.0 : FERMION_SIGN;
	}

	//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void findOperatorMatrices(MatrixType& creationMatrix,
	                          int i,
	                          int sigma,
	                          const HilbertBasisType& natBasis) const
	{
		HilbertState bra,ket;
		SizeType n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<n;ii++) {
			bra=ket=natBasis[ii];

			if (HilbertSpaceFeAsType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceFeAsType::create(bra,i,sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices: error\n");
				if (ii==SizeType(jj)) {
					std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra;
					std::cerr<<" sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("Creation op. cannot be diagonal\n");
				}

				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		transposeConjugate(creationMatrix,cm);
	}

	void findAllMatrices(VectorSparseMatrixType& vm,
	                     SizeType i,
	                     const HilbertBasisType& natBasis) const
	{
		for (SizeType sigma = 0; sigma < 2*ORBITALS; ++sigma) {
			MatrixType m;
			findOperatorMatrices(m,i,sigma,natBasis);
			vm.push_back(SparseMatrixType(m));
		}
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		VectorSizeType other(4, 0);
		SizeType offset = basis.size();
		qns.resize(offset, QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = PairType(0,0);

			SizeType naUp = HilbertSpaceFeAsType::calcNofElectrons(basis[i],
			                                                       ORBITALS*SPIN_UP);
			SizeType naDown = HilbertSpaceFeAsType::calcNofElectrons(basis[i],
			                                                         ORBITALS*SPIN_DOWN);

			SizeType flavor = 0;

			// nup
			other[1] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                        SPIN_UP);
			// ntotal
			other[0] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                        SPIN_DOWN) + other[1];

			// up ancilla
			other[2] = naUp;

			// down ancilla
			other[3] = naDown;

			bool sign = other[0] & 1;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorSparseMatrixType& cm,
	                   SizeType actualIndexOfSite) const
	{
		SizeType orbital = 0;
		SparseMatrixType nup = n(cm[orbital+SPIN_UP*ORBITALS]);
		SparseMatrixType ndown = n(cm[orbital+SPIN_DOWN*ORBITALS]);

		SizeType linSize = geometry_.numberOfSites();

		SizeType iUp = actualIndexOfSite + (orbital + 0*ORBITALS)*linSize;
		assert(iUp < modelParameters_.potentialV.size());
		hmatrix += modelParameters_.potentialV[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orbital + 1*ORBITALS)*linSize;
		assert(iDown < modelParameters_.potentialV.size());
		hmatrix += modelParameters_.potentialV[iDown] * ndown;
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

		return tmpMatrix;
	}

	//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteraction(SparseMatrixType &hmatrix,
	                    const VectorSparseMatrixType& cm,
	                    SizeType actualSite) const
	{
		SparseMatrixType tmpMatrix;

		SizeType alpha = 0; // real sites, no ancilla
		SparseMatrixType m1=cm[alpha+SPIN_UP*2];
		SparseMatrixType m2=cm[alpha+SPIN_DOWN*2];

		multiply(tmpMatrix,n(m1),n(m2));
		assert(actualSite < modelParameters_.hubbardU.size());
		hmatrix += modelParameters_.hubbardU[actualSite]*tmpMatrix;
	}

	//serializr normal modelParameters_
	ParametersHubbardAncillaType  modelParameters_;
	//serializr ref geometry_ start
	const GeometryType& geometry_;
}; //class HubbardAncilla
} // namespace Dmrg
/*@}*/
#endif

