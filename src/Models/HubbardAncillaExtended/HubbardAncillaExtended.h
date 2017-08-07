/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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
#include "../FeAsModel/HilbertSpaceFeAs.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "LinkProductHubbardAncillaExtended.h"
#include "ProgramGlobals.h"
#include "ModelCommon.h"
#include "Geometry/GeometryDca.h"
#include <cstdlib>

namespace Dmrg {
template<typename ModelBaseType>
class HubbardAncillaExtended : public ModelBaseType {

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef TargetQuantumElectrons<RealType> TargetQuantumElectronsType;
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef  HilbertSpaceFeAs<HilbertState> HilbertSpaceFeAsType;


	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef LinkProductHubbardAncillaExtended<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;

	typedef typename ModelBaseType::VectorType VectorType;
	typedef	 typename ModelBaseType::MyBasis MyBasis;
	typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef ParametersHubbardAncilla<RealType> ParametersHubbardAncillaType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;
	static SizeType const ORBITALS  = 2;

	HubbardAncillaExtended(const SolverParamsType& solverParams,
	                       InputValidatorType& io,
	                       GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry),
	      hot_(geometry_.orbitals(0,0) > 1)
	{
		LinkProductType::setHot(hot_);
		if (hot_) {
			PsimagLite::String msg("HubbardAncillaExtended: Hot ancilla mode is on");
			msg += " (EXPERIMENTAL feature)\n";
			std::cout<<msg;
			std::cerr<<msg;
		}
	}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType) const
	{
		return (1<<(2*ORBITALS));
	}

	void print(std::ostream& os) const { operator<<(os,modelParameters_); }

	//! find creation operator matrices for (i,sigma) in the natural basis,
	//! find quantum numbers and number of electrons
	//! for each state in the basis
	void setNaturalBasis(VectorOperatorType& creationMatrix,
	                     SparseMatrixType &hamiltonian,
	                     SymmetryElectronsSzType &q,
	                     const BlockType& block,
	                     const RealType& time)  const
	{
		HilbertBasisType natBasis;
		VectorSizeType qvector;
		setNaturalBasis(natBasis,qvector,block);

		setOperatorMatrices(creationMatrix,block);

		//! Set symmetry related
		setSymmetryRelated(q,natBasis,block.size());

		//! set hamiltonian
		this->calcHamiltonian(hamiltonian,creationMatrix,block,time);
	}

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         const BlockType& block) const
	{
		HilbertBasisType natBasis;
		VectorSizeType qvector;
		setNaturalBasis(natBasis,qvector,block);

		//! Set the operators c^\dagger_{i\gamma\sigma} in the natural basis
		creationMatrix.clear();
		SizeType dofs = 2*ORBITALS;
		for (SizeType i=0;i<block.size();i++) {
			for (SizeType sigma2=0;sigma2<dofs;++sigma2) {
				SizeType sigma = (hot_) ? 2*sigma2 : sigma2;
				if (sigma >= dofs) sigma -= (dofs - 1);
				if (!hot_ && (sigma & 1)) continue;
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
				creationMatrix.push_back(myOp);
			}

			setLambdaMatrices(creationMatrix,i,natBasis);
			setSplus(creationMatrix,i,0,natBasis);
			if (hot_) setSplus(creationMatrix,i,1,natBasis);
			setSz(creationMatrix,i,0,natBasis);
			if (hot_) setSz(creationMatrix,i,1,natBasis);
			setPair(creationMatrix,i,0,natBasis);
			if (hot_) setPair(creationMatrix,i,1,natBasis);
			setN(creationMatrix,i,0,natBasis);
			if (hot_) setN(creationMatrix,i,1,natBasis);
		}

		if (hot_)
			assert(creationMatrix.size() == 14);
		else
			assert(creationMatrix.size() == 8);
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		VectorOperatorType creationMatrix;
		setOperatorMatrices(creationMatrix,block);
		assert(creationMatrix.size()>0);
		SizeType nrow = creationMatrix[0].data.rows();
		PsimagLite::String what2 = what;

		if (what2 == "i" || what2=="identity") {
			SparseMatrixType tmp(nrow,nrow);
			tmp.makeDiagonal(nrow,1.0);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what2 == "0") {
			SparseMatrixType tmp(nrow,nrow);
			tmp.makeDiagonal(nrow,0.0);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what == "c") {
			if (dof >= creationMatrix.size()) {
				PsimagLite::String str("naturalOperator: dof too big ");
				str += "maximum is " + ttos(creationMatrix.size());
				str += " given is " + ttos(dof) + "\n";
				throw PsimagLite::RuntimeError(str);
			}
			return creationMatrix[dof];
		}

		if (what=="d") { // \Delta
			SizeType offset = (!hot_) ? 2 : 4;
			assert(offset < creationMatrix.size());
			return creationMatrix[offset+dof];
		}

		if (what == "splus") { // S^+
			SizeType offset = (!hot_) ? 4 : 6;
			return creationMatrix[offset+dof];
		}

		if (what == "sminus") { // S^-
			SizeType offset = (!hot_) ? 4 : 6;
			creationMatrix[offset+dof].conjugate();
			return creationMatrix[offset+dof];
		}

		if (what == "z" || what == "sz") { // S^z
			SizeType offset = (!hot_) ? 5 : 8;
			return creationMatrix[offset+dof];
		}

		if (what=="p") { // Pair
			SizeType offset = (!hot_) ? 6 : 10;
			return creationMatrix[offset+dof];
		}

		if (what=="n") { // N
			SizeType offset = (!hot_) ? 7 : 12;
			return creationMatrix[offset+dof];
		}

		PsimagLite::String str("HubbardAncillaExtended: naturalOperator: no label ");
		str += what + "\n";
		throw PsimagLite::RuntimeError(str);
	}


	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setNaturalBasis(HilbertBasisType& basis,
	                     VectorSizeType& q,
	                     const VectorSizeType& block) const
	{
		SizeType n = block.size();
		HilbertState total = hilbertSize(0);
		total = pow(total,n);

		HilbertBasisType basisTmp;
		for (HilbertState a=0;a<total;a++) basisTmp.push_back(a);

		// reorder the natural basis (needed for MULTIPLE BANDS)
		findQuantumNumbers(q,basisTmp,block.size());
		this->orderBasis(basis,q,basisTmp);
	}

	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType) const
	{
		electrons.resize(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			// nup
			SizeType nup = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                            SPIN_UP);
			// ndown
			SizeType ndown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                              SPIN_DOWN);
			electrons[i] = nup + ndown;
		}
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType&,
	                                const BlockType& block,
	                                RealType,
	                                RealType factorForDiagonals=1.0) const
	{
		SizeType n=block.size();
		HilbertBasisType natBasis;
		VectorSizeType qvector;
		setNaturalBasis(natBasis,qvector,block);
		for (SizeType i=0;i<n;i++) {
			VectorSparseMatrixType cm;
			findAllMatrices(cm,i,natBasis);
			addInteraction(hmatrix,cm,factorForDiagonals,block[i]);

			addPotentialV(hmatrix,
			              cm,
			              block[i],
			              factorForDiagonals);

		}
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

private:

	//! set creation matrices for sites in block
	void setLambdaMatrices(VectorOperatorType& cm,
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
			cm.push_back(myOp);
		}
	}

	void correctLambda(MatrixType& dlambda,
	                   SizeType spin1,
	                   VectorSparseMatrixType& vm) const
	{
		SizeType n = dlambda.n_row();
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

		SizeType n = corrector.n_row();
		SparseElementType f1 = (-1.0);
		for (SizeType i = 0; i < n; ++i)
			corrector(i,i) = std::abs(dn1(i,i) + dn2(i,i) + f1);
	}

	void setSplus(VectorOperatorType& cm,
	              SizeType i,
	              SizeType orbital,
	              const HilbertBasisType& natBasis) const
	{

		VectorSparseMatrixType vm;
		findAllMatrices(vm,i,natBasis);

		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType m;
		transposeConjugate(m,vm[2+orbital]);
		SparseMatrixType splus;
		multiply(splus,vm[0+orbital],m);

		OperatorType myOp(splus,
		                  1,
		                  typename OperatorType::PairType(0,0),
		                  -1,
		                  su2related);
		cm.push_back(myOp);
	}

	void setSz(VectorOperatorType& cm,
	           SizeType i,
	           SizeType orbital,
	           const HilbertBasisType& natBasis) const
	{

		VectorSparseMatrixType vm;
		findAllMatrices(vm,i,natBasis);

		typename OperatorType::Su2RelatedType su2related;

		SparseMatrixType cm1(vm[0+orbital]);
		SparseMatrixType cm2(vm[2+orbital]);
		SparseMatrixType n1 = n(cm1);
		SparseMatrixType n2 = n(cm2);
		MatrixType dn1;
		MatrixType dn2;
		crsMatrixToFullMatrix(dn1,n1);
		crsMatrixToFullMatrix(dn2,n2);

		SizeType n = dn1.n_row();
		MatrixType szmatrix(n,n);


		for (SizeType i = 0; i < n; ++i) {
			szmatrix(i,i) = 0.5*(dn1(i,i) - dn2(i,i));
		}

		OperatorType sz(SparseMatrixType(szmatrix),
		                1,
		                typename OperatorType::PairType(0,0),
		                1,
		                su2related);

		cm.push_back(sz);
	}

	void setPair(VectorOperatorType& cm,
	             SizeType i,
	             SizeType orbital,
	             const HilbertBasisType& natBasis) const
	{
		VectorSparseMatrixType vm;
		findAllMatrices(vm,i,natBasis);
		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType pair;
		multiply(pair,vm[0+orbital],vm[2+orbital]);

		OperatorType myOp(pair,
		                  1,
		                  typename OperatorType::PairType(0,0),
		                  1,
		                  su2related);
		cm.push_back(myOp);
	}

	void setN(VectorOperatorType& cm,
	          SizeType i,
	          SizeType orbital,
	          const HilbertBasisType& natBasis) const
	{

		VectorSparseMatrixType vm;
		findAllMatrices(vm,i,natBasis);

		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType cm1(vm[0+orbital]);
		SparseMatrixType cm2(vm[2+orbital]);
		SparseMatrixType n1 = n(cm1);
		SparseMatrixType n2 = n(cm2);
		MatrixType dn1;
		MatrixType dn2;
		crsMatrixToFullMatrix(dn1,n1);
		crsMatrixToFullMatrix(dn2,n2);
		SizeType n = dn1.n_row();
		MatrixType nmatrix(n,n);


		for (SizeType i = 0; i < n; ++i) {
			nmatrix(i,i) = 1.0*(dn1(i,i) + dn2(i,i));
		}

		OperatorType nop(SparseMatrixType(nmatrix),
		                 1,
		                 typename OperatorType::PairType(0,0),
		                 1,
		                 su2related);

		cm.push_back(nop);
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
				int jj = PsimagLite::isInVector(natBasis,bra);
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

	void findQuantumNumbers(VectorSizeType& q,const HilbertBasisType&basis,int n) const
	{
		SymmetryElectronsSzType qq;
		setSymmetryRelated(qq,basis,n);
		qq.findQuantumNumbers(q, MyBasis::useSu2Symmetry());
	}

	void setSymmetryRelated(SymmetryElectronsSzType& q,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		VectorSizeType flavors;
		VectorSizeType electrons(basis.size());
		VectorPairType jmvalues;
		VectorSizeType other(3*basis.size());
		SizeType offset = basis.size();
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair = PairType(0,0);

			jmvalues.push_back(jmpair);

			SizeType naUp = HilbertSpaceFeAsType::calcNofElectrons(basis[i],
			                                                       ORBITALS*SPIN_UP);
			SizeType naDown = HilbertSpaceFeAsType::calcNofElectrons(basis[i],
			                                                         ORBITALS*SPIN_DOWN);

			SizeType flavor = 0;

			flavors.push_back(flavor);

			// nup
			other[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                        SPIN_UP);
			// ntotal
			electrons[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                            SPIN_DOWN) +
			        other[i];

			// up ancilla
			other[i+offset] = naUp;

			// down ancilla
			other[i+2*offset] = naDown;
		}

		q.set(jmvalues,flavors,electrons,other);
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorSparseMatrixType& cm,
	                   SizeType actualIndexOfSite,
	                   RealType factorForDiagonals) const
	{
		SizeType orbital = 0;
		SparseMatrixType nup = n(cm[orbital+SPIN_UP*ORBITALS]);
		SparseMatrixType ndown = n(cm[orbital+SPIN_DOWN*ORBITALS]);

		SizeType linSize = geometry_.numberOfSites();

		SizeType iUp = actualIndexOfSite + (orbital + 0*ORBITALS)*linSize;
		assert(iUp < modelParameters_.potentialV.size());
		hmatrix += factorForDiagonals*modelParameters_.potentialV[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orbital + 1*ORBITALS)*linSize;
		assert(iDown < modelParameters_.potentialV.size());
		hmatrix += factorForDiagonals*modelParameters_.potentialV[iDown] * ndown;
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
	                    RealType factorForDiagonals,
	                    SizeType actualSite) const
	{
		SparseMatrixType tmpMatrix;
		SizeType nsites = geometry_.numberOfSites();
		SizeType factor = (hot_) ? 2 : 1;
		if (modelParameters_.hubbardU.size() != factor*nsites) {
			throw PsimagLite::RuntimeError("Number of Us is incorrect\n");
		}

		for (SizeType alpha=0;alpha<factor;++alpha) {// real sites and ancillas
			SparseMatrixType m1=cm[alpha+SPIN_UP*2];
			SparseMatrixType m2=cm[alpha+SPIN_DOWN*2];

			multiply(tmpMatrix,n(m1),n(m2));
			assert(actualSite + nsites*alpha < modelParameters_.hubbardU.size());
			hmatrix += factorForDiagonals*
			        modelParameters_.hubbardU[actualSite+nsites*alpha]*tmpMatrix;
		}
	}

	//serializr normal modelParameters_
	ParametersHubbardAncillaType  modelParameters_;
	//serializr ref geometry_ start
	const GeometryType& geometry_;
	bool hot_;
}; //class HubbardAncilla
} // namespace Dmrg
/*@}*/
#endif

