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

/*! \file HubbardHolstein.h
 *
 *  An implementation of a Hubbard Holstein model to use with the DmrgSolver
 *
 */
#ifndef DMRG_HUBBARD_HOLSTEIN_H
#define DMRG_HUBBARD_HOLSTEIN_H
#include "ModelBase.h"
#include "ParametersHubbardHolstein.h"
#include "HilbertSpaceHubbardHolstein.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "LinkProductHubbardHolstein.h"
#include "ProgramGlobals.h"
#include "Geometry/GeometryDca.h"
#include <cstdlib>

namespace Dmrg {
template<typename ModelBaseType>
class HubbardHolstein : public ModelBaseType {

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
	typedef unsigned int long WordType;
	typedef HilbertSpaceHubbardHolstein<WordType> HilbertSpaceHubbardHolsteinWordType;
	typedef HilbertSpaceHubbardHolstein<HilbertState> HilbertSpaceHubbardHolsteinType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef LinkProductHubbardHolstein<ModelHelperType, GeometryType> LinkProductType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef ParametersHubbardHolstein<RealType, QnType> ParametersHubbardHolsteinType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceHubbardHolsteinWordType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceHubbardHolsteinWordType::SPIN_DOWN;
	static SizeType const ORBITALS  = 2;

	HubbardHolstein(const SolverParamsType& solverParams,
	                InputValidatorType& io,
	                const GeometryType& geometry,
	                PsimagLite::String additional)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    new LinkProductType(additional == "SSH", io),
	                    io),
	      modelParameters_(io),
	      geometry_(geometry)
	{
		HilbertSpaceHubbardHolsteinType::setBitPhonons(modelParameters_.numberphonons);
		if (additional == "SSH") {
			PsimagLite::String warning("HubbardHolstein: ");
			warning += "SSH term in use.\n";
			std::cout<<warning;
			std::cerr<<warning;
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
		return (4*(modelParameters_.numberphonons+1));
	}

	void print(std::ostream& os) const { operator<<(os,modelParameters_); }

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         VectorQnType& qns,
	                         const BlockType& block) const
	{
		HilbertBasisType natBasis;
		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis);

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		creationMatrix.clear();
		SparseMatrixType tmpMatrix;
		for (SizeType i=0;i<block.size();i++) {
			for (int sigma=0;sigma<2;sigma++) {
				tmpMatrix = findOperatorMatrices(i,sigma,natBasis);
				int asign= 1;
				if (sigma>0) asign= 1;
				typename OperatorType::Su2RelatedType su2related;
				if (sigma==0) {
					su2related.source.push_back(i);
					su2related.source.push_back(i+1);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = 1;
				}
				OperatorType myOp(tmpMatrix,
				                  -1,
				                  typename OperatorType::PairType(1,1-sigma),
				                  asign,
				                  su2related);

				creationMatrix.push_back(myOp);
			}

			if (modelParameters_.numberphonons == 0) continue;

			tmpMatrix=findPhononadaggerMatrix(i,natBasis);

			typename OperatorType::Su2RelatedType su2related2;
			su2related2.source.push_back(i*2);
			su2related2.source.push_back(i*2+1);
			su2related2.source.push_back(i*2);
			su2related2.transpose.push_back(-1);
			su2related2.transpose.push_back(-1);
			su2related2.transpose.push_back(1);
			su2related2.offset = 1;
			OperatorType myOp2(tmpMatrix,1,PairType(2,2),-1,su2related2);
			creationMatrix.push_back(myOp2);

			if (ModelBaseType::linkProduct().terms() == 2) continue;

			// Set the operators c_(i,sigma} * x_i in the natural basis

			for (int sigma=0;sigma<2;sigma++) {
				tmpMatrix = findSSHMatrices(i,sigma,natBasis);
				int asign= 1;
				if (sigma>0) asign= 1;
				typename OperatorType::Su2RelatedType su2related3;
				if (sigma==0) {
					su2related3.source.push_back(i);
					su2related3.source.push_back(i+1);
					su2related3.transpose.push_back(-1);
					su2related3.transpose.push_back(-1);
					su2related3.offset = 1;
				}
				OperatorType myOp3(tmpMatrix,
				                   -1,
				                   typename OperatorType::PairType(1,1-sigma),
				                   asign,
				                   su2related3);

				creationMatrix.push_back(myOp3);
			}

		}
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		VectorOperatorType creationMatrix;
		VectorQnType qns;
		setOperatorMatrices(creationMatrix, qns, block);
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

		if (what=="a") { // a_dagger
			SizeType offset = 2;
			assert(offset < creationMatrix.size());
			return creationMatrix[offset+dof];
		}

		PsimagLite::String str("HubbardHolstein: naturalOperator: no label ");
		str += what + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType&,
	                                const BlockType& block,
	                                RealType,
	                                RealType factorForDiagonals=1.0) const
	{
		SizeType n=block.size();
		HilbertBasisType natBasis;
		setBasis(natBasis, block);

		for (SizeType i=0;i<n;i++) {
			VectorSparseMatrixType cm;
			findAllMatrices(cm,i,natBasis);
			addInteractionFU(hmatrix,cm,factorForDiagonals,block[i]);
			addInteractionFPhonon(hmatrix,cm,factorForDiagonals,block[i]);

			addPotentialFV(hmatrix,
			               cm,
			               block[i],
			               factorForDiagonals);

			addPotentialPhononV(hmatrix,
			                    cm,
			                    block[i],
			                    factorForDiagonals);

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

private:

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		SizeType n = block.size();
		HilbertState total = hilbertSize(0);
		total = pow(total,n);

		basis.resize(total);
		for (HilbertState a = 0; a < total; ++a) basis[a] = a;
	}

	//! Find a^+_site in the natural basis natBasis
	SparseMatrixType findPhononadaggerMatrix(SizeType site,
	                                         const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		typename HilbertSpaceHubbardHolsteinType::HilbertState bra,ket;
		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			if (HilbertSpaceHubbardHolsteinType::isNonZeroP(ket,site)) {
				SizeType nphon = SizeType(HilbertSpaceHubbardHolsteinType::getP(ket,site));
				if (modelParameters_.numberphonons<=1 || nphon >= modelParameters_.numberphonons)
					continue;
				HilbertSpaceHubbardHolsteinType::createP(bra,site);
				int jj = PsimagLite::isInVector(natBasis,bra);
				RealType x = int(HilbertSpaceHubbardHolsteinType::getP(bra,site));
				assert(x>=0);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices\n");
				cm(ii,jj) = sqrt(x);
			} else {
				HilbertSpaceHubbardHolsteinType::createP(bra,site);
				int jj = PsimagLite::isInVector(natBasis,bra);
				RealType x = 1;
				assert(x>=0);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices\n");
				cm(ii,jj) = sqrt(x);
			}
		}

		// print cm matrix
		//		std::cerr<<cm;

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceHubbardHolsteinType::HilbertState const &ket,
	              int i,
	              int sigma) const
	{
		int value=0;
		value += HilbertSpaceHubbardHolsteinType::calcNofElectrons(ket,0,i,0);
		value += HilbertSpaceHubbardHolsteinType::calcNofElectrons(ket,0,i,1);
		int tmp1 = HilbertSpaceHubbardHolsteinType::getF(ket,0) &1;
		int tmp2 = HilbertSpaceHubbardHolsteinType::getF(ket,0) &2;
		if (i>0 && tmp1>0) value++;
		if (i>0 && tmp2>0) value++;

		if (sigma==1) { // spin down
			if ((HilbertSpaceHubbardHolsteinType::getF(ket,i) &1)) value++;

		}
		if (value%2==0) return 1.0;

		return FERMION_SIGN;
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(int i,
	                                      int sigma,
	                                      const HilbertBasisType& natBasis) const
	{
		typename HilbertSpaceHubbardHolsteinType::HilbertState bra,ket;
		int n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			if (HilbertSpaceHubbardHolsteinType::isNonZeroF(ket,i,sigma)) {

			} else {
				HilbertSpaceHubbardHolsteinType::createF(bra,i,sigma);
				int jj = PsimagLite::isInVector(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices\n");
				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		// print cm matrix
		//		std::cerr<<cm;

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}


	SparseMatrixType findSSHMatrices(int site,
	                                 int sigma,
	                                 const HilbertBasisType& natBasis) const
	{
		SparseMatrixType csigma_temp = findOperatorMatrices(site,sigma,natBasis);
		SparseMatrixType a_temp = findPhononadaggerMatrix(site,natBasis);
		SparseMatrixType x_temp = displacementOp(a_temp);
		SparseMatrixType csigma_a;
		multiply(csigma_a,csigma_temp,x_temp);
		return csigma_a;
	}

	void findAllMatrices(VectorSparseMatrixType& vm,
	                     SizeType i,
	                     const HilbertBasisType& natBasis) const
	{
		for (SizeType sigma = 0; sigma < 2; ++sigma) {
			SparseMatrixType m = findOperatorMatrices(i,sigma,natBasis);
			vm.push_back(m);
		}

		if (modelParameters_.numberphonons == 0) return;
		SparseMatrixType m = findPhononadaggerMatrix(i,natBasis);
		vm.push_back(m);
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;

		qns.resize(basis.size(), ModelBaseType::QN_ZERO);
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = PairType(0,0);

			// nup
			SizeType electronsUp = HilbertSpaceHubbardHolsteinType::getNofDigits(basis[i],
			                                                                     SPIN_UP);
			// ndown
			SizeType electronsDown = HilbertSpaceHubbardHolsteinType::getNofDigits(basis[i],
			                                                                       SPIN_DOWN);

			SizeType electrons = electronsUp + electronsDown;

			qns[i] = QnType(electrons, VectorSizeType(1, electronsUp), jmpair, electrons);
		}
	}

	void addPotentialFV(SparseMatrixType &hmatrix,
	                    const VectorSparseMatrixType& cm,
	                    SizeType actualIndexOfSite,
	                    RealType factorForDiagonals) const
	{
		SparseMatrixType nup = n(cm[SPIN_UP]);
		SparseMatrixType ndown = n(cm[SPIN_DOWN]);

		SizeType linSize = geometry_.numberOfSites();
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialFV.size());
		hmatrix += factorForDiagonals*modelParameters_.potentialFV[iUp] * nup;
		SizeType iDown = actualIndexOfSite + linSize;
		assert(iDown < modelParameters_.potentialFV.size());
		hmatrix += factorForDiagonals*modelParameters_.potentialFV[iDown] * ndown;
	}

	void addPotentialPhononV(SparseMatrixType &hmatrix,
	                         const VectorSparseMatrixType& cm,
	                         SizeType actualIndexOfSite,
	                         RealType factorForDiagonals) const
	{
		if (modelParameters_.numberphonons == 0) return;
		assert(2 < cm.size());
		SparseMatrixType nphon = n(cm[2]);
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialPV.size());
		hmatrix += factorForDiagonals*modelParameters_.potentialPV[iUp] * nphon;
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

		return tmpMatrix;
	}

	SparseMatrixType displacementOp(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix = c;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		tmpMatrix+=cdagger;
		return tmpMatrix;
	}

	//! Term is U \sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteractionFU(SparseMatrixType &hmatrix,
	                      const VectorSparseMatrixType& cm,
	                      RealType factorForDiagonals,
	                      SizeType actualSite) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType m1=cm[SPIN_UP];
		SparseMatrixType m2=cm[SPIN_DOWN];

		multiply(tmpMatrix,n(m1),n(m2));
		assert(actualSite < modelParameters_.hubbardFU.size());
		hmatrix += factorForDiagonals*modelParameters_.hubbardFU[actualSite]*tmpMatrix;
	}

	//! Term is lambda\sum_{\alpha} (n_{i\alpha} -1) x_{i}
	void addInteractionFPhonon(SparseMatrixType &hmatrix,
	                           const VectorSparseMatrixType& cm,
	                           RealType factorForDiagonals,
	                           SizeType actualSite) const
	{
		if (modelParameters_.numberphonons == 0) return;
		SparseMatrixType tmpMatrix;
		SparseMatrixType m=n(cm[SPIN_UP]);
		SparseMatrixType m2=n(cm[SPIN_DOWN]);
		m+=m2;
		assert(2 < cm.size());
		SparseMatrixType x = displacementOp(cm[2]);

		multiply(tmpMatrix,m,x);
		tmpMatrix.checkValidity();
		assert(actualSite < modelParameters_.lambdaFP.size());
		hmatrix += factorForDiagonals*modelParameters_.lambdaFP[actualSite]*tmpMatrix;
	}

	ParametersHubbardHolsteinType modelParameters_;
	const GeometryType& geometry_;
}; //class HubbardHolstein
} // namespace Dmrg
/*@}*/
#endif

