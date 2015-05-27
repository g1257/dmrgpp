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

/*! \file TjAncillaC.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef DMRG_TJ_ANCILLAC_H
#define DMRG_TJ_ANCILLAC_H
#include "../Models/FeAsModel/ModelFeBasedSc.h"
#include "../Models/TjAncillaC/LinkProductTjAncillaC.h"
#include "../Models/TjAncillaC/ParametersTjAncillaC.h"
#include "ModelCommon.h"

namespace Dmrg {
//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
template<typename ModelBaseType>
class TjAncillaC : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef ModelFeBasedSc<ModelBaseType> ModelFeAsType;
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
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef LinkProductTjAncillaC<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ModelFeAsType::HilbertState HilbertStateType;
	typedef typename ModelFeAsType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelFeAsType::HilbertSpaceFeAsType HilbertSpaceType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename OperatorType::PairType PairType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename PsimagLite::Vector<HilbertStateType>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	static const int NUMBER_OF_ORBITALS = 2;
	static const int FERMION_SIGN = -1;

	enum {SPIN_UP, SPIN_DOWN};

	TjAncillaC(const SolverParamsType& solverParams,
	           InputValidatorType& io,
	           GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry)
	{
		std::cout<<"TjAncillaC: This model is EXPERIMENTAL\n";
	}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType) const
	{
		return pow(3,NUMBER_OF_ORBITALS);
	}

	//! find creation operator matrices for (i,sigma) in the natural basis,
	//! find quantum numbers and number of electrons
	//! for each state in the basis
	virtual void setNaturalBasis(VectorOperatorType& creationMatrix,
	                             SparseMatrixType &hamiltonian,
	                             SymmetryElectronsSzType& q,
	                             const BlockType& block,
	                             const RealType& time) const
	{

		HilbertBasisType natBasis;
		VectorSizeType quantumNumbs;
		setNaturalBasis(natBasis,quantumNumbs,block);

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
		VectorHilbertStateType natBasis;
		SparseMatrixType tmpMatrix;
		VectorSizeType quantumNumbs;
		setNaturalBasis(natBasis,quantumNumbs,block);

		// Set the operators c^\daggger_{i\sigma} in the natural basis
		creationMatrix.clear();
		for (SizeType i=0;i<block.size();i++) {
			for (int sigma=0;sigma<2;sigma++) {
				findOperatorMatrices(tmpMatrix,i,sigma*NUMBER_OF_ORBITALS,natBasis);
				int asign= 1;
				typename OperatorType::Su2RelatedType su2related;

				OperatorType myOp(tmpMatrix,-1,PairType(0,0),asign,su2related);

				creationMatrix.push_back(myOp);
			}

			SizeType orbital = 0;

			// Set the operators S^+_i in the natural basis
			tmpMatrix=findSplusMatrices(i,orbital,natBasis);

			typename OperatorType::Su2RelatedType su2related;

			OperatorType myOp(tmpMatrix,1,PairType(0,0),1.0,su2related);
			creationMatrix.push_back(myOp);

			// Set the operators S^z_i in the natural basis
			tmpMatrix = findSzMatrices(i,orbital,natBasis);
			typename OperatorType::Su2RelatedType su2related2;
			OperatorType myOp2(tmpMatrix,1,PairType(0,0),1.0,su2related2);
			creationMatrix.push_back(myOp2);

			// Set ni matrix:
			SparseMatrixType tmpMatrix = findNiMatrices(0,natBasis);
			RealType angularFactor= 1;
			typename OperatorType::Su2RelatedType su2related3;
			su2related3.offset = 1; //check FIXME
			OperatorType myOp3(tmpMatrix,1,PairType(0,0),angularFactor,su2related3);
			creationMatrix.push_back(myOp3);

			// Set delta_i matrix:
			tmpMatrix = findDeltaIMatrices(i,natBasis);
			typename OperatorType::Su2RelatedType su2related4;
			OperatorType myOp4(tmpMatrix,1,PairType(0,0),angularFactor,su2related4);
			creationMatrix.push_back(myOp4);

		}
	}

	/** \cppFunction{!PTEX_THISFUNCTION} returns the operator in the
	 *unmangled (natural) basis of one-site */
	MatrixType naturalOperator(const PsimagLite::String& what,
	                           SizeType site,
	                           SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatrices(creationMatrix,block);
		assert(creationMatrix.size()>0);
		SizeType nrow = creationMatrix[0].data.row();

		if (what == "i" || what=="identity") {
			MatrixType tmp(nrow,nrow);
			for (SizeType i = 0; i < tmp.n_row(); ++i) tmp(i,i) = 1.0;
			return tmp;
		}

		if (what=="+") {
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,creationMatrix[2].data);
			return tmp;
		}

		if (what=="-") {
			SparseMatrixType tmp2;
			transposeConjugate(tmp2,creationMatrix[2].data);
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,tmp2);
			return tmp;
		}

		if (what=="z") {
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,creationMatrix[3].data);
			return tmp;
		}

		if (what=="c") {
			MatrixType tmp;
			assert(dof<2);
			crsMatrixToFullMatrix(tmp,creationMatrix[dof].data);
			return tmp;
		}

		if (what=="n") {
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,creationMatrix[4].data);
			return tmp;
		}

		if (what=="d") {
			MatrixType tmp;
			crsMatrixToFullMatrix(tmp,creationMatrix[5].data);
			return tmp;
		}

		if (what=="nup") {
			MatrixType cup = naturalOperator("c",site,SPIN_UP);
			MatrixType nup = multiplyTransposeConjugate(cup,cup);
			return nup;
		}

		if (what=="ndown") {
			MatrixType cdown = naturalOperator("c",site,SPIN_DOWN);
			MatrixType ndown = multiplyTransposeConjugate(cdown,cdown);
			return ndown;
		}

		std::cerr<<"Argument: "<<what<<" "<<__FILE__<<"\n";
		throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
	}

	//! find total number of electrons for each state in the basis
	void findElectrons(typename PsimagLite::Vector<SizeType> ::Type&electrons,
	                   const VectorHilbertStateType& basis,
	                   SizeType) const
	{
		int nup,ndown;
		electrons.clear();
		for (SizeType i=0;i<basis.size();i++) {
			nup = HilbertSpaceType::getNofDigits(basis[i],0);
			ndown = HilbertSpaceType::getNofDigits(basis[i],1);
			electrons.push_back(nup+ndown);
		}
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setNaturalBasis(HilbertBasisType& basis,
	                     VectorSizeType& q,
	                     const VectorSizeType& block) const
	{
		assert(block.size()==1);
		HilbertStateType a=0;
		int sitesTimesDof=2*NUMBER_OF_ORBITALS;
		HilbertStateType total = (1<<sitesTimesDof);

		HilbertBasisType  basisTmp;
		for (a=0;a<total;a++) {
			if (!isAllowed(a)) continue;
			basisTmp.push_back(a);
		}

		assert(basisTmp.size() == pow(3,NUMBER_OF_ORBITALS));
		// reorder the natural basis (needed for MULTIPLE BANDS)
		findQuantumNumbers(q,basisTmp,block.size());
		this->orderBasis(basis,q,basisTmp);
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

	void print(std::ostream& os) const
	{
		os<<modelParameters_;
	}

private:

	bool isAllowed(HilbertStateType a) const
	{
		HilbertStateType mask1 = 1;
		mask1 <<= NUMBER_OF_ORBITALS;
		mask1--;

		HilbertStateType mask2 = mask1;
		mask2 <<= NUMBER_OF_ORBITALS;

		HilbertStateType a1 = (a & mask1);
		HilbertStateType a2 = (a & mask2);
		a2 >>= NUMBER_OF_ORBITALS;
		return ((a1 & a2) == 0);
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to
	//! basis state ket
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	RealType sign(const HilbertStateType& ket, int i,SizeType sigma) const
	{
		int value=0;
		SizeType dofs=2*NUMBER_OF_ORBITALS;
		for (SizeType alpha=0;alpha<dofs;alpha++)
			value += HilbertSpaceType::calcNofElectrons(ket,0,i,alpha);
		// add electron on site 0 if needed
		if (i>0) value += HilbertSpaceType::electrons(ket);

		//order for sign is: a up, a down, b up, b down, etc
		unsigned int x = HilbertSpaceType::get(ket,i);
		int spin = sigma/NUMBER_OF_ORBITALS;
		SizeType orb = sigma % NUMBER_OF_ORBITALS;

		for (SizeType j=0;j<orb;j++) {
			for (SizeType k=0;k<2;k++) {
				SizeType ind = j + k * NUMBER_OF_ORBITALS;
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
		int n = natBasis.size();
		PsimagLite::Matrix<SparseElementType> cm(n,n);
		findCreationDense(cm,i,sigma,natBasis);
		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,cm);
		transposeConjugate(creationMatrix,temp);

	}
	//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void findCreationDense(PsimagLite::Matrix<SparseElementType>& cm,
	                       int i,
	                       int sigma,
	                       const HilbertBasisType& natBasis) const
	{
		HilbertStateType bra,ket;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];

			if (HilbertSpaceType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceType::create(bra,i,sigma);
				int jj = PsimagLite::isInVector(natBasis,bra);
				if (jj<0) continue;
				if (ii==SizeType(jj)) {
					std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra;
					std::cerr<<" sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("Creation operator diagonal\n");
				}

				cm(ii,jj) = sign(ket,i,sigma);
			}
		}
	}

	//! Find S^+_i in the natural basis natBasis
	SparseMatrixType findSplusMatrices(int i,
	                                   SizeType orb,
	                                   const VectorHilbertStateType& natBasis) const
	{
		HilbertStateType bra,ket;
		int n = natBasis.size();
		MatrixType cm(n,n);
		SizeType sigma1 = orb + SPIN_DOWN*NUMBER_OF_ORBITALS;
		SizeType sigma2 = orb + SPIN_UP*NUMBER_OF_ORBITALS;
		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			if (HilbertSpaceType::isNonZero(ket,i,sigma1) &&
			        !HilbertSpaceType::isNonZero(ket,i,sigma2)) {
				// it is a down electron, then flip it:
				HilbertSpaceType::destroy(bra,i,sigma1);
				HilbertSpaceType::create(bra,i,sigma2);
				int jj = PsimagLite::isInVector(natBasis,bra);
				assert(jj>=0);
				cm(ii,jj)=1.0;
			}
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(int i,
	                                SizeType orb,
	                                const VectorHilbertStateType& natBasis) const
	{
		HilbertStateType ket;
		int n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			ket=natBasis[ii];
			RealType value = 0.0;
			if (HilbertSpaceType::isNonZero(ket,i,orb + SPIN_UP*NUMBER_OF_ORBITALS))
				value += 1.0;
			if (HilbertSpaceType::isNonZero(ket,i,orb + SPIN_DOWN*NUMBER_OF_ORBITALS))
				value -= 1.0;

			cm(ii,ii)=0.5*value;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findNiMatrices(int i,
	                                const VectorHilbertStateType& natBasis) const
	{
		SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			cm(ii,ii) = 0.0;
			for (SizeType sigma=0;sigma<2;sigma++)
				if (HilbertSpaceType::isNonZero(ket,i,sigma*NUMBER_OF_ORBITALS))
					cm(ii,ii) += 1.0;
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	//! Find Delta_i in the natural basis natBasis
	SparseMatrixType findDeltaIMatrices(int i,
	                                    const VectorHilbertStateType& natBasis) const
	{
		assert(NUMBER_OF_ORBITALS == 2);
		int n = natBasis.size();

		PsimagLite::Matrix<SparseElementType> cru(n,n);
		findCreationDense(cru,i,SPIN_UP*NUMBER_OF_ORBITALS,natBasis);

		PsimagLite::Matrix<SparseElementType> cad(n,n);
		findCreationDense(cad,i,1 + SPIN_DOWN*NUMBER_OF_ORBITALS,natBasis);

		PsimagLite::Matrix<SparseElementType> crd(n,n);
		findCreationDense(crd,i,SPIN_DOWN*NUMBER_OF_ORBITALS,natBasis);

		PsimagLite::Matrix<SparseElementType> cau(n,n);
		findCreationDense(cau,i,1 + SPIN_UP*NUMBER_OF_ORBITALS,natBasis);

		PsimagLite::Matrix<SparseElementType> part1 = cru*cad;
		PsimagLite::Matrix<SparseElementType> part2 = (-1.0)*crd*cau;

		RealType oneOverSqrt2 = 1.0/sqrt(2.0);
		PsimagLite::Matrix<SparseElementType> final = oneOverSqrt2*(part1 + part2);
		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,final);
		return temp;
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType&,
	                                const BlockType& block,
	                                RealType,
	                                RealType factorForDiagonals=1.0) const
	{
		SizeType n=block.size();

		SizeType linSize = geometry_.numberOfSites();
		for (SizeType i=0;i<n;i++) {
			SizeType orb = 0;
			// potentialV
			SparseMatrixType nup(naturalOperator("nup",i,orb));
			SparseMatrixType ndown(naturalOperator("ndown",i,orb));
			SparseMatrixType m = nup;
			SizeType index = block[i]+orb*linSize;
			assert(index<modelParameters_.potentialV.size());
			m *= modelParameters_.potentialV[block[i] + orb*linSize];
			m += modelParameters_.potentialV[index]*ndown;
			hmatrix += factorForDiagonals * m;

		}
	}

	void findQuantumNumbers(VectorSizeType& q,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		SymmetryElectronsSzType qq;
		setSymmetryRelated(qq,basis,n);
		MyBasis::findQuantumNumbers(q,qq);
	}

	void setSymmetryRelated(SymmetryElectronsSzType& q,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		assert(n==1);

		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		typename PsimagLite::Vector<PairType>::Type jmvalues;
		VectorSizeType flavors;
		PairType jmSaved = calcJmvalue<PairType>(basis[0]);
		jmSaved.first++;
		jmSaved.second++;

		VectorSizeType ups(basis.size());
		VectorSizeType downs(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair = calcJmvalue<PairType>(basis[i]);

			jmvalues.push_back(jmpair);
			// nup
			ups[i] = HilbertSpaceType::electronsWithGivenSpin(basis[i],SPIN_UP);
			// ndown
			downs[i] = HilbertSpaceType::electronsWithGivenSpin(basis[i],SPIN_DOWN);

			flavors.push_back(ups[i]+downs[i]);
			jmSaved = jmpair;
		}

		q.set(jmvalues,flavors,ups,downs);
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// Reinterprets 6 and 9
	template<typename PairType>
	PairType calcJmvalue(const HilbertStateType&) const
	{
		return PairType(0,0);
	}

	//serializr start class TjAncillaC
	//serializr vptr
	//serializr normal modelParameters_
	ParametersTjAncillaC<RealType>  modelParameters_;
	//serializr ref geometry_ end
	const GeometryType &geometry_;
};	//class TjAncillaC

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_ANCILLAC_H

