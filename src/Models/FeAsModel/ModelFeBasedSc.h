/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
#include "LinkProductFeAs.h"
#include "ProgramGlobals.h"
#include "ModelCommon.h"
#include "Geometry/GeometryDca.h"

namespace Dmrg {
template<typename ModelBaseType>
class ModelFeBasedSc : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
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
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef  HilbertSpaceFeAs<HilbertState> HilbertSpaceFeAsType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef LinkProductFeAs<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef	 typename ModelBaseType::MyBasis MyBasis;
	typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef ParametersModelFeAs<RealType> ParamsModelFeAsType;

	static const int FERMION_SIGN = -1;
	static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
	static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

	ModelFeBasedSc(const SolverParamsType& solverParams,
	               InputValidatorType& io,
	               GeometryType const &geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      reinterpretX_(6),
	      reinterpretY_(9),
	      modelParameters_(io),
	      geometry_(geometry),
	      geometryDca_(geometry,modelParameters_.orbitals),
	      spinSquared_(spinSquaredHelper_,
	                   modelParameters_.orbitals,
	                   2*modelParameters_.orbitals),
	      reinterpret_(true)
	{
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
			throw PsimagLite::RuntimeError(str.c_str());
		}

		LinkProductType::setOrbitals(modelParameters_.orbitals);
		HilbertSpaceFeAsType::setOrbitals(modelParameters_.orbitals);
		statesPerSite_ = (1 << (modelParameters_.orbitals*2));

		switch (modelParameters_.minElectronsPerSite) {
		case 0:
			break;
		case 1:
			statesPerSite_--;
			break;
		case 2:
			statesPerSite_ -= 11;
			break;
		case 3:
			statesPerSite_ -= 56;
			break;
		case 4:
			statesPerSite_ -= 176;
			break;
		case 5:
			statesPerSite_ -= 386;
			break;
		case 6:
			statesPerSite_ -= 638;
			break;
		case 7:
			statesPerSite_ -= 848;
			break;
		case 8:
			statesPerSite_ -= 968;
			break;
		case 9:
			statesPerSite_ -= 1013;
			break;
		default:
			throw PsimagLite::RuntimeError("ModelFeBasedSc: invalid param minElectronsPerSite");
		}
	}

	SizeType memResolv(PsimagLite::MemResolv& mres,
	                   SizeType,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "ModelFeBasedSc";

		const char* start = reinterpret_cast<const char *>(this);
		const char* end = reinterpret_cast<const char *>(&reinterpretX_);
		SizeType total = end - start;
		mres.push(PsimagLite::MemResolv::MEMORY_TEXTPTR,
		          total,
		          start,
		          msg + " ModelFeBasedSc vptr");

		start = end;
		end = reinterpret_cast<const char *>(&reinterpretY_);
		total += mres.memResolv(&reinterpretX_, end-start, str + " reinterpretX");

		start = end;
		end = reinterpret_cast<const char *>(&modelParameters_);
		total += mres.memResolv(&reinterpretY_, end-start, str + " reinterpretY");

		start = end;
		end = start + sizeof(modelParameters_);
		total += mres.memResolv(&modelParameters_, end-start, str + " modelParameters");

		start = end;
		end = reinterpret_cast<const char *>(&spinSquaredHelper_);
		mres.push(PsimagLite::MemResolv::MEMORY_HEAPPTR,
		          PsimagLite::MemResolv::SIZEOF_HEAPREF,
		          start,
		          str + " ref to geometry");
		total += (end - start);
		mres.memResolv(&geometry_, 0, str + " geometry");

		start = end;
		end = reinterpret_cast<const char *>(&spinSquared_);
		total += mres.memResolv(&spinSquaredHelper_, end-start, str + " spinSquaredHelper");

		start = end;
		end = reinterpret_cast<const char *>(&statesPerSite_);
		total += mres.memResolv(&spinSquared_, end-start, str + " spinSquared");

		assert(sizeof(*this) > total);
		total += mres.memResolv(&statesPerSite_,
		                        sizeof(*this) - total,
		                        str + " statesPerSite");

		return total;
	}

	SizeType hilbertSize(SizeType) const
	{
		return statesPerSite_;
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

		SparseMatrixType tmpMatrix2;
		tmpMatrix2.makeDiagonal(natBasis.size(),0.0);
	}

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         const BlockType& block) const
	{
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;
		VectorSizeType qvector;
		setNaturalBasis(natBasis,qvector,block);

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

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		VectorOperatorType creationMatrix;
		setOperatorMatrices(creationMatrix,block);
		SizeType orbital = dof % modelParameters_.orbitals;
		SizeType spin = dof/modelParameters_.orbitals;
		assert(creationMatrix.size()>0);
		SizeType nrow = creationMatrix[0].data.row();
		PsimagLite::String what2 = what;

		if (what2 == "i" || what2=="identity") {
			VectorSizeType allowed(1,0);
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
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
			VectorSizeType allowed(1,0);
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			SparseMatrixType tmp(nrow,nrow);
			tmp.makeDiagonal(nrow,0.0);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what2=="+") {
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < modelParameters_.orbitals; ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp(nrow,nrow);
			tmp += multiplyTc(creationMatrix[dof].data,
			                  creationMatrix[dof+modelParameters_.orbitals].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp2,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}
		if (what2=="-") {
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < modelParameters_.orbitals; ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp(nrow,nrow);
			tmp += multiplyTc(creationMatrix[dof+modelParameters_.orbitals].data,
			        creationMatrix[dof].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp2,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}
		if (what2=="z") {
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < modelParameters_.orbitals; ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp(nrow,nrow);
			MatrixType tmp2(nrow,nrow);

			tmp += multiplyTc(creationMatrix[dof].data,creationMatrix[dof].data);
			tmp2 += multiplyTc(creationMatrix[dof+modelParameters_.orbitals].data,
			        creationMatrix[dof+modelParameters_.orbitals].data);

			tmp = 0.5*(tmp-tmp2);
			SparseMatrixType tmp3(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp3,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}
		if (what2=="n") {
			VectorSizeType allowed(2*modelParameters_.orbitals,0);
			for (SizeType x = 0; x < allowed.size(); ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			MatrixType tmp =
			        multiplyTc(creationMatrix[dof].data,creationMatrix[dof].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(tmp2,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what2=="c") {
			VectorSizeType allowed(2*modelParameters_.orbitals,0);
			for (SizeType x = 0; x < allowed.size(); ++x) allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			creationMatrix[orbital + spin*modelParameters_.orbitals].conjugate();
			return creationMatrix[orbital + spin*modelParameters_.orbitals];
		}

		if (what2=="c\'") {
			VectorSizeType allowed(2*modelParameters_.orbitals,0);
			for (SizeType x = 0; x < allowed.size(); ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			return creationMatrix[orbital + spin*modelParameters_.orbitals];
		}

		if (what2=="d") { // delta = c^\dagger * c^dagger
			VectorSizeType allowed(modelParameters_.orbitals,0);
			for (SizeType x = 0; x < allowed.size(); ++x)
				allowed[x] = x;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			SparseMatrixType atmp;
			multiply(atmp,
			         creationMatrix[orbital+orbital+modelParameters_.orbitals].data,
			        creationMatrix[orbital].data);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(atmp,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		PsimagLite::String str("ModelFeBasedSc: naturalOperator: no label ");
		str += what + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setNaturalBasis(HilbertBasisType& basis,
	                     VectorSizeType& q,
	                     const VectorSizeType& block) const
	{
		int sitesTimesDof=2*modelParameters_.orbitals*block.size();
		HilbertState total = (1<<sitesTimesDof);

		HilbertBasisType basisTmp;
		for (HilbertState a=0;a<total;a++) {
			if (!isAllowedThisDof(a)) continue;
			basisTmp.push_back(a);
		}

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
	                                const VectorOperatorType& cm,
	                                const BlockType& block,
	                                RealType time,
	                                RealType factorForDiagonals=1.0) const
	{
		SizeType n=block.size();

		for (SizeType i=0;i<n;i++) {
			addInteraction(hmatrix,cm,i,factorForDiagonals,block[i]);
			addMagneticField(hmatrix,cm,i,block[i],factorForDiagonals);
			addSpinOrbit(hmatrix,cm,i,block[i],factorForDiagonals);

			if (modelParameters_.potentialT.size()==0 || time==0) {
				addPotentialV(hmatrix,
				              cm,
				              i,
				              block[i],
				              factorForDiagonals,
				              modelParameters_.potentialV);
			} else {
				addPotentialV(hmatrix,
				              cm,
				              i,
				              block[i],
				              factorForDiagonals,
				              modelParameters_.potentialT);
			}
		}
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
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

		//order for sign is: a up, a down, b up, b down, etc
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
				int jj = PsimagLite::isInVector(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices: internal error\n");
				if (ii==SizeType(jj)) {
					std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra<<" sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("Creation operator cannot be diagonal\n");
				}
				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		if (reinterpret_ && modelParameters_.orbitals==2) reinterpret(cm,natBasis);

		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,cm);
		transposeConjugate(creationMatrix,temp);
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
		typedef std::pair<SizeType,SizeType> PairType;
		typename PsimagLite::Vector<PairType>::Type jmvalues;
		VectorSizeType flavors;
		PairType jmSaved = calcJmvalue<PairType>(basis[0]);
		jmSaved.first++;
		jmSaved.second++;

		VectorSizeType electronsUp(basis.size());
		VectorSizeType electrons(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair(0,0);
			if (n == 1) jmpair = calcJmvalue<PairType>(basis[i]);

			jmvalues.push_back(jmpair);

			SizeType na = HilbertSpaceFeAsType::calcNofElectrons(basis[i],0) +
			        HilbertSpaceFeAsType::calcNofElectrons(basis[i],0+2);
			SizeType nb = HilbertSpaceFeAsType::calcNofElectrons(basis[i],1) +
			        HilbertSpaceFeAsType::calcNofElectrons(basis[i],1+2);

			SizeType flavor = na  + 3*nb;

			flavors.push_back(flavor);
			jmSaved = jmpair;

			// nup
			electronsUp[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                              SPIN_UP);
			// ndown
			SizeType electronsDown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],
			                                                                      SPIN_DOWN);
			electrons[i] = electronsDown + electronsUp[i];
			if (modelParameters_.spinOrbit.n_row() > 0)
				electronsUp[i] = 0;
		}

		q.set(jmvalues,flavors,electrons,electronsUp);
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
		int n  = cm.n_row();
		if (n!=16)
			throw PsimagLite::RuntimeError("blocks.size must be 1, and basis.size 16\n");

		MatrixType cmCopy(n,n);
		int i,j;
		int x=PsimagLite::isInVector(basis,reinterpretX_);
		int y=PsimagLite::isInVector(basis,reinterpretY_);

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
		RealType sf= -1;
		cmCopy(x,x)=0.5*(cm(x,x)+cm(x,y)+cm(y,x)+cm(y,y));
		cmCopy(x,y)=0.5*(cm(x,x)-sf*cm(x,y)+sf*cm(y,x)-cm(y,y));
		cmCopy(y,x)=0.5*(cm(x,x)+sf*cm(x,y)-sf*cm(y,x)-cm(y,y));
		cmCopy(y,y)=0.5*(cm(x,x)-cm(x,y)-cm(y,x)+cm(y,y));

		cm = cmCopy;
	}

	void addInteraction(SparseMatrixType &hmatrix,
	                    const VectorOperatorType& cm,
	                    SizeType i,
	                    RealType factorForDiagonals,
	                    SizeType actualSite) const
	{
		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_ORBITAL0) {
			return addInteractionAncilla(hmatrix,cm,i,factorForDiagonals,actualSite);
		}

		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_KSPACE) {
			return addInteractionKspace(hmatrix,cm,i,factorForDiagonals,actualSite);
		}

		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_IMPURITY) {
			return addInteractionImpurity(hmatrix,cm,i,factorForDiagonals,actualSite);
		}

		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_CODE2) {
			return addInteractionUmatrix(hmatrix,cm,i,factorForDiagonals);
		}

		addInteractionU1(hmatrix,cm,i,factorForDiagonals);
		addInteractionU2(hmatrix,cm,i,factorForDiagonals);
		if (modelParameters_.feAsMode == ParamsModelFeAsType::INT_PAPER33) {
			addInteractionJ1(hmatrix,cm,i,factorForDiagonals);
			addInteractionJ2(hmatrix,cm,i,factorForDiagonals);
		} else {
			addInteractionV(hmatrix,cm,i,factorForDiagonals);
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
	                      SizeType i,
	                      RealType factorForDiagonals) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType tmpMatrix;

		for (SizeType alpha=0;alpha<SizeType(modelParameters_.orbitals);alpha++) {
			SparseMatrixType m1=cm[alpha+SPIN_UP*modelParameters_.orbitals+i*dof].data;
			SparseMatrixType m2=cm[alpha+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;

			multiply(tmpMatrix,n(m1),n(m2));
			hmatrix += factorForDiagonals*findHubbardU(0,alpha,alpha)*tmpMatrix;
		}
	}

	//! Term is U[1] n_{i BAND0 } n_{i BAND1}
	void addInteractionU2(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      RealType factorForDiagonals) const
	{

		for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
			for (SizeType orb2=orb1+1;orb2<modelParameters_.orbitals;orb2++) {
				SparseMatrixType tmpMatrix;
				multiply(tmpMatrix, nSummedOverSpin(cm,i,orb1),nSummedOverSpin(cm,i,orb2));
				hmatrix += factorForDiagonals*findHubbardU(1,orb1,orb2)*tmpMatrix;
			}
		}
	}

	//! Term is U[2] S_{i BAND0 } S_{i BAND1}
	void addInteractionJ1(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      RealType factorForDiagonals) const
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
				hmatrix += factorForDiagonals*val*tmpMatrix;

				multiply(tmpMatrix,
				         spinOperator(cm,i,orb1,1),
				         spinOperator(cm,i,orb2,0));
				val = modelParameters_.hubbardU[2]/val2;
				// this is -2*J
				hmatrix += factorForDiagonals*val*tmpMatrix;

				multiply(tmpMatrix,
				         spinOperator(cm,i,orb1,2),
				         spinOperator(cm,i,orb2,2));
				val = modelParameters_.hubbardU[4]/val3;
				// this is -2*J
				hmatrix += factorForDiagonals*val*tmpMatrix;
			}
		}
	}

	//! Term is U[3] \sum_{\alpha}\bar{n}_{i\alpha UP} \bar{n}_{i\alpha DOWN}
	//! where \bar{n}_{i\alpha \spin} = c^\dagger_{i\alpha\spin} c_{i\bar{\alpha}\bar{spin}}
	void addInteractionJ2(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      RealType factorForDiagonals) const
	{
		SparseMatrixType tmpMatrix;

		for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
			for (SizeType orb2=0;orb2<modelParameters_.orbitals;orb2++) {
				if (orb1==orb2) continue;
				multiply(tmpMatrix,
				         nBar(cm,i,orb1,orb2,SPIN_UP),
				         nBar(cm,i,orb1,orb2,SPIN_DOWN));
				// -J
				hmatrix += factorForDiagonals*modelParameters_.hubbardU[3]*tmpMatrix;
			}
		}
	}

	void addInteractionV(SparseMatrixType& hmatrix,
	                     const VectorOperatorType& cm,
	                     SizeType site,
	                     RealType factorForDiagonals) const
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
				tmpMatrix1 = factorForDiagonals*value*tmpMatrix;
				hmatrix += tmpMatrix1;

				transposeConjugate(tmpMatrix2,tmpMatrix1);
				hmatrix += tmpMatrix2;
			}
		}
	}

	void addSpinOrbit(SparseMatrixType &hmatrix,
	                  const VectorOperatorType& cm,
	                  SizeType i,
	                  SizeType ,
	                  RealType factorForDiagonals) const
	{
		if (modelParameters_.spinOrbit.n_row()<4) return;
		SizeType orbitals = modelParameters_.orbitals;
		int dof=2*orbitals;

		for (SizeType spin1=0;spin1<2;spin1++) {
			for (SizeType spin2=0;spin2<2;spin2++) {
				for (SizeType orb1=0;orb1<orbitals;orb1++) {
					for (SizeType orb2=0;orb2<orbitals;orb2++) {
						SparseMatrixType c1 = cm[orb1+spin1*orbitals+i*dof].data;
						SparseMatrixType c1Transpose;
						transposeConjugate(c1Transpose,c1);
						SparseMatrixType c2 = cm[orb2+spin2*orbitals+i*dof].data;
						SparseMatrixType A = c1Transpose * c2;

						hmatrix += factorForDiagonals *
						        modelParameters_.spinOrbit(spin1+spin2*2,orb1+orb2*orbitals)*A;

					}
				}
			}
		}
	}

	void addMagneticField(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      SizeType actualIndexOfSite,
	                      RealType factorForDiagonals) const
	{
		if (modelParameters_.magneticField.n_row()<3) return;
		assert(actualIndexOfSite<modelParameters_.magneticField.n_col());
		for (SizeType orb=0;orb<modelParameters_.orbitals;orb++)
			addMagneticField(hmatrix,cm,i,actualIndexOfSite,orb,factorForDiagonals);
	}

	void addMagneticField(SparseMatrixType &hmatrix,
	                      const VectorOperatorType& cm,
	                      SizeType i,
	                      SizeType actualIndexOfSite,
	                      SizeType orbital,
	                      RealType factorForDiagonals) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType cup = cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data;
		SparseMatrixType cupTranspose;
		transposeConjugate(cupTranspose,cup);
		SparseMatrixType cdown = cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;
		SparseMatrixType A = cupTranspose * cdown;
		SparseMatrixType Atranspose;
		transposeConjugate(Atranspose,A);

		hmatrix += factorForDiagonals *
		        modelParameters_.magneticField(0,actualIndexOfSite) * A;

		hmatrix += factorForDiagonals *
		        modelParameters_.magneticField(1,actualIndexOfSite) * Atranspose;

		SparseMatrixType nup = n(cup);
		SparseMatrixType ndown = n(cdown);

		SparseMatrixType tmp = nup;
		tmp += (-1.0)*ndown;
		hmatrix += factorForDiagonals *
		        modelParameters_.magneticField(2,actualIndexOfSite) * tmp;

	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorOperatorType& cm,
	                   SizeType i,
	                   SizeType actualIndexOfSite,
	                   RealType factorForDiagonals,
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
				addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb,factorForDiagonals,V);
		}

		if (V.size() == v2) {
			for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {
				for (SizeType orb2=0;orb2<modelParameters_.orbitals;orb2++) {
					addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb,orb2,factorForDiagonals,V);
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
	                   RealType factorForDiagonals,
	                   const typename PsimagLite::Vector<RealType>::Type& V) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType nup = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data);
		SparseMatrixType ndown = n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data);

		SizeType linSize = geometry_.numberOfSites();

		SizeType iUp = actualIndexOfSite + (orbital + 0*modelParameters_.orbitals)*linSize;
		hmatrix += factorForDiagonals * V[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orbital + 1*modelParameters_.orbitals)*linSize;
		hmatrix += factorForDiagonals * V[iDown] * ndown;
	}

	void addPotentialV(SparseMatrixType &hmatrix,
	                   const VectorOperatorType& cm,
	                   SizeType i,
	                   SizeType actualIndexOfSite,
	                   SizeType orb,
	                   SizeType orb2,
	                   RealType factorForDiagonals,
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
		hmatrix += factorForDiagonals * V[iUp] * nup;
		SizeType iDown = actualIndexOfSite + (orb + orb2*modelParameters_.orbitals +
		                                      1*orbitalsSquared)*linSize;
		hmatrix += factorForDiagonals * V[iDown] * ndown;
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
			break;
		case 1: // S^-
			return spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_UP);
			break;
		}
		SparseMatrixType tmpMatrix=spinOperatorAux(cm,i,orbital,SPIN_UP,SPIN_UP);
		SparseMatrixType tmpMatrix2=spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_DOWN);
		tmpMatrix += (-1.0)*tmpMatrix2;
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
	                           SizeType site,
	                           RealType factorForDiagonals) const
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

						hmatrix += factorForDiagonals*U[offset]*tmpMatrix;
					}
				}
			}
		}
	}

	//! only for feAsMode == 4
	void addInteractionKspace(SparseMatrixType &hmatrix,
	                          const VectorOperatorType& cm,
	                          SizeType i,
	                          RealType factorForDiagonals,
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

					hmatrix += factorForDiagonals*modelParameters_.hubbardU[0]*tmpMatrix3;
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
	                            RealType factorForDiagonals,
	                            SizeType actualSite) const
	{
		if (actualSite > 0) return;

		addInteractionU1(hmatrix,cm,i,factorForDiagonals);
		addInteractionImp2(hmatrix,cm,i,factorForDiagonals);
		addInteractionImp3(hmatrix,cm,i,factorForDiagonals);
		addInteractionImp4(hmatrix,cm,i,factorForDiagonals);
	}

	void addInteractionImp2(SparseMatrixType &hmatrix,
	                        const VectorOperatorType& cm,
	                        SizeType site,
	                        RealType factorForDiagonals) const
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
					hmatrix += factorForDiagonals*U[1]*tmpMatrix;
				}
			}
		}
	}

	void addInteractionImp3(SparseMatrixType &hmatrix,
	                        const VectorOperatorType& cm,
	                        SizeType site,
	                        RealType factorForDiagonals) const
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
				hmatrix += factorForDiagonals*U[2]*tmpMatrix;
			}
		}
	}

	void addInteractionImp4(SparseMatrixType &hmatrix,
	                        const VectorOperatorType& cm,
	                        SizeType site,
	                        RealType factorForDiagonals) const
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
					hmatrix += factorForDiagonals*U[3]*tmpMatrix*tmpMatrix2;
				}
			}
		}
	}

	//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteractionAncilla(SparseMatrixType &hmatrix,
	                           const VectorOperatorType& cm,
	                           SizeType i,
	                           RealType factorForDiagonals,
	                           SizeType actualSite) const
	{
		int dof=2*modelParameters_.orbitals;
		SparseMatrixType tmpMatrix;

		SizeType alpha = 0; // real sites, no ancilla
		SparseMatrixType m1=cm[alpha+SPIN_UP*modelParameters_.orbitals+i*dof].data;
		SparseMatrixType m2=cm[alpha+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;

		multiply(tmpMatrix,n(m1),n(m2));
		assert(actualSite < modelParameters_.hubbardU.size());
		hmatrix += factorForDiagonals*modelParameters_.hubbardU[actualSite]*tmpMatrix;
	}

	bool isAllowedThisDof(SizeType alpha) const
	{
		SizeType electrons = HilbertSpaceFeAsType::electrons(alpha);

		return (electrons >= modelParameters_.minElectronsPerSite);
	}

	void diagTest(const SparseMatrixType& fullm,const PsimagLite::String& str) const
	{
		if (fullm.rank()!=256) return;
		MatrixType fullm2;
		crsMatrixToFullMatrix(fullm2,fullm);
		typename PsimagLite::Vector<SparseElementType>::Type eigs(fullm2.n_row());
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
}; //class ModelFeBasedSc
} // namespace Dmrg
/*@}*/
#endif

