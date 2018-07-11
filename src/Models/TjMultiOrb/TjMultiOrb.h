/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file TjMultiOrb.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef DMRG_TJ_MULTIORB_H
#define DMRG_TJ_MULTIORB_H
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "../Models/TjMultiOrb/LinkProductTjMultiOrb.h"
#include "../Models/TjMultiOrb/ParametersModelTjMultiOrb.h"
#include "../Models/FeAsModel/HilbertSpaceFeAs.h"
#include "ProgramGlobals.h"

namespace Dmrg {
//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
template<typename ModelBaseType>
class TjMultiOrb : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef ModelHubbard<ModelBaseType> ModelHubbardType;
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
	typedef LinkProductTjMultiOrb<ModelHelperType, GeometryType> LinkProductType;
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisFeAsType;
	typedef typename HilbertBasisFeAsType::value_type HilbertStateFeAs;
	typedef  HilbertSpaceFeAs<HilbertStateFeAs> HilbertSpaceFeAsType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ModelHubbardType::HilbertState HilbertStateType;
	typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename OperatorType::PairType PairType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename PsimagLite::Vector<HilbertStateType>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	static const int FERMION_SIGN = -1;

	enum {REINTERPRET_6 = 6, REINTERPRET_9 = 9};

	enum {STATE_EMPTY = 0, STATE_UP_A = 1, STATE_DOWN_A = 4};

	enum {SPIN_UP, SPIN_DOWN};

	TjMultiOrb(const SolverParamsType& solverParams,
	           InputValidatorType& io,
	           GeometryType const &geometry)
	    : ModelBaseType(solverParams, geometry, new LinkProductType(io)),
	      modelParameters_(io),
	      geometry_(geometry),
	      offset_(5*modelParameters_.orbitals), // c^\dagger_up, c^\dagger_down, S+, Sz, n
	      spinSquared_(spinSquaredHelper_,modelParameters_.orbitals,2*modelParameters_.orbitals)
	{
		if (modelParameters_.orbitals > 1) {
			PsimagLite::String str("TjMultiOrb with more than 1 orbital is EXPERIMENTAL\n");
			std::cerr<<str<<"\n";
			std::cout<<str<<"\n";
		}

		if (modelParameters_.potentialV.size() !=
		        2*geometry_.numberOfSites()*modelParameters_.orbitals)
			throw PsimagLite::RuntimeError("potentialV must be of size 2*sites*orbitals\n");

		// fill caches
		ProgramGlobals::init(modelParameters_.orbitals*geometry_.numberOfSites() + 1);
		BlockType block(1,0);
		setNaturalBasis(basis_,q_,block,true);
		setOperatorMatrices(creationMatrix_,block);
		//! Set symmetry related
		setSymmetryRelated(qq_,basis_,block.size());
	}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType) const
	{
		assert(0 < creationMatrix_.size());
		return creationMatrix_[0].data.rows();
	}

	/** \cppFunction{!PTEX_THISFUNCTION} returns the operator in the
	 *unmangled (natural) basis of one-site */
	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		assert(creationMatrix_.size()>0);
		SizeType nrow = creationMatrix_[0].data.rows();
		SizeType orbitals = modelParameters_.orbitals;

		if (what == "i" || what=="identity") {
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

		if (what == "splus") {
			VectorSizeType allowed(orbitals,0);
			for (SizeType i = 0; i < allowed.size(); ++i) allowed[i] = i;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			assert(dof < orbitals);
			return creationMatrix_[2*orbitals+dof];
		}

		if (what == "sminus") {
			VectorSizeType allowed(orbitals,0);
			for (SizeType i = 0; i < allowed.size(); ++i) allowed[i] = i;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			assert(dof < orbitals);
			OperatorType cm = creationMatrix_[2*orbitals+dof];
			cm.dagger();
			return cm;
		}

		if (what == "z" || what == "sz") { // S^z
			VectorSizeType allowed(orbitals,0);
			for (SizeType i = 0; i < allowed.size(); ++i) allowed[i] = i;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			assert(dof < orbitals);
			return creationMatrix_[3*orbitals+dof];
		}

		if (what=="c") {
			VectorSizeType allowed(2*orbitals,0);
			for (SizeType i = 0; i < allowed.size(); ++i) allowed[i] = i;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			assert(dof < 2*orbitals && creationMatrix_.size() > dof);
			return creationMatrix_[dof];
		}

		if (what=="n") {
			VectorSizeType allowed(orbitals,0);
			for (SizeType i = 0; i < allowed.size(); ++i) allowed[i] = i;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			assert(dof < orbitals);
			return creationMatrix_[4*orbitals+dof];
		}

		if (what=="nup") {
			VectorSizeType allowed(orbitals,0);
			for (SizeType i = 0; i < allowed.size(); ++i) allowed[i] = i;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			OperatorType cup = naturalOperator("c",site,dof+SPIN_UP*orbitals);
			cup.dagger();
			SparseMatrixType nup(multiplyTc(cup.data,cup.data));
			if (orbitals>1) nup = findNMatrices(dof+SPIN_UP*orbitals);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(nup,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		if (what=="ndown") {
			VectorSizeType allowed(orbitals,0);
			for (SizeType i = 0; i < allowed.size(); ++i) allowed[i] = i;
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			OperatorType cdown = naturalOperator("c",site,dof+SPIN_DOWN*orbitals);
			cdown.dagger();
			SparseMatrixType ndown(multiplyTc(cdown.data,cdown.data));
			if (orbitals>1) ndown = findNMatrices(dof+SPIN_DOWN*orbitals);
			typename OperatorType::Su2RelatedType su2Related;
			return OperatorType(ndown,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related);
		}

		PsimagLite::String str("TjMultiOrb: naturalOperator: no label ");
		str += what + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	//! find total number of electrons for each state in the basis
	void findElectrons(typename PsimagLite::Vector<SizeType> ::Type&electrons,
	                   const VectorHilbertStateType& basis,
	                   SizeType) const
	{
		SizeType numberOfDofs = 2*modelParameters_.orbitals;
		electrons.clear();
		for (SizeType i=0;i<basis.size();i++) {
			SizeType sum = 0;
			for (SizeType dof = 0; dof < numberOfDofs; ++dof) {
				HilbertStateType mask = (1<<dof);
				if (mask & basis[i]) sum++;
			}

			electrons.push_back(sum);
		}
	}

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void setBasis(HilbertBasisType  &basis,
	              SymmetryElectronsSzType& qq,
	              const VectorSizeType&) const
	{
		basis = basis_;
		qq = qq_;
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
		        io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		io.write(label + "/offset_", offset_);
		spinSquaredHelper_.write(label, io);
		spinSquared_.write(label, io);
		io.write(label + "/basis_", basis_);
		qq_.write(label, io);
		io.write(label + "/q_", q_);
		io.write(label + "/creationMatrix_", creationMatrix_);
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findCreationMatrices(int,
	                                      SizeType sigma,
	                                      const VectorHilbertStateType&) const
	{
		assert(sigma < creationMatrix_.size());
		return creationMatrix_[sigma].data;
	}

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findNiMatrices(int,
	                                SizeType orb,
	                                const VectorHilbertStateType&) const
	{
		assert(4*modelParameters_.orbitals + orb < creationMatrix_.size());
		return creationMatrix_[4*modelParameters_.orbitals + orb].data;
	}

	//! Find S^+_i in the natural basis natBasis
	SparseMatrixType findSplusMatrices(int,
	                                   SizeType orb,
	                                   const VectorHilbertStateType&) const
	{
		assert(2*modelParameters_.orbitals + orb < creationMatrix_.size());
		return creationMatrix_[2*modelParameters_.orbitals+orb].data;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(int,
	                                SizeType orb,
	                                const VectorHilbertStateType&) const
	{
		assert(3*modelParameters_.orbitals + orb< creationMatrix_.size());
		return creationMatrix_[3*modelParameters_.orbitals + orb].data;
	}

	virtual SizeType maxElectronsOneSpin() const
	{
		return modelParameters_.orbitals*geometry_.numberOfSites() + 1;
	}

private:

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findCreationMatrices(int i,
	                                      int sigma,
	                                      const VectorHilbertStateType& natBasis,
	                                      const MatrixType* rot,
	                                      const MatrixType* rotT) const
	{
		assert(i == 0);
		HilbertStateType bra,ket;
		int n = natBasis.size();
		MatrixType cm(n,n);
		SizeType orbitals = modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			HilbertStateType mask = (1<<(sigma+i*2*orbitals));
			if (ket & mask) continue;
			bra = (ket ^ mask);
			int jj = PsimagLite::isInVector(natBasis,bra);
			if (jj<0) continue;
			cm(ii,jj) = sign(ket,i,sigma);

		}

		truncateMatrix(cm,rot,rotT,natBasis);

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findNiMatrices(int i,
	                                SizeType orb,
	                                const VectorHilbertStateType& natBasis,
	                                const MatrixType* rot,
	                                const MatrixType* rotT) const
	{
		assert(i == 0);
		SizeType orbitals = modelParameters_.orbitals;
		SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);
		SizeType dofs = 2*modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			cm(ii,ii) = 0.0;
			SizeType orbitals = modelParameters_.orbitals;
			for (SizeType sigma=0;sigma<dofs;sigma+=orbitals) {
				HilbertStateType mask = (1<<(orb+sigma));
				if (ket & mask) cm(ii,ii) += 1.0;
			}
		}

		if (orbitals > 1) correctLambda(cm,natBasis);
		truncateMatrix(cm,rot,rotT,natBasis);

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	SparseMatrixType findNMatrices(SizeType sigma) const
	{
		assert(sigma<2*modelParameters_.orbitals);
		VectorHilbertStateType natBasis;
		VectorSizeType quantumNumbs;
		BlockType block(1,0);
		setNaturalBasis(natBasis,quantumNumbs,block,false);
		SizeType n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			cm(ii,ii) = 0.0;
			HilbertStateType mask = (1<<sigma);
			if (ket & mask) cm(ii,ii) += 1.0;
		}

		MatrixType rotation(n,n);
		MatrixType rotationR(n,n);
		computeRotation(rotation,rotationR,natBasis);
		truncateMatrix(cm,&rotation,&rotationR,natBasis);

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}


	//! Find S^+_i in the natural basis natBasis
	SparseMatrixType findSplusMatrices(int i,
	                                   SizeType orb,
	                                   const VectorHilbertStateType& natBasis,
	                                   const MatrixType* rot,
	                                   const MatrixType* rotT) const
	{
		assert(i == 0);
		HilbertStateType bra,ket;
		int n = natBasis.size();
		MatrixType cm(n,n);
		SizeType orbitals = modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			ket=natBasis[ii];
			SizeType l = orb+orbitals;
			bra = ket;
			HilbertStateType masklp = (1<<l);
			HilbertStateType masklm = (1<<(l-orbitals));
			if ((ket & masklp) > 0 && (ket & masklm) == 0) {
				bra = ket ^ (masklp | masklm);

				int jj = PsimagLite::isInVector(natBasis,bra);
				assert(jj>=0);
				cm(ii,jj) = 1.0;
			}
		}
		if (orbitals > 1) correctLambda(cm,natBasis);
		truncateMatrix(cm,rot,rotT,natBasis);

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(int i,
	                                SizeType orb,
	                                const VectorHilbertStateType& natBasis,
	                                const MatrixType* rot,
	                                const MatrixType* rotT) const
	{
		assert(i == 0);
		int n = natBasis.size();
		MatrixType cm(n,n);
		SizeType orbitals = modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			int counter = 0;
			SizeType l1 = orb;
			HilbertStateType masklp = (1<<l1);
			if (ket & masklp) counter++;
			SizeType l2 = orb+orbitals;
			HilbertStateType masklp1 = (1<<l2);
			if (ket & masklp1) counter--;
			cm(ii,ii) = 0.5*counter;
		}

		if (orbitals > 1) correctLambda(cm,natBasis);
		truncateMatrix(cm,rot,rotT,natBasis);

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	void correctLambda(MatrixType& cm2,
	                   const VectorHilbertStateType& natBasis) const
	{
		SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> dens(n,n);
		SizeType dofs = 2*modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			dens(ii,ii) = 0.0;
			for (SizeType sigma=0;sigma<dofs;sigma++) {
				HilbertStateType mask = (1<<sigma);
				if (ket & mask) dens(ii,ii) += 1.0;
			}
		}

		MatrixType corrector(n,n);
		SparseElementType f1 = (-1.0);
		SparseElementType f2 = 0.5;
		for (SizeType i = 0; i < n; ++i)
			corrector(i,i) = f2 * dens(i,i) * std::abs(dens(i,i) + f1);

		cm2 = cm2 * corrector;
	}

	void truncateMatrix(MatrixType& cm,
	                    const MatrixType* rot,
	                    const MatrixType* rotT,
	                    const VectorHilbertStateType& natBasis) const
	{
		if (modelParameters_.orbitals != 2 ||
		        modelParameters_.reinterpretAndTruncate == 0) return;

		if (!rot || !rotT) return;

		cm  = (*rot)*cm;
		cm = cm*(*rotT);

		VectorSizeType target;
		findIndicesToRemove(target, natBasis);
		assert(target.size() > 0);

		SizeType n = cm.rows();
		assert(n == cm.cols());
		assert(n > target.size());
		n -= target.size();

		MatrixType cm2(n,n);
		SizeType ii = 0;
		for (SizeType i = 0; i < cm.rows(); ++i) {
			VectorSizeType::const_iterator it = std::find(target.begin(),
			                                              target.end(),
			                                              i);
			if (it != target.end()) continue;
			SizeType jj = 0;
			for (SizeType j = 0; j < cm.cols(); ++j) {
				VectorSizeType::const_iterator it = std::find(target.begin(), target.end(), j);
				if (it != target.end()) continue;
				cm2(ii,jj) = cm(i,j);
				jj++;
			}

			ii++;
		}

		cm = cm2;
	}

	//! set creation matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         const BlockType& block) const
	{
		VectorHilbertStateType natBasis;
		SparseMatrixType tmpMatrix;
		VectorSizeType quantumNumbs;
		setNaturalBasis(natBasis,quantumNumbs,block,false);

		SizeType dof = 2*modelParameters_.orbitals;
		// Set the operators c^\daggger_{i\sigma} in the natural basis
		creationMatrix.clear();

		SizeType n = natBasis.size();
		MatrixType rotation(n,n);
		MatrixType rotationR(n,n);
		computeRotation(rotation,rotationR,natBasis);
		for (SizeType i=0;i<block.size();i++) {
			for (SizeType sigma=0;sigma<dof;sigma++) {
				// orbital changes first
				tmpMatrix = findCreationMatrices(i,sigma,natBasis,&rotation,&rotationR);
				int asign= 1;
				if (sigma>0) asign= 1;
				typename OperatorType::Su2RelatedType su2related;
				if (sigma==0) {
					su2related.source.push_back(i*offset_);
					su2related.source.push_back(i*offset_+1);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = modelParameters_.orbitals;
				}

				OperatorType myOp(tmpMatrix,-1,PairType(1,1-sigma),asign,su2related);

				creationMatrix.push_back(myOp);
			}

			// Set the operators S^+_i in the natural basis
			for (SizeType i=0;i<block.size();i++) {
				for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {
					tmpMatrix = findSplusMatrices(i,orb,natBasis,&rotation,&rotationR);

					typename OperatorType::Su2RelatedType su2related;
					su2related.source.push_back(i*modelParameters_.orbitals*2);
					su2related.source.push_back(i*modelParameters_.orbitals*2+modelParameters_.orbitals);
					su2related.source.push_back(i*modelParameters_.orbitals*2);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(1);
					su2related.offset = modelParameters_.orbitals;

					OperatorType myOp(tmpMatrix,1,PairType(2,2),-1,su2related);
					creationMatrix.push_back(myOp);
				}
			}
			// Set the operators S^z_i in the natural basis
			for (SizeType i=0;i<block.size();i++) {
				for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {

					tmpMatrix = findSzMatrices(i,orb,natBasis,&rotation,&rotationR);
					typename OperatorType::Su2RelatedType su2related2;
					OperatorType myOp2(tmpMatrix,1,PairType(2,1),1.0/sqrt(2.0),su2related2);
					creationMatrix.push_back(myOp2);
				}
			}
			// Set ni matrices
			for (SizeType i=0;i<block.size();i++) {
				for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {
					tmpMatrix = findNiMatrices(0,orb,natBasis,&rotation,&rotationR);
					RealType angularFactor= 1;
					typename OperatorType::Su2RelatedType su2related3;
					su2related3.offset = 1; //check FIXME
					OperatorType myOp3(tmpMatrix,1,PairType(0,0),angularFactor,su2related3);
					creationMatrix.push_back(myOp3);
				}
			}

		}
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma}
	//! to basis state ket
	RealType sign(HilbertStateType const &ket, int i, int sigma) const
	{
		SizeType orbitals = modelParameters_.orbitals;
		if (orbitals == 1) return 1;

		int value=0;
		SizeType dofs=2*modelParameters_.orbitals;
		for (SizeType alpha=0;alpha<dofs;alpha++)
			value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);

		if (i>0) value += HilbertSpaceFeAsType::electrons(ket);

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

	void computeRotation(MatrixType& u,
	                     MatrixType& uT,
	                     const VectorHilbertStateType& natBasis) const
	{
		if (modelParameters_.orbitals != 2 ||
		        modelParameters_.reinterpretAndTruncate == 0) return;

		SizeType n = natBasis.size();
		for (SizeType ii=0;ii<n;ii++) {
			for (SizeType jj=0;jj<n;jj++) {
				HilbertStateType ket = natBasis[ii];
				HilbertStateType bra = natBasis[jj];

				if (ket == bra) u(ii,jj)=1;
				if (ket == REINTERPRET_6 && bra == REINTERPRET_9) u(ii,jj)=-1/sqrt(2.0);
				if (ket == REINTERPRET_9 && bra == REINTERPRET_6) u(ii,jj)=1/sqrt(2.0);
				if (ket == REINTERPRET_6 && bra == REINTERPRET_6) u(ii,jj)=1/sqrt(2.0);
				if (ket == REINTERPRET_9 && bra == REINTERPRET_9) u(ii,jj)=1/sqrt(2.0);
			}

		}

		transposeConjugate(uT,u);
	}

	void findIndicesToRemove(VectorSizeType& indices,
	                         const VectorHilbertStateType& natBasis) const
	{
		VectorHilbertStateType target(1,REINTERPRET_6);
		if (modelParameters_.reinterpretAndTruncate > 1)
			target.push_back(STATE_EMPTY);
		if (modelParameters_.reinterpretAndTruncate > 2) {
			target.push_back(STATE_UP_A);
			target.push_back(STATE_DOWN_A);
		}

		for (SizeType i = 0; i < natBasis.size(); ++i) {
			typename VectorHilbertStateType::const_iterator it = std::find(target.begin(),
			                                                               target.end(),
			                                                               natBasis[i]);
			if (it == target.end()) continue;
			indices.push_back(i);
		}
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType&,
	                                const BlockType& block,
	                                RealType,
	                                RealType factorForDiagonals=1.0) const
	{
		SizeType n=block.size();
		SizeType orbitals = modelParameters_.orbitals;
		SizeType linSize = geometry_.numberOfSites();
		for (SizeType i=0;i<n;i++) {
			for (SizeType orb = 0; orb < orbitals; ++orb) {
				// potentialV
				SparseMatrixType nup(naturalOperator("nup",i,orb).data);
				SparseMatrixType ndown(naturalOperator("ndown",i,orb).data);
				SparseMatrixType m = nup;
				assert(block[i]+linSize*orb+linSize*orbitals<modelParameters_.potentialV.size());
				m *= modelParameters_.potentialV[block[i]+linSize*orb];
				m += modelParameters_.potentialV[block[i]+linSize*orb+linSize*orbitals]*ndown;
				hmatrix += factorForDiagonals * m;
			}
		}
	}

	void setNaturalBasis(HilbertBasisType& basis,
	                     VectorSizeType& q,
	                     const VectorSizeType& block,
	                     bool truncated) const
	{
		assert(block.size()==1);
		HilbertStateType a=0;
		HilbertStateType total = (1 << 2*modelParameters_.orbitals);

		HilbertBasisType basisTmp;
		for (a=0;a<total;a++) basisTmp.push_back(a);
		weedOutBasis(basisTmp,truncated);

		// reorder the natural basis (needed for MULTIPLE BANDS)
		SymmetryElectronsSzType qq;
		setSymmetryRelated(qq, basisTmp, 1);
		ModelBaseType::orderBasis(basis, basisTmp, qq);
	}

	void weedOutBasis(VectorHilbertStateType& basis, bool truncated) const
	{
		SizeType orbitals = modelParameters_.orbitals;
		HilbertBasisType basisTmp;
		VectorHilbertStateType electrons(orbitals,0);
		if (orbitals != 2) truncated = false;

		for (SizeType i = 0; i < basis.size(); ++i) {
			HilbertStateType ket = basis[i];
			if (truncated && modelParameters_.reinterpretAndTruncate > 0 && ket == REINTERPRET_6)
				continue;
			if (truncated && modelParameters_.reinterpretAndTruncate > 1 && ket == STATE_EMPTY)
				continue;
			bool b = (ket == STATE_UP_A || ket == STATE_DOWN_A);
			if (truncated && modelParameters_.reinterpretAndTruncate > 2 && b)
				continue;

			SizeType orb = 0;
			while (ket > 0) {
				if (ket & 1) electrons[orb]++;
				orb++;
				if (orb >= orbitals) orb = 0;
				ket >>= 1;
			}

			bool addIt = true;
			for (SizeType orb = 0; orb < electrons.size(); ++orb) {
				if (electrons[orb] > 1) addIt = false;
				electrons[orb] = 0;
			}

			if (addIt) basisTmp.push_back(basis[i]);
		}

		basis = basisTmp;
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
		SizeType orbitals = modelParameters_.orbitals;
		VectorSizeType electronsUp(basis.size());
		VectorSizeType electrons(basis.size());
		for (SizeType i=0;i<basis.size();i++) {
			PairType jmpair(0,0);

			if(orbitals==1) {
				jmpair = calcJmvalue<PairType>(basis[i]);
			}
			
			jmvalues.push_back(jmpair);
			// nup
			electronsUp[i] = 0;
			SizeType electronsDown  = 0;
			for (SizeType dof = 0; dof < orbitals; ++dof) {
				HilbertStateType mask = (1<<dof);
				if (mask & basis[i]) electronsUp[i]++;
				mask = (1<<(dof+orbitals));
				if (mask & basis[i]) electronsDown++;
			}

			electrons[i] = electronsDown + electronsUp[i];
			flavors.push_back(electrons[i]);
			jmSaved = jmpair;
		}

		q.set(jmvalues,flavors,electrons,electronsUp);
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// Reinterprets 6 and 9
	template<typename PairType>
	PairType calcJmvalue(const HilbertStateType& ket) const
	{
		return calcJmValueAux<PairType>(ket);
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// does not work for 6 or 9
	template<typename PairType>
	PairType calcJmValueAux(const HilbertStateType& ket) const
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

	//serializr start class TjMultiOrb
	//serializr vptr
	//serializr normal modelParameters_
	ParametersModelTjMultiOrb<RealType>  modelParameters_;
	//serializr ref geometry_ end
	const GeometryType &geometry_;
	//serializr normal offset_
	SizeType offset_;
	//serializr normal spinSquaredHelper_
	SpinSquaredHelper<RealType,HilbertStateType> spinSquaredHelper_;
	//serializr normal spinSquared_
	SpinSquared<SpinSquaredHelper<RealType,HilbertStateType> > spinSquared_;
	HilbertBasisType basis_;
	SymmetryElectronsSzType qq_;
	VectorSizeType q_;
	VectorOperatorType creationMatrix_;
};	//class TjMultiOrb

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_MULTIORB_H

