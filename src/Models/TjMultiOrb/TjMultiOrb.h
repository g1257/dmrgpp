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
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisFeAsType;
	typedef typename HilbertBasisFeAsType::value_type HilbertStateFeAs;
	typedef HilbertSpaceFeAs<HilbertStateFeAs> HilbertSpaceFeAsType;
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelHubbardType::HilbertState HilbertStateType;
	typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename OperatorType::PairType PairType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<HilbertStateType>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const int FERMION_SIGN = -1;

	enum {REINTERPRET_6 = 6, REINTERPRET_9 = 9};

	enum {STATE_EMPTY = 0, STATE_UP_A = 1, STATE_DOWN_A = 4};

	enum {SPIN_UP, SPIN_DOWN};

	TjMultiOrb(const SolverParamsType& solverParams,
	           InputValidatorType& io,
	           const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      superGeometry_(geometry),
	      offset_(5*modelParameters_.orbitals), // c^\dagger_up, c^\dagger_down, S+, Sz, n
	      spinSquared_(spinSquaredHelper_,modelParameters_.orbitals,2*modelParameters_.orbitals)
	{
		if (modelParameters_.orbitals > 1) {
			PsimagLite::String str("TjMultiOrb with more than 1 orbital is EXPERIMENTAL\n");
			std::cerr<<str<<"\n";
			std::cout<<str<<"\n";
		}

		if (modelParameters_.potentialV.size() !=
		        2*superGeometry_.numberOfSites()*modelParameters_.orbitals)
			throw PsimagLite::RuntimeError("potentialV must be of size 2*sites*orbitals\n");

		// fill caches
		ProgramGlobals::init(modelParameters_.orbitals*superGeometry_.numberOfSites() + 1);
		BlockType block(1,0);
		setNaturalBasis(basis_, block, true);
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findCreationMatrices(int,
	                                      SizeType sigma,
	                                      const VectorHilbertStateType&) const
	{
		assert(sigma < creationMatrix_.size());
		return creationMatrix_[sigma].getCRS();
	}

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findNiMatrices(int,
	                                SizeType orb,
	                                const VectorHilbertStateType&) const
	{
		assert(4*modelParameters_.orbitals + orb < creationMatrix_.size());
		return creationMatrix_[4*modelParameters_.orbitals + orb].getCRS();
	}

	//! Find S^+_i in the natural basis natBasis
	SparseMatrixType findSplusMatrices(int,
	                                   SizeType orb,
	                                   const VectorHilbertStateType&) const
	{
		assert(2*modelParameters_.orbitals + orb < creationMatrix_.size());
		return creationMatrix_[2*modelParameters_.orbitals+orb].getCRS();
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(int,
	                                SizeType orb,
	                                const VectorHilbertStateType&) const
	{
		assert(3*modelParameters_.orbitals + orb< creationMatrix_.size());
		return creationMatrix_[3*modelParameters_.orbitals + orb].getCRS();
	}

	SizeType maxElectronsOneSpin() const
	{
		return modelParameters_.orbitals*superGeometry_.numberOfSites() + 1;
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
		io.write(label + "/qq_", qq_);
		io.write(label + "/creationMatrix_", creationMatrix_);
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		VectorSizeType block(1, site);
		VectorHilbertStateType natBasis;
		SparseMatrixType tmpMatrix;
		setNaturalBasis(natBasis, block, false);
		setSymmetryRelated(qns, natBasis, block.size());

		SizeType dof = 2*modelParameters_.orbitals;
		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& splus = this->createOpsLabel("splus");
		OpsLabelType& sz = this->createOpsLabel("sz");
		OpsLabelType& nop = this->createOpsLabel("n");
		OpsLabelType& sminus = this->createOpsLabel("sminus");

		this->makeTrackable("c");
		this->makeTrackable("splus");
		this->makeTrackable("sz");
		this->makeTrackable("n");

		// Set the operators c^\daggger_{i\sigma} in the natural basis

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

				OperatorType myOp(tmpMatrix,
				                  ProgramGlobals::FermionOrBosonEnum::FERMION,
				                  PairType(1, 1 - sigma),
				                  asign,
				                  su2related);

				c.push(myOp);
			}

			// Set the operators S^+_i in the natural basis
			for (SizeType i=0;i<block.size();i++) {
				for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {
					tmpMatrix = findSplusMatrices(i,orb,natBasis,&rotation,&rotationR);

					typename OperatorType::Su2RelatedType su2related;
					su2related.source.push_back(i*modelParameters_.orbitals*2);
					su2related.source.push_back(i*modelParameters_.orbitals*2 +
					                            modelParameters_.orbitals);
					su2related.source.push_back(i*modelParameters_.orbitals*2);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(1);
					su2related.offset = modelParameters_.orbitals;

					OperatorType myOp(tmpMatrix,
					                  ProgramGlobals::FermionOrBosonEnum::BOSON,
					                  PairType(2, 2),
					                  -1,
					                  su2related);
					splus.push(myOp);
					myOp.dagger();
					sminus.push(myOp);
				}
			}
			// Set the operators S^z_i in the natural basis
			for (SizeType i=0;i<block.size();i++) {
				for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {

					tmpMatrix = findSzMatrices(i,orb,natBasis,&rotation,&rotationR);
					typename OperatorType::Su2RelatedType su2related2;
					OperatorType myOp2(tmpMatrix,
					                   ProgramGlobals::FermionOrBosonEnum::BOSON,
					                   PairType(2, 1),
					                   1.0/sqrt(2.0),
					                   su2related2);
					sz.push(myOp2);
				}
			}
			// Set ni matrices
			for (SizeType i=0;i<block.size();i++) {
				for (SizeType orb=0;orb<modelParameters_.orbitals;orb++) {
					tmpMatrix = findNiMatrices(0,orb,natBasis,&rotation,&rotationR);
					RealType angularFactor= 1;
					typename OperatorType::Su2RelatedType su2related3;
					su2related3.offset = 1; //check FIXME
					OperatorType myOp3(tmpMatrix,
					                   ProgramGlobals::FermionOrBosonEnum::BOSON,
					                   PairType(0, 0),
					                   angularFactor,
					                   su2related3);
					nop.push(myOp3);
				}
			}
		}

		OpsLabelType& nupop = this->createOpsLabel("nup");
		for (SizeType dof = 0; dof < modelParameters_.orbitals; ++dof) {
			SizeType x = dof + SPIN_UP*modelParameters_.orbitals;
			OperatorType cup = this->naturalOperator("c", site, x);
			cup.dagger();
			SparseMatrixType nup(multiplyTc(cup.getCRS(),cup.getCRS()));
			if (modelParameters_.orbitals > 1)
				nup = findNMatrices(dof + SPIN_UP*modelParameters_.orbitals);
			typename OperatorType::Su2RelatedType su2Related;
			nupop.push(OperatorType(nup,
			                        ProgramGlobals::FermionOrBosonEnum::BOSON,
			                        typename OperatorType::PairType(0,0),
			                        1.0,
			                        su2Related));
		}

		OpsLabelType& ndownpop = this->createOpsLabel("ndown");
		for (SizeType dof = 0; dof < modelParameters_.orbitals; ++dof) {
			SizeType x = dof + SPIN_DOWN*modelParameters_.orbitals;
			OperatorType cdown = this->naturalOperator("c", site, x);
			cdown.dagger();
			SparseMatrixType ndown(multiplyTc(cdown.getCRS(),cdown.getCRS()));
			if (modelParameters_.orbitals > 1)
				ndown = findNMatrices(dof + SPIN_DOWN*modelParameters_.orbitals);
			typename OperatorType::Su2RelatedType su2Related;
			ndownpop.push(OperatorType(ndown,
			                           ProgramGlobals::FermionOrBosonEnum::BOSON,
			                           typename OperatorType::PairType(0,0),
			                           1.0,
			                           su2Related));
		}
	}

	void fillModelLinks()
	{
		//! There are orbitals*orbitals different orbitals
		//! and 2 spins. Spin is diagonal so we end up with 2*orbitals*orbitals possiblities
		//! a up a up, a up b up, b up a up, b up, b up, etc
		//! and similarly for spin down.

		const SizeType orbitals = modelParameters_.orbitals;
		ModelTermType& hop = ModelBaseType::createTerm("hopping");
		bool isSu2 = BasisType::useSu2Symmetry();
		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");
		ModelTermType& szsz = ModelBaseType::createTerm("szsz");
		ModelTermType& ninj = ModelBaseType::createTerm("ninj");

		for (SizeType spin = 0; spin < 2; ++spin) {
			for (SizeType orb1 = 0; orb1 < orbitals; ++orb1) {
				OpForLinkType c1("c", orb1 + spin*orbitals, orb1);
				for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
					OpForLinkType c2("c", orb2 + spin*orbitals, orb2);

					hop.push(c1,
					         'N',
					         c2,
					         'C',
					         typename ModelTermType::Su2Properties(1, (spin == 1) ? -1 : 1, spin));
				}
			}

		}

		auto valueModiferTerm0 = [isSu2](ComplexOrRealType& value)
		{ value *= (isSu2) ? -0.5 : 0.5;};
		auto valueModifierTermOther = [isSu2](ComplexOrRealType& value)
		{ if (isSu2) value = -value;};

		for (SizeType orb1 = 0; orb1 < orbitals; ++orb1) {
			OpForLinkType splus1("splus", orb1, orb1);
			OpForLinkType sz1("sz", orb1, orb1);
			OpForLinkType n1("n", orb1, orb1);

			for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
				OpForLinkType splus2("splus", orb2, orb2);
				OpForLinkType sz2("sz", orb2, orb2);
				OpForLinkType n2("n", orb2, orb2);

				spsm.push(splus1,
				          'N',
				          splus2,
				          'C',
				          valueModiferTerm0,
				          typename ModelTermType::Su2Properties(2, -1, 2));

				if (!isSu2)
					szsz.push(sz1, 'N', sz2, 'N', typename ModelTermType::Su2Properties(2, 0.5));
				else
					spsm.push(splus1,
					          'N',
					          splus2,
					          'C',
					          valueModifierTermOther,
					          typename ModelTermType::Su2Properties(2, -1, 2));

				ninj.push(n1, 'N', n2, 'N');
			}
		}
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
			int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
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
		BlockType block(1,0);
		setNaturalBasis(natBasis, block, false);
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

				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
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
		ComplexOrRealType f1 = (-1.0);
		ComplexOrRealType f2 = 0.5;
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
	                                const BlockType& block,
	                                RealType) const
	{
		SizeType n=block.size();
		SizeType orbitals = modelParameters_.orbitals;
		SizeType linSize = superGeometry_.numberOfSites();
		for (SizeType i=0;i<n;i++) {
			for (SizeType orb = 0; orb < orbitals; ++orb) {
				// potentialV
				SparseMatrixType nup(this->naturalOperator("nup",i,orb).getCRS());
				SparseMatrixType ndown(this->naturalOperator("ndown",i,orb).getCRS());
				SparseMatrixType m = nup;
				assert(block[i]+linSize*orb+linSize*orbitals<modelParameters_.potentialV.size());
				m *= modelParameters_.potentialV[block[i]+linSize*orb];
				m += modelParameters_.potentialV[block[i]+linSize*orb+linSize*orbitals]*ndown;
				hmatrix += m;
			}
		}
	}

	void setNaturalBasis(HilbertBasisType& basis,
	                     const VectorSizeType& block,
	                     bool truncated) const
	{
		assert(block.size()==1);
		HilbertStateType total = (1 << 2*modelParameters_.orbitals);

		basis.resize(total);
		for (SizeType a = 0; a< total; ++a) basis[a] = a;
		weedOutBasis(basis, truncated);
		if (modelParameters_.orbitals == 1 && basis.size() == 3) {
			basis[0] = 0;
			basis[1] = 2;
			basis[2] = 1;
		}
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

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		assert(n==1);

		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		SizeType orbitals = modelParameters_.orbitals;
        VectorSizeType other(2, 0);
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0,0);

			if(orbitals==1)
				jmpair = calcJmvalue<PairType>(basis[i]);

			// nup
			SizeType electronsUp = 0;
			SizeType electronsDown  = 0;
			for (SizeType dof = 0; dof < orbitals; ++dof) {
				HilbertStateType mask = (1<<dof);
				if (mask & basis[i]) electronsUp++;
				mask = (1<<(dof+orbitals));
				if (mask & basis[i]) electronsDown++;
			}

			SizeType electrons = electronsDown + electronsUp;
			other[0] = electrons;
            other[1] = electronsUp;
			bool sign = electrons & 1;
			qns[i] = QnType(sign, other, jmpair, electrons);
		}
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


	ParametersModelTjMultiOrb<RealType, QnType>  modelParameters_;
	const SuperGeometryType& superGeometry_;
	SizeType offset_;
	SpinSquaredHelper<RealType,HilbertStateType> spinSquaredHelper_;
	SpinSquared<SpinSquaredHelper<RealType,HilbertStateType> > spinSquared_;
	HilbertBasisType basis_;
	VectorQnType qq_;
	VectorOperatorType creationMatrix_;
};	//class TjMultiOrb

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_MULTIORB_H

