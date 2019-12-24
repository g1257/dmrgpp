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

/*! \file TjAnisotropic.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef DMRG_TJ_Anisotropic_H
#define DMRG_TJ_Anisotropic_H
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "../Models/TjAnisotropic/ParametersModelTjAnisotropic.h"
#include "../Models/FeAsModel/HilbertSpaceFeAs.h"
#include "ProgramGlobals.h"
#include "Complex.h"

namespace Dmrg {
//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
template<typename ModelBaseType>
class TjAnisotropic : public ModelBaseType {

	enum InternalDir {DIR_X, DIR_Y, DIR_Z, DIR_NUP, DIR_NDOWN, DIR_N};

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

	TjAnisotropic(const SolverParamsType& solverParams,
	              InputValidatorType& io,
	              const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      superGeometry_(geometry),
	      offset_(5*modelParameters_.orbitals), // c^\dagger_up, c^\dagger_down, S+, Sz, n
	      spinSquared_(spinSquaredHelper_,modelParameters_.orbitals,2*modelParameters_.orbitals)
	{
		if (modelParameters_.orbitals != 1)
			throw PsimagLite::RuntimeError("TjAnisotropic: must use Orbital=1 \n");

		SizeType n = superGeometry_.numberOfSites();
		SizeType mx = modelParameters_.magneticFieldX.size();
		SizeType my = modelParameters_.magneticFieldY.size();
		SizeType mz = modelParameters_.magneticFieldZ.size();
		SizeType m = mz;

		if (modelParameters_.potentialV.size() != 2*n)
			throw PsimagLite::RuntimeError("potentialV must be of size 2*sites\n");

		if (mx != my || my != mz || mz != mx ) {
			PsimagLite::String msg("tJKitaev: If provided, ");
			msg += " MagneticField must be a vector of " + ttos(n) + " entries.\n";
			msg += " MagneticFieldX, MagneticFieldY, MagneticFieldZ must be ";
			msg += "provided in all 3 (x,y,z) directions.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (m > 0 && m != n) {
			PsimagLite::String msg("Kitaev: If provided, ");
			msg += " MagneticField must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (BasisType::useSu2Symmetry())
			err("Kitaev does not have SU(2) symmetry\n");

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
		typename MatrixType::value_type dummy = 0.0;
		setNaturalBasis(natBasis, block, false);
		setSymmetryRelated(qns, natBasis, block.size());

		SizeType dof = 2*modelParameters_.orbitals;
		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& sx = this->createOpsLabel("sx");
		OpsLabelType& sy = this->createOpsLabel("sy");
		OpsLabelType& sz = this->createOpsLabel("sz");
		OpsLabelType& nop = this->createOpsLabel("n");
		OpsLabelType& nupop = this->createOpsLabel("nup");
		OpsLabelType& ndownpop = this->createOpsLabel("ndown");

		this->makeTrackable("c");
		this->makeTrackable("sx");
		this->makeTrackable("sy");
		this->makeTrackable("sz");
		this->makeTrackable("nup");
		this->makeTrackable("ndown");
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
		}

		// Set local spin and charge matrices
		for (SizeType i=0;i<block.size();i++) {

			typename OperatorType::Su2RelatedType su2related;

			// Sx
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_X, dummy);
			OperatorType myOp1(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			sx.push(myOp1);

			// Sy
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_Y, dummy);
			OperatorType myOp2(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			sy.push(myOp2);

			// Sz
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_Z, dummy);
			OperatorType myOp3(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			sz.push(myOp3);

			// N_up
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_NUP, dummy);
			OperatorType myOp4(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			nupop.push(myOp4);

			// N_down
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_NDOWN, dummy);
			OperatorType myOp5(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			ndownpop.push(myOp5);

			// total N
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_N, dummy);
			OperatorType myOp6(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			nop.push(myOp6);
		}
	}

	void fillModelLinks()
	{
		OpForLinkType sx("sx");
		OpForLinkType sy("sy");
		OpForLinkType sz("sz");
		ModelTermType& hop = ModelBaseType::createTerm("hopping");

		for (SizeType spin = 0; spin < 2; ++spin) {
			OpForLinkType c1("c", spin, 0);
			OpForLinkType c2("c", spin, 0);
			hop.push(c1,
			         'N',
			         c2,
			         'C',
			         typename ModelTermType::Su2Properties(1, (spin == 1) ? -1 : 1, spin));
		}

		ModelBaseType::createTerm("sxsx").push(sx, 'N', sx, 'N');
		ModelBaseType::createTerm("sysy").push(sy, 'N', sy, 'N');
		ModelBaseType::createTerm("szsz").push(sz, 'N', sz, 'N');
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

	SparseMatrixType findSdirMatrices(SizeType,
	                                  const HilbertBasisType&,
	                                  InternalDir,
	                                  RealType) const
	{
		err("Kitaev needs useComplex in SolverOptions in the input file\n");
		throw PsimagLite::RuntimeError("FATAL\n");
	}

	SparseMatrixType findSdirMatrices(SizeType,// site,
	                                  const HilbertBasisType& natBasis,
	                                  InternalDir dir,
	                                  std::complex<RealType>) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		cm.setTo(0.0);
		assert(total == 3);
		// basis labels --> empty, down, up

		if (dir == DIR_X) {
			cm(1,2) = cm(2,1) = 0.5;
		} else if (dir == DIR_Y) {
			cm(1, 2) = std::complex<RealType>(0.0, 0.5);
			cm(2, 1) = std::complex<RealType>(0.0, -0.5);
		} else if (dir == DIR_Z) {
			cm(1, 1) = -0.5;
			cm(2, 2) = 0.5;
		} else if (dir == DIR_NUP) {
			cm(2, 2) = 1.0;
		} else if (dir == DIR_NDOWN) {
			cm(1, 1) = 1.0;
		} else if (dir == DIR_N) {
			cm(1, 1) = 1.0;
			cm(2, 2) = 1.0;
		} else {
			assert(false);
		}

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

		for (SizeType i = 0; i < n; ++i) {
			SizeType site = block[i];

			// magnetic field x
			const OperatorType& sx = ModelBaseType::naturalOperator("sx", site, 0);
			RealType tmp = modelParameters_.magneticFieldX[block[0]];
			hmatrix += tmp*sx.getCRS();

			// magnetic field y
			const OperatorType& sy = ModelBaseType::naturalOperator("sy", site, 0);
			tmp = modelParameters_.magneticFieldY[block[0]];
			hmatrix += tmp*sy.getCRS();

			// magnetic field z
			const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);
			tmp = modelParameters_.magneticFieldZ[block[0]];
			hmatrix += tmp*sz.getCRS();
		}
	}

	void setNaturalBasis(HilbertBasisType& basis,
	                     const VectorSizeType& block,
	                     bool truncated) const
	{
		assert(block.size()==1);
		assert(modelParameters_.orbitals == 1);
		HilbertStateType total = (1 << 2*modelParameters_.orbitals);

		basis.resize(total);
		for (SizeType a = 0; a< total; ++a) basis[a] = a;
		weedOutBasis(basis, truncated);
		if (modelParameters_.orbitals == 1 && basis.size() == 3) {
			basis[0] = 0; // empty
			basis[1] = 2; // down
			basis[2] = 1; // up
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
		VectorSizeType other(1, 0); // <----- conserve only 1 number, not two
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
			// other[1] = electronsUp; <--- conserve only other[0]
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


	ParametersModelTjAnisotropic<RealType, QnType>  modelParameters_;
	const SuperGeometryType& superGeometry_;
	SizeType offset_;
	SpinSquaredHelper<RealType,HilbertStateType> spinSquaredHelper_;
	SpinSquared<SpinSquaredHelper<RealType,HilbertStateType> > spinSquared_;
	HilbertBasisType basis_;
	VectorQnType qq_;
	VectorOperatorType creationMatrix_;
};	//class TjAnisotropic

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_Anisotropic_H

