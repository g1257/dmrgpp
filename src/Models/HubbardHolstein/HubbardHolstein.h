/*
Copyright (c) 2009-2015-2020, UT-Battelle, LLC
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
#include "CrsMatrix.h"
#include "Geometry/GeometryDca.h"
#include "HilbertSpaceHubbardHolstein.h"
#include "ModelBase.h"
#include "ParametersHubbardHolstein.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "SpinSquared.h"
#include "SpinSquaredHelper.h"
#include "VerySparseMatrix.h"
#include <cstdlib>
#include <numeric>

namespace Dmrg {
template <typename ModelBaseType> class HubbardHolstein : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
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
	typedef typename ModelHelperType::SparseElementType SparseElementType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelBaseType::HilbertBasisType HilbertBasisType;
	typedef typename HilbertBasisType::value_type HilbertState;
	typedef unsigned int long WordType;
	typedef HilbertSpaceHubbardHolstein<WordType> HilbertSpaceHubbardHolsteinWordType;
	typedef HilbertSpaceHubbardHolstein<HilbertState> HilbertSpaceHubbardHolsteinType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::MyBasis BasisType;
	typedef typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef ParametersHubbardHolstein<RealType, QnType> ParametersHubbardHolsteinType;
	typedef std::pair<SizeType, SizeType> PairType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const int FERMION_SIGN = -1;
	static SizeType const ORBITALS = 2;

	HubbardHolstein(const SolverParamsType& solverParams,
	                InputValidatorType& io,
	                const SuperGeometryType& geometry,
	                PsimagLite::String additional,
	                PsimagLite::String hdf5fileIfAny)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , isSsh_(additional == "SSH")
	    , isLrh_(additional == "LRH")
	    , oStruncActive_(false)
	    , wantsOneSiteTruncation_(0)
	{
		if (isSsh_) {
			PsimagLite::String warning("HubbardHolstein: ");
			warning += "SSH term in use.\n";
			std::cout << warning;
			std::cerr << warning;
		}
		if (isLrh_) {
			PsimagLite::String warning("HubbardHolstein: ");
			warning += "Long Range Holstein term in use.\n";
			std::cout << warning;
			std::cerr << warning;
		}

		restartHook(hdf5fileIfAny);
	}

	void print(std::ostream& os) const { operator<<(os, modelParameters_); }

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType& block,
	                                RealType time) const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		SizeType phonons = (oStruncActive_ || U_.rows() > 0)
		    ? modelParameters_.oStruncPhonons
		    : modelParameters_.numberphonons;
		if (phonons == 0)
			err("addDiagonalsInNaturalBasis fatal error when OSTRUNC active\n");

		SizeType n = block.size();
		HilbertBasisType natBasis;
		setBasis(natBasis, block, phonons);

		VectorSparseMatrixType cm;
		findAllMatrices(cm, natBasis, phonons);

		for (SizeType i = 0; i < n; ++i) {

			addInteractionFU(hmatrix, cm, block[i]);

			addInteractionFPhonon(hmatrix, cm, block[i], phonons);

			addPotentialFV(hmatrix, cm, block[i]);

			addPotentialPhononV(hmatrix, cm, block[i], phonons);
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);

		if (U_.rows() > 0) {
			U_.write("OneSiteTruncationU", io);
			io.write("OsTruncPhonons", modelParameters_.oStruncPhonons);
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);

		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& a = this->createOpsLabel("a");
		OpsLabelType& cx = this->createOpsLabel("cx");
		this->makeTrackable("c");
		const SizeType phonons = (U_.rows() > 0) ? modelParameters_.oStruncPhonons
		                                         : modelParameters_.numberphonons;
		if (phonons > 0) {
			this->makeTrackable("a");
			if (isSsh_)
				this->makeTrackable("cx");
		}

		HilbertBasisType natBasis;
		setBasis(natBasis, block, phonons);
		setSymmetryRelated(qns, natBasis);

		//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
		SparseMatrixType nmatrix;
		SparseMatrixType tmpMatrix;
		for (SizeType sigma = 0; sigma < 2; ++sigma) {
			tmpMatrix = findOperatorMatrices(sigma, natBasis);

			if (sigma == 0)
				nmatrix = n(tmpMatrix);
			else
				nmatrix += n(tmpMatrix);

			int asign = 1;
			if (sigma > 0)
				asign = 1;
			typename OperatorType::Su2RelatedType su2related;
			if (sigma == 0) {
				su2related.source.push_back(site);
				su2related.source.push_back(site + 1);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(-1);
				su2related.offset = 1;
			}

			transformByU(tmpMatrix);

			SparseMatrixType tmpMatrix2;
			transposeConjugate(tmpMatrix2, tmpMatrix);

			OperatorType myOp(tmpMatrix2,
			                  ProgramGlobals::FermionOrBosonEnum::FERMION,
			                  typename OperatorType::PairType(1, 1 - sigma),
			                  asign,
			                  su2related);

			c.push(myOp, (sigma == 0) ? "up" : "down");
		}

		OpsLabelType& n = this->createOpsLabel("n");
		typename OperatorType::Su2RelatedType su2relatedA;

		transformByU(nmatrix);

		OperatorType myOp(nmatrix,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  typename OperatorType::PairType(0, 0),
		                  1,
		                  su2relatedA);

		n.push(myOp);

		SizeType iup = 0;
		SizeType idown = 1;
		{
			OpsLabelType& splus = this->createOpsLabel("splus");
			OpsLabelType& sminus = this->createOpsLabel("sminus");
			SparseMatrixType t1MatrixUp = findOperatorMatrices(iup, natBasis);
			SparseMatrixType t1MatrixDown = findOperatorMatrices(idown, natBasis);
			PsimagLite::Matrix<typename SparseMatrixType::value_type> t1
			    = multiplyTc(t1MatrixUp, t1MatrixDown);
			typename OperatorType::Su2RelatedType su2Related;
			SparseMatrixType spmatrix(t1);
			splus.push(OperatorType(spmatrix,
			                        ProgramGlobals::FermionOrBosonEnum::BOSON,
			                        typename OperatorType::PairType(0, 0),
			                        1.0,
			                        su2Related));
			SparseMatrixType smmatrix;
			transposeConjugate(smmatrix, spmatrix);
			sminus.push(OperatorType(smmatrix,
			                         ProgramGlobals::FermionOrBosonEnum::BOSON,
			                         typename OperatorType::PairType(0, 0),
			                         1.0,
			                         su2Related));
		}
		{
			OpsLabelType& sz = this->createOpsLabel("sz");
			SparseMatrixType t1MatrixUp = findOperatorMatrices(iup, natBasis);
			SparseMatrixType t1MatrixDown = findOperatorMatrices(idown, natBasis);
			PsimagLite::Matrix<SparseElementType> t1
			    = multiplyTc(t1MatrixUp, t1MatrixUp);
			PsimagLite::Matrix<SparseElementType> t2
			    = multiplyTc(t1MatrixDown, t1MatrixDown);
			t1 = 0.5 * (t1 - t2);
			SparseMatrixType szmatrix(t1);
			typename OperatorType::Su2RelatedType su2Related;
			sz.push(OperatorType(szmatrix,
			                     ProgramGlobals::FermionOrBosonEnum::BOSON,
			                     typename OperatorType::PairType(0, 0),
			                     1.0,
			                     su2Related));
		}

		{
			OpsLabelType& doubleOcc = this->createOpsLabel("double");
			OpsLabelType& doubleOccM = this->createOpsLabel("doubleM");
			SparseMatrixType t1MatrixUp = findOperatorMatrices(iup, natBasis);
			SparseMatrixType t1MatrixDown = findOperatorMatrices(idown, natBasis);
			PsimagLite::Matrix<SparseElementType> t1
			    = multiplyTc(t1MatrixUp, t1MatrixUp);
			PsimagLite::Matrix<SparseElementType> t2
			    = multiplyTc(t1MatrixDown, t1MatrixDown);
			SparseMatrixType nMatrixUp(t1);
			t1 = t1 + t2;
			SparseMatrixType nMatrixDown(t2);
			PsimagLite::Matrix<SparseElementType> t3
			    = multiplyTc(nMatrixDown, nMatrixUp);

			t1 = t1 + (-2.0) * t3;

			SparseMatrixType tmp4(t1);
			typename OperatorType::Su2RelatedType su2Related;
			doubleOccM.push(OperatorType(tmp4,
			                             ProgramGlobals::FermionOrBosonEnum::BOSON,
			                             typename OperatorType::PairType(0, 0),
			                             1.0,
			                             su2Related));

			SparseMatrixType tmp5(t3);
			doubleOcc.push(OperatorType(tmp5,
			                            ProgramGlobals::FermionOrBosonEnum::BOSON,
			                            typename OperatorType::PairType(0, 0),
			                            1.0,
			                            su2Related));
		}

		if (phonons == 0)
			return; //<<--- EARLY EXIT

		tmpMatrix = findPhononadaggerMatrix(natBasis, phonons);

		typename OperatorType::Su2RelatedType su2related2;
		su2related2.source.push_back(site * 2);
		su2related2.source.push_back(site * 2 + 1);
		su2related2.source.push_back(site * 2);
		su2related2.transpose.push_back(-1);
		su2related2.transpose.push_back(-1);
		su2related2.transpose.push_back(1);
		su2related2.offset = 1;

		transformByU(tmpMatrix);

		SparseMatrixType tmpMatrix2;
		transposeConjugate(tmpMatrix2, tmpMatrix);

		OperatorType myOp2(tmpMatrix2,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(2, 2),
		                   -1,
		                   su2related2);
		a.push(myOp2);

		{
			OpsLabelType& disp = this->createOpsLabel("disp");
			SparseMatrixType tmp1 = findPhononadaggerMatrix(natBasis, phonons);
			SparseMatrixType tmp2;
			transposeConjugate(tmp2, tmp1);
			tmp2 += tmp1;
			typename OperatorType::Su2RelatedType su2Related2;
			disp.push(OperatorType(tmp2,
			                       ProgramGlobals::FermionOrBosonEnum::BOSON,
			                       typename OperatorType::PairType(0, 0),
			                       1.0,
			                       su2Related2));
		}

		if (!isSsh_ && !isLrh_)
			return; //<<--- EARLY EXIT

		if (isSsh_) {

			// Set the operators c_(i,sigma} * x_i in the natural basis

			for (SizeType sigma = 0; sigma < 2; ++sigma) {
				tmpMatrix = findSSHMatrices(sigma, natBasis, phonons);
				int asign = 1;
				if (sigma > 0)
					asign = 1;
				typename OperatorType::Su2RelatedType su2related3;
				if (sigma == 0) {
					su2related3.source.push_back(site);
					su2related3.source.push_back(site + 1);
					su2related3.transpose.push_back(-1);
					su2related3.transpose.push_back(-1);
					su2related3.offset = 1;
				}

				transformByU(tmpMatrix);

				SparseMatrixType tmpMatrix2;
				transposeConjugate(tmpMatrix2, tmpMatrix);

				OperatorType myOp3(tmpMatrix2,
				                   ProgramGlobals::FermionOrBosonEnum::FERMION,
				                   typename OperatorType::PairType(1, 1 - sigma),
				                   asign,
				                   su2related3);

				cx.push(myOp3, (sigma == 0) ? "up" : "down");
			}
		} else if (isLrh_) {
			this->makeTrackable("disp");
		}
	}

	void fillModelLinks()
	{
		ModelTermType& hopf = ModelBaseType::createTerm("HoppingFermionic");

		OpForLinkType cup("c", 0);
		hopf.push(cup, 'C', cup, 'N');

		OpForLinkType cdown("c", 1);
		hopf.push(cdown, 'C', cdown, 'N');

		const SizeType phonons = modelParameters_.numberphonons;
		if (phonons > 0) {
			ModelTermType& hopb = ModelBaseType::createTerm("HoppingBosonic");

			OpForLinkType a("a");
			hopb.push(a, 'C', a, 'N');
		}

		if (isSsh_) {
			ModelTermType& hopSsh = ModelBaseType::createTerm("HoppingSSH");

			OpForLinkType cx0("cx", 0);
			OpForLinkType cx1("cx", 1);

			auto modifier = [](ComplexOrRealType& value) { value *= (-1.0); };

			hopSsh.push(cup, 'C', cx0, 'N');

			hopSsh.push(cx0, 'C', cup, 'N', modifier);

			hopSsh.push(cdown, 'C', cx1, 'N');

			hopSsh.push(cx1, 'C', cdown, 'N', modifier);

		} else if (isLrh_) {
			ModelTermType& hopLrh = ModelBaseType::createTerm("LongRangeH");

			OpForLinkType x("disp");

			OpForLinkType n("n");

			hopLrh.push(x, 'N', n, 'N');

		} else {
			return;
		}
	}

	void announce(PsimagLite::String str) const
	{
		const PsimagLite::String msg("finite loop");
		const SizeType l = msg.length();

		if (str.substr(0, l) != msg)
			return;

		PsimagLite::Vector<PsimagLite::String>::Type tokens;
		PsimagLite::split(tokens, str, ";");
		if (tokens.size() != 2)
			err("Model::announce()\n");

		wantsOneSiteTruncation_ = PsimagLite::atoi(tokens[1]);
	}

	void oneSiteTruncationUpdate(OutputFileOrNot& ioOut, const MatrixType& U, SizeType start)
	{
		bool firstCall = (U_.rows() == 0);

		notReallySortU(U, start);
		ModelBaseType::oneSiteTruncationUpdate(ioOut, U, start);

		static const bool verbose = true;
		if (verbose) {
			std::cout << "U UPDATED\n";
			std::cout << U_;
			std::cout << "--------\n";
		}

		if (firstCall) {
			ioOut.write(U_, "OneSiteTruncationU");
			ioOut.write(modelParameters_.oStruncPhonons, "OsTruncPhonons");
		} else {
			ioOut.overwrite(U_, "OneSiteTruncationU");
			ioOut.write(modelParameters_.oStruncPhonons,
			            "OsTruncPhonons",
			            PsimagLite::IoNgSerializer::ALLOW_OVERWRITE);
		}
	}

	// virtual override
	SizeType
	setOperatorMatrices(VectorOperatorType& ops, VectorQnType& qm, const BlockType& block) const
	{
		oStruncActive_ = false;

		const bool b1 = (modelParameters_.oStruncPhonons == 0);

		assert(block.size() == 1);

		const bool b2 = (modelParameters_.oStruncSite != block[0]);

		const bool b3 = (wantsOneSiteTruncation_ == 0);

		if (b1 || b2 || b3)
			return ModelBaseType::setOperatorMatrices(ops, qm, block);

		oStruncActive_ = true;
		HilbertBasisType natBasis;
		setBasis(natBasis, block, modelParameters_.oStruncPhonons);
		setSymmetryRelated(qm, natBasis);

		const SizeType oneSiteTruncSize = natBasis.size();

		const SizeType ind = 0;
		for (SizeType sigma = 0; sigma < 2; ++sigma) {
			SparseMatrixType tmpMatrix = findOperatorMatrices(sigma, natBasis);
			int asign = 1;
			if (sigma > 0)
				asign = 1;
			typename OperatorType::Su2RelatedType su2related;
			if (sigma == 0) {
				su2related.source.push_back(ind);
				su2related.source.push_back(ind + 1);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(-1);
				su2related.offset = 1;
			}

			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::FERMION,
			                  typename OperatorType::PairType(1, 1 - sigma),
			                  asign,
			                  su2related);

			ops.push_back(myOp);
		}

		// operator n should NOT be pushed because it isn't tracked

		if (modelParameters_.oStruncPhonons == 0)
			return oneSiteTruncSize;

		SparseMatrixType tmpMatrix
		    = findPhononadaggerMatrix(natBasis, modelParameters_.oStruncPhonons);

		typename OperatorType::Su2RelatedType su2related2;
		su2related2.source.push_back(ind * 2);
		su2related2.source.push_back(ind * 2 + 1);
		su2related2.source.push_back(ind * 2);
		su2related2.transpose.push_back(-1);
		su2related2.transpose.push_back(-1);
		su2related2.transpose.push_back(1);
		su2related2.offset = 1;
		OperatorType myOp2(tmpMatrix,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(2, 2),
		                   -1,
		                   su2related2);
		ops.push_back(myOp2);

		if (!isSsh_ && !isLrh_)
			return oneSiteTruncSize; //<<--- EARLY EXIT

		// Set the operators c_(i,sigma} * x_i in the natural basis

		for (SizeType sigma = 0; sigma < 2; ++sigma) {
			tmpMatrix
			    = findSSHMatrices(sigma, natBasis, modelParameters_.oStruncPhonons);
			int asign = 1;
			if (sigma > 0)
				asign = 1;
			typename OperatorType::Su2RelatedType su2related3;
			if (sigma == 0) {
				su2related3.source.push_back(ind);
				su2related3.source.push_back(ind + 1);
				su2related3.transpose.push_back(-1);
				su2related3.transpose.push_back(-1);
				su2related3.offset = 1;
			}

			transformByU(tmpMatrix);

			OperatorType myOp3(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::FERMION,
			                   typename OperatorType::PairType(1, 1 - sigma),
			                   asign,
			                   su2related3);

			ops.push_back(myOp3);
		}

		return oneSiteTruncSize;
	}

private:

	//! find all states in the natural basis for a block of n sites
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	static void setBasis(HilbertBasisType& basis, const VectorSizeType& block, SizeType phonons)
	{
		SizeType n = block.size();
		HilbertState total = 4 * (phonons + 1);
		total = pow(total, n);

		basis.resize(total);

		SizeType counter = 0;
		SizeType npPlusOne = phonons + 1;
		for (SizeType i = 0; i < 4; ++i) {
			for (SizeType b = 0; b < npPlusOne; ++b) {
				basis[counter++] = b * 4 + i;
			}
		}

		assert(counter == total);

		SizeType sum
		    = std::accumulate(basis.begin(), basis.end(), static_cast<SizeType>(0));
		if (sum != total * (total - 1) / 2)
			err("Could not set up basis\n");
	}

	//! Find a^+ in the natural basis natBasis
	SparseMatrixType findPhononadaggerMatrix(const HilbertBasisType& natBasis,
	                                         SizeType phonons) const
	{
		const SizeType total = natBasis.size();

		MatrixType cm(total, total);

		for (SizeType ii = 0; ii < total; ++ii) {

			typename HilbertSpaceHubbardHolsteinType::HilbertState ket = natBasis[ii];

			SizeType nphon = SizeType(HilbertSpaceHubbardHolsteinType::getP(ket));
			assert(nphon <= phonons);
			if (nphon == phonons)
				continue;

			typename HilbertSpaceHubbardHolsteinType::HilbertState bra = ket;
			HilbertSpaceHubbardHolsteinType::createP(bra);

			const int jj = PsimagLite::indexOrMinusOne(natBasis, bra);
			if (jj < 0)
				err("findOperatorMatrices\n");

			const RealType x = HilbertSpaceHubbardHolsteinType::getP(bra);

			assert(x >= 0);

			cm(ii, jj) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp, cm);
		transposeConjugate(operatorMatrix, temp);
		return operatorMatrix;
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceHubbardHolsteinType::HilbertState ket) const
	{
		const SizeType ndown = HilbertSpaceHubbardHolsteinType::electronsWithGivenSpin(
		    ket, HilbertSpaceHubbardHolsteinType::SpinEnum::SPIN_DOWN);

		return (ndown % 2 == 0) ? 1 : FERMION_SIGN;
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(SizeType sigma,
	                                      const HilbertBasisType& natBasis) const
	{
		const SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n, n);

		for (SizeType ii = 0; ii < n; ++ii) {

			typename HilbertSpaceHubbardHolsteinType::HilbertState ket = natBasis[ii];

			const typename HilbertSpaceHubbardHolsteinType::SpinEnum spin = (sigma == 0)
			    ? HilbertSpaceHubbardHolsteinType::SpinEnum::SPIN_UP
			    : HilbertSpaceHubbardHolsteinType::SpinEnum::SPIN_DOWN;

			SizeType neSpin
			    = HilbertSpaceHubbardHolsteinType::electronsWithGivenSpin(ket, spin);
			if (neSpin == 0) {
				typename HilbertSpaceHubbardHolsteinType::HilbertState bra = ket;
				HilbertSpaceHubbardHolsteinType::createF(bra, sigma);
				const int jj = PsimagLite::indexOrMinusOne(natBasis, bra);
				if (jj < 0)
					err("findOperatorMatrices\n");
				cm(ii, jj) = sign(ket);
			}
		}

		SparseMatrixType creationMatrix(cm);
		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp, cm);
		transposeConjugate(creationMatrix, temp);

		return creationMatrix;
	}

	SparseMatrixType
	findSSHMatrices(SizeType sigma, const HilbertBasisType& natBasis, SizeType phonons) const
	{
		SparseMatrixType csigma_temp = findOperatorMatrices(sigma, natBasis);
		SparseMatrixType a_temp = findPhononadaggerMatrix(natBasis, phonons);
		SparseMatrixType x_temp = displacementOp(a_temp);
		SparseMatrixType csigma_a;
		multiply(csigma_a, csigma_temp, x_temp);
		return csigma_a;
	}

	void findAllMatrices(VectorSparseMatrixType& vm,
	                     const HilbertBasisType& natBasis,
	                     SizeType phonons) const
	{
		for (SizeType sigma = 0; sigma < 2; ++sigma) {
			SparseMatrixType m = findOperatorMatrices(sigma, natBasis);
			transformByU(m);
			vm.push_back(m);
		}

		if (phonons == 0)
			return;
		SparseMatrixType m = findPhononadaggerMatrix(natBasis, phonons);
		transformByU(m);
		vm.push_back(m);
	}

	void setSymmetryRelated(VectorQnType& qns, const HilbertBasisType& basis) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType, SizeType> PairType;
		VectorSizeType other(2, 0);
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = PairType(0, 0);

			// nup
			SizeType electronsUp
			    = HilbertSpaceHubbardHolsteinType::electronsWithGivenSpin(
			        basis[i], HilbertSpaceHubbardHolsteinType::SpinEnum::SPIN_UP);
			// ndown
			SizeType electronsDown
			    = HilbertSpaceHubbardHolsteinType::electronsWithGivenSpin(
			        basis[i], HilbertSpaceHubbardHolsteinType::SpinEnum::SPIN_DOWN);

			SizeType electrons = electronsUp + electronsDown;
			other[0] = electrons;
			other[1] = electronsUp;
			bool sign = electrons & 1;
			qns[i] = QnType(sign, other, jmpair, electrons);
		}
	}

	void addPotentialFV(SparseMatrixType& hmatrix,
	                    const VectorSparseMatrixType& cm,
	                    SizeType actualIndexOfSite) const
	{
		SparseMatrixType nup = n(cm[0]); // spin up
		SparseMatrixType ndown = n(cm[1]); // spin down

		SizeType linSize = ModelBaseType::superGeometry().numberOfSites();
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialFV.size());
		hmatrix += modelParameters_.potentialFV[iUp] * nup;
		SizeType iDown = actualIndexOfSite + linSize;
		assert(iDown < modelParameters_.potentialFV.size());
		hmatrix += modelParameters_.potentialFV[iDown] * ndown;
	}

	void addPotentialPhononV(SparseMatrixType& hmatrix,
	                         const VectorSparseMatrixType& cm,
	                         SizeType actualIndexOfSite,
	                         SizeType phonons) const
	{
		if (phonons == 0)
			return;
		assert(2 < cm.size());
		SparseMatrixType nphon = n(cm[2]);
		SizeType iUp = actualIndexOfSite;
		assert(iUp < modelParameters_.potentialPV.size());
		hmatrix += modelParameters_.potentialPV[iUp] * nphon;
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger, c);
		multiply(tmpMatrix, c, cdagger);

		return tmpMatrix;
	}

	SparseMatrixType displacementOp(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix = c;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger, c);
		tmpMatrix += cdagger;
		return tmpMatrix;
	}

	//! Term is U \sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
	void addInteractionFU(SparseMatrixType& hmatrix,
	                      const VectorSparseMatrixType& cm,
	                      SizeType actualSite) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType m1 = cm[0]; // spin up
		SparseMatrixType m2 = cm[1]; // spin down

		multiply(tmpMatrix, n(m1), n(m2));
		assert(actualSite < modelParameters_.hubbardFU.size());
		hmatrix += modelParameters_.hubbardFU[actualSite] * tmpMatrix;
	}

	//! Term is lambda\sum_{\alpha} (n_{i\alpha} -1) x_{i}
	void addInteractionFPhonon(SparseMatrixType& hmatrix,
	                           const VectorSparseMatrixType& cm,
	                           SizeType actualSite,
	                           SizeType phonons) const
	{
		if (phonons == 0)
			return;
		SparseMatrixType tmpMatrix;
		SparseMatrixType m = n(cm[0]); // spin up
		SparseMatrixType m2 = n(cm[1]); // spin down
		m += m2;
		assert(2 < cm.size());
		SparseMatrixType x = displacementOp(cm[2]);

		multiply(tmpMatrix, m, x);
		tmpMatrix.checkValidity();
		assert(actualSite < modelParameters_.lambdaFP.size());
		hmatrix += modelParameters_.lambdaFP[actualSite] * tmpMatrix;
	}

	void transformByU(SparseMatrixType& m) const
	{
		if (U_.rows() == 0)
			return;

		if (U_.rows() != m.rows())
			err("HubbardHolstein::transformByU(): wrong sizes\n");

		MatrixType mdense;
		crsMatrixToFullMatrix(mdense, m);
		rotate(mdense, U_);
		fullMatrixToCrsMatrix(m, mdense);
	}

	void notReallySortU(const MatrixType& U, SizeType start)
	{
		const SizeType site = 0;
		BlockType block(1, site);
		HilbertBasisType natBasis;
		const SizeType phonons = modelParameters_.oStruncPhonons;
		setBasis(natBasis, block, phonons);
		VectorQnType qns;
		setSymmetryRelated(qns, natBasis);
		BasisType::notReallySortU(U_, U, qns, start);
	}

	void restartHook(PsimagLite::String restartFilename)
	{
		if (restartFilename == "")
			return;

		PsimagLite::IoSelector::In io(restartFilename);
		try {
			io.read(U_, "OneSiteTruncationU");
			std::cout << "OneSiteTruncationU = " << U_.rows() << " x " << U_.cols();
			std::cout << " read from " << restartFilename << "\n";
		} catch (...) { }

		if (U_.rows() > 0) {
			io.read(modelParameters_.oStruncPhonons, "OsTruncPhonons");
			std::cout << "OsTruncPhonons set to " << modelParameters_.oStruncPhonons
			          << "\n";
		}
	}

	ParametersHubbardHolsteinType modelParameters_;
	bool isSsh_;
	bool isLrh_;
	mutable bool oStruncActive_;
	mutable SizeType wantsOneSiteTruncation_;
	mutable MatrixType U_;
}; // class HubbardHolstein
} // namespace Dmrg
/*@}*/
#endif
