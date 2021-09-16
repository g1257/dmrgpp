/*
Copyright (c) 2009-2017-2018 UT-Battelle, LLC
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

/*! \file Kitaev.h
 *
 *  An implementation of the Kitaev model (started March 2018)
 *
 */

#ifndef DMRG_KITAEV_H
#define DMRG_KITAEV_H

#include <algorithm>
#include "ModelBase.h"
#include "ParametersKitaev.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramGlobals.h"
#include "Utils.h"
#include "Complex.h"

namespace Dmrg {

template<typename ModelBaseType>
class Kitaev : public ModelBaseType {

	static const int TWICE_THE_SPIN = 1;

	enum class InternalDir {DIR_X, DIR_Y, DIR_Z, DIR_PLUS, DIR_MINUS};

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef	typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef unsigned int long WordType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef typename PsimagLite::Vector<unsigned int long>::Type HilbertBasisType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	Kitaev(const SolverParamsType& solverParams,
	       InputValidatorType& io,
	       const SuperGeometryType& geometry,
	       PsimagLite::String additional)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      extended_(additional.length() > 7 && additional.substr(0, 8) == "Extended"),
	      withGammas_(additional.length() > 9 && additional.substr(0, 10) == "WithGammas"),
	      withCharge_(additional.length() > 9 && (additional.substr(8, 10) == "WithCharge" ||
	                                              additional.substr(10, 10) == "WithCharge" ||
	                                              additional.substr(0, 10) == "WithCharge"))
	{
		if (withCharge_ and TWICE_THE_SPIN != 1)
			err("Kitaev: Charged model only for s=1/2\n");

		SizeType n = geometry.numberOfSites();
		SizeType mx = modelParameters_.magneticFieldX.size();
		SizeType my = modelParameters_.magneticFieldY.size();
		SizeType mz = modelParameters_.magneticFieldZ.size();

		if (mx > 0 && mx != n) {
			PsimagLite::String msg("Kitaev: If provided, ");
			msg += " MagneticFieldX must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}
		if (mz > 0 && mz != n) {
			PsimagLite::String msg("Kitaev: If provided, ");
			msg += " MagneticFieldZ must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}
		if (my > 0 && my != n) {
			PsimagLite::String msg("Kitaev: If provided, ");
			msg += " MagneticFieldY must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}
		if (BasisType::useSu2Symmetry())
			err("Kitaev does not have SU(2) symmetry\n");
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType)  const
	{
		SizeType linSize = ModelBaseType::superGeometry().numberOfSites();
		bool hasX = (modelParameters_.magneticFieldX.size() == linSize);
		bool hasY = (modelParameters_.magneticFieldY.size() == linSize);
		bool hasZ = (modelParameters_.magneticFieldZ.size() == linSize);

		SizeType n=block.size();
		for (SizeType i = 0; i < n; ++i) {
			SizeType site = block[i];

			if (hasX) {
				// magnetic field x
				const OperatorType& sx = ModelBaseType::naturalOperator("sx", site, 0);
				RealType tmp = modelParameters_.magneticFieldX[block[0]];
				hmatrix += tmp*sx.getCRS();
			}

			if (hasY) {
				// magnetic field y
				const OperatorType& sy = ModelBaseType::naturalOperator("sy", site, 0);
				RealType tmp = modelParameters_.magneticFieldY[block[0]];
				hmatrix += tmp*sy.getCRS();
			}

			if (hasZ) {
				// magnetic field z
				const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);
				RealType tmp = modelParameters_.magneticFieldZ[block[0]];
				hmatrix += tmp*sz.getCRS();
			}
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;

		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis, block.size());

		OpsLabelType& sx = this->createOpsLabel("sx");
		OpsLabelType& sybar = this->createOpsLabel("sybar");
		OpsLabelType& sz = this->createOpsLabel("sz");
		OpsLabelType& splus = this->createOpsLabel("splus");
		OpsLabelType& sminus = this->createOpsLabel("sminus");

		this->makeTrackable("sx");
		this->makeTrackable("sybar"); // Sybar = iSy
		this->makeTrackable("sz");

		for (SizeType i=0;i<block.size();i++) {

			typename OperatorType::Su2RelatedType su2related;

			// Set the operators S^x_i in the natural basis
			tmpMatrix = findSdirMatrices(i, natBasis, InternalDir::DIR_X);
			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(0, 0),
			                  1.0,
			                  su2related);
			sx.push(myOp);

			// Set the operators S^y_i in the natural basis
			tmpMatrix = findSdirMatrices(i, natBasis, InternalDir::DIR_Y);
			OperatorType myOp2(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			sybar.push(myOp2); // Sybar = iSy

			// Set the operators S^z_i in the natural basis
			tmpMatrix = findSdirMatrices(i, natBasis, InternalDir::DIR_Z);
			OperatorType myOp3(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			sz.push(myOp3);

			tmpMatrix = findSdirMatrices(i, natBasis, InternalDir::DIR_PLUS);
			OperatorType myOp4(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);

			splus.push(myOp4);

			OperatorType myOp5 = myOp4;
			myOp5.conjugate();
			sminus.push(myOp5);
		}

		if (!withCharge_) return; // <<=== EARLY EXIT HERE

		OpsLabelType& destructionC = this->createOpsLabel("c");
		this->makeTrackable("c");
		for (SizeType spin = 0; spin < 2; ++spin) {
			typename OperatorType::Su2RelatedType su2related;

			tmpMatrix = findDestructionC(natBasis, spin);
			OperatorType myOpC(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::FERMION,
			                   PairType(0, 0),
			                   1.0,
			                   su2related);
			destructionC.push(myOpC);
		}
	}

	void fillModelLinks()
	{

		if (withCharge_) {
			ModelTermType& hop = ModelBaseType::createTerm("hopping");//(A)

			OpForLinkType cup("c", 0); // (B)
			hop.push(cup, 'N', cup, 'C', typename ModelTermType::Su2Properties(1, 1, 0)); // (C)

			OpForLinkType cdown("c", 1); // (D)
			hop.push(cdown, 'N', cdown, 'C', typename ModelTermType::Su2Properties(1, -1, 1));
		}

		VectorStringType labels = {"sx", "sybar", "sz"};

		for (SizeType i = 0; i < labels.size(); ++i) {
			OpForLinkType smu(labels[i]);
			auto modif = [i] (ComplexOrRealType& value)
			{
				ComplexOrRealType x = (i == 1) ? -1 : 1; // sybar * sybar = -sy*sy
				value *= x;
			};

			ModelBaseType::createTerm(labels[i] + labels[i], false).push(smu, 'N', smu, 'N', modif);
		}

		OpForLinkType sx("sx");
		OpForLinkType sybar("sybar");

		const bool wantsHermit = true;

		typename MatrixType::value_type dummy = 0.0;

		if (extended_) {
			createTermSxSy(sx, sybar, wantsHermit, dummy);
			return; // <<---- EARLY EXIT HERE
		}

		if (!withGammas_) return; // <<---- EARLY EXIT HERE

		OpForLinkType sz("sz");

		createTermSySz(sybar, sz, wantsHermit, dummy);

		ModelBaseType::createTerm("sxsz").push(sx, 'N', sz, 'N');
		ModelTermType& szsx = ModelBaseType::createTerm("szsx", wantsHermit, "sxsz");
		szsx.push(sz, 'N', sx, 'N');

		createTermSxSy(sx, sybar, wantsHermit, dummy);
	}

private:

	void createTermSySz(const OpForLinkType&,
	                    const OpForLinkType&,
	                    bool,
	                    RealType)
	{
		PsimagLite::String str = "needs useComplex in SolverOptions in the input\n";
		err("FATAL: createTermSySz(): This Kitaev variant needs useComplex " + str);
	}

	void createTermSySz(const OpForLinkType& sybar,
	                    const OpForLinkType& sz,
	                    bool wantsHermit,
	                    std::complex<RealType>)
	{
		auto modifMinusSqrtMinusOne = [](ComplexOrRealType& value)
		{
			value *= std::complex<RealType>(0, -1);
		};

		ModelBaseType::createTerm("sysz").push(sybar, 'N', sz, 'N', modifMinusSqrtMinusOne);
		ModelTermType& szsy = ModelBaseType::createTerm("szsy", wantsHermit, "sysz");
		szsy.push(sz, 'N', sybar, 'N', modifMinusSqrtMinusOne);
	}

	void createTermSxSy(const OpForLinkType&,
	                    const OpForLinkType&,
	                    bool,
	                    RealType) const
	{
		PsimagLite::String str = "needs useComplex in SolverOptions in the input\n";
		err("createTermSxSy(): This Kitaev variant " + str);
	}

	void createTermSxSy(const OpForLinkType& sx,
	                    const OpForLinkType& sybar,
	                    bool wantsHermit,
	                    std::complex<RealType>) const
	{
		auto modifMinusSqrtMinusOne = [](ComplexOrRealType& value)
		{
			value *= std::complex<RealType>(0, -1);
		};

		ModelBaseType::createTerm("sxsy").push(sx, 'N', sybar, 'N', modifMinusSqrtMinusOne);
		ModelTermType& sysx = ModelBaseType::createTerm("sysx", wantsHermit, "sxsy");
		sysx.push(sybar, 'N', sx, 'N', modifMinusSqrtMinusOne);

		SizeType site = 0;
		BlockType block(1, site);
		HilbertBasisType natBasis;

		setBasis(natBasis, block);

		OpsLabelType& sy = this->createOpsLabel("sy");
		SparseMatrixType tmpMatrix = findSdirMatrices(0, natBasis, InternalDir::DIR_Y);
		tmpMatrix *= std::complex<RealType>(0, -1); // Sy = -iSybar
		typename OperatorType::Su2RelatedType su2related;
		OperatorType myOp2(tmpMatrix,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(0, 0),
		                   1.0,
		                   su2related);

		sy.push(myOp2); // Sybar = iSy
	}

	//! find all states in the natural basis for a block of n sites
	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		const SizeType total1 = (withCharge_) ? 3
		                                      : TWICE_THE_SPIN + 1;

		SizeType total = utils::powUint(total1, block.size());

		basis.resize(total);
		for (SizeType i = 0; i < total; ++i) basis[i] = i;
	}

	SizeType logBase2(SizeType x) const
	{
		SizeType counter = 0;
		while (x > 0) {
			x >>= 1;
			counter++;
		}

		return (counter == 0) ? counter : counter - 1;
	}

	SparseMatrixType findSdirMatrices(SizeType,// site,
	                                  const HilbertBasisType& natBasis,
	                                  InternalDir dir) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total, total);
		SizeType offset = (withCharge_) ? 1 : 0;

		if (withCharge_) {
			if (total != 3)
				err("findSdirMatrices: with charge Hilbert space should be 3\n");
		} else {
			if (total != TWICE_THE_SPIN + 1 || TWICE_THE_SPIN != 1)
				err("findSdirMatrices: only for spin 1/2 AND block of one site\n");
		}

		if (dir == InternalDir::DIR_X) {
			cm(0 + offset, 1 + offset) = cm(1 + offset, 0 + offset) = 0.5;
		} else if (dir == InternalDir::DIR_Y) {
			cm(0 + offset, 1 + offset) = 0.5;
			cm(1 + offset, 0 + offset) = -0.5;
		} else if (dir == InternalDir::DIR_Z) {
			cm(0 + offset, 0 + offset) = 0.5;
			cm(1 + offset, 1 + offset) = -0.5;
		} else if (dir == InternalDir::DIR_PLUS) {
			cm(0 + offset, 1 + offset) = 1.0;
		} else if (dir == InternalDir::DIR_MINUS) {
			cm(1 + offset, 0 + offset) = 1.0;
		} else {
			assert(false);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	// destruction
	SparseMatrixType findDestructionC(const HilbertBasisType&, SizeType spin) const
	{
		assert(withCharge_);
		assert(spin == 0 || spin == 1);

		MatrixType m(3, 3);

		if (spin == 0)
			m(0, 1) = 1;
		else
			m(0, 2) = 1;
		SparseMatrixType operatorMatrix(m);
		return operatorMatrix;
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		assert(n == 1);

		SizeType nbasis = basis.size();

		qns.resize(nbasis, QnType::zero());

		if (!withCharge_) return; // <<---- EARLY EXIT HERE

		setSymmetryRelatedWithCharge(qns, basis);
	}

	void setSymmetryRelatedWithCharge(VectorQnType& qns,
	                                  const HilbertBasisType& basis) const
	{
		typedef std::pair<SizeType,SizeType> PairType;

		VectorSizeType other(1);
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0, 0);

			assert(basis[i] == 0 || basis[i] == 1 || basis[i] == 2);

			other[0] = (basis[i] == 0) ? 0 : 1;

			bool sign = other[0] & 1;
			qns[i] = QnType(sign, other, jmpair, other[0]);
		}
	}

	ParametersKitaev<RealType, QnType>  modelParameters_;
	const bool extended_;
	const bool withGammas_;
	const bool withCharge_;
}; // class Kitaev

} // namespace Dmrg
/*@}*/
#endif //DMRG_KITAEV_H

