/*
Copyright (c) 2009, 2017-2019, UT-Battelle, LLC
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

/*! \file GaugeSpin.h
 *
 *  An implementation of the Quantum Heisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_MODEL_GAUGESPIN_H
#define DMRG_MODEL_GAUGESPIN_H

#include <algorithm>
#include "ParametersGaugeSpin.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include "Utils.h"
#include "SuperOpHelperPlaquette.h"

namespace Dmrg {

template<typename ModelBaseType>
class GaugeSpin : public ModelBaseType {

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::RealType RealType;
	typedef	typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename ModelBaseType::VectorQnType VectorQnType;
	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef unsigned int long WordType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef typename PsimagLite::Vector<SizeType>::Type HilbertBasisType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;

	static const SizeType TWICE_THE_SPIN = 2;

	GaugeSpin(const SolverParamsType& solverParams,
	          InputValidatorType& io,
	          const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      superOpHelper_(nullptr)
	{}

	~GaugeSpin()
	{
		delete superOpHelper_;
		superOpHelper_ = nullptr;
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
	                                RealType) const
	{
		SizeType linSize = ModelBaseType::superGeometry().numberOfSites();
		assert(block.size() == 1);

		SizeType site = block[0];

		const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);

		if (modelParameters_.magneticFieldV.size() == linSize) {

			RealType tmp = modelParameters_.magneticFieldV[site];
			hmatrix += tmp*sz.getCRS();
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		const SizeType total = TWICE_THE_SPIN + 1;
		HilbertBasisType natBasis(total);
		for (SizeType i = 0; i < total; ++i) natBasis[i] = i;

		setSymmetryRelated(qns, natBasis);

		// Set the operators S^+_i in the natural basis
		SparseMatrixType tmpMatrix = findSplusMatrices(0, natBasis);

		typename OperatorType::Su2RelatedType su2related;

		OperatorType myOp(tmpMatrix,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  PairType(2, 2),
		                  -1,
		                  su2related);
		this->createOpsLabel("splus").push(myOp);
		// this->makeTrackable("splus");

		myOp.dagger();
		this->createOpsLabel("sminus").push(myOp);

		// Set the operators S^z_i in the natural basis
		tmpMatrix = findSzMatrices(0, natBasis);
		typename OperatorType::Su2RelatedType su2related2;
		OperatorType myOp2(tmpMatrix,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(2, 1),
		                   1.0/sqrt(2.0),
		                   su2related2);
		this->createOpsLabel("sz").push(myOp2);
		// this->makeTrackable("sz");

		// Set the operators S^x_i in the natural basis
		tmpMatrix = findSxMatrices(0, natBasis);
		typename OperatorType::Su2RelatedType su2related3;
		OperatorType myOp3(tmpMatrix,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   PairType(2, 1),
		                   1.0/sqrt(2.0),
		                   su2related3);
		this->createOpsLabel("sx").push(myOp3);
		this->makeTrackable("sx");

	}

	void fillModelLinks()
	{
		if (BasisType::useSu2Symmetry())
			err("SU(2): no longer supported\n");

		//auto lambda = []
		ModelTermType& plaquetteX = ModelBaseType::createTerm("PlaquetteX");

		OpForLinkType sx("sx");

		plaquetteX.push4(sx, 'N', sx, 'N', sx, 'N', sx, 'N');
	}

	const SuperOpHelperPlaquette& superOpHelper(const VectorSizeType& bigBlock,
	                                            const VectorSizeType& smallBlock,
	                                            ProgramGlobals::DirectionEnum dir) const
	{
		if (superOpHelper_) {
			delete superOpHelper_;
			superOpHelper_ = nullptr;
		}

		superOpHelper_ = new SuperOpHelperPlaquette(bigBlock, smallBlock, dir);
		return *superOpHelper_;
	}


private:

	//! Find S^+_site in the natural basis natBasis
	SparseMatrixType findSplusMatrices(SizeType site,
	                                   const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5*TWICE_THE_SPIN;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(TWICE_THE_SPIN);
		SizeType bits = 1 + ProgramGlobals::logBase2(TWICE_THE_SPIN);
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site*bitsForOneSite);

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);

			assert(ketsite == ket);
			SizeType brasite = ketsite + 1;
			if (brasite >= TWICE_THE_SPIN+1) continue;

			SizeType bra = ket & (~mask);
			assert(bra == 0);
			brasite <<= (site*bitsForOneSite);
			bra |= brasite;
			assert(bra == brasite);

			RealType m = ketsite - j;
			RealType x = j*(j+1)-m*(m+1);
			assert(x>=0);

			cm(ket,bra) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(SizeType site,
	                                const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5*TWICE_THE_SPIN;
		SizeType bitsForOneSite = utils::bitSizeOfInteger(TWICE_THE_SPIN);
		SizeType bits = ProgramGlobals::logBase2(TWICE_THE_SPIN) + 1;
		SizeType mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site*bitsForOneSite);

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);
			assert(ketsite == ket);
			RealType m = ketsite - j;
			cm(ket,ket) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findSxMatrices(SizeType site,
	                                const HilbertBasisType& natBasis) const
	{
		SparseMatrixType Splus_temp=findSplusMatrices(site,natBasis);
		SparseMatrixType Sminus_temp,Sx;
		transposeConjugate(Sminus_temp,Splus_temp);
		RealType tmp=0.5;

		Sx = tmp*Splus_temp;
		Sx += tmp*Sminus_temp;

		return Sx;
	}

	void setSymmetryRelated(VectorQnType& qns, const HilbertBasisType& basis) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;

		VectorSizeType other;
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0, 0);
			SizeType flavor = 1;
			qns[i] = QnType(false, other, jmpair, flavor);
		}
	}

	ParametersGaugeSpin<RealType, QnType> modelParameters_;
	mutable SuperOpHelperPlaquette* superOpHelper_;
}; // class GaugeSpin

} // namespace Dmrg
/*@}*/
#endif //DMRG_MODEL_GAUGESPIN_H

