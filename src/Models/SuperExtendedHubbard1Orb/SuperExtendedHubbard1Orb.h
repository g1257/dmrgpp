/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file SuperExtendedHubbard1Orb.h
 *
 *  Hubbard + V^n_{ij} n_i n_j + V^S_{ij} S_i S_j + V^P_{ij} P^dag_i P_j
 *  + interaction-assisted hopping
 *
 */
#ifndef EXTENDED_SUPER_HUBBARD_1ORB_H
#define EXTENDED_SUPER_HUBBARD_1ORB_H
#include "../Models/ExtendedHubbard1Orb/ExtendedHubbard1Orb.h"

namespace Dmrg {
//! Super Extended Hubbard for DMRG solver, uses ModelHubbard by containment
template<typename ModelBaseType>
class ExtendedSuperHubbard1Orb : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef ExtendedHubbard1Orb<ModelBaseType> ExtendedHubbard1OrbType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ExtendedHubbard1OrbType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename ExtendedHubbard1OrbType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
	typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename PsimagLite::Vector<HilbertState>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	ExtendedSuperHubbard1Orb(const SolverParamsType& solverParams,
	                         InputValidatorType& io,
	                         const SuperGeometryType& superGeometry,
	                         PsimagLite::String extension)
	    : ModelBaseType(solverParams,
	                    superGeometry,
	                    io),
	      modelParameters_(io),
	      superGeometry_(superGeometry),
	      extension_(extension),
	      extendedHubbard_(solverParams, io, superGeometry, extension)
	{}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		extendedHubbard_.write(label, io);
	}

	virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                        const BlockType& block,
	                                        RealType time)  const
	{
		extendedHubbard_.addDiagonalsInNaturalBasis(hmatrix, block, time);
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		const SizeType site = 0;
		BlockType block(1, site);
		extendedHubbard_.fillLabeledOperators(qns);
		VectorOperatorType cm(3);
		for (SizeType dof = 0; dof < 2; ++dof)
			cm[dof] = this->naturalOperator("c", site, dof);

		cm[2] = this->naturalOperator("n", site, 0);

		// BEGIN define ph operator
		MatrixType tmp1 = multiplyTc(cm[1].getCRS(), cm[1].getCRS());
		MatrixType tmp2 = multiplyTc(cm[0].getCRS(), cm[0].getCRS());
		MatrixType tmp3 = tmp1*tmp2;
		SparseMatrixType tmp3crs;
		fullMatrixToCrsMatrix(tmp3crs, tmp3);
		OpsLabelType& ph = this->createOpsLabel("ph");
		typename OperatorType::Su2RelatedType su2related;
		OperatorType phOp(tmp3crs,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  typename OperatorType::PairType(0,0),
		                  1,
		                  su2related);
		ph.push(phOp);
		// END define ph operator

		OpsLabelType& p = this->createOpsLabel("pair");

		setPairi(p, cm);

		this->makeTrackable("splus");
		this->makeTrackable("sz");
		this->makeTrackable("pair");

		if (extension_ == "InteractionAssistedHopping") {
			// cdn_i == c^\dagger_{i \bar{sigma} n_{i sigma}
			OpsLabelType& cdn = this->createOpsLabel("cdn");
			for (SizeType i = 0; i < 2; ++i) {
				setCdn(cdn, cm, i);
			}
		}
	}

	void fillModelLinks()
	{
		extendedHubbard_.fillModelLinks();

		const bool spinOrbit = (superGeometry_.orbitals(0, 0) == 2);
		if (spinOrbit) {
			ModelTermType& so = ModelBaseType::createTerm("SpinOrbit");
			// spin dependence of the hopping parameter (spin orbit)
			OpForLinkType cup("c", 0, 0);
			OpForLinkType cdown("c", 1, 1);

			// spin1 = 1 spin2 = 0
			so.push(cdown, 'N', cup, 'C', typename ModelTermType::Su2Properties(1, -1, 1));

			// spin1 = 0 spin2 = 1
			so.push(cup, 'N', cdown, 'C', typename ModelTermType::Su2Properties(1, 1, 0));
		}

		const bool isSu2 = BasisType::useSu2Symmetry();

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");
		OpForLinkType splus("splus");

		spsm.push(splus,
		          'N',
		          splus,
		          'C',
		          [isSu2](SparseElementType& value) { value *= (isSu2) ? -0.5 : 0.5;},
		typename ModelTermType::Su2Properties(2, -1, 2));

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");
		OpForLinkType sz("sz");

		szsz.push(sz,
		          'N',
		          sz,
		          'N',
		          [isSu2](SparseElementType& value) { if (isSu2) value = -value; },
		typename ModelTermType::Su2Properties(2, 0.5, 1));

		ModelTermType& pp = ModelBaseType::createTerm("PairPair");
		OpForLinkType pair("pair");

		pp.push(pair, 'N', pair, 'C', typename ModelTermType::Su2Properties(2, 1, 2));

		if (extension_ == "InteractionAssistedHopping") {
			addInteractionAssistedHopping();
		}
	}

private:

	// c^\dagger_{i \bar{sigma} n_{i sigma} c_{j \bar{sigma}}
	void addInteractionAssistedHopping()
	{
		assert(extension_ == "InteractionAssistedHopping");
		ModelTermType& iah = ModelBaseType::createTerm("InteractionAssistedHopping");

		OpForLinkType cdn0("cdn", 0);
		OpForLinkType cdn1("cdn", 1);
		OpForLinkType cup("c", 0);
		OpForLinkType cdown("c", 1);

		iah.push(cdn0, 'N', cup, 'N');
		iah.push(cdn1, 'N', cdown, 'N');
	}

	// cdn_i == c^\dagger_{i \bar{sigma} n_{i sigma}
	void setCdn(OpsLabelType& p,
	            const VectorOperatorType& cm,
	            SizeType sigma) const
	{
		typename OperatorType::Su2RelatedType su2related;



		std::string nSigmaStr = (sigma == 0) ? "nup" : "ndown";
		SparseMatrixType nSigmaOp = this->naturalOperator(nSigmaStr, 0, 0).getCRS();
		OperatorType cDaggerBarSigma = cm[1 - sigma];
		cDaggerBarSigma.dagger();

		SparseMatrixType cdn;
		multiply(cdn, cDaggerBarSigma.getCRS(), nSigmaOp);

		OperatorType myOp(cdn,
		                  ProgramGlobals::FermionOrBosonEnum::FERMION,
		                  typename OperatorType::PairType(0,0),
		                  1,
		                  su2related);
		p.push(myOp);
	}


	void setPairi(OpsLabelType& p, const VectorOperatorType& cm) const
	{
		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType pair;
		multiply(pair,cm[0].getCRS(),cm[1].getCRS());

		OperatorType myOp(pair,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  typename OperatorType::PairType(0,0),
		                  1,
		                  su2related);
		p.push(myOp);
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceHubbardType::HilbertState const &ket,
	              int i,
	              int sigma) const
	{
		int value=0;
		value += HilbertSpaceHubbardType::calcNofElectrons(ket,0,i,0);
		value += HilbertSpaceHubbardType::calcNofElectrons(ket,0,i,1);
		int tmp1 = HilbertSpaceHubbardType::get(ket,0) &1;
		int tmp2 = HilbertSpaceHubbardType::get(ket,0) &2;
		if (i>0 && tmp1>0) value++;
		if (i>0 && tmp2>0) value++;

		if (sigma==1) { // spin down
			if ((HilbertSpaceHubbardType::get(ket,i) &1)) value++;

		}
		if (value%2==0) return 1.0;

		return -1;
	}

	ParametersModelHubbard<RealType, QnType>  modelParameters_;
	const SuperGeometryType& superGeometry_;
	PsimagLite::String extension_;
	ExtendedHubbard1OrbType extendedHubbard_;
};	//class ExtendedSuperHubbard1Orb

} // namespace Dmrg
/*@}*/
#endif // EXTENDED_SUPER_HUBBARD_1ORB_H

