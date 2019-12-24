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

/*! \file HeisenbergAncillaC.h
 *
 *  An implementation of the Quantum Heisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_HEISENBERG_ANCILLAC_H
#define DMRG_HEISENBERG_ANCILLAC_H

#include <algorithm>
#include "ModelBase.h"
#include "ParametersHeisenbergAncillaC.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramGlobals.h"

namespace Dmrg {

/* Order of the basis
  For S=1/2 and NUMBER_OF_ORBITALS=2 it is
  down A, down B; down A, up B; up A, down B; up A, up B;
  */
template<typename ModelBaseType>
class HeisenbergAncillaC : public ModelBaseType {

	static const int NUMBER_OF_ORBITALS=2;

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::RealType RealType;
	typedef	typename ModelBaseType::VectorType VectorType;
	typedef	typename std::pair<SizeType,SizeType> PairSizeType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef unsigned int long WordType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef typename PsimagLite::Vector<unsigned int long>::Type HilbertBasisType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	HeisenbergAncillaC(const SolverParamsType& solverParams,
	                   InputValidatorType& io,
	                   const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      hot_(geometry.orbitals(0,0) > 1)
	{
		SizeType n = geometry.numberOfSites();
		SizeType m = modelParameters_.magneticField.size();
		if (m > 0 && m != n) {
			PsimagLite::String msg("HeisenbergAncillaC: If provided, ");
			msg += " MagneticField must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (BasisType::useSu2Symmetry()) {
			PsimagLite::String msg("HeisenbergAncillaC: SU(2) symmetry");
			msg += " is not implemented yet.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (modelParameters_.twiceTheSpin != 1) {
			PsimagLite::String msg("HeisenbergAncillaC: spin > 1/2");
			msg += " HIGHLY EXPERIMENTAL!\n";
			std::cout<<msg;
		}

		if (hot_) {
			PsimagLite::String msg("HeisenbergAncillaC: Hot ancilla mode is on");
			msg += " (EXPERIMENTAL feature)\n";
			std::cout<<msg;
			std::cerr<<msg;
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		io.write(label + "/hot_", hot_);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType)  const
	{
		if (modelParameters_.magneticField.size() !=
		        ModelBaseType::superGeometry().numberOfSites())
			return;

		SizeType n = block.size();

		for (SizeType i = 0; i < n; ++i) {
			SizeType site = block[i];
			const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);
			// magnetic field
			RealType tmp = modelParameters_.magneticField[site];
			hmatrix += tmp * sz.getCRS();
		}
	}

	virtual PsimagLite::String oracle() const
	{
		const RealType nup = ModelBaseType::targetQuantum().qn(0).other[0];
		const RealType n = ModelBaseType::superGeometry().numberOfSites();
		RealType energy = -nup*(n - nup);
		return ModelBaseType::oracle(energy, "-Nup*Ndown");
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

		OpsLabelType& splus = this->createOpsLabel("splus");
		OpsLabelType& sminus = this->createOpsLabel("sminus");
		OpsLabelType& sz = this->createOpsLabel("sz");
		OpsLabelType& d = this->createOpsLabel("d");

		this->makeTrackable("splus");
		this->makeTrackable("sz");
		this->makeTrackable("d");

		for (SizeType i=0;i<block.size();i++) {
			// Set the operators S^+_i for orbital a in the natural basis
			tmpMatrix=findSplusMatrices(i,0,natBasis);

			typename OperatorType::Su2RelatedType su2related;
			su2related.source.push_back(i*2);
			su2related.source.push_back(i*2+NUMBER_OF_ORBITALS);
			su2related.source.push_back(i*2);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(1);
			su2related.offset = NUMBER_OF_ORBITALS;

			OperatorType myOp(tmpMatrix,
			                  ProgramGlobals::FermionOrBosonEnum::BOSON,
			                  PairType(2, 2),
			                  -1,
			                  su2related);
			splus.push(myOp);

			myOp.dagger();
			sminus.push(myOp);

			if (hot_) {
				// Set the operators S^+_i for orbital b in the natural basis
				tmpMatrix=findSplusMatrices(i,1,natBasis);

				typename OperatorType::Su2RelatedType su2related;
				su2related.source.push_back(i*2);
				su2related.source.push_back(i*2+NUMBER_OF_ORBITALS);
				su2related.source.push_back(i*2);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(1);
				su2related.offset = NUMBER_OF_ORBITALS;

				OperatorType myOp(tmpMatrix,
				                  ProgramGlobals::FermionOrBosonEnum::BOSON,
				                  PairType(2, 2),
				                  -1,
				                  su2related);
				splus.push(myOp);
			}

			// Set the operators S^z_i orbital a in the natural basis
			tmpMatrix = findSzMatrices(i,0,natBasis);
			typename OperatorType::Su2RelatedType su2related2;
			OperatorType myOp2(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(2, 1),
			                   1.0/sqrt(2.0),
			                   su2related2);
			sz.push(myOp2);
			if (hot_) {
				// Set the operators S^z_i orbital b in the natural basis
				tmpMatrix = findSzMatrices(i,1,natBasis);
				typename OperatorType::Su2RelatedType su2related2;
				OperatorType myOp2(tmpMatrix,
				                   ProgramGlobals::FermionOrBosonEnum::BOSON,
				                   PairType(2, 1),
				                   1.0/sqrt(2.0),
				                   su2related2);
				sz.push(myOp2);
			}

			// Set the operators \Delta_i in the natural basis
			tmpMatrix = findDeltaMatrices(i,natBasis);
			typename OperatorType::Su2RelatedType su2related3;
			OperatorType myOp3(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(0, 0),
			                   1.0,
			                   su2related3);
			d.push(myOp3);
		}
	}

	void fillModelLinks()
	{
		const bool isSu2 = BasisType::useSu2Symmetry();

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");
		OpForLinkType splus0("splus");

		auto modifierTerm0 = [isSu2](SparseElementType& value) { value *= (isSu2) ? -0.5 : 0.5;};
		spsm.push(splus0,
		          'N',
		          splus0,
		          'C',
		          modifierTerm0,
		          typename ModelTermType::Su2Properties(2, -1, 2));

		if (hot_) {
			OpForLinkType splus1("splus", 1, 1);
			spsm.push(splus1,
			          'N',
			          splus1,
			          'C',
			          modifierTerm0,
			          typename ModelTermType::Su2Properties(2, -1, 2));
		}

		auto modifierTerm1 = [isSu2](SparseElementType& value) {if (isSu2) value = -value;};

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");
		OpForLinkType sz0("sz");

		szsz.push(sz0,
		          'N',
		          sz0,
		          'C',
		          modifierTerm1,
		          typename ModelTermType::Su2Properties(2, 0.5, 1));

		if (hot_) {
			OpForLinkType splus1("splus", 1, 1);
			OpForLinkType sz1("sz", 1, 1);
			if (isSu2)
				szsz.push(splus1,
				          'N',
				          splus1,
				          'C',
				          modifierTerm1,
				          typename ModelTermType::Su2Properties(2, 0.5, 1));
			else
				szsz.push(sz1, 'N', sz1, 'C', typename ModelTermType::Su2Properties(2, 0.5, 1));
		}

		ModelTermType& ancilla = ModelBaseType::createTerm("ancilla");
		OpForLinkType d("d");
		ancilla.push(d, 'N', d, 'C', typename ModelTermType::Su2Properties(2, 1, 0));
	}

private:

	//! find all states in the natural basis for a block of n sites
	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		assert(block.size()==1);
		SizeType total = pow(modelParameters_.twiceTheSpin + 1, NUMBER_OF_ORBITALS);

		basis.resize(total);
		for (SizeType i = 0; i < total; ++i) basis[i] = i;
	}

	//! Find S^+_i a in the natural basis natBasis
	SparseMatrixType findSplusMatrices(int,
	                                   SizeType orbital,
	                                   const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5*modelParameters_.twiceTheSpin;
		SizeType total1 = modelParameters_.twiceTheSpin + 1;
		for (SizeType ii=0;ii<total;ii++) {
			PairSizeType ket = getOneOrbital(natBasis[ii]);

			SizeType ket1 = (orbital == 0) ? ket.first : ket.second;
			SizeType ket2 = (orbital == 0) ? ket.second : ket.first;
			SizeType bra1 = ket1 + 1;
			if (bra1 >= total1) continue;
			PairSizeType bra = (orbital == 0) ? PairSizeType(bra1,ket2) : PairSizeType(ket2,bra1);
			SizeType jj = getFullIndex(bra);
			RealType m = ket1 - j;

			RealType x = j*(j+1)-m*(m+1);
			assert(x>=0);
			cm(ii,jj) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(int,
	                                SizeType orbital,
	                                const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5*modelParameters_.twiceTheSpin;

		for (SizeType ii=0;ii<total;ii++) {
			PairSizeType ket = getOneOrbital(natBasis[ii]);
			SizeType ket1 = (orbital == 0) ? ket.first : ket.second;
			RealType m = ket1 - j;
			cm(ii,ii) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findDeltaMatrices(int,const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5*modelParameters_.twiceTheSpin;
		SizeType total1 = modelParameters_.twiceTheSpin + 1;
		for (SizeType ii=0;ii<total;ii++) {
			PairSizeType ket = getOneOrbital(natBasis[ii]);

			SizeType bra1 = ket.first;
			if (bra1 >= total1) continue;

			SizeType bra2 = ket.second;
			if (bra2 >= total1) continue;
			if (bra2 >= bra1) continue;

			PairSizeType bra(bra2,bra1);
			SizeType jj = getFullIndex(bra);
			RealType m = bra.first - j;
			RealType x1 = j*(j+1)-m*(m+1);
			assert(x1>=0);
			m = ket.second - j;
			RealType x2 = j*(j+1)-m*(m+1);
			assert(x2>=0);
			cm(ii,jj) = sqrt(x1)*sqrt(x2);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		if (n!=1)
			PsimagLite::RuntimeError("HeisenbergAncillaC::setSymmetryRelated n=1 only\n");

		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		VectorSizeType other(2, 0);
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair;
			jmpair.first = modelParameters_.twiceTheSpin;
			jmpair.second = basis[i];
			PairSizeType ket = getOneOrbital(basis[i]);
			other[1] = ket.first + ket.second;
			other[0] = ket.first;
			SizeType flavor = 1;
			bool sign = other[0] & 1;
			qns[i] = QnType(sign, other, jmpair, flavor);
		}
	}

	PairSizeType getOneOrbital(SizeType state) const
	{
		SizeType total1 = modelParameters_.twiceTheSpin + 1;
		assert(NUMBER_OF_ORBITALS == 2);
		SizeType first = static_cast<SizeType>(state/total1);
		SizeType second = state % total1;
		return PairSizeType(first, second);
	}

	SizeType getFullIndex(PairSizeType p) const
	{
		return p.first*(modelParameters_.twiceTheSpin + 1) + p.second;
	}

	ParametersHeisenbergAncillaC<RealType, QnType> modelParameters_;
	bool hot_;
}; // class HeisenbergAncillaC

} // namespace Dmrg
/*@}*/
#endif //DMRG_HEISENBERG_ANCILLAC_H

