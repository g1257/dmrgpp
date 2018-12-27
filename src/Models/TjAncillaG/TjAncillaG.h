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

/*! \file TjAncillaG.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef DMRG_TJ_ANCILLAG_H
#define DMRG_TJ_ANCILLAG_H
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "../Models/TjMultiOrb/TjMultiOrb.h"
#include "../Models/TjMultiOrb/ParametersModelTjMultiOrb.h"

namespace Dmrg {
//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
template<typename ModelBaseType>
class TjAncillaG : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef ModelHubbard<ModelBaseType> ModelHubbardType;
	typedef TjMultiOrb<ModelBaseType> TjMultiOrbType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
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
	typedef	typename ModelBaseType::MyBasis BasisType;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelHubbardType::HilbertState HilbertStateType;
	typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelHubbardType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename OperatorType::PairType PairType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<HilbertStateType>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const int DEGREES_OF_FREEDOM=2;
	static const int NUMBER_OF_ORBITALS = 1;
	static const int FERMION_SIGN = -1;

	enum {SPIN_UP, SPIN_DOWN};

	TjAncillaG(const SolverParamsType& solverParams,
	           InputValidatorType& io,
	           GeometryType const &geometry)
	    : ModelBaseType(solverParams, geometry, io),
	      modelParameters_(io),
	      geometry_(geometry),
	      TjMultiOrb_(solverParams,io,geometry),
	      offset_(DEGREES_OF_FREEDOM+3), // c^\dagger_up, c^\dagger_down, S+, Sz, n
	      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM)
	{
		if (BasisType::useSu2Symmetry())
			err("TjAncillaG: SU(2) not supported\n");
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType& block,
	                                RealType) const
	{
		SizeType n=block.size();

		SizeType linSize = geometry_.numberOfSites();
		for (SizeType i=0;i<n;i++) {
			// potentialV
			SparseMatrixType nup = this->naturalOperator("nup",i,0).data;
			SparseMatrixType ndown = this->naturalOperator("ndown",i,0).data;
			SparseMatrixType m = nup;
			assert(block[i]+linSize<modelParameters_.potentialV.size());
			m *= modelParameters_.potentialV[block[i]];
			m += modelParameters_.potentialV[block[i]+linSize]*ndown;
			hmatrix += m;
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		TjMultiOrb_.write(label, io);
		io.write(label + "/offset_", offset_);
		spinSquaredHelper_.write(label, io);
		spinSquared_.write(label, io);
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		VectorHilbertStateType natBasis;
		SparseMatrixType tmpMatrix;
		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis, block.size());

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
		for (SizeType i=0;i<block.size();i++) {
			for (int sigma=0;sigma<DEGREES_OF_FREEDOM;sigma++) {
				tmpMatrix = TjMultiOrb_.findCreationMatrices(i,sigma,natBasis);
				int asign= 1;
				if (sigma>0) asign= 1;
				typename OperatorType::Su2RelatedType su2related;
				if (sigma==0) {
					su2related.source.push_back(i*offset_);
					su2related.source.push_back(i*offset_+1);
					su2related.transpose.push_back(-1);
					su2related.transpose.push_back(-1);
					su2related.offset = NUMBER_OF_ORBITALS;
				}

				OperatorType myOp(tmpMatrix,-1,PairType(1,1-sigma),asign,su2related);

				c.push(myOp);
			}

			// Set the operators S^+_i in the natural basis
			tmpMatrix=findSplusMatrices(i,natBasis);

			typename OperatorType::Su2RelatedType su2related;
			su2related.source.push_back(i*DEGREES_OF_FREEDOM);
			su2related.source.push_back(i*DEGREES_OF_FREEDOM+NUMBER_OF_ORBITALS);
			su2related.source.push_back(i*DEGREES_OF_FREEDOM);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(-1);
			su2related.transpose.push_back(1);
			su2related.offset = NUMBER_OF_ORBITALS;

			OperatorType myOp(tmpMatrix,1,PairType(2,2),-1,su2related);
			splus.push(myOp);
			myOp.dagger();
			sminus.push(myOp);

			// Set the operators S^z_i in the natural basis
			tmpMatrix = findSzMatrices(i,natBasis);
			typename OperatorType::Su2RelatedType su2related2;
			OperatorType myOp2(tmpMatrix,1,PairType(2,1),1.0/sqrt(2.0),su2related2);
			sz.push(myOp2);

			// Set ni matrices:
			SparseMatrixType tmpMatrix = findNiMatrices(0,natBasis);
			RealType angularFactor= 1;
			typename OperatorType::Su2RelatedType su2related3;
			su2related3.offset = 1; //check FIXME
			OperatorType myOp3(tmpMatrix,1,PairType(0,0),angularFactor,su2related3);

			nop.push(myOp3);
		}

		{
			OpsLabelType& nupop = this->createOpsLabel("nup");
			OperatorType tmp = this->naturalOperator("c", site, 0);
			tmp.dagger();
			SparseMatrixType c = tmp.data;
			SparseMatrixType tmp3(multiplyTc(c,c));
			typename OperatorType::Su2RelatedType su2Related;
			nupop.push(OperatorType(tmp3,
			                                              1.0,
			                                              typename OperatorType::PairType(0,0),
			                                              1.0,
			                                              su2Related));
		}

		{
			OperatorType tmp = this->naturalOperator("c", site, 1);
			tmp.dagger();
			SparseMatrixType c = tmp.data;
			SparseMatrixType tmp3(multiplyTc(c,c));
			typename OperatorType::Su2RelatedType su2Related;
			this->createOpsLabel("ndown").push(OperatorType(tmp3,
			                                                1.0,
			                                                typename OperatorType::PairType(0,0),
			                                                1.0,
			                                                su2Related));
		}
	}

	void fillModelLinks()
	{
		const bool isSu2 = BasisType::useSu2Symmetry();

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");
		OpForLinkType splus0("splus");

		auto modifierTerm0 = [isSu2](ComplexOrRealType& value) { value *= (isSu2) ? -0.5 : 0.5;};
		spsm.push(splus0, 'N', splus0, 'C', 2, -1, 2, modifierTerm0);

		auto modifierTerm1 = [isSu2](ComplexOrRealType& value) {if (isSu2) value = -value;};

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");
		OpForLinkType sz0("sz");

		szsz.push(sz0, 'N', sz0, 'C', 2, 0.5, 1, modifierTerm1);

		ModelTermType& ancilla = ModelBaseType::createTerm("ancilla");
		OpForLinkType d("d");
		ancilla.push(d, 'N', d, 'C', 2, 1, 0);
	}

private:

	//! find all states in the natural basis for a block of n sites
	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		assert(block.size() == 1);
		int sitesTimesDof = DEGREES_OF_FREEDOM;
		HilbertStateType total = (1<<sitesTimesDof);
		--total;

		basis.resize(total);
		for (HilbertStateType a = 0; a < total; ++a) basis[a] = a;
	}

	SparseMatrixType findSplusMatrices(int i,
	                                   const VectorHilbertStateType& natBasis) const
	{
		assert(i == 0);
		HilbertStateType bra,ket;
		int n = natBasis.size();
		MatrixType cm(n,n);
		SizeType orbitals = modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			ket=natBasis[ii];
			for (SizeType l = orbitals; l < 2*orbitals; ++l) {
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
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findSzMatrices(int i,
	                                const VectorHilbertStateType& natBasis) const
	{
		assert(i == 0);
		int n = natBasis.size();
		MatrixType cm(n,n);
		SizeType orbitals = modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			int counter = 0;
			for (SizeType l = 0; l < orbitals; ++l) {
				HilbertStateType masklp = (1<<l);
				if (ket & masklp) counter++;
			}

			for (SizeType l = orbitals; l < 2*orbitals; ++l) {
				HilbertStateType masklp = (1<<l);
				if (ket & masklp) counter--;
			}

			cm(ii,ii) = 0.5*counter;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findNiMatrices(int,
	                                const VectorHilbertStateType& natBasis) const
	{
		SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);
		SizeType dofs = 2*modelParameters_.orbitals;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			cm(ii,ii) = 0.0;
			for (SizeType sigma=0;sigma<dofs;sigma++) {
				HilbertStateType mask = (1<<sigma);
				if (ket & mask) cm(ii,ii) += 1.0;
			}
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
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

		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = calcJmvalue<PairType>(basis[i]);
			SizeType szPlusConst = (basis[i] == 2) ? 0 : basis[i] + 1;

			qns[i] = QnType(false, VectorSizeType(1, szPlusConst), jmpair, 0);
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

	//serializr start class TjAncillaG
	//serializr vptr
	//serializr normal modelParameters_
	ParametersModelTjMultiOrb<RealType, QnType>  modelParameters_;
	//serializr ref geometry_ end
	const GeometryType &geometry_;
	TjMultiOrbType TjMultiOrb_;
	//serializr normal offset_
	SizeType offset_;
	//serializr normal spinSquaredHelper_
	SpinSquaredHelper<RealType,HilbertStateType> spinSquaredHelper_;
	//serializr normal spinSquared_
	SpinSquared<SpinSquaredHelper<RealType,HilbertStateType> > spinSquared_;

};	//class TjAncillaG

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_ANCILLAG_H

