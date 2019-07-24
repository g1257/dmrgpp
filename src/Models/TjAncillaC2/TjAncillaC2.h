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

/*! \file TjAncillaC2.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef DMRG_TJ_ANCILLAC2_H
#define DMRG_TJ_ANCILLAC2_H
#include "../Models/FeAsModel/ModelFeBasedSc.h"
#include "../Models/TjAncillaC/ParametersTjAncillaC.h"

namespace Dmrg {
//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
template<typename ModelBaseType>
class TjAncillaC2 : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorRealType VectorRealType;
	typedef ModelFeBasedSc<ModelBaseType> ModelFeAsType;
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
	typedef typename ModelFeAsType::HilbertState HilbertStateType;
	typedef typename ModelFeAsType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelFeAsType::HilbertSpaceFeAsType HilbertSpaceType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename OperatorType::PairType PairType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<HilbertStateType>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;

	static const int NUMBER_OF_ORBITALS = 2;
	static const int FERMION_SIGN = -1;

	enum {SPIN_UP, SPIN_DOWN};

	TjAncillaC2(const SolverParamsType& solverParams,
	            InputValidatorType& io,
	            GeometryType const &geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      geometry_(geometry),
	      hot_(geometry_.orbitals(0,0) > 1)
	{}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		io.write(label + "/hot_", hot_);
	}

	virtual PsimagLite::String oracle() const
	{
		const RealType ne = ModelBaseType::targetQuantum().qn(0).other[0];
		const RealType n = ModelBaseType::geometry().numberOfSites();
		RealType energy = -ne*(n - ne);
		return ModelBaseType::oracle(energy, "-Ne*(L-Ne)");
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);

		SizeType orbitals = (hot_) ? 2 : 1;
		VectorHilbertStateType natBasis;
		SparseMatrixType tmpMatrix;
		setBasis(natBasis, block);
		setSymmetryRelated(qns, natBasis, block.size());

		OpsLabelType& c = this->createOpsLabel("c");
		OpsLabelType& splus = this->createOpsLabel("splus");
		OpsLabelType& sminus = this->createOpsLabel("sminus");
		OpsLabelType& sz = this->createOpsLabel("sz");
		OpsLabelType& nop = this->createOpsLabel("n");
		OpsLabelType& d = this->createOpsLabel("d");

		this->makeTrackable("c");
		this->makeTrackable("splus");
		this->makeTrackable("sz");
		this->makeTrackable("n");
		this->makeTrackable("d");

		// Set the operators c^\daggger_{i\sigma} in the natural basis
		for (SizeType i=0;i<block.size();i++) {
			for (int sigma=0;sigma<2;sigma++) {
				for (SizeType orb = 0; orb < orbitals; ++orb) {
					findOperatorMatrices(tmpMatrix,i,orb+sigma*NUMBER_OF_ORBITALS,natBasis);
					int asign= 1;
					typename OperatorType::Su2RelatedType su2related;

					OperatorType myOp(tmpMatrix,
					                  ProgramGlobals::FermionOrBosonEnum::FERMION,
					                  PairType(0, 0),
					                  asign,
					                  su2related);

					c.push(myOp);
				}
			}

			for (SizeType orbital = 0; orbital < orbitals; ++orbital) {
				// Set the operators S^+_i in the natural basis
				tmpMatrix=findSplusMatrices(i,orbital,natBasis);

				typename OperatorType::Su2RelatedType su2related;

				OperatorType myOp(tmpMatrix,
				                  ProgramGlobals::FermionOrBosonEnum::BOSON,
				                  PairType(0, 0),
				                  1.0,
				                  su2related);
				splus.push(myOp);
				myOp.dagger();
				sminus.push(myOp);
			}

			for (SizeType orbital = 0; orbital < orbitals; ++orbital) {
				// Set the operators S^z_i in the natural basis
				tmpMatrix = findSzMatrices(i,orbital,natBasis);
				typename OperatorType::Su2RelatedType su2related2;
				OperatorType myOp2(tmpMatrix,
				                   ProgramGlobals::FermionOrBosonEnum::BOSON,
				                   PairType(0, 0),
				                   1.0,
				                   su2related2);
				sz.push(myOp2);
			}

			for (SizeType orbital = 0; orbital < orbitals; ++orbital) {
				// Set ni matrix:
				SparseMatrixType tmpMatrix = findNiMatrices(i,orbital,natBasis);
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

			// Set delta_{i, sigma} matrix:
			for (int sigma=0;sigma<2;sigma++) {
				tmpMatrix = findDeltaIMatrices(i,natBasis,sigma);
				RealType angularFactor= 1;
				typename OperatorType::Su2RelatedType su2related4;
				OperatorType myOp4(tmpMatrix,
				                   ProgramGlobals::FermionOrBosonEnum::BOSON,
				                   PairType(0, 0),
				                   angularFactor,
				                   su2related4);
				d.push(myOp4);
			}
		}

		{
			OperatorType tmp = this->naturalOperator("c", site, 0);
			tmp.dagger();
			SparseMatrixType c = tmp.getCRS();
			SparseMatrixType tmp3(multiplyTc(c,c));
			typename OperatorType::Su2RelatedType su2Related;
			this->createOpsLabel("nup").push(OperatorType(tmp3,
			                                              ProgramGlobals::FermionOrBosonEnum::BOSON,
			                                              typename OperatorType::PairType(0,0),
			                                              1.0,
			                                              su2Related));
		}

		{
			OperatorType tmp = this->naturalOperator("c", site, 1);
			tmp.dagger();
			SparseMatrixType c = tmp.getCRS();
			SparseMatrixType tmp3(multiplyTc(c,c));
			typename OperatorType::Su2RelatedType su2Related;
			this->createOpsLabel("ndown").push(OperatorType(tmp3,
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

		const SizeType orbitals = (hot_) ? 2 : 1;
		const bool isSu2 = BasisType::useSu2Symmetry();
		ModelTermType& hop = ModelBaseType::createTerm("hopping");
		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");
		ModelTermType& szsz = ModelBaseType::createTerm("szsz");
		ModelTermType& ninj = ModelBaseType::createTerm("ninj");
		ModelTermType& dd = ModelBaseType::createTerm("didj");

		for (SizeType spin = 0; spin < 2; ++spin) {
			for (SizeType orb1 = 0; orb1 < orbitals; ++orb1) {
				OpForLinkType c1("c", orb1 + spin*orbitals, orb1);
				for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
					OpForLinkType c2("c", orb2 + spin*orbitals, orb2);

					hop.push(c1, 'N', c2, 'C', 1, (spin == 1) ? -1 : 1, spin);
				}
			}

			OpForLinkType d("d", spin);
			dd.push(d, 'N', d, 'C');
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

				spsm.push(splus1, 'N', splus2, 'C', 2, -1, 2, valueModiferTerm0);

				if (!isSu2)
					szsz.push(sz1, 'N', sz2, 'N', 2, 0.5);
				else
					spsm.push(splus1, 'N', splus2, 'C', 2, -1, 2, valueModifierTermOther);

				ninj.push(n1, 'N', n2, 'N');
			}
		}
	}

private:

	//! find all states in the natural basis for a block of n sites
	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		assert(block.size()==1);
		int sitesTimesDof = 2*NUMBER_OF_ORBITALS;
		HilbertStateType total = (1<<sitesTimesDof);

		basis.clear();
		for (HilbertStateType a = 0; a < total; ++a) {
			if (!isAllowed(a)) continue;
			basis.push_back(a);
		}

		assert(basis.size() == pow(3,NUMBER_OF_ORBITALS));
	}

	bool isAllowed(HilbertStateType a) const
	{
		HilbertStateType mask1 = 1;
		mask1 <<= NUMBER_OF_ORBITALS;
		mask1--;

		HilbertStateType mask2 = mask1;
		mask2 <<= NUMBER_OF_ORBITALS;

		HilbertStateType a1 = (a & mask1);
		HilbertStateType a2 = (a & mask2);
		a2 >>= NUMBER_OF_ORBITALS;
		return ((a1 & a2) == 0);
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to
	//! basis state ket
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	RealType sign(const HilbertStateType& ket, int i,SizeType sigma) const
	{
		int value=0;
		SizeType dofs=2*NUMBER_OF_ORBITALS;
		for (SizeType alpha=0;alpha<dofs;alpha++)
			value += HilbertSpaceType::calcNofElectrons(ket,0,i,alpha);
		// add electron on site 0 if needed
		if (i>0) value += HilbertSpaceType::electrons(ket);

		//order for sign is: a up, a down, b up, b down, etc
		unsigned int x = HilbertSpaceType::get(ket,i);
		int spin = sigma/NUMBER_OF_ORBITALS;
		SizeType orb = sigma % NUMBER_OF_ORBITALS;

		for (SizeType j=0;j<orb;j++) {
			for (SizeType k=0;k<2;k++) {
				SizeType ind = j + k * NUMBER_OF_ORBITALS;
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
		int n = natBasis.size();
		PsimagLite::Matrix<ComplexOrRealType> cm(n,n);
		findCreationDense(cm,i,sigma,natBasis);
		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,cm);
		transposeConjugate(creationMatrix,temp);

	}
	//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
	//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
	void findCreationDense(PsimagLite::Matrix<ComplexOrRealType>& cm,
	                       int i,
	                       int sigma,
	                       const HilbertBasisType& natBasis) const
	{
		HilbertStateType bra,ket;

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];

			if (HilbertSpaceType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceType::create(bra,i,sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				if (jj<0) continue;
				if (ii==SizeType(jj)) {
					std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra;
					std::cerr<<" sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("Creation operator diagonal\n");
				}

				cm(ii,jj) = sign(ket,i,sigma);
			}
		}
	}

	//! Find S^+_i in the natural basis natBasis
	SparseMatrixType findSplusMatrices(int i,
	                                   SizeType orb,
	                                   const VectorHilbertStateType& natBasis) const
	{
		HilbertStateType bra,ket;
		int n = natBasis.size();
		MatrixType cm(n,n);
		SizeType sigma1 = orb + SPIN_DOWN*NUMBER_OF_ORBITALS;
		SizeType sigma2 = orb + SPIN_UP*NUMBER_OF_ORBITALS;
		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			if (HilbertSpaceType::isNonZero(ket,i,sigma1) &&
			        !HilbertSpaceType::isNonZero(ket,i,sigma2)) {
				// it is a down electron, then flip it:
				HilbertSpaceType::destroy(bra,i,sigma1);
				HilbertSpaceType::create(bra,i,sigma2);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				assert(jj>=0);
				cm(ii,jj)=1.0;
			}
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(int i,
	                                SizeType orb,
	                                const VectorHilbertStateType& natBasis) const
	{
		HilbertStateType ket;
		int n = natBasis.size();
		MatrixType cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			ket=natBasis[ii];
			RealType value = 0.0;
			if (HilbertSpaceType::isNonZero(ket,i,orb + SPIN_UP*NUMBER_OF_ORBITALS))
				value += 1.0;
			if (HilbertSpaceType::isNonZero(ket,i,orb + SPIN_DOWN*NUMBER_OF_ORBITALS))
				value -= 1.0;

			cm(ii,ii)=0.5*value;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findNiMatrices(int i,
	                                SizeType orb,
	                                const VectorHilbertStateType& natBasis) const
	{
		SizeType n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			HilbertStateType ket=natBasis[ii];
			cm(ii,ii) = 0.0;
			for (SizeType sigma=0;sigma<2;sigma++)
				if (HilbertSpaceType::isNonZero(ket,i,orb + sigma*NUMBER_OF_ORBITALS))
					cm(ii,ii) += 1.0;
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	//! Find Delta_i in the natural basis natBasis
	SparseMatrixType findDeltaIMatrices(int i,
	                                    const VectorHilbertStateType& natBasis,
	                                    SizeType sigma) const
	{
		assert(NUMBER_OF_ORBITALS == 2);
		int n = natBasis.size();

		PsimagLite::Matrix<ComplexOrRealType> cr(n,n);
		findCreationDense(cr,i,sigma*NUMBER_OF_ORBITALS,natBasis);

		PsimagLite::Matrix<ComplexOrRealType> ca(n,n);
		SizeType sigmaBar = 1 - sigma;
		findCreationDense(ca,i,1 + sigmaBar*NUMBER_OF_ORBITALS,natBasis);

		PsimagLite::Matrix<ComplexOrRealType> part1 = cr*ca;

		SparseMatrixType temp;
		fullMatrixToCrsMatrix(temp,part1);
		return temp;
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType) const
	{
		SizeType n=block.size();

		SizeType linSize = geometry_.numberOfSites();
		for (SizeType i=0;i<n;i++) {
			SizeType orb = 0;
			// potentialV
			SparseMatrixType nup(this->naturalOperator("nup",i,orb).getCRS());
			SparseMatrixType ndown(this->naturalOperator("ndown",i,orb).getCRS());
			SparseMatrixType m = nup;
			SizeType index = block[i]+orb*linSize;
			assert(index<modelParameters_.potentialV.size());
			m *= modelParameters_.potentialV[block[i] + orb*linSize];
			m += modelParameters_.potentialV[index]*ndown;
			hmatrix += m;

		}
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

		SizeType basisSize = basis.size();
		VectorSizeType other(4, 0);
		qns.resize(basisSize, QnType::zero());
		for (SizeType i = 0; i < basisSize; ++i) {
			PairType jmpair = calcJmvalue<PairType>(basis[i]);

			// nup
			other[1] = HilbertSpaceType::electronsWithGivenSpin(basis[i],SPIN_UP);
			// ndown
			SizeType ndown = HilbertSpaceType::electronsWithGivenSpin(basis[i],SPIN_DOWN);

			other[2] =  HilbertSpaceType::calcNofElectrons(basis[i],
			                                               0);
			other[3] = HilbertSpaceType::calcNofElectrons(basis[i],
			                                              NUMBER_OF_ORBITALS);

			other[0] = other[1] + ndown;
			bool sign = other[0] & 1;
			qns[i] = QnType(sign, other, jmpair, other[0]);
		}
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// Reinterprets 6 and 9
	template<typename PairType>
	PairType calcJmvalue(const HilbertStateType&) const
	{
		return PairType(0,0);
	}

	ParametersTjAncillaC<RealType, QnType>  modelParameters_;
	const GeometryType &geometry_;
	bool hot_;
};	//class TjAncillaC2

} // namespace Dmrg
/*@}*/
#endif // DMRG_TJ_ANCILLAC2_H

