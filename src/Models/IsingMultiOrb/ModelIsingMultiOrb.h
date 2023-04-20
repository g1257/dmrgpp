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

/*! \file ModelIsingMultiOrb.h
 *
 *  An implementation of the Quantum Ising MultiOrb Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_MODEL_ISINGMULTIORB_HEADER_H
#define DMRG_MODEL_ISINGMULTIORB_HEADER_H

#include "HilbertSpaceIsingMultiOrb.h"
#include <algorithm>
#include "ParametersModelIsingMultiOrb.h"
#include "CrsMatrix.h"
#include "../../Engine/VerySparseMatrix.h"
#include "../../Engine/ProgramGlobals.h"
#include "../../Engine/Utils.h"

namespace Dmrg {

template<typename ModelBaseType>
class ModelIsingMultiOrb : public ModelBaseType {

	static const int DEGREES_OF_FREEDOM=2; // spin up and down

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
        typedef typename HilbertBasisType::value_type HilbertState;	
        typedef  HilbertSpaceIsingMultiOrb<HilbertState> HilbertSpaceIsingMultiOrbType;	
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	using ModelParametersType = ParametersModelIsingMultiOrb<RealType, QnType>;

        static const int SPIN_UP=HilbertSpaceIsingMultiOrbType::SPIN_UP;
        static const int SPIN_DOWN=HilbertSpaceIsingMultiOrbType::SPIN_DOWN;
	
	ModelIsingMultiOrb(const SolverParamsType& solverParams,
	                InputValidatorType& io,
	                const SuperGeometryType& geometry,
	                PsimagLite::String additional)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      superGeometry_(geometry),
	      additional_(additional)
	{
                ProgramGlobals::init(modelParameters_.orbitals*superGeometry_.numberOfSites() + 1);

		HilbertSpaceIsingMultiOrbType::setOrbitals(modelParameters_.orbitals);
		
		SizeType n = superGeometry_.numberOfSites();
		SizeType orbs = modelParameters_.orbitals;
		SizeType md_cols = modelParameters_.onsitelinksSzSz.cols();
                SizeType md_rows = modelParameters_.onsitelinksSzSz.rows();
		SizeType md = md_rows*md_cols;
		SizeType nn = ModelParametersType::combinations(orbs,2);

		ModelParametersType::checkMagneticField(modelParameters_.magneticFieldX.cols(), 'X', n, orbs);

		ModelParametersType::checkMagneticField(modelParameters_.magneticFieldZ.cols(), 'Z', n, orbs);

		if (md > 0 && md != nn) {
			PsimagLite::String msg("ModelIsingMultiOrb: If provided, ");
			msg += " onsitelinksSzSz must be a vector of " + ttos(nn) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

	/* PSIDOC IsingMultiOrb::addDiagonalsInNaturalBasis
	 We describe only the addition of a Zeeman term to the Ising model here; note
	 that this function is more complicated.
	 Please look carefully at the following C++ lines:
	 PSIDOCCOPY $FirstFunctionBelow::MagneticField
	 */
	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType time) const
	{
//		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		SizeType n = block.size();
                SizeType orbs = modelParameters_.orbitals;
                SizeType orbs1 = ModelParametersType::combinations(orbs,2);		

		for (SizeType i = 0; i < n; ++i) {

			// PSIDOCMARK_BEGIN MagneticField
			SizeType site = block[i];
			for (SizeType orb = 0; orb < orbs; ++orb) {
				addMagneticField(hmatrix, 'X', site, orb);
				addMagneticField(hmatrix, 'Z', site, orb);
			}
                        for (SizeType orb1 = 0; orb1 < orbs1; ++orb1) {
                        // onsitelinksSzSz			
				//addonsitelinksSzSz(hmatrix,site,orb1);
			// PSIDOCMARK_END
			}
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns)
	{
		SizeType site = 0;
		BlockType block(1, site);
		const SizeType orbs = modelParameters_.orbitals;
		SizeType total = utils::powUint(2, orbs*block.size());
		HilbertBasisType natBasis(total);
		for (SizeType i = 0; i < total; ++i) natBasis[i] = i;

		setSymmetryRelated(qns, natBasis, block.size());

		for (SizeType i=0;i<block.size();i++) {

			OpsLabelType& splus = this->createOpsLabel("splus");
			OpsLabelType& sminus = this->createOpsLabel("sminus");

                        for (SizeType orb=0;orb<orbs;orb++) {
                        // Set the operators S^+_i in the natural basis			
                        	SparseMatrixType tmpMatrix = findSplusMatrices(i,orb,natBasis);
                        	typename OperatorType::Su2RelatedType su2related;
                        	su2related.source.push_back(i*DEGREES_OF_FREEDOM);
                        	su2related.source.push_back(i*DEGREES_OF_FREEDOM+orbs);
                        	su2related.source.push_back(i*DEGREES_OF_FREEDOM);
                        	su2related.transpose.push_back(-1);
                        	su2related.transpose.push_back(-1);
                        	su2related.transpose.push_back(1);
                        	su2related.offset = orbs;

                        	splus.push(OperatorType(tmpMatrix,
                                                ProgramGlobals::FermionOrBosonEnum::BOSON,
                                                typename OperatorType::PairType(0,0),
                                                1.0,
                                                su2related));
                        	SparseMatrixType tmp;
                        	transposeConjugate(tmp, tmpMatrix);
                        	sminus.push(OperatorType(tmp,
                                                 ProgramGlobals::FermionOrBosonEnum::BOSON,
                                                 typename OperatorType::PairType(0,0),
                                                 1.0,
                                                 su2related));
			}
			
	                OpsLabelType& sz = this->createOpsLabel("sz");
                        this->makeTrackable("sz");			
                        for (SizeType orb=0;orb<orbs;orb++) {			
			// Set the operators S^z_i in the natural basis
				SparseMatrixType tmpMatrix = findSzMatrices(i,orb,natBasis);
				typename OperatorType::Su2RelatedType su2related2;
				sz.push(OperatorType(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(2, 1),
			                   1.0/sqrt(2.0),
			                   su2related2));
			}

                        OpsLabelType& sx = this->createOpsLabel("sx");
                        OpsLabelType& sybar = this->createOpsLabel("sybar");
			for (SizeType orb=0;orb<orbs;orb++) {
			// Set the operators S^x_i in the natural basis
				SparseMatrixType tmpMatrix = findSxOrSyBarMatrices(i,orb,natBasis, "sx");
				typename OperatorType::Su2RelatedType su2related3;
				sx.push(OperatorType(tmpMatrix,
			                   ProgramGlobals::FermionOrBosonEnum::BOSON,
			                   PairType(2, 1),
			                   1.0/sqrt(2.0),
			                   su2related3));

			// Set the operators S^ybar_i in the natural basis
				tmpMatrix = findSxOrSyBarMatrices(i, orb, natBasis, "sybar");
				sybar.push(OperatorType(tmpMatrix,
			                     ProgramGlobals::FermionOrBosonEnum::BOSON,
			                     PairType(0, 0),
			                     1,
			                     su2related3));
			}
		}
	}

	void fillModelLinks()
	{
		if (BasisType::useSu2Symmetry())
			err("SU(2) no longer supported\n");

		connectionSzSz();
	}

private:

	void addMagneticField(SparseMatrixType& hmatrix,
	                      char c,
	                      SizeType site,
			      SizeType orb) const
	{
		assert(c == 'X' || c == 'Z');

		const SizeType linSize = superGeometry_.numberOfSites();
                const SizeType orbs = modelParameters_.orbitals;
		const PsimagLite::Matrix<RealType>& v = (c == 'X') ? modelParameters_.magneticFieldX
		                                     : modelParameters_.magneticFieldZ;

		if (v.cols() != linSize || v.rows() != orbs) return;

		assert(site < v.cols());

		RealType tmp = v(orb,site);

		const OperatorType& sminus = ModelBaseType::naturalOperator("sminus", site, orb);

		const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, orb);

		if (c == 'Z') {
			hmatrix += tmp*sz.getCRS();
			return;
		}

		assert(c == 'X');
		const OperatorType& splus = ModelBaseType::naturalOperator("splus", site, orb);
		constexpr RealType zeroPointFive = 0.5;
		hmatrix += zeroPointFive*tmp*splus.getCRS();
		hmatrix += zeroPointFive*tmp*sminus.getCRS();
	}

        void connectionSzSz()
        {
                const SizeType orbitals = modelParameters_.orbitals;
                ModelTermType& szsz = ModelBaseType::createTerm("szsz");
                for (SizeType orb1 = 0; orb1 < orbitals; ++orb1) {
                        OpForLinkType c1("sz", orb1, orb1); 
                        for (SizeType orb2 = 0; orb2 < orbitals; ++orb2) {
                                OpForLinkType c2("sz", orb2, orb2);
                                szsz.push(c1,'N',c2,'N',typename ModelTermType::Su2Properties(2,0.5));
                        }
                }
        }
	
        //! Find S^+_i in the natural basis natBasis
        SparseMatrixType findSplusMatrices(SizeType site,
					   SizeType orb,
                                           const HilbertBasisType& natBasis) const
        {
                SizeType total = natBasis.size();
                MatrixType cm(total,total);
		SizeType mask = 1;
		mask<<=orb;
                SizeType ket = natBasis[0];
                SizeType bra = ket | mask;
                cm(ket,bra) = 1.0;
		std::cout<<cm;
                SparseMatrixType operatorMatrix(cm);
                return operatorMatrix;
        }
	
	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(SizeType site,
					SizeType orb,
	                                const HilbertBasisType& natBasis) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		RealType j = 0.5;
		SizeType mask = 1;
		mask <<= orb; // mask = 2^orb

		for (SizeType ii=0;ii<total;ii++) {
			SizeType ket = natBasis[ii];
			SizeType ketsite = ket & mask;
			if (ketsite>1) ketsite/=mask;
			RealType m = ketsite - j;
			cm(ket,ket) = m;
		}
		std::cout<<cm;
		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findSxOrSyBarMatrices(SizeType site,
					       SizeType orb,
	                                       const HilbertBasisType& natBasis,
	                                       PsimagLite::String what) const
	{
		if (what != "sx" && what != "sybar")
			err("findSxOrSyBarMatrices: don't know how to calculate " + what + "\n");

		const RealType sign = (what == "sx") ? 1.0 : -1.0;
		SparseMatrixType Splus_temp=findSplusMatrices(site,orb,natBasis);
		SparseMatrixType Sminus_temp,Sx;
		transposeConjugate(Sminus_temp,Splus_temp);
		const RealType tmp = 0.5;

		Sx = tmp*Splus_temp;
		Sx += sign*tmp*Sminus_temp;

		return Sx;
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;

		bool isCanonical = (ModelBaseType::targetQuantum().sizeOfOther() == 1);

		if (isCanonical && modelParameters_.magneticFieldX.cols()>0)
			err(PsimagLite::String(__FILE__) +
			    ": MagneticFieldX CANNOT be canonical. Please " +
			    "delete the TargetSzPlusConst= from the input file\n");

		VectorSizeType other;
		if (isCanonical) other.resize(1, 0);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(1, basis[i]);
			if (isCanonical)
				other[0] = getSzPlusConst(basis[i], n);
			SizeType flavor = 1;
			qns[i] = QnType(false, other, jmpair, flavor);
		}
	}

	SizeType getSzPlusConst(SizeType ket, SizeType n) const
	{
		if (n == 1) return ket;

		SizeType bitsForOneSite = 1;

		SizeType sum = 0;
		for (SizeType site = 0; site < n; ++site) {
			SizeType mask = 1;
			mask <<= (site*bitsForOneSite);
			SizeType ketsite = ket & mask;
			ketsite >>= (site*bitsForOneSite);
			sum += ketsite;
		}

		return sum;
	}

	ModelParametersType  modelParameters_;
	const SuperGeometryType& superGeometry_;
	PsimagLite::String additional_;

}; // class ModelIsingMultiOrb

} // namespace Dmrg
/*@}*/
#endif //DMRG_MODEL_ISINGMULTIORB_HEADER_H

