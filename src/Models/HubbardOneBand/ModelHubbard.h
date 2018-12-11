/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file ModelHubbard.h
 *
 *  An implementation of the Hubbard Model to use with the DmrgSolver
 *
 */
#ifndef MODEL_HUBBARD_DMRG
#define MODEL_HUBBARD_DMRG
#include <cassert>
#include "Sort.h" // in PsimagLite
#include "ParametersModelHubbard.h"
#include "HilbertSpaceHubbard.h"
#include "LinkProductHubbardOneBand.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"
#include "MemResolv.h"

namespace Dmrg {
template<typename ModelBaseType> class ExtendedHubbard1Orb;

//! Model Hubbard for DMRG solver, inherits from ModelBase and implements its interface:
template<typename ModelBaseType>
class ModelHubbard : public ModelBaseType {

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::SparseElementType SparseElementType;
	typedef unsigned int long WordType;
	typedef  HilbertSpaceHubbard<WordType> HilbertSpaceHubbardType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;

private:

	static const int FERMION_SIGN = -1;
	static const int DEGREES_OF_FREEDOM=2;
	static const int NUMBER_OF_ORBITALS=1;

	enum {SPIN_UP = HilbertSpaceHubbardType::SPIN_UP,
		  SPIN_DOWN = HilbertSpaceHubbardType::SPIN_DOWN};

	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;

public:

	typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
	typedef LinkProductHubbardOneBand<ModelHelperType, GeometryType> LinkProductType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename PsimagLite::Vector<HilbertState>::Type HilbertBasisType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;

	ModelHubbard(const SolverParamsType& solverParams,
	             InputValidatorType& io,
	             GeometryType const &geometry,
	             PsimagLite::String terms)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    new LinkProductType(io, terms),
	                    io),
	      modelParameters_(io),
	      geometry_(geometry),
	      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM)
	{
		SizeType usize = modelParameters_.hubbardU.size();
		SizeType vsize = modelParameters_.potentialV.size();
		SizeType totalSites = geometry_.numberOfSites();

		if (usize != totalSites) {
			PsimagLite::String msg("ModelHubbard: hubbardU expecting ");
			msg += ttos(totalSites) + " entries, got " + ttos(usize) + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (vsize != 2*totalSites) {
			PsimagLite::String msg("ModelHubbard: potentialV expecting ");
			msg += ttos(2*totalSites) + " entries, got " + ttos(vsize) + "\n";
			throw PsimagLite::RuntimeError(msg);
		}
	}

	/** \cppFunction{!PTEX_THISFUNCTION} returns the size of the one-site Hilbert space. */
	SizeType hilbertSize(SizeType) const
	{
		return (SizeType)pow(2,2*NUMBER_OF_ORBITALS);
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

		OpsLabelType& c = this->createOpsLabel(OpsLabelType::TRACKABLE_YES, "c");
		VectorOperatorType creationMatrix(DEGREES_OF_FREEDOM);
		//! Set the operators c^\daggger_{i\sigma} in the natural basis
		SizeType ind = 0;
		for (int sigma=0;sigma<DEGREES_OF_FREEDOM;sigma++) {
			tmpMatrix = findOperatorMatrices(ind,sigma,natBasis);
			int asign= 1;
			if (sigma>0) asign= 1;
			typename OperatorType::Su2RelatedType su2related;
			if (sigma==0) {
				su2related.source.push_back(ind*DEGREES_OF_FREEDOM);
				su2related.source.push_back(ind*DEGREES_OF_FREEDOM+1);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(-1);
				su2related.offset = NUMBER_OF_ORBITALS;
			}

			OperatorType myOp(tmpMatrix,
			                  -1,
			                  typename OperatorType::PairType(1,1-sigma),
			                  asign,
			                  su2related);

			c.push(myOp);
			creationMatrix[sigma] = myOp;
		}

		SizeType iup = SPIN_UP;
		SizeType idown = SPIN_DOWN;

		{
			OpsLabelType& splus = this->createOpsLabel(OpsLabelType::TRACKABLE_NO, "splus");
			OpsLabelType& sminus = this->createOpsLabel(OpsLabelType::TRACKABLE_NO, "sminus");

			PsimagLite::Matrix<SparseElementType> tmp =
			        multiplyTc(creationMatrix[iup].data,creationMatrix[idown].data);
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			splus.push(OperatorType(tmp2,
			                        1.0,
			                        typename OperatorType::PairType(0,0),
			                        1.0,
			                        su2Related));
			SparseMatrixType tmp3;
			transposeConjugate(tmp3, tmp2);
			sminus.push(OperatorType(tmp3,
			                         1.0,
			                         typename OperatorType::PairType(0,0),
			                         1.0,
			                         su2Related));
		}

		{
			OpsLabelType& sz = this->createOpsLabel(OpsLabelType::TRACKABLE_NO, "sz");
			PsimagLite::Matrix<SparseElementType> tmp =
			        multiplyTc(creationMatrix[iup].data,creationMatrix[iup].data);
			PsimagLite::Matrix<SparseElementType> tmp2 =
			        multiplyTc(creationMatrix[idown].data,creationMatrix[idown].data);
			tmp = tmp-tmp2;
			SparseMatrixType tmp3(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			sz.push(OperatorType(tmp3,
			                     1.0,
			                     typename OperatorType::PairType(0,0),
			                     1.0,
			                     su2Related));
		}


		PsimagLite::Matrix<SparseElementType> dense1;
		{
			OpsLabelType& nupop = this->createOpsLabel(OpsLabelType::TRACKABLE_NO, "nup");
			OperatorType cup = creationMatrix[SPIN_UP];
			cup.dagger();
			SparseMatrixType tmp3(multiplyTc(cup.data,cup.data));
			typename OperatorType::Su2RelatedType su2Related;
			nupop.push(OperatorType(tmp3,
			                        1.0,
			                        typename OperatorType::PairType(0,0),
			                        1.0,
			                        su2Related));
			crsMatrixToFullMatrix(dense1, tmp3);
		}


		PsimagLite::Matrix<SparseElementType> dense2;
		{
			OpsLabelType& ndownop = this->createOpsLabel(OpsLabelType::TRACKABLE_NO, "ndown");

			OperatorType cdown = creationMatrix[SPIN_DOWN];
			cdown.dagger();
			SparseMatrixType tmp3(multiplyTc(cdown.data,cdown.data));
			typename OperatorType::Su2RelatedType su2Related;
			ndownop.push(OperatorType(tmp3,
			                          1.0,
			                          typename OperatorType::PairType(0,0),
			                          1.0,
			                          su2Related));
			crsMatrixToFullMatrix(dense2, tmp3);

		}

		{
			OpsLabelType& nop = this->createOpsLabel(OpsLabelType::TRACKABLE_NO, "n");
			dense1 += dense2;
			SparseMatrixType tmp(dense1);
			typename OperatorType::Su2RelatedType su2Related;
			nop.push(OperatorType(tmp,
			                      1.0,
			                      typename OperatorType::PairType(0,0),
			                      1.0,
			                      su2Related));
		}

		{
			OpsLabelType& d = this->createOpsLabel(OpsLabelType::TRACKABLE_NO, "d");
			PsimagLite::Matrix<SparseElementType> cup;
			crsMatrixToFullMatrix(cup,creationMatrix[SPIN_UP].data);
			PsimagLite::Matrix<SparseElementType> cdown;
			crsMatrixToFullMatrix(cdown,creationMatrix[SPIN_DOWN].data);
			cup = (cup*cdown);
			SparseMatrixType tmp3(cup);
			typename OperatorType::Su2RelatedType su2Related;
			d.push(OperatorType(tmp3,
			                    1.0,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related));
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		spinSquaredHelper_.write(label, io);
		spinSquared_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType& cm,
	                                const BlockType& block,
	                                RealType time)  const
	{
		SizeType n=block.size();
		SparseMatrixType tmpMatrix,niup,nidown,Szsquare,Szi;
		SizeType linSize = geometry_.numberOfSites();

		for (SizeType i=0;i<n;i++) {
			// onsite U hubbard
			//n_i up
			SizeType sigma =0; // up sector
			transposeConjugate(tmpMatrix,cm[sigma+i*DEGREES_OF_FREEDOM].data);
			multiply(niup,tmpMatrix,cm[sigma+i*DEGREES_OF_FREEDOM].data);
			//n_i down
			sigma =1; // down sector
			transposeConjugate(tmpMatrix,cm[sigma+i*DEGREES_OF_FREEDOM].data);
			multiply(nidown,tmpMatrix,cm[sigma+i*DEGREES_OF_FREEDOM].data);

			multiply(tmpMatrix,niup,nidown);
			RealType tmp = modelParameters_.hubbardU[block[i]];
			hmatrix += tmp*tmpMatrix;

			// V_iup term
			tmp = modelParameters_.potentialV[block[i]+0*linSize];
			hmatrix += tmp*niup;

			// V_idown term
			tmp = modelParameters_.potentialV[block[i]+1*linSize];
			hmatrix += tmp*nidown;

			// anisotropy
			if (modelParameters_.anisotropy.size() == linSize) {
				RealType mult1 =  1.0;
				RealType mult2 = -1.0;
				operatorPlus(Szi,niup,mult1,nidown,mult2);
				multiply(Szsquare,Szi,Szi);
				assert(i*2 < block.size());
				assert(block[i*2] < modelParameters_.anisotropy.size());
				RealType tmp = modelParameters_.anisotropy[block[i*2]]*0.25;
				hmatrix += tmp*Szsquare;
			}

			if (modelParameters_.potentialT.size()==0) continue;
			RealType cosarg = cos(time*modelParameters_.omega +
			                      modelParameters_.phase);
			// VT_iup term
			tmp = modelParameters_.potentialT[block[i]];
			tmp *= cosarg;
			hmatrix += tmp*niup;

			// VT_idown term
			tmp = modelParameters_.potentialT[block[i]];
			tmp *= cosarg;
			hmatrix += tmp*nidown;
		}
	}

	friend class ExtendedHubbard1Orb<ModelBaseType>;

private:

	void setBasis(HilbertBasisType& basis,
	              const VectorSizeType& block) const
	{
		int sitesTimesDof = DEGREES_OF_FREEDOM*block.size();
		HilbertState total = (1<<sitesTimesDof);

		basis.resize(total);
		for (HilbertState a = 0; a < total; ++a) basis[a] = a;
		if (basis.size() == 4) {
			basis[0] = 0;
			basis[1] = 2;
			basis[2] = 1;
			basis[3] = 3;
		}
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

		return FERMION_SIGN;
	}

	//! Find c^\dagger_isigma in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(int i,
	                                      int sigma,
	                                      const HilbertBasisType& natBasis) const
	{
		typename HilbertSpaceHubbardType::HilbertState bra,ket;
		int n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

		for (SizeType ii=0;ii<natBasis.size();ii++) {
			bra=ket=natBasis[ii];
			if (HilbertSpaceHubbardType::isNonZero(ket,i,sigma)) {

			} else {
				HilbertSpaceHubbardType::create(bra,i,sigma);
				int jj = PsimagLite::indexOrMinusOne(natBasis,bra);
				if (jj<0)
					throw PsimagLite::RuntimeError("findOperatorMatrices\n");
				cm(ii,jj) =sign(ket,i,sigma);
			}
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	void setSymmetryRelated(VectorQnType& qns,
	                        const HilbertBasisType& basis,
	                        int) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;

		bool isCanonical = (ModelBaseType::targetQuantum().isCanonical);

		qns.resize(basis.size(), QnType::zero());
		VectorSizeType other((isCanonical) ? 2 : 1, 0);
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair = calcJmValue<PairType>(basis[i]);
			// nup
			SizeType electronsUp = HilbertSpaceHubbardType::getNofDigits(basis[i],SPIN_UP);
			// ndown
			SizeType electronsDown = HilbertSpaceHubbardType::getNofDigits(basis[i],SPIN_DOWN);

			other[0] = electronsUp + electronsDown;

			if (isCanonical) other[1] = electronsUp;

			bool sign = other[0] & 1;
			qns[i] = QnType(sign, other, jmpair, other[0]);
		}
	}

	// note: we use 2j instead of j
	// note: we use m+j instead of m
	// This assures us that both j and m are SizeType
	// does not work for 6 or 9
	template<typename PairType>
	PairType calcJmValue(const HilbertState& ket) const
	{
		if (!MyBasis::useSu2Symmetry()) return PairType(0,0);
		SizeType site0=0;
		SizeType site1=0;

		spinSquared_.doOnePairOfSitesA(ket,site0,site1);
		spinSquared_.doOnePairOfSitesB(ket,site0,site1);
		spinSquared_.doDiagonal(ket,site0,site1);

		RealType sz = spinSquared_.spinZ(ket,site0);
		PairType jm= spinSquaredHelper_.getJmPair(sz);

		return jm;
	}

	ParametersModelHubbard<RealType, QnType>  modelParameters_;
	const GeometryType &geometry_;
	SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
	SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
};	//class ModelHubbard

} // namespace Dmrg
/*@}*/
#endif

