/*
Copyright (c) 2009-2013-2021, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 6.]
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
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"

namespace Dmrg {
template<typename ModelBaseType> class ExtendedHubbard1Orb;

//! Model Hubbard for DMRG solver, inherits from ModelBase and implements its interface:
template<typename ModelBaseType>
class ModelHubbard : public ModelBaseType {

	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;

	static const int FERMION_SIGN = -1;
	static const int DEGREES_OF_FREEDOM=2;
	static const int NUMBER_OF_ORBITALS=1;

public:

	typedef ModelHubbard<ModelBaseType> ThisType;
	typedef unsigned int long WordType;
	typedef  HilbertSpaceHubbard<WordType> HilbertSpaceHubbardType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::SparseElementType SparseElementType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef typename ModelBaseType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename PsimagLite::Vector<HilbertState>::Type HilbertBasisType;
	typedef typename ModelBaseType::ModelTermType ModelTermType;
	typedef typename ModelBaseType::OpForLinkType OpForLinkType;
	typedef typename ModelBaseType::OpsLabelType OpsLabelType;

	enum {SPIN_UP = HilbertSpaceHubbardType::SPIN_UP,
		  SPIN_DOWN = HilbertSpaceHubbardType::SPIN_DOWN};

	ModelHubbard(const SolverParamsType& solverParams,
	             InputValidatorType& io,
	             const SuperGeometryType& geometry,
	             PsimagLite::String extension)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    io),
	      modelParameters_(io),
	      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM),
	      extension_(extension)
	{
		ModelBaseType::onSiteHLegacyFix(modelParameters_.onSiteHaddLegacy);

		SizeType usize = modelParameters_.hubbardU.size();
		SizeType vsize = modelParameters_.potentialV.size();
		SizeType totalSites = geometry.numberOfSites();

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

		if (modelParameters_.magneticX.size()>0 && ModelBaseType::targetQuantum().sizeOfOther() == 2) {
			PsimagLite::String msg("ModelHubbard: Sz not conserved: You should remove ");
			msg += "TargeteElectronsDown or TargetSzPlusConst from the input file.\n";
			msg += "TargeteElectronsUp is then interpreted as total number of electrons\n";
			throw PsimagLite::RuntimeError(msg);
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

		OpsLabelType& c = this->createOpsLabel("c");
		this->makeTrackable("c");

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
			                  ProgramGlobals::FermionOrBosonEnum::FERMION,
			                  typename OperatorType::PairType(1,1-sigma),
			                  asign,
			                  su2related);

			c.push(myOp, (sigma == 0) ? "up" : "down");
			creationMatrix[sigma] = myOp;
		}

		SizeType iup = SPIN_UP;
		SizeType idown = SPIN_DOWN;

		{
			OpsLabelType& splus = this->createOpsLabel("splus");
			OpsLabelType& sminus = this->createOpsLabel("sminus");

			PsimagLite::Matrix<SparseElementType> tmp =
			        multiplyTc(creationMatrix[iup].getCRS(),creationMatrix[idown].getCRS());
			SparseMatrixType tmp2(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			splus.push(OperatorType(tmp2,
			                        ProgramGlobals::FermionOrBosonEnum::BOSON,
			                        typename OperatorType::PairType(0,0),
			                        1.0,
			                        su2Related));
			SparseMatrixType tmp3;
			transposeConjugate(tmp3, tmp2);
			sminus.push(OperatorType(tmp3,
			                         ProgramGlobals::FermionOrBosonEnum::BOSON,
			                         typename OperatorType::PairType(0,0),
			                         1.0,
			                         su2Related));
		}

		{
			OpsLabelType& sz = this->createOpsLabel("sz");
			PsimagLite::Matrix<SparseElementType> tmp =
			        multiplyTc(creationMatrix[iup].getCRS(),creationMatrix[iup].getCRS());
			PsimagLite::Matrix<SparseElementType> tmp2 =
			        multiplyTc(creationMatrix[idown].getCRS(),creationMatrix[idown].getCRS());
			tmp = 0.5*(tmp - tmp2);
			SparseMatrixType tmp3(tmp);
			typename OperatorType::Su2RelatedType su2Related;
			sz.push(OperatorType(tmp3,
			                     ProgramGlobals::FermionOrBosonEnum::BOSON,
			                     typename OperatorType::PairType(0,0),
			                     1.0,
			                     su2Related));
		}


		PsimagLite::Matrix<SparseElementType> dense1;
		{
			OpsLabelType& nupop = this->createOpsLabel("nup");
			OperatorType cup = creationMatrix[SPIN_UP];
			cup.dagger();
			SparseMatrixType tmp3(multiplyTc(cup.getCRS(),cup.getCRS()));
			typename OperatorType::Su2RelatedType su2Related;
			nupop.push(OperatorType(tmp3,
			                        ProgramGlobals::FermionOrBosonEnum::BOSON,
			                        typename OperatorType::PairType(0,0),
			                        1.0,
			                        su2Related));
			crsMatrixToFullMatrix(dense1, tmp3);
		}


		PsimagLite::Matrix<SparseElementType> dense2;
		{
			OpsLabelType& ndownop = this->createOpsLabel("ndown");

			OperatorType cdown = creationMatrix[SPIN_DOWN];
			cdown.dagger();
			SparseMatrixType tmp3(multiplyTc(cdown.getCRS(),cdown.getCRS()));
			typename OperatorType::Su2RelatedType su2Related;
			ndownop.push(OperatorType(tmp3,
			                          ProgramGlobals::FermionOrBosonEnum::BOSON,
			                          typename OperatorType::PairType(0,0),
			                          1.0,
			                          su2Related));
			crsMatrixToFullMatrix(dense2, tmp3);

		}

		{
			OpsLabelType& nop = this->createOpsLabel("n");
			dense1 += dense2;
			SparseMatrixType tmp(dense1);
			typename OperatorType::Su2RelatedType su2Related;
			nop.push(OperatorType(tmp,
			                      ProgramGlobals::FermionOrBosonEnum::BOSON,
			                      typename OperatorType::PairType(0,0),
			                      1.0,
			                      su2Related));
		}

		{
			OpsLabelType& d = this->createOpsLabel("d");
			PsimagLite::Matrix<SparseElementType> cup;
			crsMatrixToFullMatrix(cup,creationMatrix[SPIN_UP].getCRS());
			PsimagLite::Matrix<SparseElementType> cdown;
			crsMatrixToFullMatrix(cdown,creationMatrix[SPIN_DOWN].getCRS());
			cup = (cup*cdown);
			SparseMatrixType tmp3(cup);
			typename OperatorType::Su2RelatedType su2Related;
			d.push(OperatorType(tmp3,
			                    ProgramGlobals::FermionOrBosonEnum::BOSON,
			                    typename OperatorType::PairType(0,0),
			                    1.0,
			                    su2Related));
		}
	}

	/* PSIDOC Hubbard::fillModelLinks
	The Hubbard model has one one term, the hopping,
	and within this term it has two connections: up-up and down-down.
	First, define the term, and name it whatever you like; see line (A) below.
	In this case the variable is called \verb!hop! and the name is ``hopping''.
	Create then a special type of opaque object, an OpForLinkType, to represent
	the operators you need to connect. In this case, we need to connect $c^\dagger_\uparrow$
	with $c_\uparrow$. Because $c_\uparrow$ is not tracked in this model, we
	must use $c^\dagger_\uparrow$ only, and transpose conjugate it as needed.
	In line (B) below, we create a C++ variable \verb!cup! and the label name is ``c'' with
	degree of freedom 0 indicating spin up. Note that (``c'', 0) must be
	a trackable operator; see \verb!fillLabeledOperators! above.
	We now push to the variable \verb!hop! the connection as seen in (C) below.
	Here we must supply
	\emph{two} operators: in this case they are cup with 'N' and cup with 'C' to
	indicate that we want $c^\dagger c$. The last three numbers are SU(2) related; just
	use \verb!1, 1, 0! if unsure.
	Note that the operators are the same, but the second one is transpose conjugated.
	Likewise, we create cdown in (D) and add the connection to \verb!hop! in (E)
	PSIDOCCOPY $FirstFunctionBelow
	*/
	void fillModelLinks()
	{
		ModelTermType& hop = ModelBaseType::createTerm("hopping");//(A)

		OpForLinkType cup("c", 0); // (B)
		hop.push(cup, 'N', cup, 'C', typename ModelTermType::Su2Properties(1, 1, 0)); // (C)

		OpForLinkType cdown("c", 1); // (D)
		hop.push(cdown, 'N', cdown, 'C', typename ModelTermType::Su2Properties(1, -1, 1)); // (E)

		if (extension_ != "RashbaSOC") return;

		assert(extension_ == "RashbaSOC");

		ModelTermType& rashbaSOC = ModelBaseType::createTerm("RashbaSOC");

		// up-down
		rashbaSOC.push(cup, 'N', cdown, 'C');

		// down-up
		auto valueModifer = [](SparseElementType& value)
		{ value = -PsimagLite::conj(value);};

		rashbaSOC.push(cdown,
		               'N',
		               cup,
		               'C',
		               valueModifer);
	}

	/* PSIDOC Hubbard::write
	 We are passed a label in \verb!label1! (the first argument),
	and something like a file handle in \verb!io! (the second argument).
	We must interpret the label as a directory within the file.
	We then first create a top directory if it doesn't exist in lines (A) and (B) below.
	Oh, and then create another directory with this model's name in (C).
	Hmmm, maybe this two-step process ought to be simplified\ldots.
	Anyway, we go ahead and start writing things into directory \verb!label! being
	\verb!label1/model!,
	where label1 is the label passed as first argument to this function,
	and model is this model's name. What things \emph{need to} be written here?
	Probably nothing needs to be written (yeah, big FIXME TODO on the ``Probably'', sorry).
	Typically though, you'd like to write all data members for this class.
	For most models, it's the model parameters like in (D) below,
	and that's that. Here we write also some SU(2) auxiliary members in the last two lines.
	PSIDOCCOPY $FirstFunctionBelow
	 */
	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))    // (A)
			io.createGroup(label1);          // (B)

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);             // (C)
		modelParameters_.write(label, io); // (D)
		spinSquaredHelper_.write(label, io);
		spinSquared_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const BlockType& block,
	                                RealType time)  const
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		SizeType n=block.size();
		SparseMatrixType tmpMatrix,niup,nidown,Szsquare,Szi,Splusi,Sminusi;
		SizeType linSize = ModelBaseType::superGeometry().numberOfSites();

		for (SizeType i = 0; i < n; ++i) {
			const SizeType site = block[i];

			// onsite U hubbard
			//n_i up
			const OperatorType& cup = ModelBaseType::naturalOperator("c", site, 0);
			transposeConjugate(tmpMatrix, cup.getCRS());
			multiply(niup,tmpMatrix, cup.getCRS());
			//n_i down
			const OperatorType& cdown = ModelBaseType::naturalOperator("c", site, 1);
			transposeConjugate(tmpMatrix, cdown.getCRS());
			multiply(Sminusi, tmpMatrix, cup.getCRS());
			transposeConjugate(Splusi, Sminusi);
			multiply(nidown,tmpMatrix, cdown.getCRS());

			multiply(tmpMatrix, niup, nidown);
			RealType tmp = modelParameters_.hubbardU[block[i]];
			hmatrix += tmp*tmpMatrix;

			// V_iup term
			tmp = modelParameters_.potentialV[site + 0*linSize];
			hmatrix += tmp*niup;

			// V_idown term
			tmp = modelParameters_.potentialV[site + 1*linSize];
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

			// magneticX
			if (modelParameters_.magneticX.size() == linSize) {
				RealType tmp = modelParameters_.magneticX[block[i]]*0.5;
				hmatrix += tmp*Sminusi;
				hmatrix += tmp*Splusi;
			}
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

		const bool isCanonical = (ModelBaseType::targetQuantum().sizeOfOther() == 2);
		if (isCanonical && extension_ == "RashbaSOC")
			err(PsimagLite::String(__FILE__) +
			    ": RashbaSOC sub-model must be canonical. Please " +
			    "delete the TargetSzPlusConst= from the input file\n");

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
	SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
	SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
	PsimagLite::String extension_;
};	//class ModelHubbard

} // namespace Dmrg
/*@}*/
#endif

