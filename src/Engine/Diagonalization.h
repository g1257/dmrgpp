/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

/*! \file Diagonalization.h
 *
 *  FIXME needs doc
 */
#ifndef DIAGONALIZATION_HEADER_H
#define DIAGONALIZATION_HEADER_H
#include "ProgressIndicator.h"
#include "VectorWithOffset.h" // includes the std::norm functions
#include "VectorWithOffsets.h" // includes the std::norm functions
#include "ProgramGlobals.h"
#include "LanczosSolver.h"
#include "DavidsonSolver.h"
#include "ParametersForSolver.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename ParametersType, typename TargettingType>
class Diagonalization {

public:

	typedef typename TargettingType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename TargettingType::ModelType ModelType;
	typedef typename TargettingType::BasisType BasisType;
	typedef typename TargettingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename TargettingType::BlockType BlockType;
	typedef typename TargettingType::TargetVectorType TargetVectorType;
	typedef typename TargettingType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename  OperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef typename TargettingType::MatrixVectorType MatrixVectorType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	Diagonalization(const ParametersType& parameters,
	                const ModelType& model,
	                const bool& verbose,
	                ReflectionSymmetryType& reflectionOperator,
	                InputValidatorType& io,
	                const SizeType& quantumSector,
	                WaveFunctionTransfType& waveFunctionTransformation)
	    : parameters_(parameters),
	      model_(model),
	      verbose_(verbose),
	      reflectionOperator_(reflectionOperator),
	      io_(io),
	      progress_("Diag."),
	      quantumSector_(quantumSector),
	      wft_(waveFunctionTransformation),
	      oldEnergy_(0)
	{}

	//!PTEX_LABEL{Diagonalization}
	RealType operator()(TargettingType& target,
	                    SizeType direction,
	                    const BlockType& blockLeft,
	                    const BlockType& blockRight)
	{
		if (direction!=WaveFunctionTransfType::INFINITE) throw PsimagLite::RuntimeError(
		            "Diagonalization::operator(): expecting INFINITE direction\n");
		SizeType loopIndex = 0;
		  VectorSizeType sectors;
		targetedSymmetrySectors(sectors,target.leftRightSuper());
		reflectionOperator_.update(sectors);
		RealType gsEnergy = internalMain_(target,direction,loopIndex,false,blockLeft);
		//  targetting:
		target.evolve(gsEnergy,direction,blockLeft,blockRight,loopIndex);
		wft_.triggerOff(target.leftRightSuper()); //,m);
		return gsEnergy;
	}

	RealType operator()(TargettingType& target,
	                    SizeType direction,
	                    const BlockType& block,
	                    SizeType loopIndex,
	                    bool)
	{
		assert(direction!=WaveFunctionTransfType::INFINITE);

		RealType gsEnergy = internalMain_(target,direction,loopIndex,false,block);
		//  targetting:
		target.evolve(gsEnergy,direction,block,block,loopIndex);
		wft_.triggerOff(target.leftRightSuper()); //,m);
		return gsEnergy;
	}

private:

	void targetedSymmetrySectors(VectorSizeType& mVector,
	                             const LeftRightSuperType& lrs) const
	{
		SizeType total = lrs.super().partition()-1;
		for (SizeType i=0;i<total;i++) {
			//SizeType bs = lrs.super().partition(i+1)-lrs.super().partition(i);
			if (lrs.super().pseudoEffectiveNumber(lrs.super().partition(i))!=quantumSector_ )
				continue;
			mVector.push_back(i);
		}
	}

	RealType internalMain_(TargettingType& target,
	                       SizeType direction,
	                       SizeType loopIndex,
	                       bool,
	                       const VectorSizeType& block)

	{
		const LeftRightSuperType& lrs= target.leftRightSuper();
		wft_.triggerOn(lrs);

		RealType gsEnergy = 0;

		bool onlyWft = false;
		if (direction != WaveFunctionTransfType::INFINITE)
			onlyWft = ((parameters_.finiteLoop[loopIndex].saveOption & 2)>0) ? true : false;

		if (parameters_.options.find("MettsTargetting")!=PsimagLite::String::npos ||
		        parameters_.options.find("TargetingAncilla")!=PsimagLite::String::npos)
			return gsEnergy;

		PsimagLite::OstringStream msg0;
		msg0<<"Setting up Hamiltonian basis of size="<<lrs.super().size();
		progress_.printline(msg0,std::cout);

		typename PsimagLite::Vector<TargetVectorType>::Type vecSaved;
		typename PsimagLite::Vector<RealType>::Type energySaved;

		SizeType total = lrs.super().partition()-1;

		energySaved.resize(total);
		vecSaved.resize(total);
		  VectorSizeType weights(total);

		SizeType counter=0;
		for (SizeType i=0;i<total;i++) {
			SizeType bs = lrs.super().partition(i+1)-lrs.super().partition(i);
			if (verbose_) {
				SizeType j = lrs.super().qn(lrs.super().partition(i));
				  VectorSizeType qns = BasisType::decodeQuantumNumber(j);
				//std::cerr<<"partition "<<i<<" of size="<<bs<<" has qns=";
				for (SizeType k=0;k<qns.size();k++) std::cerr<<qns[k]<<" ";
				//std::cerr<<"\n";
			}

			weights[i]=bs;

			// Do only one sector unless doing su(2) with j>0, then do all m's
			if (lrs.super().pseudoEffectiveNumber(lrs.super().partition(i))!=quantumSector_ )
				weights[i]=0;

			counter+=bs;
			vecSaved[i].resize(weights[i]);
		}

		typedef typename TargettingType::VectorWithOffsetType
		        VectorWithOffsetType;
		VectorWithOffsetType initialVector(weights,lrs.super());

		target.initialGuess(initialVector,block);

		for (SizeType i=0;i<total;i++) {
			if (weights[i]==0) continue;
			PsimagLite::OstringStream msg;
			msg<<"About to diag. sector with quantum numbs. ";
			SizeType j = lrs.super().qn(lrs.super().partition(i));
			  VectorSizeType qns = BasisType::decodeQuantumNumber(j);
			for (SizeType k=0;k<qns.size();k++) msg<<qns[k]<<" ";
			msg<<" pseudo="<<lrs.super().pseudoEffectiveNumber(
			         lrs.super().partition(i));
			msg<<" quantumSector="<<quantumSector_;

			if (verbose_ && PsimagLite::Concurrency::root()) {
				msg<<" diagonaliseOneBlock, i="<<i;
				msg<<" and weight="<<weights[i];
			}
			progress_.printline(msg,std::cout);
			TargetVectorType initialVectorBySector(weights[i]);
			initialVector.extract(initialVectorBySector,i);
			RealType norma = PsimagLite::norm(initialVectorBySector);
			if (fabs(norma)<1e-12) {
				PsimagLite::String str("Norm of initial vector is zero\n");
				if (onlyWft) throw PsimagLite::RuntimeError("FATAL " + str);
			} else {
				initialVectorBySector /= norma;
			}

			if (onlyWft) {
				vecSaved[i]=initialVectorBySector;
				gsEnergy = oldEnergy_;
			} else {
				diagonaliseOneBlock(i,vecSaved[i],gsEnergy,lrs,initialVectorBySector);
			}

			energySaved[i]=gsEnergy;
		}

		// calc gs energy
		if (verbose_ && PsimagLite::Concurrency::root())
			std::cerr<<"About to calc gs energy\n";
		gsEnergy=1e6;
		for (SizeType i=0;i<total;i++) {
			if (weights[i]==0) continue;
			if (energySaved[i]<gsEnergy) gsEnergy=energySaved[i];
		}

		if (verbose_ && PsimagLite::Concurrency::root())
			std::cerr<<"About to calc gs vector\n";

		counter=0;
		for (SizeType i=0;i<lrs.super().partition()-1;i++) {
			if (weights[i]==0) continue;

			SizeType j = lrs.super().qn(lrs.super().partition(i));
			  VectorSizeType qns = BasisType::decodeQuantumNumber(j);
			PsimagLite::OstringStream msg;
			msg<<"Found targetted symmetry sector in partition "<<i;
			msg<<" of size="<<vecSaved[i].size();
			progress_.printline(msg,std::cout);

			PsimagLite::OstringStream msg2;
			msg2<<"Norm of vector is "<<PsimagLite::norm(vecSaved[i]);
			msg2<<" and quantum numbers are ";
			for (SizeType k=0;k<qns.size();k++) msg2<<qns[k]<<" ";
			progress_.printline(msg2,std::cout);
			counter++;
		}

		target.setGs(vecSaved,lrs.super());

		if (PsimagLite::Concurrency::root()) {
			oldEnergy_=gsEnergy;
		}

		return gsEnergy;
	}

	/** Diagonalise the i-th block of the matrix, return its eigenvectors
			in tmpVec and its eigenvalues in energyTmp
		!PTEX_LABEL{diagonaliseOneBlock} */
	void diagonaliseOneBlock(int i,
	                         TargetVectorType &tmpVec,
	                         RealType &energyTmp,
	                         const LeftRightSuperType& lrs,
	                         const TargetVectorType& initialVector)
	{
		typename PsimagLite::Vector<RealType>::Type tmpVec1,tmpVec2;
		//srand48(7123443);

		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(i,lrs,threadId);

		if (parameters_.options.find("debugmatrix")!=PsimagLite::String::npos) {
			SparseMatrixType fullm;

			model_.fullHamiltonian(fullm,modelHelper);

			PsimagLite::Matrix<typename SparseMatrixType::value_type> fullm2;
			crsMatrixToFullMatrix(fullm2,fullm);
			if (PsimagLite::isZero(fullm2)) std::cerr<<"Matrix is zero\n";
			if (fullm.row()>40) {
				printNonZero(fullm2,std::cerr);
			} else {
				printFullMatrix(fullm,"matrix",1);
			}

			if (!isHermitian(fullm,true))
				throw PsimagLite::RuntimeError("Not hermitian matrix block\n");

			typename PsimagLite::Vector<RealType>::Type eigs(fullm2.n_row());
			PsimagLite::diag(fullm2,eigs,'V');
			std::cerr<<"eigs[0]="<<eigs[0]<<"\n";
			if (parameters_.options.find("test")!=PsimagLite::String::npos)
				throw std::logic_error("Exiting due to option test in the input file\n");
		}
		PsimagLite::OstringStream msg;
		msg<<"I will now diagonalize a matrix of size="<<modelHelper.size();
		progress_.printline(msg,std::cout);
		diagonaliseOneBlock(i,tmpVec,energyTmp,modelHelper,initialVector);
	}

	void diagonaliseOneBlock(int i,
	                         TargetVectorType &tmpVec,
	                         RealType &energyTmp,
	                         typename ModelType::ModelHelperType& modelHelper,
	                         const TargetVectorType& initialVector)
	{
		int n = modelHelper.size();
		if (verbose_)
			std::cerr<<"Lanczos: About to do block number="<<i<<" of size="<<n<<"\n";

		ReflectionSymmetryType *rs = 0;
		if (reflectionOperator_.isEnabled()) rs = &reflectionOperator_;
		typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
		typedef PsimagLite::LanczosOrDavidsonBase<ParametersForSolverType,
		        MatrixVectorType,
		        TargetVectorType> LanczosOrDavidsonBaseType;
		typedef PsimagLite::DavidsonSolver<ParametersForSolverType,
		        MatrixVectorType,
		        TargetVectorType> DavidsonSolverType;
		typedef PsimagLite::LanczosSolver<ParametersForSolverType,
		        MatrixVectorType,
		        TargetVectorType> LanczosSolverType;

		typename LanczosOrDavidsonBaseType::MatrixType lanczosHelper(&model_,
		                                                             &modelHelper,
		                                                             rs);

		ParametersForSolverType params(io_,"Lanczos");
		LanczosOrDavidsonBaseType* lanczosOrDavidson = 0;

		bool useDavidson = (parameters_.options.find("useDavidson")!=PsimagLite::String::npos);
		if (useDavidson) {
			lanczosOrDavidson = new DavidsonSolverType(lanczosHelper,params);
		} else {
			lanczosOrDavidson = new LanczosSolverType(lanczosHelper,params);
		}

		if (lanczosHelper.rank()==0) {
			energyTmp=10000;
			if (lanczosOrDavidson) delete lanczosOrDavidson;
			return;
		}

		RealType norma = PsimagLite::norm(initialVector);

		if (!reflectionOperator_.isEnabled()) {
			tmpVec.resize(lanczosHelper.rank());
			try {
				if (fabs(norma)<1e-12) {
					PsimagLite::OstringStream msg;
					msg<<"WARNING: diagonaliseOneBlock: Norm of guess vector is zero, ";
					msg<<"ignoring guess\n";
					progress_.printline(msg,std::cout);
					lanczosOrDavidson->computeGroundState(energyTmp,
					                                      tmpVec);
				} else {
					lanczosOrDavidson->computeGroundState(energyTmp,
					                                      tmpVec,
					                                      initialVector);
				}
			} catch (std::exception& e) {
				PsimagLite::OstringStream msg0;
				msg0<<"Lanczos or Davidson solver failed, ";
				msg0<<"trying with exact diagonalization...";
				progress_.printline(msg0,std::cout);

				VectorRealType eigs(lanczosHelper.rank());
				PsimagLite::Matrix<ComplexOrRealType> fm;
				lanczosHelper.fullDiag(eigs,fm);
				for (SizeType j = 0; j < eigs.size(); ++j)
					tmpVec[j] = fm(j,0);
				energyTmp = eigs[0];

				PsimagLite::OstringStream msg1;
				msg1<<"Found lowest eigenvalue= "<<energyTmp<<" ";
				progress_.printline(msg1,std::cout);
			}

			if (lanczosOrDavidson) delete lanczosOrDavidson;
			return;
		}

		TargetVectorType initialVector1,initialVector2;
		reflectionOperator_.setInitState(initialVector,initialVector1,initialVector2);
		tmpVec.resize(initialVector1.size());
		lanczosOrDavidson->computeGroundState(energyTmp,tmpVec,initialVector1);

		RealType gsEnergy1 = energyTmp;
		TargetVectorType gsVector1 = tmpVec;

		lanczosHelper.reflectionSector(1);
		TargetVectorType gsVector2(initialVector2.size());
		RealType gsEnergy2 = 0;
		lanczosOrDavidson->computeGroundState(gsEnergy2,gsVector2,initialVector2);

		energyTmp=reflectionOperator_.setGroundState(tmpVec,
		                                             gsEnergy1,
		                                             gsVector1,
		                                             gsEnergy2,
		                                             gsVector2);

		if (lanczosOrDavidson) delete lanczosOrDavidson;
	}

	const ParametersType& parameters_;
	const ModelType& model_;
	const bool& verbose_;
	ReflectionSymmetryType& reflectionOperator_;
	InputValidatorType& io_;
	PsimagLite::ProgressIndicator progress_;
	// quantumSector_ needs to be a reference since DmrgSolver will change it
	const SizeType& quantumSector_;
	WaveFunctionTransfType& wft_;
	RealType oldEnergy_;
}; // class Diagonalization
} // namespace Dmrg

/*@}*/
#endif

