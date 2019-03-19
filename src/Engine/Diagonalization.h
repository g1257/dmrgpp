/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file Diagonalization.h
 *
 *  FIXME needs doc
 */
#ifndef DIAGONALIZATION_HEADER_H
#define DIAGONALIZATION_HEADER_H
#include "ProgressIndicator.h"
#include "VectorWithOffset.h" // includes the PsimagLite::norm functions
#include "VectorWithOffsets.h" // includes the PsimagLite::norm functions
#include "ProgramGlobals.h"
#include "LanczosSolver.h"
#include "DavidsonSolver.h"
#include "ParametersForSolver.h"
#include "Concurrency.h"
#include "Profiling.h"

namespace Dmrg {

template<typename ParametersType, typename TargetingType>
class Diagonalization {

public:

	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename TargetingType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename TargetingType::ModelType ModelType;
	typedef typename TargetingType::BasisType BasisType;
	typedef typename TargetingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename TargetingType::BlockType BlockType;
	typedef typename TargetingType::TargetVectorType TargetVectorType;
	typedef typename TargetingType::RealType RealType;
	typedef typename TargetingType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename BasisType::QnType QnType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::HamiltonianConnectionType HamiltonianConnectionType;
	typedef typename OperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::ParamsForKroneckerDumperType ParamsForKroneckerDumperType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef typename TargetingType::MatrixVectorType MatrixVectorType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
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

	Diagonalization(const ParametersType& parameters,
	                const ModelType& model,
	                const bool& verbose,
	                ReflectionSymmetryType& reflectionOperator,
	                InputValidatorType& io,
	                const QnType& quantumSector,
	                WaveFunctionTransfType& waveFunctionTransformation,
	                RealType oldEnergy)
	    : parameters_(parameters),
	      model_(model),
	      verbose_(verbose),
	      reflectionOperator_(reflectionOperator),
	      io_(io),
	      progress_("Diag."),
	      quantumSector_(quantumSector),
	      wft_(waveFunctionTransformation),
	      oldEnergy_(oldEnergy)
	{}

	//!PTEX_LABEL{Diagonalization}
	RealType operator()(TargetingType& target,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType& blockLeft,
	                    const BlockType& blockRight)
	{
		PsimagLite::Profiling profiling("Diagonalization", std::cout);
		assert(direction == ProgramGlobals::INFINITE);
		SizeType loopIndex = 0;
		VectorSizeType sectors;
		targetedSymmetrySectors(sectors,target.lrs());
		reflectionOperator_.update(sectors);
		RealType gsEnergy = internalMain_(target,direction,loopIndex,false,blockLeft);
		//  targeting:
		target.evolve(gsEnergy,direction,blockLeft,blockRight,loopIndex);
		wft_.triggerOff(target.lrs());
		return gsEnergy;
	}

	RealType operator()(TargetingType& target,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType& block,
	                    SizeType loopIndex)
	{
		PsimagLite::Profiling profiling("Diagonalization", std::cout);
		assert(direction != ProgramGlobals::INFINITE);

		RealType gsEnergy = internalMain_(target,direction,loopIndex,false,block);
		//  targeting:
		target.evolve(gsEnergy,direction,block,block,loopIndex);
		wft_.triggerOff(target.lrs());
		return gsEnergy;
	}

private:

	void targetedSymmetrySectors(VectorSizeType& mVector,
	                             const LeftRightSuperType& lrs) const
	{
		SizeType total = lrs.super().partition()-1;
		for (SizeType i=0;i<total;i++) {
			if (lrs.super().pseudoQn(i) != quantumSector_) continue;
			mVector.push_back(i);
		}
	}

	RealType internalMain_(TargetingType& target,
	                       ProgramGlobals::DirectionEnum direction,
	                       SizeType loopIndex,
	                       bool,
	                       const VectorSizeType& block)

	{
		PsimagLite::String options = parameters_.options;
		bool findSymmetrySector = (options.find("findSymmetrySector") != PsimagLite::String::npos);
		const LeftRightSuperType& lrs= target.lrs();
		wft_.triggerOn();

		RealType gsEnergy = 0;
		SizeType saveOption = parameters_.finiteLoop[loopIndex].saveOption;
		checkSaveOption(saveOption);

		bool onlyWft = false;
		if (direction != ProgramGlobals::INFINITE)
			onlyWft = ((saveOption & 2)>0);

		bool noguess = ((saveOption & 8) > 0); // bit 3 set means guess is random vector

		if (parameters_.options.find("MettsTargeting")!=PsimagLite::String::npos)
			return gsEnergy;

		if (parameters_.options.find("TargetingAncilla")!=PsimagLite::String::npos)
			onlyWft = true;

		PsimagLite::OstringStream msg0;
		msg0<<"Setting up Hamiltonian basis of size="<<lrs.super().size();
		progress_.printline(msg0,std::cout);

		SizeType total = lrs.super().partition()-1;
		SizeType weightsTotal = 0;
		VectorSizeType weights(total);
		VectorSizeType sectors;
		for (SizeType i = 0; i < total; ++i) {
			SizeType bs = lrs.super().partition(i+1)-lrs.super().partition(i);
			if (verbose_)
				std::cerr<<lrs.super().qnEx(i);

			// Do only one sector unless doing su(2) with j>0, then do all m's
			if (lrs.super().pseudoQn(i) != quantumSector_ && !findSymmetrySector)
				bs = 0;
			else
				sectors.push_back(i);

			weightsTotal += bs;
			weights[i] = bs;
		}

		if (weightsTotal == 0) {
			PsimagLite::String msg("Diagonalization: ");
			msg += "No symmetry sectors found. Perhaps there are too many particles?\n";
			throw PsimagLite::RuntimeError(msg);
		}

		SizeType totalSectors = sectors.size();
		VectorWithOffsetType initialVector(weights, lrs.super());

		target.initialGuess(initialVector, block, noguess);

		typename PsimagLite::Vector<RealType>::Type energySaved(totalSectors);
		typename PsimagLite::Vector<TargetVectorType>::Type vecSaved(totalSectors);

		for (SizeType j = 0; j < totalSectors; ++j) {
			SizeType i = sectors[j];
			PsimagLite::OstringStream msg;
			msg<<"About to diag. sector with";
			msg<<" quantumSector="<<quantumSector_;
			progress_.printline(msg,std::cout);
			TargetVectorType initialVectorBySector(weights[i]);
			initialVector.extract(initialVectorBySector,i);
			RealType norma = PsimagLite::norm(initialVectorBySector);

			if (fabs(norma) < 1e-12) {
				if (onlyWft)
					err("FATAL Norm of initial vector is zero\n");
			} else {
				initialVectorBySector /= norma;
			}

			if (onlyWft) {
				vecSaved[j] = initialVectorBySector;
				gsEnergy = oldEnergy_;
				PsimagLite::OstringStream msg;
				msg<<"Early exit due to user requesting (fast) WFT only, ";
				msg<<"(non updated) energy= "<<gsEnergy;
				progress_.printline(msg,std::cout);
			} else {
				vecSaved[j].resize(initialVectorBySector.size());
				diagonaliseOneBlock(i,
				                    vecSaved[j],
				                    gsEnergy,
				                    lrs,
				                    target.time(),
				                    initialVectorBySector,
				                    saveOption);
			}

			energySaved[j] = gsEnergy;
		}

		// calc gs energy
		if (verbose_ && PsimagLite::Concurrency::root())
			std::cerr<<"About to calc gs energy\n";
		gsEnergy=1e6;
		for (SizeType i = 0; i < totalSectors; ++i)
			if (energySaved[i] < gsEnergy) gsEnergy = energySaved[i];

		PsimagLite::OstringStream msg3;
		msg3<<"Ground state energy= "<<gsEnergy;
		progress_.printline(msg3,std::cout);

		if (verbose_ && PsimagLite::Concurrency::root())
			std::cerr<<"About to calc gs vector\n";

		SizeType counter = 0;
		RealType degeneracyMax = parameters_.degeneracyMax;
		for (SizeType j = 0; j < totalSectors; ++j) {
			if (findSymmetrySector && fabs(energySaved[j] - gsEnergy) > degeneracyMax)
				continue;

			SizeType i = sectors[j];
			SizeType bs = lrs.super().partition(i + 1) - lrs.super().partition(i);

			const QnType& q = lrs.super().qnEx(i);
			PsimagLite::OstringStream msg;
			msg<<"Found targetted symmetry sector in partition "<<i;
			msg<<" of size="<<bs;
			progress_.printline(msg, std::cout);

			PsimagLite::OstringStream msg2;
			msg2<<"Norm of vector is "<<PsimagLite::norm(vecSaved[j]);
			msg2<<" and quantum numbers are ";
			msg2<<q;
			progress_.printline(msg2, std::cout);
			++counter;
		}

		PsimagLite::OstringStream msg4;
		msg4<<"Number of Sectors found "<<counter;
		progress_.printline(msg4, std::cout);

		target.gs().set(vecSaved, sectors, lrs.super());

		if (PsimagLite::Concurrency::root())
			oldEnergy_ = gsEnergy;

		return gsEnergy;
	}

	/** Diagonalise the i-th block of the matrix, return its eigenvectors
			in tmpVec and its eigenvalues in energyTmp
		!PTEX_LABEL{diagonaliseOneBlock} */
	void diagonaliseOneBlock(SizeType partitionIndex,
	                         TargetVectorType& tmpVec,
	                         RealType& energyTmp,
	                         const LeftRightSuperType& lrs,
	                         RealType targetTime,
	                         const TargetVectorType& initialVector,
	                         SizeType saveOption)
	{
		PsimagLite::String options = parameters_.options;
		bool dumperEnabled = (options.find("KroneckerDumper") != PsimagLite::String::npos);
		ParamsForKroneckerDumperType paramsKrDumper(dumperEnabled,
		                                            parameters_.dumperBegin,
		                                            parameters_.dumperEnd,
		                                            parameters_.precision);
		ParamsForKroneckerDumperType* paramsKrDumperPtr = 0;
		if (lrs.super().block().size() == model_.geometry().numberOfSites())
			paramsKrDumperPtr = &paramsKrDumper;

		HamiltonianConnectionType hc(partitionIndex,
		                             lrs,
		                             model_.geometry(),
		                             ModelType::modelLinks(),
		                             targetTime,
		                             paramsKrDumperPtr);

		if (options.find("debugmatrix")!=PsimagLite::String::npos && !(saveOption & 4) ) {
			SparseMatrixType fullm;

			model_.fullHamiltonian(fullm, hc);

			PsimagLite::Matrix<typename SparseMatrixType::value_type> fullm2;
			crsMatrixToFullMatrix(fullm2,fullm);
			if (PsimagLite::isZero(fullm2)) std::cerr<<"Matrix is zero\n";
			if (options.find("printmatrix")!=PsimagLite::String::npos)
				printFullMatrix(fullm,"matrix",1);

			if (!isHermitian(fullm,true))
				throw PsimagLite::RuntimeError("Not hermitian matrix block\n");

			typename PsimagLite::Vector<RealType>::Type eigs(fullm2.rows());
			PsimagLite::diag(fullm2,eigs,'V');
			std::cerr<<"eigs[0]="<<eigs[0]<<"\n";
			if (options.find("test")!=PsimagLite::String::npos)
				throw std::logic_error("Exiting due to option test in the input\n");

			if (options.find("exactdiag")!=PsimagLite::String::npos &&
			        (saveOption & 4) == 0) {
				energyTmp = eigs[0];
				for (SizeType i = 0; i < tmpVec.size(); ++i)
					tmpVec[i] = fullm2(i,0);
				PsimagLite::OstringStream msg;
				msg<<"Uses exact due to user request. ";
				msg<<"Found lowest eigenvalue= "<<energyTmp;
				progress_.printline(msg,std::cout);
				return;
			}
		}

		PsimagLite::OstringStream msg;
		msg<<"I will now diagonalize a matrix of size="<<hc.modelHelper().size();
		progress_.printline(msg,std::cout);
		diagonaliseOneBlock(partitionIndex,
		                    tmpVec,
		                    energyTmp,
		                    hc,
		                    initialVector,
		                    saveOption);
	}

	void diagonaliseOneBlock(SizeType partitionIndex,
	                         TargetVectorType& tmpVec,
	                         RealType &energyTmp,
	                         HamiltonianConnectionType& hc,
	                         const TargetVectorType& initialVector,
	                         SizeType saveOption)
	{
		ReflectionSymmetryType *rs = 0;
		if (reflectionOperator_.isEnabled()) rs = &reflectionOperator_;

		typename LanczosOrDavidsonBaseType::MatrixType lanczosHelper(model_,
		                                                             hc,
		                                                             rs);

		if ((saveOption & 4)>0) {
			energyTmp = slowWft(lanczosHelper, tmpVec, initialVector);
			PsimagLite::OstringStream msg;
			msg<<"Early exit due to user requesting (slow) WFT, energy= "<<energyTmp;
			progress_.printline(msg,std::cout);
			return;
		}

		ParametersForSolverType params(io_,"Lanczos");
		LanczosOrDavidsonBaseType* lanczosOrDavidson = 0;

		bool useDavidson = (parameters_.options.find("useDavidson") !=
		        PsimagLite::String::npos);
		if (useDavidson) {
			lanczosOrDavidson = new DavidsonSolverType(lanczosHelper, params);
		} else {
			lanczosOrDavidson = new LanczosSolverType(lanczosHelper, params);
		}

		if (lanczosHelper.rows()==0) {
			energyTmp=10000;
			PsimagLite::OstringStream msg;
			msg<<"Early exit due to matrix rank being zero.";
			msg<<" BOGUS energy= "<<energyTmp;
			progress_.printline(msg,std::cout);
			if (lanczosOrDavidson) delete lanczosOrDavidson;
			return;
		}

		if (reflectionOperator_.isEnabled())
			err("ReflectionOperator enabled is not longer supported\n");


		try {
			energyTmp = computeLevel(*lanczosOrDavidson,tmpVec,initialVector);
		} catch (std::exception& e) {
			PsimagLite::OstringStream msg0;
			msg0<<e.what()<<"\n";
			msg0<<"Lanczos or Davidson solver failed, ";
			msg0<<"trying with exact diagonalization...";
			progress_.printline(msg0,std::cout);
			progress_.printline(msg0,std::cerr);

			VectorRealType eigs(lanczosHelper.rows());
			PsimagLite::Matrix<ComplexOrRealType> fm;
			lanczosHelper.fullDiag(eigs,fm);
			for (SizeType j = 0; j < eigs.size(); ++j)
				tmpVec[j] = fm(j, 0);
			energyTmp = eigs[0];

			PsimagLite::OstringStream msg1;
			msg1<<"Found lowest eigenvalue= "<<energyTmp<<" ";
			progress_.printline(msg1,std::cout);
		}

		if (lanczosOrDavidson) delete lanczosOrDavidson;
	}

	RealType computeLevel(LanczosOrDavidsonBaseType& object,
	                      TargetVectorType& gsVector,
	                      const TargetVectorType& initialVector) const
	{
		SizeType excited = parameters_.excited;
		RealType norma = PsimagLite::norm(initialVector);
		RealType gsEnergy = 0;
		if (fabs(norma) < 1e-12) {
			PsimagLite::OstringStream msg;
			msg<<"WARNING: diagonaliseOneBlock: Norm of guess vector is zero, ";
			msg<<"ignoring guess\n";
			progress_.printline(msg, std::cout);
			TargetVectorType init(initialVector.size());
			PsimagLite::fillRandom(init);
			object.computeOneState(gsEnergy, gsVector, init, excited);
		} else {
			object.computeOneState(gsEnergy, gsVector, initialVector, excited);
		}

		return gsEnergy;
	}

	RealType slowWft(const typename LanczosOrDavidsonBaseType::MatrixType& object,
	                 TargetVectorType& gsVector,
	                 const TargetVectorType& initialVector) const
	{
		SizeType excited = parameters_.excited;
		if (excited > 0) {
			PsimagLite::OstringStream msg;
			msg<<"FATAL: slowWft: Not possible when excited > 0\n";
			throw PsimagLite::RuntimeError(msg.str());
		}

		RealType norma = PsimagLite::norm(initialVector);
		if (fabs(norma) < 1e-12)
			err("FATAL: slowWft: Norm of guess vector is zero\n");

		object.matrixVectorProduct(gsVector, initialVector);
		RealType gsEnergy = PsimagLite::real(initialVector*gsVector);
		gsVector = initialVector;
		PsimagLite::String options = parameters_.options;
		bool debugmatrix = (options.find("debugmatrix") != PsimagLite::String::npos);

		if (debugmatrix)
			std::cout<<"VECTOR: Printing of BlockOffDiagMatrix not supported\n";

		return gsEnergy;
	}

	void checkSaveOption(SizeType saveOption) const
	{
		bool bit1 = (saveOption & 2);
		bool bit2 = (saveOption & 4);
		if (!bit1 || !bit2) return;
		PsimagLite::OstringStream msg;
		msg<<"FATAL: Third number of triplet cannot have both bits 1 and 2 set\n";
		throw PsimagLite::RuntimeError(msg.str());
	}

	const ParametersType& parameters_;
	const ModelType& model_;
	const bool& verbose_;
	ReflectionSymmetryType& reflectionOperator_;
	InputValidatorType& io_;
	PsimagLite::ProgressIndicator progress_;
	// quantumSector_ needs to be a reference since DmrgSolver will change it
	const QnType& quantumSector_;
	WaveFunctionTransfType& wft_;
	RealType oldEnergy_;
}; // class Diagonalization
} // namespace Dmrg

/*@}*/
#endif

