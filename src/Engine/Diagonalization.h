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
#include "FiniteLoop.h"
#include "PackIndices.h"
#include "PredicateAwesome.h"

namespace Dmrg {

template<typename ParametersType, typename TargetingType>
class Diagonalization {

public:

	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename ParametersType::OptionsType OptionsType;
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
	typedef typename TargetingType::MatrixVectorType MatrixVectorType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
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
	typedef typename PsimagLite::Vector<TargetVectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Vector<VectorVectorType>::Type VectorVectorVectorType;
	typedef FiniteLoop FiniteLoopType;

	Diagonalization(const ParametersType& parameters,
	                const ModelType& model,
	                const bool& verbose,
	                InputValidatorType& io,
	                const typename QnType::VectorQnType& quantumSector,
	                WaveFunctionTransfType& waveFunctionTransformation,
	                const VectorVectorRealType& oldEnergy)
	    : parameters_(parameters),
	      model_(model),
	      verbose_(verbose),
	      io_(io),
	      progress_("Diag."),
	      quantumSector_(quantumSector),
	      wft_(waveFunctionTransformation),
	      oldEnergy_(oldEnergy)
	{}

	//!PTEX_LABEL{Diagonalization}
	void operator()(TargetingType& target,
	                VectorVectorRealType& energies,
	                ProgramGlobals::DirectionEnum direction,
	                const BlockType& blockLeft,
	                const BlockType& blockRight)
	{
		PsimagLite::Profiling profiling("Diagonalization", std::cout);
		assert(direction == ProgramGlobals::DirectionEnum::INFINITE);
		SizeType loopIndex = 0;
		internalMain_(target, energies, direction, loopIndex, blockLeft);
		//  targeting:
		target.evolve(energies[0], direction, blockLeft, blockRight, loopIndex);
		wft_.triggerOff(target.lrs());
	}

	void operator()(TargetingType& target,
	                VectorVectorRealType& energies,
	                ProgramGlobals::DirectionEnum direction,
	                const BlockType& block,
	                SizeType loopIndex)
	{
		PsimagLite::Profiling profiling("Diagonalization", std::cout);
		assert(direction != ProgramGlobals::DirectionEnum::INFINITE);

		internalMain_(target, energies, direction, loopIndex, block);
		//  targeting:
		target.evolve(energies[0], direction, block, block, loopIndex);
		wft_.triggerOff(target.lrs());
	}

private:

	SizeType targetedSymmetrySectors(VectorSizeType& mVector,
	                                 VectorSizeType& compactedWeights,
	                                 const LeftRightSuperType& lrs) const
	{
		const SizeType total = lrs.super().partition() - 1;
		SizeType sum = 0;
		const bool findSymmetrySector = (parameters_.findSymmetrySector != "");

		for (SizeType j = 0; j < quantumSector_.size(); ++j) {
			for (SizeType i = 0; i < total; ++i) {

				if (j == 0 && verbose_)
					std::cerr<<lrs.super().qnEx(i);

				const QnType& qn = lrs.super().pseudoQn(i);
				const bool b1 =  (qn != quantumSector_[j] ||
				        std::find(mVector.begin(), mVector.end(), i) != mVector.end());

				if (findSymmetrySector) {
					if (!passSymmetryConstraints(qn, lrs.super().block().size()))
						continue;
				} else {
					if (b1) continue;
				}

				SizeType bs = lrs.super().partition(i + 1) - lrs.super().partition(i);
				mVector.push_back(i);
				compactedWeights.push_back(bs);
				sum += bs;
			}
		}

		return sum;
	}

	void internalMain_(TargetingType& target,
	                   VectorVectorRealType& energies,
	                   ProgramGlobals::DirectionEnum direction,
	                   SizeType loopIndex,
	                   const VectorSizeType& block)

	{
		const LeftRightSuperType& lrs= target.lrs();
		wft_.triggerOn();

		SizeType numberOfExcited = parameters_.numberOfExcited;
		const FiniteLoopType& finiteLoop = parameters_.finiteLoop[loopIndex];

		bool onlyWft = false;
		if (direction != ProgramGlobals::DirectionEnum::INFINITE)
			onlyWft = finiteLoop.wantsOnlyFastWft();

		bool noguess = finiteLoop.wantsRandomGuess();

		if (parameters_.options.isSet("MettsTargeting"))
			return;

		if (parameters_.options.isSet("TargetingAncilla"))
			onlyWft = true;

		PsimagLite::OstringStream msgg0(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg0 = msgg0();
		msg0<<"Setting up Hamiltonian basis of size="<<lrs.super().size();
		progress_.printline(msgg0, std::cout);

		VectorSizeType sectors;
		VectorSizeType compactedWeights;
		const SizeType weightsTotal = targetedSymmetrySectors(sectors,
		                                                      compactedWeights,
		                                                      lrs);

		if (weightsTotal == 0) {
			PsimagLite::String msg("Diagonalization: ");
			msg += "No symmetry sectors found. Perhaps there are too many particles?\n";
			throw PsimagLite::RuntimeError(msg);
		}

		SizeType totalSectors = sectors.size();
		VectorVectorRealType energySaved(totalSectors);
		VectorVectorVectorType vecSaved(totalSectors);

		target.initPsi(totalSectors, numberOfExcited);

		if (oldEnergy_.size() != totalSectors)
			oldEnergy_.resize(totalSectors);

		for (SizeType j = 0; j < totalSectors; ++j) {

			energySaved[j].resize(numberOfExcited);

			if (oldEnergy_[j].size() != numberOfExcited)
				oldEnergy_[j].resize(numberOfExcited);

			vecSaved[j].resize(numberOfExcited);
		}

		const bool isVwoS = (VectorWithOffsetType::name() == "vectorwithoffsets");
		VectorVectorType onlyForVwoS;
		if (isVwoS)
			target.initialGuess(onlyForVwoS,
			                    block,
			                    noguess,
			                    compactedWeights,
			                    sectors,
			                    lrs.super());

		for (SizeType j = 0; j < totalSectors; ++j) {

			SizeType i = sectors[j];

			if (parameters_.options.isSet("entangler")) {
				assert(j < vecSaved.size());
				assert(vecSaved[j].size() == 1);
				doEntangler(vecSaved[j][0], lrs, block[0], direction, i);
				target.set(vecSaved, sectors, lrs.super());
				energies = energySaved;
				return;
			}

			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"About to diag. sector with";
			msg<<" quantumSector="<<lrs.super().qnEx(i);
			msg<<" and numberOfExcited="<<parameters_.numberOfExcited;
			progress_.printline(msgg, std::cout);

			TargetVectorType* initialBySector = (isVwoS) ? &onlyForVwoS[j]
			                                               : new TargetVectorType;

			if (onlyWft) {
				for (SizeType excitedIndex = 0; excitedIndex < numberOfExcited; ++excitedIndex) {
					if (!isVwoS)
						target.initialGuess(*initialBySector,
						                    block,
						                    noguess,
						                    compactedWeights,
						                    sectors,
						                    j,
						                    excitedIndex,
						                    lrs.super());
					RealType norma = PsimagLite::norm(*initialBySector);

					if (fabs(norma) < 1e-12) {
						err("FATAL Norm of initial vector is zero\n");
					} else {
						*initialBySector /= norma;
					}

					vecSaved[j][excitedIndex] = *initialBySector;
					energySaved[j][excitedIndex] = oldEnergy_[j][excitedIndex];
					PsimagLite::OstringStream msgg(std::cout.precision());
					PsimagLite::OstringStream::OstringStreamType& msg = msgg();
					msg<<"Early exit due to user requesting (fast) WFT only, ";
					msg<<"(non updated) energy= "<<energySaved[j][excitedIndex];
					progress_.printline(msgg, std::cout);
				}
			} else {
				if (!isVwoS)
					target.initialGuess(*initialBySector,
					                    block,
					                    noguess,
					                    compactedWeights,
					                    sectors,
					                    j,
					                    numberOfExcited, // sum all excited
					                    lrs.super());
				RealType norma = PsimagLite::norm(*initialBySector);

				if (fabs(norma) >= 1e-12)
					*initialBySector /= norma;

				for (SizeType excitedIndex = 0; excitedIndex < numberOfExcited; ++excitedIndex)
					vecSaved[j][excitedIndex].resize(initialBySector->size());

				VectorRealType myEnergy;
				diagonaliseOneBlock(myEnergy,
				                    vecSaved[j],
				                    i,
				                    lrs,
				                    target.time(),
				                    *initialBySector,
				                    loopIndex);

				for (SizeType excitedIndex = 0; excitedIndex < numberOfExcited; ++excitedIndex) {
					energySaved[j][excitedIndex] = myEnergy[excitedIndex];
					oldEnergy_[j][excitedIndex] = myEnergy[excitedIndex];
				}
			} // end if

			if (!isVwoS) {
				delete initialBySector;
				initialBySector = 0;
			}

		} // end sectors

		// calc gs energy
		if (verbose_ && PsimagLite::Concurrency::root())
			std::cerr<<"About to calc gs energy\n";
		if (totalSectors == 0)
			err("FATAL: No sectors\n");

		assert(vecSaved.size() == totalSectors);
		SizeType counter = 0;
		for (SizeType j = 0; j < totalSectors; ++j) {
			const SizeType sector = sectors[j];
			SizeType nexcited = energySaved[j].size();
			assert(vecSaved[j].size() == nexcited);
			for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {
				if (vecSaved[j][excitedIndex].size() == 0)
					continue;

				PsimagLite::OstringStream msgg4(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg4 = msgg4();
				msg4<<" Sector["<<j<<"]="<<sectors[j];
				msg4<<" excited="<<excitedIndex;
				msg4<<" sector energy = "<<energySaved[j][excitedIndex];
				progress_.printline(msgg4, std::cout);

				const QnType& q = lrs.super().qnEx(sector);
				const SizeType bs = lrs.super().partition(sector + 1) -
				        lrs.super().partition(sector);
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg<<"Found targetted symmetry sector in partition "<<j;
				msg<<" of size="<<bs;
				progress_.printline(msgg, std::cout);

				PsimagLite::OstringStream msgg2(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
				msg2<<"Norm of vector is "<<PsimagLite::norm(vecSaved[j][excitedIndex]);
				msg2<<" and quantum numbers are ";
				msg2<<q;
				progress_.printline(msgg2, std::cout);
				++counter;
			}

			PsimagLite::OstringStream msgg4(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg4 = msgg4();
			msg4<<"Number of Sectors found "<<counter;
			msg4<<" for numberOfExcited="<<parameters_.numberOfExcited;
			progress_.printline(msgg4, std::cout);
		}

		target.set(vecSaved, sectors, lrs.super());
		energies = energySaved;
	}

	/** Diagonalise the i-th block of the matrix, return its eigenvectors
			in tmpVec and its eigenvalues in energyTmp
		!PTEX_LABEL{diagonaliseOneBlock} */
	void diagonaliseOneBlock(VectorRealType& energyTmp,
	                         VectorVectorType& tmpVec,
	                         SizeType partitionIndex,
	                         const LeftRightSuperType& lrs,
	                         RealType targetTime,
	                         const TargetVectorType& initialVector,
	                         SizeType loopIndex)
	{
		const OptionsType& options = parameters_.options;

		HamiltonianConnectionType hc(lrs,
		                             ModelType::modelLinks(),
		                             targetTime,
		                             model_.superOpHelper());

		const FiniteLoopType& finiteLoop = parameters_.finiteLoop[loopIndex];
		typename ModelHelperType::Aux aux(partitionIndex, lrs);

		if (options.isSet("debugmatrix") && !finiteLoop.wantsOnlySlowWft()) {
			SparseMatrixType fullm;

			model_.fullHamiltonian(fullm, hc, aux);

			PsimagLite::Matrix<typename SparseMatrixType::value_type> fullm2;
			crsMatrixToFullMatrix(fullm2,fullm);
			if (PsimagLite::isZero(fullm2)) std::cerr<<"Matrix is zero\n";
			if (options.isSet("printmatrix"))
				printFullMatrix(fullm,"matrix",1);

			if (!isHermitian(fullm,true))
				throw PsimagLite::RuntimeError("Not hermitian matrix block\n");

			typename PsimagLite::Vector<RealType>::Type eigs(fullm2.rows());
			std::cerr<<"Diagonalizing full matrix of rank "<<fullm2.rows()<<"\n";
			PsimagLite::diag(fullm2,eigs,'V');
			std::cerr<<"eigs[0]="<<eigs[0]<<"\n";
			if (options.isSet("printmatrix")) {
				for (SizeType i = 0; i < eigs.size(); ++i)
					std::cout<<eigs[i]<<" ";
				std::cout<<"\n";
			}

			if (options.isSet("test"))
				throw std::logic_error("Exiting due to option test in the input\n");

			if (options.isSet("exactdiag") && !finiteLoop.wantsOnlySlowWft()) {
				SizeType nexcited = std::min(energyTmp.size(), eigs.size());
				for (SizeType excited = 0; excited < nexcited; ++excited) {
					energyTmp[excited] = eigs[excited];
					for (SizeType i = 0; i < tmpVec.size(); ++i)
						tmpVec[excited][i] = fullm2(i, excited);
					PsimagLite::OstringStream msgg(std::cout.precision());
					PsimagLite::OstringStream::OstringStreamType& msg = msgg();
					msg<<"Uses exact due to user request. ";
					msg<<"Found lowest eigenvalue= "<<energyTmp[0];
					progress_.printline(msgg, std::cout);
				}
				return;
			}
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"I will now diagonalize a matrix of size="<<hc.modelHelper().size(partitionIndex);
		progress_.printline(msgg, std::cout);
		diagonaliseOneBlock(energyTmp,
		                    tmpVec,
		                    hc,
		                    initialVector,
		                    loopIndex,
		                    aux);
	}

	void diagonaliseOneBlock(VectorRealType& energyTmp,
	                         VectorVectorType& tmpVec,
	                         HamiltonianConnectionType& hc,
	                         const TargetVectorType& initialVector,
	                         SizeType loopIndex,
	                         const typename ModelHelperType::Aux& aux)
	{
		const SizeType nexcited = energyTmp.size();
		typename LanczosOrDavidsonBaseType::MatrixType lanczosHelper(model_,
		                                                             hc,
		                                                             aux);

		if (parameters_.finiteLoop[loopIndex].wantsOnlySlowWft()) {
			slowWft(energyTmp, tmpVec, lanczosHelper, initialVector);
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Early exit due to user requesting (slow) WFT, energy= ";
			for (SizeType i = 0; i < nexcited; ++i)
				msg<<energyTmp[i]<<" ";
			progress_.printline(msgg, std::cout);
			return;
		}

		ParametersForSolverType params(io_, "Lanczos", loopIndex);
		LanczosOrDavidsonBaseType* lanczosOrDavidson = 0;

		const bool useDavidson = parameters_.options.isSet("useDavidson");

		if (useDavidson) {
			lanczosOrDavidson = new DavidsonSolverType(lanczosHelper, params);
		} else {
			lanczosOrDavidson = new LanczosSolverType(lanczosHelper, params);
		}

		if (lanczosHelper.rows()==0) {
			static const RealType val = 10000;
			std::fill(energyTmp.begin(), energyTmp.end(), val);
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Early exit due to matrix rank being zero.";
			msg<<" BOGUS energy= "<<val;
			progress_.printline(msgg, std::cout);
			if (lanczosOrDavidson) delete lanczosOrDavidson;
			return;
		}

		try {
			computeAllLevelsBelow(energyTmp, tmpVec, *lanczosOrDavidson, initialVector);
		} catch (std::exception& e) {
			PsimagLite::OstringStream msgg0(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg0 = msgg0();
			msg0<<e.what()<<"\n";
			msg0<<"Lanczos or Davidson solver failed, ";
			msg0<<"trying with exact diagonalization...";
			progress_.printline(msgg0, std::cout);
			progress_.printline(msgg0, std::cerr);

			VectorRealType eigs(lanczosHelper.rows());
			PsimagLite::Matrix<ComplexOrRealType> fm;
			lanczosHelper.fullDiag(eigs,fm);
			SizeType minExcited = std::min(energyTmp.size(), eigs.size());
			for (SizeType excited = 0; excited < minExcited; ++excited) {
				for (SizeType j = 0; j < eigs.size(); ++j)
					tmpVec[excited][j] = fm(j, excited);
				energyTmp[excited] = eigs[excited];
				PsimagLite::OstringStream msgg2(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
				msg2<<"Found eigenvalue["<<excited<<"]= "<<energyTmp[excited];
				progress_.printline(msgg2, std::cout);
			}
		}

		PsimagLite::OstringStream msgg1(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg1 = msgg1();
		msg1<<"Found lowest eigenvalue= "<<energyTmp[0];
		progress_.printline(msgg1, std::cout);

		if (lanczosOrDavidson) delete lanczosOrDavidson;
	}

	void computeAllLevelsBelow(VectorRealType& energyTmp,
	                           VectorVectorType& gsVector,
	                           LanczosOrDavidsonBaseType& object,
	                           const TargetVectorType& initialVector) const
	{
		const SizeType nexcited = gsVector.size();
		RealType norma = PsimagLite::norm(initialVector);
		if (fabs(norma) < 1e-12) {
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"WARNING: diagonaliseOneBlock: Norm of guess vector is zero, ";
			msg<<"ignoring guess\n";
			progress_.printline(msgg, std::cout);
			TargetVectorType init(initialVector.size());
			PsimagLite::fillRandom(init);
			object.computeAllStatesBelow(energyTmp, gsVector,init, nexcited);
		} else {
			object.computeAllStatesBelow(energyTmp, gsVector, initialVector, nexcited);
		}
	}

	void slowWft(VectorRealType& energyTmp,
	             VectorVectorType& gsVector,
	             const typename LanczosOrDavidsonBaseType::MatrixType& object,
	             const TargetVectorType& initialVector) const
	{
		const SizeType nexcited = gsVector.size();

		RealType norma = PsimagLite::norm(initialVector);
		if (fabs(norma) < 1e-12)
			err("FATAL: slowWft: Norm of guess vector is zero\n");

		for (SizeType i = 0; i < nexcited; ++i) {
			object.matrixVectorProduct(gsVector[i], initialVector);
			energyTmp[i] = PsimagLite::real(initialVector*gsVector[i]);
			gsVector[i] = initialVector;
		}

		if (parameters_.options.isSet("debugmatrix"))
			std::cout<<"VECTOR: Printing of BlockOffDiagMatrix not supported\n";
	}

	void doEntangler(TargetVectorType& v,
	                 const LeftRightSuperType& lrs,
	                 SizeType site,
	                 ProgramGlobals::DirectionEnum direction,
	                 SizeType p)
	{
		const SizeType offset = lrs.super().partition(p);
		const SizeType bs = lrs.super().partition(p + 1) - offset;
		PsimagLite::PackIndices packSuper(lrs.left().size());
		const SizeType h = model_.hilbertSize(site);
		const SizeType leftLeftSize = lrs.left().size()/h;
		PsimagLite::PackIndices packLeft(leftLeftSize);
		PsimagLite::PackIndices packRight(h);
		RealType sum = 0;
		bool breakLeft = (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM ||
		                  direction == ProgramGlobals::DirectionEnum::INFINITE);
		bool breakRight = (direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON ||
		                   direction == ProgramGlobals::DirectionEnum::INFINITE);
		bool firstCall = (direction == ProgramGlobals::DirectionEnum::INFINITE
		                  && leftLeftSize == h && lrs.right().size() == h*h);

		v.resize(bs);
		for (SizeType i = 0; i < bs; ++i) {
			SizeType la = 0;
			SizeType br = 0;
			packSuper.unpack(la, br, lrs.super().permutation(i + offset));

			SizeType a = 0;
			SizeType b = 0;

			if (breakLeft) {
				SizeType l = 0;
				packLeft.unpack(l, a, lrs.left().permutation(la));
				if (firstCall && !model_.isCorrectlyPaired(l)) continue;
			}

			if (breakRight) {
				SizeType r = 0;
				packRight.unpack(b, r, lrs.right().permutation(br));
				if (firstCall && !model_.isCorrectlyPaired(r)) continue;
			}

			if (breakLeft && !model_.isCorrectlyPaired(a)) continue;
			if (breakRight && !model_.isCorrectlyPaired(b)) continue;

			v[i] = 1.0;
			sum += 1.0;
		}

		assert(sum > 0);
		const RealType factor = 1.0/sqrt(sum);
		for (SizeType i = 0; i < bs; ++i)
			v[i] *= factor;
	}

	// FIXME TODO: Write two versions: an unscaled one and an scaled one
	bool passSymmetryConstraints(const QnType& qn, SizeType superSize) const
	{
        //const SizeType latticeSize = model_.superGeometry().numberOfSites();
		PsimagLite::String predicate = parameters_.findSymmetrySector;

		const SizeType total = qn.other.size();

		for (SizeType ind = 0; ind < total; ++ind) {
			const SizeType x = qn.other[ind];
			PsimagLite::PredicateAwesome<>::replaceAll(predicate, "%" + ttos(ind), ttos(x));
		}

		PsimagLite::PredicateAwesome<> pAwesome(predicate);
		const RealType superSizeReal = superSize; // conversion from SizeType to RealType
		return pAwesome.isTrue("n", superSizeReal);
	}

	const ParametersType& parameters_;
	const ModelType& model_;
	const bool& verbose_;
	InputValidatorType& io_;
	PsimagLite::ProgressIndicator progress_;
	// quantumSector_ needs to be a reference since DmrgSolver will change it
	const typename QnType::VectorQnType& quantumSector_;
	WaveFunctionTransfType& wft_;
	VectorVectorRealType oldEnergy_;
}; // class Diagonalization
} // namespace Dmrg

/*@}*/
#endif

