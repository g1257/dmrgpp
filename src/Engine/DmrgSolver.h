/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file DmrgSolver.h
 *
 *  A class to represent a generic solver for the Dmrg method
 *
 */
#ifndef SOLVER_DMRG_HEADER_H
#define SOLVER_DMRG_HEADER_H

#include "ApplicationInfo.h"
#include "ParametersDmrgSolver.h"
#include "Diagonalization.h"
#include "ProgressIndicator.h"
#include "DmrgSerializer.h"
#include "Recovery.h"
#include "Truncation.h"
#include "ObservablesInSitu.h"
#include "TargetingGroundState.h"
#include "TargetingTimeStep.h"
#include "TargetingDynamic.h"
#include "TargetingAdaptiveDynamic.h"
#include "TargetingCorrection.h"
#include "TargetingCorrectionVector.h"
#include "TargetingMetts.h"
#include "TargetingCorrelations.h"
#include "TargetingInSitu.h"
#include "TargetingRixsStatic.h"
#include "TargetingRixsDynamic.h"
#include "PsiBase64.h"
#include "PrinterInDetail.h"
#include "Io/IoSelector.h"

namespace Dmrg {

//  A class to represent a generic solver for the Dmrg method
template<typename SolverType, typename VectorWithOffsetType_>
class DmrgSolver {

	typedef TargetingBase<SolverType,VectorWithOffsetType_> TargetingType;
	typedef typename TargetingType::ModelType ModelType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef ObservablesInSitu<typename TargetingType::TargetVectorType>
	ObservablesInSituType;

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename  OperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename ModelType::MyBasis MyBasis;
	typedef typename MyBasis::RealType RealType;
	typedef typename MyBasis::BlockType BlockType;
	typedef typename ModelType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelType::ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename TargetingType::TargetVectorType TargetVectorType;
	typedef typename TargetVectorType::value_type ComplexOrRealType;
	typedef typename TargetingType::TargetParamsType TargetParamsType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename ModelType::SolverParamsType ParametersType;
	typedef Diagonalization<ParametersType,TargetingType> DiagonalizationType;
	typedef typename TargetingType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef Truncation<ParametersType,TargetingType> TruncationType;
	typedef DmrgSerializer<LeftRightSuperType,VectorWithOffsetType> DmrgSerializerType;
	typedef typename ModelType::GeometryType GeometryType;
	typedef Checkpoint<ParametersType, TargetingType> CheckpointType;
	typedef Recovery<ParametersType, CheckpointType> RecoveryType;
	typedef typename DmrgSerializerType::FermionSignType FermionSignType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef typename PsimagLite::Vector<BlockType>::Type VectorBlockType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename TargetingType::LanczosSolverType LanczosSolverType;
	typedef TargetingGroundState<LanczosSolverType,VectorWithOffsetType> TargetingGroundStateType;
	typedef TargetingTimeStep<LanczosSolverType,VectorWithOffsetType> TargetingTimeStepType;
	typedef TargetingDynamic<LanczosSolverType,VectorWithOffsetType> TargetingDynamicType;
	typedef TargetingAdaptiveDynamic<LanczosSolverType,VectorWithOffsetType>
	TargetingAdaptiveDynamicType;
	typedef TargetingCorrectionVector<LanczosSolverType,VectorWithOffsetType>
	TargetingCorrectionVectorType;
	typedef TargetingCorrection<LanczosSolverType,VectorWithOffsetType> TargetingCorrectionType;
	typedef TargetingMetts<LanczosSolverType,VectorWithOffsetType> TargetingMettsType;
	typedef TargetingCorrelations<LanczosSolverType,VectorWithOffsetType> TargetingCorrelationsType;
	typedef TargetingInSitu<LanczosSolverType,VectorWithOffsetType> TargetingInSituType;
	typedef TargetingRixsStatic<LanczosSolverType,VectorWithOffsetType> TargetingRixsStaticType;
	typedef TargetingRixsDynamic<LanczosSolverType,VectorWithOffsetType> TargetingRixsDynamicType;
	typedef PrinterInDetail<LeftRightSuperType> PrinterInDetailType;
	typedef typename DiagonalizationType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BasisWithOperatorsType::QnType QnType;
	typedef typename QnType::PairSizeType PairSizeType;

	enum {SAVE_ALL=MyBasis::SAVE_ALL, SAVE_PARTIAL=MyBasis::SAVE_PARTIAL};

	DmrgSolver(ModelType const &model,
	           InputValidatorType& ioIn)
	    : model_(model),
	      parameters_(model_.params()),
	      ioIn_(ioIn),
	      appInfo_("DmrgSolver:"),
	      verbose_(false),
	      lrs_("pSprime", "pEprime", "pSE"),
	      ioOut_(parameters_.filename, PsimagLite::IoSelector::ACC_TRUNC),
	      progress_("DmrgSolver"),
	      quantumSector_(QnType::zero()),
	      stepCurrent_(0),
	      checkpoint_(parameters_, ioIn, model, false),
	      wft_(parameters_),
	      reflectionOperator_(lrs_,
	                          model_.hilbertSize(0),
	                          parameters_.useReflectionSymmetry,
	                          ProgramGlobals::EXPAND_SYSTEM),
	      diagonalization_(parameters_,
	                       model,
	                       verbose_,
	                       reflectionOperator_,
	                       ioIn,
	                       quantumSector_,
	                       wft_,
	                       checkpoint_.energy()),
	      truncate_(reflectionOperator_,
	                wft_,
	                parameters_,
	                model.geometry(),
	                ioOut_),
	      energy_(0.0),
	      saveData_(parameters_.options.find("noSaveData") == PsimagLite::String::npos)
	{
		std::cout<<appInfo_;
		PsimagLite::OstringStream msg;
		msg<<"Turning the engine on";
		progress_.printline(msg,std::cout);
		ioOut_.write(appInfo_, "ApplicationInfo");

		PsimagLite::PsiBase64::Encode base64encode(ioIn.data());
		ioOut_.write(base64encode, "InputBase64Encoded");
		ioOut_.write(parameters_, "PARAMETERS");
		ioOut_.write(model_, "MODEL");
		ioOut_.write(PsimagLite::IsComplexNumber<ComplexOrRealType>::True, "IsComplex");
		if (parameters_.options.find("verbose")!=PsimagLite::String::npos)
			verbose_=true;
	}

	~DmrgSolver()
	{
		SizeType site = 0; // FIXME FOR IMMM
		typename BasisWithOperatorsType::VectorBoolType oddElectrons;
		model_.findOddElectronsOfOneSite(oddElectrons, site);
		ioOut_.write(oddElectrons, "OddElectronsOneSite");

		appInfo_.finalize();
		ioOut_.write(appInfo_, "ApplicationInfo");
		ioOut_.close();

		PsimagLite::OstringStream msg2;
		msg2<<"Turning off the engine.";
		progress_.printline(msg2,std::cout);
	}

	void main(const GeometryType& geometry, PsimagLite::String targeting)
	{
		ioOut_.write(geometry, "GEOMETRY");

		BlockType S,E;
		VectorBlockType X,Y;

		bool allInSystem = (parameters_.options.find("geometryallinsystem")!=
		        PsimagLite::String::npos);

		geometry.split(parameters_.sitesPerBlock,S,X,Y,E,allInSystem);
		for (SizeType i=0;i<X.size();i++)
			sitesIndices_.push_back(X[i]);
		for (SizeType i=0;i<Y.size();i++) sitesIndices_.push_back(Y[Y.size()-i-1]);

		TargetingType* psi = 0;

		if (targeting=="TimeStepTargeting" || targeting == "TargetingAncilla") {
			psi = new TargetingTimeStepType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="DynamicTargeting") {
			psi = new TargetingDynamicType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="AdaptiveDynamicTargeting") {
			psi = new TargetingAdaptiveDynamicType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="CorrectionVectorTargeting") {
			psi = new TargetingCorrectionVectorType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting=="CorrectionTargeting") {
			psi = new TargetingCorrectionType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "GroundStateTargeting") {
			psi = new TargetingGroundStateType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "MettsTargeting") {
			psi = new TargetingMettsType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "TargetingCorrelations") {
			psi = new TargetingCorrelationsType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "TargetingInSitu") {
			psi = new TargetingInSituType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "TargetingRixsStatic") {
			psi = new TargetingRixsStaticType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else if (targeting == "TargetingRixsDynamic") {
			psi = new TargetingRixsDynamicType(lrs_,model_,wft_,quantumSector_,ioIn_);
		} else {
			throw PsimagLite::RuntimeError("Unknown targeting " + targeting + "\n");
		}

		ioIn_.printUnused(std::cerr);

		MyBasisWithOperators pS("BasisWithOperators.System");
		MyBasisWithOperators pE("BasisWithOperators.Environ");

		if (checkpoint_.isRestart()) {
			checkpoint_.read(pS, pE, *psi, false, "FinalPsi");
		} else { // move this block elsewhere:

			RealType time = 0;
			pE.setVarious(E, model_, time);
			pS.setVarious(S, model_, time);

			infiniteDmrgLoop(X,Y,E,pS,pE,*psi);
		}

		RecoveryType recovery(sitesIndices_, ioOut_, checkpoint_, wft_, pS, pE);
		finiteDmrgLoops(pS, pE, *psi, recovery);

		inSitu_.init(*psi,geometry.numberOfSites());

		delete psi;
		psi = 0;
	}

	const ComplexOrRealType& inSitu(SizeType i) const
	{
		return inSitu_(i);
	}

	RealType energy() const { return energy_; }

private:

	/* PSIDOC DmrgSolverInfiniteDmrgLoop
		I shall give a procedural description of the DMRG method in the following.
		We start with an initial block $S$ (the initial system) and $E$ (the initial environment).
		Consider two sets of blocks $X$ and $Y$.
		We will be adding blocks from $X$ to $S$, one at a time, and from $Y$ to $E$, one at a time.
		Again, note that $X$ and $Y$ are sets of blocks whereas $S$ and $E$ are blocks.
This is shown schematically in Fig.~\ref{fig:sxye}.
		All sites in $S$, $X$, $Y$ and $E$ are numbered as shown in the figure.
		\begin{figure}
		\caption{Labeling of blocks for the DMRG procedure. Blocks from  vector of blocks X are
 added one at a time to block $S$ to form
		the system and blocks from  vector of blocks
		Y are added one at a time to E to form the environment. Blocks are vectors of integers.
The integers (numbers at the top of the figure)
		label all sites in a fixed and unique way.\label{fig:sxye}}
		\end{figure}

		Now we start a loop for the DMRG ``infinite'' algorithm
		 by setting $step=0$ and $\mathcal{V}_R(S)\equiv\mathcal{V}(S)$ and
$\mathcal{V}_R(E)\equiv\mathcal{V}(E)$.

		The system is grown by adding the sites in $X_{step}$ to it, and let
		$S'=S\cup X_{step}$, i.e. the $step$-th block of $X$ to $S$ is added to
form the block $S'$; likewise, let $E'=E\cup Y_{step}$.
		Let us form the following product Hilbert spaces:
		$\mathcal{V}(S')=\mathcal{V}_R(S)\otimes \mathcal{V}(X_{step})$ and
		$\mathcal{V}(E')=\mathcal{V}_R(E)\otimes \mathcal{V}(Y_{step})$ and their union
$\mathcal{V}(S')\otimes\mathcal{V}(E')$ which is disjoint.

		Consider $\hat{H}_{S'\cup E'}$, the Hamiltonian operator,
acting on $\mathcal{V}(S')\otimes\mathcal{V}(E')$.
		Using Lanczos\ref{sec:lanczos},
		we  diagonalize $\hat{H}_{S'\cup E'}$ to obtain its lowest eigenvector:
		\begin{equation}
		|\psi\rangle = \sum_{\alpha\in \mathcal{V}(S'),
\beta\in\mathcal{V}(E')}\psi_{\alpha,\beta}|\alpha\rangle\otimes|\beta\rangle,
		\label{eq:psi}
		\end{equation}
		where $\{|\alpha\rangle\}$ is a basis of $\mathcal{V}(S')$ and
$\{|\beta\rangle\}$ is a basis of $\mathcal{V}(E')$.

		We proceed in the same way for the environment,  diagonalize $\hat{\rho}_E$ to
obtain ordered
		eigenvectors $w^E$, and define $(H^{ E' {\rm new\,\,basis}})_{\alpha,\alpha'}$.
		Now we set $S\leftarrow S'$, $\mathcal{V}_R(S)\leftarrow\mathcal{V}_R(S')$,
		$H_{S'}\leftarrow H_{S}$,
		and similarly for the environment, increase step by one,
		and continue with the growth phase of the algorithm.
		*/
	void infiniteDmrgLoop(const VectorBlockType& X,
	                      const VectorBlockType& Y,
	                      BlockType const &E,
	                      MyBasisWithOperators &pS,
	                      MyBasisWithOperators &pE,
	                      TargetingType& psi)
	{
		bool twoSiteDmrg = (parameters_.options.find("twositedmrg")!=
		        PsimagLite::String::npos);
		bool extendedPrint = (parameters_.options.find("extendedPrint") !=
		        PsimagLite::String::npos);
		PrinterInDetailType printerInDetail(lrs_, extendedPrint);

		lrs_.left(pS);
		lrs_.right(pE);
		checkpoint_.push(pS,pE);

		RealType time = 0; // no time advancement possible in the infiniteDmrgLoop
		for (SizeType step=0;step<X.size();step++) {
			PsimagLite::OstringStream msg;
			msg<<"Infinite-loop: step="<<step<<" ( of "<<X.size()<<"), ";
			msg<<" size of blk. added="<<X[step].size();
			progress_.printline(msg,std::cout);
			printerInDetail.print(std::cout, "infinite");

			lrs_.growLeftBlock(model_,pS,X[step],time); // grow system
			bool needsRightPush = false;
			if (step < Y.size()) {
				lrs_.growRightBlock(model_,pE,Y[step],time); // grow environment
				needsRightPush = true;
			}

			progress_.print("Growth done.\n",std::cout);
			lrs_.printSizes("Infinite",std::cout);

			updateQuantumSector(lrs_.sites(),ProgramGlobals::INFINITE,step);

			lrs_.setToProduct(quantumSector_);

			const BlockType& ystep = findRightBlock(Y,step,E);
			energy_ = diagonalization_(psi,ProgramGlobals::INFINITE,X[step],ystep);
			printEnergy(energy_);

			truncate_.changeBasisInfinite(pS, pE, psi, parameters_.keptStatesInfinite);

			if (needsRightPush) {
				if (!twoSiteDmrg) checkpoint_.push(pS,pE);
				else checkpoint_.push(lrs_.left(),lrs_.right());
			} else {
				checkpoint_.push((twoSiteDmrg) ? lrs_.left() : pS,
				                 ProgramGlobals::SYSTEM);
			}

			progress_.printMemoryUsage();
		}
		progress_.print("Infinite dmrg loop has been done!\n",std::cout);
	}

	/* PSIDOC DmrgSolverFiniteDmrgLoops
		In the infinite algorithm, the  number of sites in the
		system and environment grows as more steps are performed.
		After this infinite algorithm, a finite algorithm is applied where the
		environment is shrunk at the expense of the system, and the system is grown
		at the expense of the environment. During the finite algorithm
		(\todo{Section to be written}) phase the total number of sites remains
		constant allowing for a formulation
		of DMRG as a variational method in a basis of matrix product states.
  */
	void finiteDmrgLoops(MyBasisWithOperators &pS,
	                     MyBasisWithOperators &pE,
	                     TargetingType& psi,
	                     RecoveryType& recovery)
	{
		if (parameters_.options.find("nofiniteloops") != PsimagLite::String::npos)
			return;

		SizeType loopsTotal = parameters_.finiteLoop.size();
		if (loopsTotal == 0) {
			PsimagLite::String msg("finiteDmrgLoops(...): there are no finite loops!");
			throw PsimagLite::RuntimeError(msg + " (and nofiniteloops is not set)\n");
		}

		SizeType indexOfFirstFiniteLoop = recovery.indexOfFirstFiniteLoop();

		// let us set the initial direction first:
		assert(indexOfFirstFiniteLoop < parameters_.finiteLoop.size());
		ProgramGlobals::DirectionEnum direction =
		        (parameters_.finiteLoop[indexOfFirstFiniteLoop].stepLength < 0) ?
		            ProgramGlobals::EXPAND_ENVIRON :  ProgramGlobals::EXPAND_SYSTEM;

		int lastSign = 1;

		stepCurrent_ = recovery.stepCurrent(direction);

		for (SizeType i = indexOfFirstFiniteLoop; i < loopsTotal; ++i)  {

			lastSign = (parameters_.finiteLoop[i].stepLength < 0) ? -1 : 1;
			PsimagLite::OstringStream msg;
			msg<<"Finite loop number "<<i;
			msg<<" with l="<<parameters_.finiteLoop[i].stepLength;
			msg<<" keptStates="<<parameters_.finiteLoop[i].keptStates;
			msg<<". "<<(parameters_.finiteLoop.size()-i)<<" more loops to go.";
			progress_.printline(msg,std::cout);

			if (i > 0) {
				int sign = parameters_.finiteLoop[i].stepLength*
				        parameters_.finiteLoop[i-1].stepLength;
				if (sign>0) {
					if (parameters_.finiteLoop[i].stepLength>0) stepCurrent_++;
					if (parameters_.finiteLoop[i].stepLength<0) stepCurrent_--;
				}
			}

			finiteStep(pS, pE, i, psi, recovery);

			if (psi.end()) break;

			if (recovery.byLoop(i))
				recovery.write(psi, i + 1, stepCurrent_, lastSign, ioOut_);
		}

		if (!saveData_) return;

		checkpoint_.write(pS, pE, ioOut_);

		ioOut_.createGroup("FinalPsi");
		psi.write(sitesIndices_[stepCurrent_], ioOut_, "FinalPsi");
		ioOut_.write(lastSign, "LastLoopSign");
	}

	void finiteStep(MyBasisWithOperators &pS,
	                MyBasisWithOperators &pE,
	                SizeType loopIndex,
	                TargetingType& target,
	                RecoveryType& recovery)
	{
		bool extendedPrint = (parameters_.options.find("extendedPrint") !=
		        PsimagLite::String::npos);
		PrinterInDetailType printerInDetail(lrs_, extendedPrint);
		int stepLength = parameters_.finiteLoop[loopIndex].stepLength;
		SizeType keptStates = parameters_.finiteLoop[loopIndex].keptStates;
		int saveOption = parameters_.finiteLoop[loopIndex].saveOption;

		ProgramGlobals::DirectionEnum direction = (stepLength < 0) ?
		            ProgramGlobals::EXPAND_ENVIRON : ProgramGlobals::EXPAND_SYSTEM;

		wft_.setStage(direction);

		SizeType sitesPerBlock = parameters_.sitesPerBlock;
		int stepLengthCorrected = int((stepLength+1-sitesPerBlock)/sitesPerBlock);
		if (stepLength<0)
			stepLengthCorrected = int((stepLength+sitesPerBlock-1)/sitesPerBlock);
		int stepFinal = stepCurrent_ + stepLengthCorrected;

		while (true) {

			if (static_cast<SizeType>(stepCurrent_) >= sitesIndices_.size())
				err("stepCurrent_ too large!\n");

			RealType time = target.time();
			printerInDetail.print(std::cout, "finite");
			if (direction == ProgramGlobals::EXPAND_SYSTEM) {
				lrs_.growLeftBlock(model_,pS,sitesIndices_[stepCurrent_],time);
				lrs_.right(checkpoint_.shrink(ProgramGlobals::ENVIRON,target));
			} else {
				lrs_.growRightBlock(model_,pE,sitesIndices_[stepCurrent_],time);
				lrs_.left(checkpoint_.shrink(ProgramGlobals::SYSTEM,target));
			}

			lrs_.printSizes("finite",std::cout);
			if (verbose_) {
				PsimagLite::OstringStream msg;
				msg<<" stackS="<<checkpoint_.stackSize(ProgramGlobals::SYSTEM);
				msg<<" stackE="<<checkpoint_.stackSize(ProgramGlobals::ENVIRON);
				msg<< " step="<<stepCurrent_;
				msg<<" loopIndex="<<loopIndex<<" length="<<stepLength;
				msg<<" StepFinal="<<stepFinal;
				progress_.printline(msg,std::cout);
			}

			updateQuantumSector(lrs_.sites(),direction,stepCurrent_);

			lrs_.setToProduct(quantumSector_);

			bool needsPrinting = (saveOption & 1);
			energy_ = diagonalization_(target,
			                           direction,
			                           sitesIndices_[stepCurrent_],
			                           loopIndex,
			                           needsPrinting);
			printEnergy(energy_);

			changeTruncateAndSerialize(pS,pE,target,keptStates,direction,loopIndex);

			if (finalStep(stepLength, stepFinal)) break;

			if (stepCurrent_ < 0)
				err("DmrgSolver::finiteStep() currentStep_ is negative\n");

			progress_.printMemoryUsage();

			if (target.end()) break;
			if (recovery.byTime()) {
				int lastSign = (parameters_.finiteLoop[loopIndex].stepLength < 0) ? -1 : 1;
				recovery.write(target, loopIndex, stepCurrent_, lastSign, ioOut_);
			}
		}

		if (direction == ProgramGlobals::EXPAND_SYSTEM)
			pE = lrs_.right();
		else
			pS = lrs_.left();
	}

	void changeTruncateAndSerialize(MyBasisWithOperators& pS,
	                                MyBasisWithOperators& pE,
	                                const TargetingType& target,
	                                SizeType keptStates,
	                                ProgramGlobals::DirectionEnum direction,
	                                SizeType loopIndex)
	{
		bool twoSiteDmrg = (parameters_.options.find("twositedmrg") != PsimagLite::String::npos);
		FermionSignType fsS(pS.signs());

		FermionSignType fsE(pE.signs());

		truncate_.changeBasisFinite(pS, pE, target, keptStates, direction);

		if (direction == ProgramGlobals::EXPAND_SYSTEM)
			checkpoint_.push((twoSiteDmrg) ? lrs_.left() : pS, ProgramGlobals::SYSTEM);
		else
			checkpoint_.push((twoSiteDmrg) ? lrs_.right() : pE, ProgramGlobals::ENVIRON);

		write(fsS,fsE,target,direction,loopIndex);
	}

	void write(const FermionSignType& fsS,
	           const FermionSignType& fsE,
	           const TargetingType& target,
	           ProgramGlobals::DirectionEnum direction,
	           SizeType loopIndex)
	{
		int saveOption = parameters_.finiteLoop[loopIndex].saveOption;

		if (!(saveOption & 1)) return;
		if (!saveData_) return;

		const BlockDiagonalMatrixType& transform = truncate_.transform(direction);
		DmrgSerializerType ds(fsS,fsE,lrs_,target.gs(),transform, direction);

		SizeType saveOption2 = (saveOption & 4) ? SAVE_ALL : SAVE_PARTIAL;
		SizeType numberOfSites = model_.geometry().numberOfSites();

		static SizeType counter = 0;
		PsimagLite::String prefix("Serializer");
		ds.write(ioOut_, prefix, saveOption2, numberOfSites, counter);
		PsimagLite::String prefixForTarget = TargetingType::buildPrefix(ioOut_, counter);
		target.write(sitesIndices_[stepCurrent_], ioOut_, prefixForTarget);
		++counter;
	}

	bool finalStep(int stepLength,int stepFinal)
	{
		if (stepLength<0) {
			stepCurrent_--;
			if (stepCurrent_<=stepFinal) {
				stepCurrent_++; // revert
				return true;
			}
			return false;
		}
		stepCurrent_++;
		if (stepCurrent_>=stepFinal) {
			stepCurrent_--; //revert
			return true;

		}
		return false;
	}

	PsimagLite::String getDirection(ProgramGlobals::DirectionEnum dir) const
	{
		if (dir == ProgramGlobals::INFINITE) return  "INFINITE";
		return (dir == ProgramGlobals::EXPAND_ENVIRON) ?
		            "EXPAND_ENVIRON" : "EXPAND_SYSTEM";
	}

	void updateQuantumSector(SizeType sites,
	                         ProgramGlobals::DirectionEnum direction,
	                         SizeType step)
	{
		SizeType maxSites = model_.geometry().numberOfSites();

		if (direction == ProgramGlobals::INFINITE &&
		        sites < maxSites &&
		        parameters_.adjustQuantumNumbers.size() > step) {
			quantumSector_ = parameters_.adjustQuantumNumbers[step];
			return;
		} else {
			quantumSector_ = model_.targetQuantum().qn;
		}


		quantumSector_.scale(sites,
		                     model_.geometry().numberOfSites(),
		                     direction,
		                     MyBasis::useSu2Symmetry());
	}

	void printEnergy(RealType energy)
	{
		if (!saveData_) return;
		static SizeType counter = 0;
		if (counter == 0) {
			try {
				PsimagLite::IoSelector::In ioIn(ioOut_.filename());
				SizeType x = 0;
				ioIn.read(x, "Energy/Size");
				ioIn.close();
				counter = x;
			} catch (...) {}
		}

		ioOut_.writeVectorEntry(energy, "Energy", counter++);
	}

	const BlockType& findRightBlock(const VectorBlockType& y,
	                                SizeType step,
	                                const BlockType& E) const
	{
		if (step < y.size()) return y[step];

		return E;
	}

	const ModelType& model_;
	const ParametersType& parameters_;
	InputValidatorType& ioIn_;
	PsimagLite::ApplicationInfo appInfo_;
	bool verbose_;
	LeftRightSuperType lrs_;
	PsimagLite::IoSelector::Out ioOut_;
	PsimagLite::ProgressIndicator progress_;
	QnType quantumSector_;
	int stepCurrent_;
	CheckpointType checkpoint_;
	WaveFunctionTransfType wft_;
	VectorBlockType sitesIndices_;
	ReflectionSymmetryType reflectionOperator_;
	DiagonalizationType diagonalization_;
	TruncationType truncate_;
	ObservablesInSituType inSitu_;
	RealType energy_;
	bool saveData_;
}; //class DmrgSolver
} // namespace Dmrg

/*@}*/
#endif

