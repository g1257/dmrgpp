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
#include "ObservablesOnePointInSitu.h"
#include "TargetSelector.h"
#include "PsiBase64.h"
#include "PrinterInDetail.h"
#include "OutputFileOrNot.h"
#include "TargetingBase.h"
#include "OneSiteTruncation.h"

namespace Dmrg {

//  A class to represent a generic solver for the Dmrg method
template<typename SolverType, typename VectorWithOffsetType_>
class DmrgSolver {

public:

	typedef TargetingBase<SolverType,VectorWithOffsetType_> TargetingType;
	typedef typename TargetingType::ModelType ModelType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef ObservablesOnePointInSitu<typename TargetingType::TargetVectorType>
	ObservablesOnePointInSituType;
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
	typedef typename ModelType::SuperGeometryType SuperGeometryType;
	typedef Checkpoint<ModelType, WaveFunctionTransfType> CheckpointType;
	typedef Recovery<CheckpointType, TargetingType> RecoveryType;
	typedef typename DmrgSerializerType::FermionSignType FermionSignType;
	typedef typename PsimagLite::Vector<BlockType>::Type VectorBlockType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename TargetingType::LanczosSolverType LanczosSolverType;
	typedef PrinterInDetail<LeftRightSuperType> PrinterInDetailType;
	typedef typename DiagonalizationType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BasisWithOperatorsType::QnType QnType;
	typedef typename QnType::PairSizeType PairSizeType;
	typedef typename DiagonalizationType::VectorRealType VectorRealType;
	typedef typename DiagonalizationType::VectorVectorRealType VectorVectorRealType;
	typedef typename TargetingType::VectorVectorVectorWithOffsetType
	VectorVectorVectorWithOffsetType;
	typedef OneSiteTruncation<ModelType, VectorWithOffsetType> OneSiteTruncationType;

	DmrgSolver(ModelType& model, InputValidatorType& ioIn)
	    : model_(model),
	      parameters_(model_.params()),
	      ioIn_(ioIn),
	      appInfo_("DmrgSolver:"),
	      verbose_(parameters_.options.isSet("verbose")),
	      lrs_("pSprime", "pEprime", "pSE"),
	      ioOut_(parameters_.filename,
	             PsimagLite::IoSelector::ACC_TRUNC,
	             parameters_.options.isSet("minimizeDisk")),
	      progress_("DmrgSolver"),
	      stepCurrent_(0),
	      checkpoint_(parameters_, ioIn, model, quantumSector_.size(), false),
	      wft_(parameters_, model.hilbertSize(0)),
	      diagonalization_(parameters_,
	                       model,
	                       verbose_,
	                       ioIn,
	                       quantumSector_,
	                       wft_,
	                       checkpoint_.energies()),
	      truncate_(lrs_, wft_, parameters_, model.superGeometry(), ioOut_),
	      saveData_(!parameters_.options.isSet("noSaveData")),
	      oneSiteTruncation_(lrs_,
	                         model,
	                         ioIn,
	                         ioOut_)
	{
		firstCall_ = true;
		counter_ = 0;

		std::cout<<appInfo_;
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Turning the engine on";
		progress_.printline(msgg, std::cout);
		ioOut_.write(appInfo_, "ApplicationInfo");

		printFiniteLoops(std::cout);

		PsimagLite::PsiBase64::Encode base64encode(ioIn.data());
		ioOut_.write(base64encode, "InputBase64Encoded");
		ioOut_.write(parameters_, "PARAMETERS");
		ioOut_.write(model_, "MODEL");
		ioOut_.write(PsimagLite::IsComplexNumber<ComplexOrRealType>::True, "IsComplex");

		const SizeType n = model_.targetQuantum().size();
		for (SizeType i = 0; i < n; ++i)
			quantumSector_.push_back(model_.targetQuantum().qn(i));		
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

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Turning off the engine.";
		progress_.printline(msgg, std::cout);
	}

	void main(const SuperGeometryType& geometry)
	{
		ioOut_.write(geometry, "GEOMETRY");

		BlockType S,E;
		VectorBlockType X,Y;

		const bool allInSystem = (parameters_.options.isSet("geometryallinsystem"));

		geometry.split(parameters_.sitesPerBlock,S,X,Y,E,allInSystem);
		for (SizeType i=0;i<X.size();i++)
			sitesIndices_.push_back(X[i]);
		for (SizeType i=0;i<Y.size();i++) sitesIndices_.push_back(Y[Y.size()-i-1]);

		TargetSelector<TargetingType> targetSelector(lrs_,
		                                             model_,
		                                             wft_,
		                                             quantumSector_,
		                                             ioIn_);
		TargetingType& psi = targetSelector();

		ioIn_.printUnused(std::cerr);

		MyBasisWithOperators pS("BasisWithOperators.System");
		MyBasisWithOperators pE("BasisWithOperators.Environ");

		if (checkpoint_.isRestart()) {
			PsimagLite::IoSelector::In io(parameters_.checkpoint.filename());

			checkpoint_.read(pS, pE, false);
			psi.read(io, "FinalPsi");
		} else { // move this block elsewhere:

			RealType time = 0;
			pE.setOneSite(E, model_, time);
			pS.setOneSite(S, model_, time);

			infiniteDmrgLoop(X,Y,E,pS,pE,psi);
		}

		RecoveryType recovery(sitesIndices_, checkpoint_, wft_, pS, pE);
		finiteDmrgLoops(pS, pE, psi, recovery);

		inSitu_.init(psi,geometry.numberOfSites());
	}

	const ComplexOrRealType& inSitu(SizeType i) const
	{
		return inSitu_(i);
	}

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
		bool twoSiteDmrg = parameters_.options.isSet("twositedmrg");
		bool extendedPrint = parameters_.options.isSet("extendedPrint");
		PrinterInDetailType printerInDetail(lrs_, extendedPrint);

		lrs_.left(pS);
		lrs_.right(pE);
		checkpoint_.push(pS,pE);

		const SizeType ten = 10;
		const SizeType initialSizeOfHashTable = std::max(ten, parameters_.keptStatesInfinite);

		RealType time = 0; // no time advancement possible in the infiniteDmrgLoop
		VectorVectorRealType energies;
		for (SizeType step=0;step<X.size();step++) {
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Infinite-loop: step="<<step<<" ( of "<<X.size()<<"), ";
			msg<<" size of blk. added="<<X[step].size();
			progress_.printline(msgg, std::cout);
			printerInDetail.print(std::cout, "infinite");

			lrs_.growLeftBlock(model_, pS, X[step], time); // grow system
			bool needsRightPush = false;
			if (step < Y.size()) {
				lrs_.growRightBlock(model_, pE, Y[step], time); // grow environment
				needsRightPush = true;
			}

			progress_.print("Growth done.\n",std::cout);
			lrs_.printSizes("Infinite",std::cout);

			model_.targetQuantum().updateQuantumSector(quantumSector_,
			                                           lrs_.sites(),
			                                           ProgramGlobals::DirectionEnum::INFINITE,
			                                           step,
			                                           parameters_.adjustQuantumNumbers);

			assert(0 < quantumSector_.size()); // used only for SU(2)
			lrs_.setToProduct(initialSizeOfHashTable);

			const BlockType& ystep = findRightBlock(Y,step,E);

			diagonalization_(psi,
			                 energies,
			                 ProgramGlobals::DirectionEnum::INFINITE,
			                 X[step],
			                 ystep);

			truncate_.changeBasisInfinite(pS, pE, psi, parameters_.keptStatesInfinite);

			if (needsRightPush) {
				if (!twoSiteDmrg) checkpoint_.push(pS,pE);
				else checkpoint_.push(lrs_.left(),lrs_.right());
			} else {
				checkpoint_.push((twoSiteDmrg) ? lrs_.left() : pS,
				                 ProgramGlobals::SysOrEnvEnum::SYSTEM);
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
		if (parameters_.options.isSet("nofiniteloops"))
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
		        (parameters_.finiteLoop[indexOfFirstFiniteLoop].stepLength() < 0) ?
		            ProgramGlobals::DirectionEnum::EXPAND_ENVIRON :
		            ProgramGlobals::DirectionEnum::EXPAND_SYSTEM;

		int lastSign = 1;

		stepCurrent_ = recovery.stepCurrent(direction);

		for (SizeType i = indexOfFirstFiniteLoop; i < loopsTotal; ++i)  {

			lastSign = (parameters_.finiteLoop[i].stepLength() < 0) ? -1 : 1;
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Finite loop number "<<i;
			msg<<" with l="<<parameters_.finiteLoop[i].stepLength();
			msg<<" keptStates="<<parameters_.finiteLoop[i].keptStates();
			msg<<". "<<(parameters_.finiteLoop.size()-i)<<" more loops to go.";
			progress_.printline(msgg, std::cout);

			if (i > 0) {
				int signPrev = parameters_.finiteLoop[i - 1].stepLength();
				int signThis = parameters_.finiteLoop[i].stepLength();
				int sign = signPrev*signThis;
				if (sign > 0) {
					if (parameters_.finiteLoop[i].stepLength() > 0) stepCurrent_++;
					if (parameters_.finiteLoop[i].stepLength() < 0) stepCurrent_--;
				} else { // has bounced
					checkForWft((signThis > 0) ? ProgramGlobals::SysOrEnvEnum::SYSTEM :
					                             ProgramGlobals::SysOrEnvEnum::ENVIRON,
					            pS);
				}
			}

			const SizeType wantsOneSiteTrunc = (parameters_.finiteLoop[i].wantsOneSiteTruncation())
			        ? 1 : 0;
			model_.announce("finite loop;" + ttos(wantsOneSiteTrunc));

			finiteStep(pS, pE, i, psi);

			if (psi.end()) break;

			if (recovery.byLoop(i))
				recovery.write(psi,
				               i + 1,
				               stepCurrent_,
				               lastSign,
				               ioOut_,
				               ioIn_.data());
		}

		if (!saveData_) return; //<<== EARLY EXIT HERE

		if (!ioOut_.nonNull()) return; //<<== EARLY EXIT HERE

		checkpoint_.write(pS, pE, ioOut_.handle());

		ioOut_.createGroup("FinalPsi");
		psi.write(sitesIndices_[stepCurrent_], ioOut_.handle(), "FinalPsi");
		ioOut_.write(lastSign, "LastLoopSign");
	}

	void checkForWft(ProgramGlobals::SysOrEnvEnum what,
	                 const MyBasis& pS) const
	{
		SizeType numberOfSites = model_.superGeometry().numberOfSites();
		assert(numberOfSites > 2);
		SizeType last = pS.block().size();
		assert(last > 0);
		SizeType lastSiteOfSystem = pS.block()[--last];
		PsimagLite::String lOrR = "";

		if (what == ProgramGlobals::SysOrEnvEnum::SYSTEM &&
		        lastSiteOfSystem != 0)
			lOrR = "right";

		if (what == ProgramGlobals::SysOrEnvEnum::ENVIRON &&
		        lastSiteOfSystem != numberOfSites - 2)
			lOrR = "left";

		if (lOrR == "") return;

		PsimagLite::String msg = "WARNING: To-the-" + lOrR;
		msg += " movement from middle of lattice not supported\n";
		std::cout<<msg;
		std::cerr<<msg;
	}

	void finiteStep(MyBasisWithOperators &pS,
	                MyBasisWithOperators &pE,
	                SizeType loopIndex,
	                TargetingType& target)
	{
		const bool extendedPrint = parameters_.options.isSet("extendedPrint");
		PrinterInDetailType printerInDetail(lrs_, extendedPrint);
		int stepLength = parameters_.finiteLoop[loopIndex].stepLength();
		SizeType keptStates = parameters_.finiteLoop[loopIndex].keptStates();

		const SizeType ten = 10;
		const SizeType initialSizeOfHashTable = std::max(ten, keptStates);

		ProgramGlobals::DirectionEnum direction = (stepLength < 0) ?
		            ProgramGlobals::DirectionEnum::EXPAND_ENVIRON :
		            ProgramGlobals::DirectionEnum::EXPAND_SYSTEM;

		wft_.setStage(direction);

		SizeType sitesPerBlock = parameters_.sitesPerBlock;
		int stepLengthCorrected = int((stepLength+1-sitesPerBlock)/sitesPerBlock);
		if (stepLength<0)
			stepLengthCorrected = int((stepLength+sitesPerBlock-1)/sitesPerBlock);
		int stepFinal = stepCurrent_ + stepLengthCorrected;

		VectorVectorRealType energies;
		while (true) {

			if (static_cast<SizeType>(stepCurrent_) >= sitesIndices_.size())
				err("stepCurrent_ too large!\n");

			RealType time = target.time();
			SizeType oneSiteTruncSize = 0;

			printerInDetail.print(std::cout, "finite");

			if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				oneSiteTruncSize = lrs_.growLeftBlock(model_,
				                                      pS,
				                                      sitesIndices_[stepCurrent_],
				                                      time);
				BasisWithOperatorsType* dummyBwo =
				        &checkpoint_.shrink(ProgramGlobals::SysOrEnvEnum::ENVIRON);
				target.updateOnSiteForCorners(*dummyBwo); // <-- only updates extreme sites
				lrs_.right(*dummyBwo);
			} else {
				oneSiteTruncSize = lrs_.growRightBlock(model_,
				                                       pE,
				                                       sitesIndices_[stepCurrent_],
				                                       time);
				BasisWithOperatorsType* dummyBwo =
				        &checkpoint_.shrink(ProgramGlobals::SysOrEnvEnum::SYSTEM);
				target.updateOnSiteForCorners(*dummyBwo); // <-- only updates extreme sites
				lrs_.left(*dummyBwo);
			}


			lrs_.printSizes("finite",std::cout);

			model_.targetQuantum().updateQuantumSector(quantumSector_,
			                                           lrs_.sites(),
			                                           direction,
			                                           stepCurrent_,
			                                           parameters_.adjustQuantumNumbers);

			assert(0 < quantumSector_.size()); // used only for SU(2)
			lrs_.setToProduct(initialSizeOfHashTable);

			diagonalization_(target,
			                 energies,
			                 direction,
			                 sitesIndices_[stepCurrent_],
			                 loopIndex);
			printEnergies(energies);

			assert(target.psiConst().size() > 0 || oneSiteTruncSize == 0);

			// METTS may not have a g.s.
			if (target.psiConst().size() > 0) {
				assert(target.psiConst()[0].size() > 0);
				oneSiteTruncation_.update(oneSiteTruncSize, *(target.psiConst()[0][0]), direction);
			}

			changeTruncateAndSerialize(pS,pE,target,keptStates,direction,loopIndex);

			if (finalStep(stepLength, stepFinal)) break;

			if (stepCurrent_ < 0)
				err("DmrgSolver::finiteStep() currentStep_ is negative\n");

			progress_.printMemoryUsage();

			if (target.end()) break;
		}

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
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
		const bool twoSiteDmrg = parameters_.options.isSet("twositedmrg");
		FermionSignType fsS(pS.signs());

		FermionSignType fsE(pE.signs());

		truncate_.changeBasisFinite(pS, pE, target, keptStates, direction);

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			checkpoint_.push((twoSiteDmrg) ? lrs_.left() : pS,
			                 ProgramGlobals::SysOrEnvEnum::SYSTEM);
		else
			checkpoint_.push((twoSiteDmrg) ? lrs_.right() : pE,
			                 ProgramGlobals::SysOrEnvEnum::ENVIRON);

		write(fsS,fsE,target,direction,loopIndex);
	}

	void write(const FermionSignType& fsS,
	           const FermionSignType& fsE,
	           const TargetingType& target,
	           ProgramGlobals::DirectionEnum direction,
	           SizeType loopIndex)
	{
		if (!saveData_) return;

		const FiniteLoop& finiteLoop = parameters_.finiteLoop[loopIndex];

		if (!finiteLoop.wantsSave() && !finiteLoop.wantsMultiSitePush()) return;

		const BlockDiagonalMatrixType& transform = truncate_.transform(direction);
		// FIXME: Serializer will for now save only one psi target
		const VectorVectorVectorWithOffsetType& psi = target.psiConst();
		if (psi.size() == 0)
			err("No psi targets?\n");

		if (psi[0].size() == 0)
			err("No psi[0] targets?\n");

		DmrgSerializerType* ds = new DmrgSerializerType(fsS,
		                                                fsE,
		                                                lrs_,
		                                                psi,
		                                                transform,
		                                                direction);
		if (finiteLoop.wantsMultiSitePush()) {
			target.multiSitePush(ds);
			return;
		}

		typename BasisWithOperatorsType::SaveEnum saveOption2 = (finiteLoop.wantsOnlySlowWft())
		        ? BasisWithOperatorsType::SaveEnum::ALL
		        : BasisWithOperatorsType::SaveEnum::PARTIAL;
		SizeType numberOfSites = model_.superGeometry().numberOfSites();
		PsimagLite::String prefix("Serializer");
		if (ioOut_.nonNull()) {
			ds->write(ioOut_.handle(), prefix, saveOption2, numberOfSites, counter_);
			PsimagLite::String prefixForTarget = TargetingType::buildPrefix(ioOut_.handle(),
			                                                                counter_);
			target.write(sitesIndices_[stepCurrent_], ioOut_.handle(), prefixForTarget);
		}

		++counter_;
		delete ds;
		ds = 0;
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

	void printEnergies(const VectorVectorRealType& energies)
	{
		if (!saveData_ || !ioOut_.nonNull()) return;

		CheckpointType::writeEnergies(firstCall_, "Energies", energies, ioOut_.handle());

		firstCall_ = false;
	}

	const BlockType& findRightBlock(const VectorBlockType& y,
	                                SizeType step,
	                                const BlockType& E) const
	{
		if (step < y.size()) return y[step];

		return E;
	}

	void printFiniteLoops(std::ostream& os) const
	{
		const SizeType n = parameters_.finiteLoop.size();
		os<<"FiniteLoops printing starts\n";
		for (SizeType i = 0; i < n; ++i) {
			parameters_.finiteLoop[i].print(os);
			std::cout<<"\n";
		}

		os<<"FiniteLoops printing ends\n";
	}

	const ModelType& model_;
	const ParametersType& parameters_;
	InputValidatorType& ioIn_;
	PsimagLite::ApplicationInfo appInfo_;
	bool verbose_;
	LeftRightSuperType lrs_;
	OutputFileOrNot ioOut_;
	PsimagLite::ProgressIndicator progress_;
	typename QnType::VectorQnType quantumSector_;
	int stepCurrent_;
	CheckpointType checkpoint_;
	WaveFunctionTransfType wft_;
	VectorBlockType sitesIndices_;
	DiagonalizationType diagonalization_;
	TruncationType truncate_;
	ObservablesOnePointInSituType inSitu_;
	bool saveData_;
	OneSiteTruncationType oneSiteTruncation_;
	static bool firstCall_;
	static SizeType counter_;
}; //class DmrgSolver

template<typename T1, typename T2>
bool DmrgSolver<T1, T2>::firstCall_ = true;

template<typename T1, typename T2>
SizeType DmrgSolver<T1, T2>::counter_ = 0;
} // namespace Dmrg

/*@}*/
#endif
