/*
Copyright (c) 2009, UT-Battelle, LLC
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
#include "Checkpoint.h"
#include "WaveFunctionTransfFactory.h"
#include "Truncation.h"
#include "ObservablesInSitu.h"

namespace Dmrg {

//  A class to represent a generic solver for the Dmrg method
template<typename TargettingType>
class DmrgSolver {

	typedef typename TargettingType::ModelType ModelType;
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef ObservablesInSitu<typename TargettingType::TargetVectorType> ObservablesInSituType;

public:

	typedef typename  OperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename ModelType::MyBasis MyBasis;
	typedef typename MyBasis::RealType RealType;
	typedef typename MyBasis::BlockType BlockType;
	typedef typename ModelType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename ModelType::ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename MyBasis::BasisDataType BasisDataType;
	typedef typename TargettingType::TargetVectorType TargetVectorType;
	typedef typename TargetVectorType::value_type DensityMatrixElementType;
	typedef typename TargettingType::TargettingParamsType TargettingParamsType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename ModelType::SolverParamsType ParametersType;
	typedef Diagonalization<ParametersType,TargettingType> DiagonalizationType;
	typedef typename TargettingType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef Truncation<LeftRightSuperType,ParametersType,TargettingType> TruncationType;
	typedef DmrgSerializer<LeftRightSuperType,VectorWithOffsetType> DmrgSerializerType;
	typedef typename ModelType::GeometryType GeometryType;
	typedef Checkpoint<ParametersType,TargettingType> CheckpointType;
	typedef typename DmrgSerializerType::FermionSignType FermionSignType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef typename PsimagLite::Vector<BlockType>::Type VectorBlockType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
	      EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
	      INFINITE=WaveFunctionTransfType::INFINITE};

	enum {SAVE_ALL=MyBasis::SAVE_ALL, SAVE_PARTIAL=MyBasis::SAVE_PARTIAL};

	DmrgSolver(ModelType const &model,
	           TargettingParamsType& targetStruct,
	           InputValidatorType& ioIn)
	    : model_(model),
	      parameters_(model_.params()),
	      targetStruct_(targetStruct),
	      ioIn_(ioIn),
	      appInfo_("DmrgSolver:"),
	      verbose_(false),
	      lrs_("pSprime","pEprime","pSE"),
	      ioOut_(parameters_.filename),
	      progress_("DmrgSolver"),
	      quantumSector_(0),
	      stepCurrent_(0),
	      checkpoint_(parameters_),
	      wft_(parameters_),
	      reflectionOperator_(lrs_,
	                          model_.hilbertSize(0),
	                          parameters_.useReflectionSymmetry,
	                          EXPAND_SYSTEM),
	      diagonalization_(parameters_,model,verbose_,
	                       reflectionOperator_,ioIn,quantumSector_,wft_),
	      truncate_(reflectionOperator_,wft_,parameters_,
	                model_.geometry().maxConnections(),verbose_),
	      energy_(0.0)
	{
		ioOut_.print(appInfo_);
		ioOut_.print("PARAMETERS",parameters_);
		ioOut_.print("TARGETSTRUCT",targetStruct_);
		if (parameters_.options.find("verbose")!=PsimagLite::String::npos) verbose_=true;
	}

	~DmrgSolver()
	{
		ioOut_.print(appInfo_);
	}

	void main(const GeometryType& geometry)
	{
		ioOut_.print("GEOMETRY",geometry);
		bool allInSystem = (parameters_.options.find("geometryallinsystem")!=
		        PsimagLite::String::npos);

		if (checkpoint_()) {
			std::cerr<<"WARNING: Will not check finite loops for ";
			std::cerr<<"consistency while checkpoint is in use\n";
		} else {
			if (parameters_.options.find("nofiniteloops")==PsimagLite::String::npos)
				checkFiniteLoops(parameters_.finiteLoop,
				                 geometry.numberOfSites(),
				                 allInSystem);
		}

		PsimagLite::OstringStream msg;
		msg<<"Turning the engine on";
		progress_.printline(msg,std::cout);

		ioOut_.print("MODEL",model_);
		BlockType S,E;
		VectorBlockType X,Y;

		geometry.split(parameters_.sitesPerBlock,S,X,Y,E,allInSystem);
		for (SizeType i=0;i<X.size();i++)
			sitesIndices_.push_back(X[i]);
		for (SizeType i=0;i<Y.size();i++) sitesIndices_.push_back(Y[Y.size()-i-1]);

		//wft_.init();
		//if (parameters_.options.find("nowft")!=PsimagLite::String::npos) wft_.disable();

		TargettingType psi(lrs_,model_,targetStruct_,wft_,quantumSector_,ioIn_);
		ioOut_.print("PSI",psi);

		MyBasisWithOperators pS("pS");
		MyBasisWithOperators pE("pE");

		if (checkpoint_()) {
			checkpoint_.load(pS,pE,psi);
		} else { // move this block elsewhere:
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			SparseMatrixType hmatrix;
			BasisDataType q;

			RealType time = 0;
			model_.setNaturalBasis(creationMatrix,hmatrix,q,E,time);
			pE.setVarious(E,hmatrix,q,creationMatrix);

			model_.setNaturalBasis(creationMatrix,hmatrix,q,S,time);
			pS.setVarious(S,hmatrix,q,creationMatrix);
			infiniteDmrgLoop(S,X,Y,E,pS,pE,psi);
		}

		finiteDmrgLoops(S,E,pS,pE,psi);

		inSitu_.init(psi,geometry.numberOfSites());

		PsimagLite::OstringStream msg2;
		msg2<<"Turning off the engine.";
		progress_.printline(msg2,std::cout);
	}

	const DensityMatrixElementType& inSitu(SizeType i) const
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
		Again, note that $X$ and $Y$ are sets of blocks whereas $S$ and $E$ are blocks. This is shown schematically in Fig.~\ref{fig:sxye}.
		All sites in $S$, $X$, $Y$ and $E$ are numbered as shown in the figure.
		\begin{figure}
		\centering{
		\includegraphics[width=8cm]{dmrg_sxye}}
		\caption{Labeling of blocks for the DMRG procedure. Blocks from  vector of blocks X are added one at a time to block $S$ to form
		the system and blocks from  vector of blocks
		Y are added one at a time to E to form the environment. Blocks are vectors of integers. The integers (numbers at the top of the figure)
		label all sites in a fixed and unique way.\label{fig:sxye}}
		\end{figure}

		Now we start a loop for the DMRG ``infinite'' algorithm
		 by setting $step=0$ and $\mathcal{V}_R(S)\equiv\mathcal{V}(S)$ and $\mathcal{V}_R(E)\equiv\mathcal{V}(E)$.

		The system is grown by adding the sites in $X_{step}$ to it, and let
		$S'=S\cup X_{step}$, i.e. the $step$-th block of $X$ to $S$ is added to form the block $S'$; likewise, let $E'=E\cup Y_{step}$.
		Let us form the following product Hilbert spaces:
		$\mathcal{V}(S')=\mathcal{V}_R(S)\otimes \mathcal{V}(X_{step})$ and
		$\mathcal{V}(E')=\mathcal{V}_R(E)\otimes \mathcal{V}(Y_{step})$ and their union $\mathcal{V}(S')\otimes\mathcal{V}(E')$ which is disjoint.

		Consider $\hat{H}_{S'\cup E'}$, the Hamiltonian operator, acting on $\mathcal{V}(S')\otimes\mathcal{V}(E')$.
		Using Lanczos\ref{sec:lanczos},
		we  diagonalize $\hat{H}_{S'\cup E'}$ to obtain its lowest eigenvector:
		\begin{equation}
		|\psi\rangle = \sum_{\alpha\in \mathcal{V}(S'), \beta\in\mathcal{V}(E')}\psi_{\alpha,\beta}|\alpha\rangle\otimes|\beta\rangle,
		\label{eq:psi}
		\end{equation}
		where $\{|\alpha\rangle\}$ is a basis of $\mathcal{V}(S')$ and $\{|\beta\rangle\}$ is a basis of $\mathcal{V}(E')$.

		We proceed in the same way for the environment,  diagonalize $\hat{\rho}_E$ to obtain ordered
		eigenvectors $w^E$, and define $(H^{ E' {\rm new\,\,basis}})_{\alpha,\alpha'}$.
		Now we set $S\leftarrow S'$, $\mathcal{V}_R(S)\leftarrow\mathcal{V}_R(S')$,
		$H_{S'}\leftarrow H_{S}$,
		and similarly for the environment, increase step by one,
		and continue with the growth phase of the algorithm.
		*/
	void infiniteDmrgLoop(
	        BlockType const &,
	        const VectorBlockType& X,
	        const VectorBlockType& Y,
	        BlockType const &E,
	        MyBasisWithOperators &pS,
	        MyBasisWithOperators &pE,
	        TargettingType& psi)
	{
		bool twoSiteDmrg = (parameters_.options.find("twositedmrg")!=
		        PsimagLite::String::npos);

		lrs_.left(pS);
		lrs_.right(pE);
		checkpoint_.push(pS,pE);

		RealType time = 0; // no time advancement possible in the infiniteDmrgLoop
		for (SizeType step=0;step<X.size();step++) {
			PsimagLite::OstringStream msg;
			msg<<"Infinite-loop: step="<<step<<" ( of "<<X.size()<<"), ";
			msg<<" size of blk. added="<<X[step].size();
			progress_.printline(msg,std::cout);

			lrs_.growLeftBlock(model_,pS,X[step],time); // grow system
			bool needsRightPush = false;
			if (step < Y.size()) {
				lrs_.growRightBlock(model_,pE,Y[step],time); // grow environment
				needsRightPush = true;
			}

			progress_.print("Growth done.\n",std::cout);
			lrs_.printSizes("Infinite",std::cout);

			updateQuantumSector(lrs_.sites(),INFINITE,step);

			lrs_.setToProduct(quantumSector_);

			const BlockType& ystep = findRightBlock(Y,step,E);
			energy_ = diagonalization_(psi,INFINITE,X[step],ystep);
			printEnergy(energy_);

			truncate_.changeBasis(pS,pE,psi,parameters_.keptStatesInfinite);

			if (needsRightPush) {
				if (!twoSiteDmrg) checkpoint_.push(pS,pE);
				else checkpoint_.push(lrs_.left(),lrs_.right());
			} else {
				checkpoint_.push((twoSiteDmrg) ? lrs_.left() : pS,ProgramGlobals::SYSTEM);
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
	void finiteDmrgLoops(BlockType const &S,
	                     BlockType const &E,
	                     MyBasisWithOperators &pS,
	                     MyBasisWithOperators &pE,
	                     TargettingType& psi)
	{
		if (parameters_.options.find("nofiniteloops")!=PsimagLite::String::npos) return;

		if (parameters_.finiteLoop.size()==0) {
			PsimagLite::String msg("finiteDmrgLoops(...): there are no finite loops!");
			throw PsimagLite::RuntimeError(msg + " (and nofiniteloops is not set)\n");
		}

		// set initial site to add to either system or environment:
		// this is a bit tricky and has been a source of endless bugs
		// basically we have pS on the left and pE on the right,
		// and we need to determine which site is to be added
		// let us set the initial direction first:
		SizeType direction = EXPAND_SYSTEM;
		if (parameters_.finiteLoop[0].stepLength<0) direction=EXPAND_ENVIRON;
		// all right, now we can get the actual site to add:

		SizeType sitesPerBlock = parameters_.sitesPerBlock;
		VectorSizeType siteToAdd(sitesPerBlock);
		// left-most site of pE
		for (SizeType j=0;j<siteToAdd.size();j++)
			siteToAdd[j] = pE.block()[j];

		if (direction==EXPAND_ENVIRON) {
			// right-most site of pS
			for (SizeType j=0;j<siteToAdd.size();j++)
				siteToAdd[j] = pS.block()[pS.block().size()-1-j];
		}
		// now stepCurrent_ is such that sitesIndices_[stepCurrent_] = siteToAdd
		// so:
		int sc = PsimagLite::isInVector(sitesIndices_,siteToAdd);
		if (sc<0) {
			PsimagLite::String msg("finiteDmrgLoops(...): internal error: ");
			throw PsimagLite::RuntimeError(msg + "siteIndices_\n");
		}

		stepCurrent_ = sc; // phew!!, that's all folks, now bugs, go away!!

		for (SizeType i=0;i<parameters_.finiteLoop.size();i++)  {
			PsimagLite::OstringStream msg;
			msg<<"Finite loop number "<<i;
			msg<<" with l="<<parameters_.finiteLoop[i].stepLength;
			msg<<" keptStates="<<parameters_.finiteLoop[i].keptStates;
			msg<<". "<<(parameters_.finiteLoop.size()-i)<<" more loops to go.";
			progress_.printline(msg,std::cout);

			if (i>0) {
				int sign = parameters_.finiteLoop[i].stepLength*
				        parameters_.finiteLoop[i-1].stepLength;
				if (sign>0) {
					if (parameters_.finiteLoop[i].stepLength>0) stepCurrent_++;
					if (parameters_.finiteLoop[i].stepLength<0) stepCurrent_--;
				}
			}
			finiteStep(S,E,pS,pE,i,psi);
			if (psi.end()) break;
		}

		checkpoint_.save(pS,pE,ioOut_);
		psi.save(sitesIndices_[stepCurrent_],ioOut_);
	}

	void finiteStep(
	        BlockType const &,
	        BlockType const &,
	        MyBasisWithOperators &pS,
	        MyBasisWithOperators &pE,
	        SizeType loopIndex,
	        TargettingType& target)
	{
		int stepLength = parameters_.finiteLoop[loopIndex].stepLength;
		SizeType keptStates = parameters_.finiteLoop[loopIndex].keptStates;
		int saveOption = parameters_.finiteLoop[loopIndex].saveOption;

		SizeType direction=EXPAND_SYSTEM;
		if (stepLength<0) direction=EXPAND_ENVIRON;

		wft_.setStage(direction);

		SizeType sitesPerBlock = parameters_.sitesPerBlock;
		int stepLengthCorrected = int((stepLength+1-sitesPerBlock)/sitesPerBlock);
		if (stepLength<0)
			stepLengthCorrected = int((stepLength+sitesPerBlock-1)/sitesPerBlock);
		int stepFinal = stepCurrent_+stepLengthCorrected;

		while (true) {
			if (SizeType(stepCurrent_)>=sitesIndices_.size())
				throw PsimagLite::RuntimeError("stepCurrent_ too large!\n");

			RealType time = target.time();
			if (direction==EXPAND_SYSTEM) {
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

			changeTruncateAndSerialize(pS,pE,target,keptStates,direction,saveOption);

			if (finalStep(stepLength,stepFinal)) break;
			if (stepCurrent_<0) {
				PsimagLite::String msg("DmrgSolver::finiteStep()");
				throw PsimagLite::RuntimeError(msg + " currentStep_ is negative\n");
			}

			progress_.printMemoryUsage();

			if (target.end()) break;
		}

		if (direction==EXPAND_SYSTEM) {
			pE = lrs_.right();

		} else {
			pS = lrs_.left();
		}
		if (saveOption & 1) {
			PsimagLite::String s="#WAVEFUNCTION_ENERGY="+ttos(energy_);
			ioOut_.printline(s);
		}
	}

	void changeTruncateAndSerialize(MyBasisWithOperators& pS,
	                                MyBasisWithOperators& pE,
	                                const TargettingType& target,
	                                SizeType keptStates,
	                                SizeType direction,
	                                int saveOption)
	{
		bool twoSiteDmrg = (parameters_.options.find("twositedmrg")!=
		        PsimagLite::String::npos);
		const VectorSizeType& eS = pS.electronsVector();
		FermionSignType fsS(eS);

		const VectorSizeType& eE = pS.electronsVector();
		FermionSignType fsE(eE);

		truncate_(pS,pE,target,keptStates,direction);
		PsimagLite::OstringStream msg2;
		msg2<<"#Error="<<truncate_.error();
		ioOut_.printline(msg2);

		if (direction==EXPAND_SYSTEM) {
			checkpoint_.push((twoSiteDmrg) ? lrs_.left() : pS,ProgramGlobals::SYSTEM);
		} else {
			checkpoint_.push((twoSiteDmrg) ? lrs_.right() : pE,ProgramGlobals::ENVIRON);
		}
		serialize(fsS,fsE,target,truncate_.transform(),direction,saveOption);
	}

	void serialize(const FermionSignType& fsS,
	               const FermionSignType& fsE,
	               const TargettingType& target,
	               const SparseMatrixType& transform,
	               SizeType direction,
	               int saveOption)
	{
		if (!(saveOption & 1)) return;

		DmrgSerializerType ds(fsS,fsE,lrs_,target.gs(),transform,direction);

		SizeType saveOption2 = (saveOption & 4) ? SAVE_ALL : SAVE_PARTIAL;
		ds.save(ioOut_,saveOption2);

		target.save(sitesIndices_[stepCurrent_],ioOut_);
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

	PsimagLite::String getDirection(SizeType dir) const
	{
		if (dir==INFINITE) return  "INFINITE";
		if (dir==EXPAND_ENVIRON) return "EXPAND_ENVIRON";
		return "EXPAND_SYSTEM";
	}

	void updateQuantumSector(SizeType sites,SizeType direction,SizeType step)
	{
		if (direction==INFINITE && parameters_.adjustQuantumNumbers.size()>0) {
			VectorSizeType targetQuantumNumbers(2,0);
			if (2*step+1 >= parameters_.adjustQuantumNumbers.size()) {
				PsimagLite::String msg("adjustQuantumNumbers must be a vector");
				msg += " of size N-2, where N is the TotalNumberOfSites\n";
				throw PsimagLite::RuntimeError(msg);
			}

			targetQuantumNumbers[0]=parameters_.adjustQuantumNumbers[2*step];
			targetQuantumNumbers[1]=parameters_.adjustQuantumNumbers[2*step+1];
			setQuantumSector(targetQuantumNumbers,direction);
			return;
		}
		if (parameters_.targetQuantumNumbers.size()>0) {
			updateQuantumSectorT(sites,direction);
		} else {
			updateQuantumSectorUd(sites,direction);
		}

	}

	void updateQuantumSectorT(SizeType sites,SizeType direction)
	{
		VectorSizeType targetQuantumNumbers(parameters_.targetQuantumNumbers.size());
		for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
			targetQuantumNumbers[ii]=
			        SizeType(round(parameters_.targetQuantumNumbers[ii]*sites));
		if (MyBasis::useSu2Symmetry()) {
			SizeType ne = targetQuantumNumbers[0]+targetQuantumNumbers[1];
			if (ne%2==0) {
				if (targetQuantumNumbers[2]%2!=0) targetQuantumNumbers[2]++;
			} else {
				if (targetQuantumNumbers[2]%2==0) targetQuantumNumbers[2]++;
			}
		}
		setQuantumSector(targetQuantumNumbers,direction);
	}

	void updateQuantumSectorUd(SizeType sites,SizeType direction)
	{
		assert(!MyBasis::useSu2Symmetry());
		VectorSizeType targetQuantumNumbers(2);

		if (direction==INFINITE) {
			SizeType totalSites = model_.geometry().numberOfSites();
			targetQuantumNumbers[0]=
			        SizeType(round(parameters_.electronsUp*sites/totalSites));
			targetQuantumNumbers[1]=
			        SizeType(round(parameters_.electronsDown*sites/totalSites));
		} else {
			targetQuantumNumbers[0]=parameters_.electronsUp;
			targetQuantumNumbers[1]=parameters_.electronsDown;
		}

		setQuantumSector(targetQuantumNumbers,direction);
	}

	void setQuantumSector(const VectorSizeType& targetQuantumNumbers,SizeType direction)
	{
		PsimagLite::OstringStream msg;
		msg<<"Integer target quantum numbers are: ";
		for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
			msg<<targetQuantumNumbers[ii]<<" ";
		progress_.printline(msg,std::cout);
		if (direction==INFINITE)
			ioOut_.printVector(targetQuantumNumbers,"TargetedQuantumNumbers");
		quantumSector_=MyBasis::pseudoQuantumNumber(targetQuantumNumbers);
	}

	void printEnergy(RealType energy)
	{
		PsimagLite::OstringStream msg;
		msg.precision(8);
		msg<<"#Energy="<<energy;
		ioOut_.printline(msg);
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
	const TargettingParamsType& targetStruct_;
	InputValidatorType& ioIn_;
	PsimagLite::ApplicationInfo appInfo_;
	bool verbose_;
	LeftRightSuperType lrs_;
	PsimagLite::IoSimple::Out ioOut_;
	PsimagLite::ProgressIndicator progress_;
	SizeType quantumSector_;
	int stepCurrent_;
	CheckpointType checkpoint_;
	WaveFunctionTransfType wft_;
	VectorBlockType sitesIndices_;
	ReflectionSymmetryType reflectionOperator_;
	DiagonalizationType diagonalization_;
	TruncationType truncate_;
	ObservablesInSituType inSitu_;
	RealType energy_;
}; //class DmrgSolver
} // namespace Dmrg

/*@}*/
#endif

