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
#include "MemoryUsage.h"

namespace Dmrg {

	//!  A class to represent a generic solver for the Dmrg method
	template<template<typename,typename> class InternalProductTemplate,typename TargettingType>
	class DmrgSolver {

		typedef typename TargettingType::ModelType ModelType;
		typedef typename TargettingType::IoType IoType;
		typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
		typedef typename ModelType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;

	public:
		typedef typename  OperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename ModelType::MyBasis MyBasis;
		typedef typename MyBasis::RealType RealType;
		typedef typename MyBasis::BlockType BlockType;
		typedef typename ModelType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename ModelType::ModelHelperType::LeftRightSuperType LeftRightSuperType;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
		typedef typename TargettingType::TargettingParamsType TargettingParamsType;
		typedef typename ModelType::InputValidatorType InputValidatorType;
		typedef ParametersDmrgSolver<RealType,InputValidatorType> ParametersType;
		typedef Diagonalization<ParametersType,TargettingType,InternalProductTemplate> DiagonalizationType;
		typedef typename TargettingType::WaveFunctionTransfType WaveFunctionTransfType;
		typedef Truncation<LeftRightSuperType,ParametersType,TargettingType> TruncationType;
		typedef DmrgSerializer<LeftRightSuperType,VectorWithOffsetType> DmrgSerializerType;
		typedef typename ModelType::GeometryType GeometryType;
		typedef Checkpoint<ParametersType,TargettingType> CheckpointType;
		typedef typename DmrgSerializerType::FermionSignType FermionSignType;
		typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;

		enum {SAVE_TO_DISK=1,DO_NOT_SAVE=0};
		enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};

		DmrgSolver(ParametersDmrgSolver<RealType,InputValidatorType> const &parameters,
		           ModelType const &model,
		           TargettingParamsType& targetStruct)
		: parameters_(parameters),
		  model_(model),
		  targetStruct_(targetStruct),
		  appInfo_("DmrgSolver:"),
		  verbose_(false),
		  lrs_("pSprime","pEprime","pSE"),
		  io_(parameters_.filename),
		  ioIn_(parameters_.filename),
		  progress_("DmrgSolver"),
		  quantumSector_(0),
		  stepCurrent_(0),
		  checkpoint_(parameters_),
		  wft_(parameters_),
		  reflectionOperator_(lrs_,model_.hilbertSize(0),parameters_.useReflectionSymmetry,EXPAND_SYSTEM),
		  diagonalization_(parameters,model,verbose_,
				   reflectionOperator_,io_,quantumSector_,wft_),
		  truncate_(reflectionOperator_,wft_,parameters_,
			    model_.geometry().maxConnections(),verbose_)
		{
			io_.print(appInfo_);
			io_.print("PARAMETERS",parameters_);
			io_.print("TARGETSTRUCT",targetStruct_);
			if (parameters_.options.find("verbose")!=PsimagLite::String::npos) verbose_=true;
		}

		~DmrgSolver()
		{
			io_.print(appInfo_);
		}

		void main(const GeometryType& geometry)
		{
			io_.print("GEOMETRY",geometry);
			if (checkpoint_())
				std::cerr<<"WARNING: Will not check finite loops for consistency while checkpoint is in use\n";
			 else 
				if (parameters_.options.find("nofiniteloops")==PsimagLite::String::npos)
					checkFiniteLoops(parameters_.finiteLoop,geometry.numberOfSites());

			PsimagLite::OstringStream msg;
			msg<<"Turning the engine on";
			progress_.printline(msg,std::cout);

			io_.print("MODEL",model_);
			BlockType S,E;
			typename PsimagLite::Vector<BlockType>::Type X,Y;
			geometry.split(parameters_.sitesPerBlock,S,X,Y,E);
			for (SizeType i=0;i<X.size();i++) 
				sitesIndices_.push_back(X[i]);
			for (SizeType i=0;i<Y.size();i++) sitesIndices_.push_back(Y[Y.size()-i-1]);

			//wft_.init();
			//if (parameters_.options.find("nowft")!=PsimagLite::String::npos) wft_.disable();

			TargettingType psi(lrs_,model_,targetStruct_,wft_,quantumSector_);
			io_.print("PSI",psi);

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

			PsimagLite::OstringStream msg2;
			msg2<<"Turning off the engine.";
			progress_.printline(msg2,std::cout);
		}

	private:

		/** I shall give a procedural description of the DMRG method in the following.
		We start with an initial block $S$ (the initial system) and $E$ (the initial environment). 
		Consider two sets of blocks $X$ and $Y$. 
		We will be adding blocks from $X$ to $S$, one at a time, and from $Y$ to $E$, one at a time. 
		Again, note that $X$ and $Y$ are sets of blocks whereas $S$ and $E$ are blocks. This is shown schematically in Fig.~\\ref{fig:sxye}.
		All sites in $S$, $X$, $Y$ and $E$ are numbered as shown in the figure.
		\\begin{figure}
		\\centering{
		\\includegraphics[width=8cm]{dmrg_sxye}}
		\\caption{Labeling of blocks for the DMRG procedure. Blocks from  vector of blocks X are added one at a time to block $S$ to form 
		the system and blocks from  vector of blocks
		Y are added one at a time to E to form the environment. Blocks are vectors of integers. The integers (numbers at the top of the figure)
		label all sites in a fixed and unique way.\\label{fig:sxye}}
		\\end{figure}
		
		Now we start a loop for the DMRG ``infinite'' algorithm
		 by setting $step=0$ and $\\mathcal{V}_R(S)\\equiv\\mathcal{V}(S)$ and $\\mathcal{V}_R(E)\\equiv\\mathcal{V}(E)$.

		The system is grown by adding the sites in $X_{step}$ to it, and let
		$S'=S\\cup X_{step}$, i.e. the $step$-th block of $X$ to $S$ is added to form the block $S'$; likewise, let $E'=E\\cup Y_{step}$. 
		Let us form the following product Hilbert spaces:
		$\\mathcal{V}(S')=\\mathcal{V}_R(S)\\otimes \\mathcal{V}(X_{step})$ and 
		$\\mathcal{V}(E')=\\mathcal{V}_R(E)\\otimes \\mathcal{V}(Y_{step})$ and their union $\\mathcal{V}(S')\\otimes\\mathcal{V}(E')$ which is disjoint.
		
		Consider $\\hat{H}_{S'\\cup E'}$, the Hamiltonian operator, acting on $\\mathcal{V}(S')\\otimes\\mathcal{V}(E')$.
		Using Lanczos\\ref{sec:lanczos},
		we  diagonalize $\\hat{H}_{S'\\cup E'}$ to obtain its lowest eigenvector:
		\\begin{equation}
		|\\psi\\rangle = \\sum_{\\alpha\\in \\mathcal{V}(S'), \\beta\\in\\mathcal{V}(E')}\\psi_{\\alpha,\\beta}|\\alpha\\rangle\\otimes|\\beta\\rangle,
		\\label{eq:psi}
		\\end{equation}
		where $\\{|\\alpha\\rangle\\}$ is a basis of $\\mathcal{V}(S')$ and $\\{|\\beta\\rangle\\}$ is a basis of $\\mathcal{V}(E')$.
		
		We proceed in the same way for the environment,  diagonalize $\\hat{\\rho}_E$ to obtain ordered
		eigenvectors $w^E$, and define $(H^{ E' {\\rm new\\,\\,basis}})_{\\alpha,\\alpha'}$.
		Now we set $S\\leftarrow S'$, $\\mathcal{V}_R(S)\\leftarrow\\mathcal{V}_R(S')$, 
		$H_{S'}\\leftarrow H_{S}$,
		and similarly for the environment, increase step by one,
		and continue with the growth phase of the algorithm.
		*/
		void infiniteDmrgLoop(
				BlockType const &S,
				const typename PsimagLite::Vector<BlockType>::Type& X,
				const typename PsimagLite::Vector<BlockType>::Type& Y,
				BlockType const &E,
				MyBasisWithOperators &pS,
				MyBasisWithOperators &pE,
				TargettingType& psi)
		{
			bool twoSiteDmrg = (parameters_.options.find("twositedmrg")!=PsimagLite::String::npos);

			checkpoint_.push(pS,pE);

			RealType time = 0; // no time advancement possible in the infiniteDmrgLoop
			for (SizeType step=0;step<X.size();step++) {
				PsimagLite::OstringStream msg;
				msg<<"Infinite-loop: step="<<step<<" ( of "<<Y.size()<<"), size of blk. added="<<Y[step].size();
				progress_.printline(msg,std::cout);

				lrs_.growLeftBlock(model_,pS,X[step],time); // grow system
				lrs_.growRightBlock(model_,pE,Y[step],time); // grow environment

				progress_.print("Growth done.\n",std::cout);
				lrs_.printSizes("Infinite",std::cout);

				updateQuantumSector(lrs_.sites(),INFINITE,step);

				lrs_.setToProduct(quantumSector_);

				diagonalization_(psi,INFINITE,X[step],Y[step]);

				truncate_.changeBasis(pS,pE,psi,parameters_.keptStatesInfinite);

				if (!twoSiteDmrg) checkpoint_.push(pS,pE);
				else checkpoint_.push(lrs_.left(),lrs_.right());

				printMemoryUsage();
			}
			progress_.print("Infinite dmrg loop has been done!\n",std::cout);
		}

		/** 
		In the infinite algorithm, the  number of sites in the
		system and environment grows as more steps are performed.
		After this infinite algorithm, a finite algorithm is applied where the 
		environment is shrunk at the expense of the system, and the system is grown
		at the expense of the environment. During the finite algorithm 
		(\\todo{Section to be written}) phase the total number of sites remains 
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
			if (parameters_.finiteLoop.size()==0)
				throw PsimagLite::RuntimeError("finiteDmrgLoops(...): there are no finite loops! (and nofiniteloops is not set)\n");
			
			// set initial site to add to either system or environment:
			// this is a bit tricky and has been a source of endless bugs
			// basically we have pS on the left and pE on the right, 
			// and we need to determine which site is to be added
			// let us set the initial direction first:
			SizeType direction = EXPAND_SYSTEM;
			if (parameters_.finiteLoop[0].stepLength<0) direction=EXPAND_ENVIRON;
			// all right, now we can get the actual site to add:

			SizeType sitesPerBlock = parameters_.sitesPerBlock;
			typename PsimagLite::Vector<SizeType>::Type siteToAdd(sitesPerBlock);
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
			if (sc<0) throw PsimagLite::RuntimeError("finiteDmrgLoops(...): internal error: siteIndices_\n");
			stepCurrent_ = sc; // phew!!, that's all folks, now bugs, go away!!
			
			for (SizeType i=0;i<parameters_.finiteLoop.size();i++)  {
				PsimagLite::OstringStream msg;
				msg<<"Finite loop number "<<i;
				msg<<" with l="<<parameters_.finiteLoop[i].stepLength;
				msg<<" keptStates="<<parameters_.finiteLoop[i].keptStates;
				msg<<". "<<(parameters_.finiteLoop.size()-i)<<" more loops to go.";
				progress_.printline(msg,std::cout);
				
				if (i>0) {
					int sign = parameters_.finiteLoop[i].stepLength*parameters_.finiteLoop[i-1].stepLength;
					if (sign>0) {
						if (parameters_.finiteLoop[i].stepLength>0) stepCurrent_++;
						if (parameters_.finiteLoop[i].stepLength<0) stepCurrent_--;
					}
				}
				finiteStep(S,E,pS,pE,i,psi);
			}
			checkpoint_.save(pS,pE,io_);
			psi.save(sitesIndices_[stepCurrent_],io_);
		}

		void finiteStep(
				BlockType const &S,
				BlockType const &E,
				MyBasisWithOperators &pS,
				MyBasisWithOperators &pE,
				SizeType loopIndex,
    				TargettingType& target)
		{
			int stepLength = parameters_.finiteLoop[loopIndex].stepLength;
			SizeType keptStates = parameters_.finiteLoop[loopIndex].keptStates;
			int saveOption = (parameters_.finiteLoop[loopIndex].saveOption & 1);
			RealType gsEnergy=0;
			
			SizeType direction=EXPAND_SYSTEM;
			if (stepLength<0) direction=EXPAND_ENVIRON;

			wft_.setStage(direction);

			SizeType sitesPerBlock = parameters_.sitesPerBlock;
			int stepLengthCorrected = int((stepLength+1-sitesPerBlock)/sitesPerBlock);
			if (stepLength<0)
				stepLengthCorrected = int((stepLength+sitesPerBlock-1)/sitesPerBlock);
			int stepFinal = stepCurrent_+stepLengthCorrected;
			
			while(true) {
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
					msg<<" stackE="<<checkpoint_.stackSize(ProgramGlobals::ENVIRON)<< " step="<<stepCurrent_;
					msg<<" loopIndex="<<loopIndex<<" length="<<stepLength<<" StepFinal="<<stepFinal;
					progress_.printline(msg,std::cout);
				}

				updateQuantumSector(lrs_.sites(),direction,stepCurrent_);

				lrs_.setToProduct(quantumSector_);

				bool needsPrinting = (saveOption==SAVE_TO_DISK);
				gsEnergy =diagonalization_(target,direction,sitesIndices_[stepCurrent_],loopIndex,needsPrinting);

				changeTruncateAndSerialize(pS,pE,target,keptStates,direction,saveOption);

				if (finalStep(stepLength,stepFinal)) break;
				if (stepCurrent_<0) throw PsimagLite::RuntimeError("DmrgSolver::finiteStep() currentStep_ is negative\n");

				printMemoryUsage();
				
			}
			if (direction==EXPAND_SYSTEM) {
				pE = lrs_.right();
				
			} else {
				pS = lrs_.left();
			}
			if (saveOption==SAVE_TO_DISK) {
				PsimagLite::String s="#WAVEFUNCTION_ENERGY="+ttos(gsEnergy);
				io_.printline(s);
//				io_.print("#WAVEFUNCTION_ENERGY=",gsEnergy);
			}
		}

		void changeTruncateAndSerialize(MyBasisWithOperators& pS,
						MyBasisWithOperators& pE,
						const TargettingType& target,
						SizeType keptStates,
						SizeType direction,
						SizeType saveOption)
		{
			bool twoSiteDmrg = (parameters_.options.find("twositedmrg")!=PsimagLite::String::npos);
			const typename PsimagLite::Vector<SizeType>::Type& eS = pS.electronsVector();
			FermionSignType fsS(eS);

			const typename PsimagLite::Vector<SizeType>::Type& eE = pS.electronsVector();
			FermionSignType fsE(eE);

			truncate_(pS,pE,target,keptStates,direction);
			PsimagLite::OstringStream msg2;
			msg2<<"#Error="<<truncate_.error();
			io_.printline(msg2);

			if (direction==EXPAND_SYSTEM) {
				checkpoint_.push((twoSiteDmrg) ? lrs_.left() : pS,ProgramGlobals::SYSTEM);
			} else {
				checkpoint_.push((twoSiteDmrg) ? lrs_.right() : pE,ProgramGlobals::ENVIRON);
			}
			if (saveOption==SAVE_TO_DISK)
				serialize(fsS,fsE,target,truncate_.transform(),direction);
		}

		void serialize(const FermionSignType& fsS,
		               const FermionSignType& fsE,
		               const TargettingType& target,
			       const SparseMatrixType& transform,
		               SizeType direction)
		{
			DmrgSerializerType ds(fsS,fsE,lrs_,target.gs(),transform,direction);

			ds.save(io_);

			target.save(sitesIndices_[stepCurrent_],io_);
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
				typename PsimagLite::Vector<SizeType>::Type targetQuantumNumbers(2,0);
				assert(step<parameters_.adjustQuantumNumbers.size());
				targetQuantumNumbers[0]=parameters_.adjustQuantumNumbers[step];
				targetQuantumNumbers[1]=parameters_.adjustQuantumNumbers[step];
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
			typename PsimagLite::Vector<SizeType>::Type targetQuantumNumbers(parameters_.targetQuantumNumbers.size());
			for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++) 
				targetQuantumNumbers[ii]=SizeType(round(parameters_.targetQuantumNumbers[ii]*sites));
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
			typename PsimagLite::Vector<SizeType>::Type targetQuantumNumbers(2);

			if (direction==INFINITE) {
				SizeType totalSites = model_.geometry().numberOfSites();
				targetQuantumNumbers[0]=SizeType(round(parameters_.electronsUp*sites/totalSites));
				targetQuantumNumbers[1]=SizeType(round(parameters_.electronsDown*sites/totalSites));
			} else {
				targetQuantumNumbers[0]=parameters_.electronsUp;
				targetQuantumNumbers[1]=parameters_.electronsDown;
			}

			setQuantumSector(targetQuantumNumbers,direction);
		}

		void setQuantumSector(const typename PsimagLite::Vector<SizeType>::Type& targetQuantumNumbers,SizeType direction)
		{
			PsimagLite::OstringStream msg;
			msg<<"Integer target quantum numbers are: ";
			for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
				msg<<targetQuantumNumbers[ii]<<" ";
			progress_.printline(msg,std::cout);
			if (direction==INFINITE) io_.printVector(targetQuantumNumbers,"TargetedQuantumNumbers");
			quantumSector_=MyBasis::pseudoQuantumNumber(targetQuantumNumbers);
		}

		void printMemoryUsage()
		{
			musage_.update();
			PsimagLite::String vmPeak = musage_.findEntry("VmPeak:");
			PsimagLite::String vmSize = musage_.findEntry("VmSize:");
			PsimagLite::OstringStream msg;
			msg<<"Current virtual memory is "<<vmSize<<" maximum was "<<vmPeak;
			progress_.printline(msg,std::cout);
			PsimagLite::OstringStream msg2;
			msg2<<"Engine clock: "<<musage_.time()<<" seconds";
			progress_.printline(msg2,std::cout);
		}

		PsimagLite::MemoryUsage musage_;
		ParametersDmrgSolver<RealType,InputValidatorType> parameters_;
		const ModelType& model_;
		const TargettingParamsType& targetStruct_;
		PsimagLite::ApplicationInfo appInfo_;
		bool verbose_;
		LeftRightSuperType lrs_;
		typename IoType::Out io_;
		typename IoType::In ioIn_;
		PsimagLite::ProgressIndicator progress_;
		SizeType quantumSector_;
		int stepCurrent_;
		CheckpointType checkpoint_;
		WaveFunctionTransfType wft_;
		typename PsimagLite::Vector<BlockType>::Type sitesIndices_;
		ReflectionSymmetryType reflectionOperator_;
		DiagonalizationType diagonalization_;
		TruncationType truncate_;

	}; //class DmrgSolver
} // namespace Dmrg

/*@}*/
#endif
