// BEGIN LICENSE BLOCK
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file DmrgSolver.h
 *
 *  A class to represent a generic solver for the Dmrg method
 *
 */
#ifndef SOLVER_DMRG_HEADER_H
#define SOLVER_DMRG_HEADER_H

#include "BasisWithOperators.h"
#include "ParametersDmrgSolver.h"
#include "LanczosSolver.h"
#include "Diagonalization.h"
#include "ProgressIndicator.h"
#include "FermionSign.h"
#include "DmrgSerializer.h"
#include "Checkpoint.h"
#include "WaveFunctionTransfFactory.h"

namespace Dmrg {

	//!  A class to represent a generic solver for the Dmrg method
	template<
		template<typename,typename> class InternalProductTemplate,
		template<typename,typename,typename,typename> class DensityMatrixTemplate,
		class ModelType,
		class ConcurrencyType,
		class IoType,
  		template<template<typename,typename,typename> class,
  			template<typename,typename> class,
  			template<typename,typename> class,
  			typename,typename,typename,
  			template<typename> class> class TargettingTemplate,
	 	template<typename> class VectorWithOffsetTemplate>
	class DmrgSolver {
			
		typedef typename ModelType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		
	public:
		typedef typename  OperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename ModelType::MyBasis MyBasis;
		typedef typename MyBasis::RealType RealType;
		typedef typename MyBasis::BlockType BlockType;
		typedef typename ModelType::MyBasisWithOperators MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;

		typedef TargettingTemplate<LanczosSolver,InternalProductTemplate,WaveFunctionTransfFactory,
  				ModelType,ConcurrencyType,IoType,VectorWithOffsetTemplate> TargettingType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
		typedef typename TargettingType::TargettingParamsType TargettingParamsType;
		typedef ParametersDmrgSolver<RealType> ParametersType;
		typedef Diagonalization<ParametersType,TargettingType,InternalProductTemplate> DiagonalizationType;
		typedef DensityMatrixTemplate<RealType,MyBasis,MyBasisWithOperators,TargettingType> DensityMatrixType;
		typedef typename DensityMatrixType::BuildingBlockType TransformType;
		typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
		typedef typename TargettingType::WaveFunctionTransfType WaveFunctionTransfType;

		typedef DmrgSerializer<RealType,VectorWithOffsetType,TransformType,MyBasis,FermionSign> DmrgSerializerType;
		typedef typename ModelType::GeometryType GeometryType;
		typedef Checkpoint<ParametersType,TargettingType> CheckpointType;
				
		enum {SAVE_TO_DISK=1,DO_NOT_SAVE=0};
		enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};
		
		DmrgSolver(
				ParametersDmrgSolver<RealType> const &parameters,
				ModelType const &model,
				ConcurrencyType &concurrency,
			  	TargettingParamsType& targetStruct) :
				parameters_(parameters),
				model_(model),
				concurrency_(concurrency),
				targetStruct_(targetStruct),
				verbose_(false),
				useReflection_(false),
				pSprime_("pSprime"),
				pEprime_("pEprime"),
				pSE_("pSE"),
				io_(parameters_.filename,concurrency.rank()),
				ioIn_(parameters_.filename),
				progress_("DmrgSolver",concurrency.rank()),
				quantumSector_(0),
				stepCurrent_(0),
				checkpoint_(parameters_,concurrency.rank()),
				diagonalization_(parameters,model,concurrency,verbose_,
					useReflection_,io_,quantumSector_,waveFunctionTransformation_)
		{
			io_.print(parameters_);
			io_.print(targetStruct_);
			std::string s = utils::getTimeDate();
			io_.print(s);
			if (parameters_.options.find("verbose")!=std::string::npos) verbose_=true;
			if (parameters_.options.find("useReflection")!=std::string::npos)
				useReflection_=true;
			ModelType::SharedMemoryType::setThreads(parameters_.nthreads);
		}

		~DmrgSolver()
		{
			std::string s = utils::getTimeDate();
			io_.print(s);
		}

		void main(const GeometryType& geometry)
		{
			io_.print(geometry);
			if (checkpoint_())
				std::cerr<<"WARNING: Will not check finite loops for consistency while checkpoint is in use\n";
			 else 
				checkFiniteLoops(parameters_.finiteLoop,geometry.numberOfSites());
			std::ostringstream msg;
			msg<<"Turning the engine on";
			progress_.printline(msg,std::cout);

			io_.print(model_);
			BlockType S,E;
			std::vector<BlockType> X,Y;
			geometry.split(S,X,Y,E);
			for (size_t i=0;i<X.size();i++) 
				sitesIndices_.push_back(X[i]);
			for (size_t i=0;i<Y.size();i++) sitesIndices_.push_back(Y[Y.size()-i-1]);

			waveFunctionTransformation_.init(model_.hilbertSize());
			if (parameters_.options.find("nowft")!=std::string::npos) waveFunctionTransformation_.disable();

			
			TargettingType psi(pSprime_,pEprime_,pSE_,model_,targetStruct_,waveFunctionTransformation_);

			MyBasisWithOperators pS("pS");
			MyBasisWithOperators pE("pE");

			if (checkpoint_()) {	
				checkpoint_.load(pS,pE,psi);
			} else { // move this block elsewhere:
				std::vector<OperatorType> creationMatrix;
				SparseMatrixType hmatrix;
				BasisDataType q;
				
				model_.setNaturalBasis(creationMatrix,hmatrix,q,E);
				pE.setVarious(E,hmatrix,q,creationMatrix);
			
				model_.setNaturalBasis(creationMatrix,hmatrix,q,S);
				pS.setVarious(S,hmatrix,q,creationMatrix);
				infiniteDmrgLoop(S,X,Y,E,pS,pE,psi);
			}

			finiteDmrgLoops(S,E,pS,pE,psi);
			
			std::ostringstream msg2;
			msg2<<"Turning off the engine.";
			progress_.printline(msg2,std::cout);
		}

	private:
		ParametersDmrgSolver<RealType> parameters_;
		const ModelType& model_;
		ConcurrencyType& concurrency_;
		const TargettingParamsType& targetStruct_;
		bool verbose_,useReflection_;
		MyBasisWithOperators pSprime_,pEprime_;
		MyBasis pSE_;
		typename IoType::Out io_;
		typename IoType::In ioIn_;
		PsimagLite::ProgressIndicator progress_;
		size_t quantumSector_;
		int stepCurrent_;
		CheckpointType checkpoint_;
		WaveFunctionTransfType waveFunctionTransformation_;
		std::vector<BlockType> sitesIndices_;
		DiagonalizationType diagonalization_;

		void infiniteDmrgLoop(
				BlockType const &S,
				std::vector<BlockType> const &X,
				std::vector<BlockType> const &Y,
				BlockType const &E,
				MyBasisWithOperators &pS,
				MyBasisWithOperators &pE,
				TargettingType& psi)
		{
			int ns,ne;
			checkpoint_.push(pS,pE);
			
			for (size_t step=0;step<X.size();step++) {
				std::ostringstream msg;
				msg<<"Infinite-loop: step="<<step<<" ( of "<<Y.size()<<"), size of blk. added="<<Y[step].size();
				progress_.printline(msg,std::cout);
				
				grow(pSprime_,pS,X[step],GROW_RIGHT); // grow system
				grow(pEprime_,pE,Y[step],GROW_LEFT); // grow environment
					
				
				progress_.print("Growth done.\n",std::cout);
				msg<<"Infinite-loop: sys-env. size="<<pSprime_.size()<<"x"<<pEprime_.size();
				msg<<" and block="<<pSprime_.block().size()<<"+"<<pEprime_.block().size()<<"*";
				progress_.printline(msg,std::cout);
				
				updateQuantumSector(pSprime_.block().size()+pEprime_.block().size());
				
				pSE_.setToProduct(pSprime_,pEprime_,quantumSector_);
				
				diagonalization_(psi,INFINITE,X[step]);

				progress_.print("Truncating basis now...\n",std::cout);
				
				ns=pSprime_.size();
				ne=pEprime_.size();
				
				TransformType transform;
				changeAndTruncateBasis(pS,psi,pSprime_,pEprime_,transform,parameters_.keptStatesInfinite,EXPAND_SYSTEM);
				changeAndTruncateBasis(pE,psi,pEprime_,pSprime_,transform,parameters_.keptStatesInfinite,EXPAND_ENVIRON);
				
				checkpoint_.push(pS,pE);
			}
			progress_.print("Infinite dmrg loop has been done!\n",std::cout);
		}

		void finiteDmrgLoops(
					BlockType const &S,
     					BlockType const &E,
					MyBasisWithOperators &pS,
     					MyBasisWithOperators &pE,
	  				//int l,
       					TargettingType& psi)
		{
			useReflection_=false; // disable reflection symmetry for finite loop if it was enabled:
			if (parameters_.options.find("nofiniteloops")!=std::string::npos) return;
			if (parameters_.finiteLoop.size()==0)
				throw std::runtime_error("finiteDmrgLoops(...): there are no finite loops! (and nofiniteloops is not set)\n");
			
			// set initial site to add to either system or environment:
			// this is a bit tricky and has been a source of endless bugs
			// basically we have pS on the left and pE on the right, 
			// and we need to determine which site is to be added
			// let us set the initial direction first:
			size_t direction = EXPAND_SYSTEM;
			if (parameters_.finiteLoop[0].stepLength<0) direction=EXPAND_ENVIRON;
			// all right, now we can get the actual site to add:

			std::vector<size_t> siteToAdd(1,pE.block()[0]); // left-most site of pE
			if (direction==EXPAND_ENVIRON) {
				siteToAdd[0] = pS.block()[pS.block().size()-1]; // right-most site of pS
			}
			// now stepCurrent_ is such that sitesIndices_[stepCurrent_] = siteToAdd
			// so:
			int sc = utils::isInVector(sitesIndices_,siteToAdd);
			if (sc<0) throw std::runtime_error("finiteDmrgLoops(...): internal error: siteIndices_\n");
			stepCurrent_ = sc; // phew!!, that's all folks, now bugs, go away!!
			
			for (size_t i=0;i<parameters_.finiteLoop.size();i++)  {
				std::ostringstream msg;
				msg<<"Finite loop number "<<i;
				msg<<" with l="<<parameters_.finiteLoop[i].stepLength;
				msg<<" keptStates="<<parameters_.finiteLoop[i].keptStates;
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
		}

		void finiteStep(
				BlockType const &S,
				BlockType const &E,
				MyBasisWithOperators &pS,
				MyBasisWithOperators &pE,
				size_t loopIndex,
    				TargettingType& target)
		{
			static size_t prevDirection = 0;
			
			int stepLength = parameters_.finiteLoop[loopIndex].stepLength;
			size_t keptStates = parameters_.finiteLoop[loopIndex].keptStates;
			int saveOption = (parameters_.finiteLoop[loopIndex].saveOption & 1);
			RealType gsEnergy=0;
			
			size_t direction=EXPAND_SYSTEM;
			if (stepLength<0) direction=EXPAND_ENVIRON;
			//std::cerr<<"PUSHING DIRECTION="<<getDirection(direction)<<"\n";
			int resetCounter = WaveFunctionTransfType::RESET_COUNTER;
			if (prevDirection ==  direction)
				resetCounter = WaveFunctionTransfType::DO_NOT_RESET_COUNTER;
			prevDirection = direction;

			waveFunctionTransformation_.setStage(direction,resetCounter); 

			int stepFinal = stepCurrent_+stepLength;
			
			while(true) {
				
				std::ostringstream msg;
				if (size_t(stepCurrent_)>=sitesIndices_.size()) throw std::runtime_error("stepCurrent_ too large!\n");
				if (direction==EXPAND_SYSTEM) {
					grow(pSprime_,pS,sitesIndices_[stepCurrent_],GROW_RIGHT);             //grow system
					checkpoint_.shrink(pEprime_,CheckpointType::ENVIRON); //shrink env
				} else {
					grow(pEprime_,pE,sitesIndices_[stepCurrent_],GROW_LEFT);   // grow env.
					checkpoint_.shrink(pSprime_,CheckpointType::SYSTEM); // shrink system
				}
				
				msg<<"finite (dir="<<direction<<"): sys-env: "<<pSprime_.size()<<"x"<<pEprime_.size();
				msg<<" and block="<<pSprime_.block().size()<<"+"<<pEprime_.block().size();
				if (verbose_) {
					msg<<" stackS="<<checkpoint_.stackSize(CheckpointType::SYSTEM);
					msg<<" stackE="<<checkpoint_.stackSize(CheckpointType::ENVIRON)<< " step="<<stepCurrent_;
					msg<<" loopIndex="<<loopIndex<<" length="<<stepLength<<" StepFinal="<<stepFinal;
				}
				progress_.printline(msg,std::cout);
				
				updateQuantumSector(pSprime_.block().size()+pEprime_.block().size());
				
				pSE_.setToProduct(pSprime_,pEprime_,quantumSector_);
				
				bool needsPrinting = (saveOption==SAVE_TO_DISK);
				gsEnergy =diagonalization_(target,direction,sitesIndices_[stepCurrent_],loopIndex,needsPrinting);
				
				changeTruncateAndSerialize(pS,pE,target,keptStates,direction,saveOption);
				
				if (finalStep(stepLength,stepFinal)) break;
				if (stepCurrent_<0) throw std::runtime_error("DmrgSolver::finiteStep() currentStep_ is negative\n");
				
			}
			if (direction==EXPAND_SYSTEM) {
				pE = pEprime_;
				
			} else {
				pS = pSprime_;
			}
			if (saveOption==SAVE_TO_DISK) {
				std::string s="#WAVEFUNCTION_ENERGY="+utils::ttos(gsEnergy);
				io_.printline(s);
			}
			checkpoint_.save(pS,pE,loopIndex,io_);
		}

		void changeTruncateAndSerialize(MyBasisWithOperators& pS,MyBasisWithOperators& pE,
			    const TargettingType& target,size_t keptStates,size_t direction,size_t saveOption)
		{
			std::vector<size_t> electronsVector = pS.electronsVector();

			progress_.print("Truncating (env) basis now...\n",std::cout);
			TransformType transform;
			if (direction==EXPAND_SYSTEM) {
				changeAndTruncateBasis(pS,target,pSprime_,pEprime_,transform,keptStates,direction);
				//systemStack_.push(pS);
				checkpoint_.push(pS,CheckpointType::SYSTEM);
			} else {
				changeAndTruncateBasis(pE,target,pEprime_,pSprime_,transform,keptStates,direction);
				//envStack_.push(pE);
				checkpoint_.push(pE,CheckpointType::ENVIRON);
			}
			
			if (saveOption==SAVE_TO_DISK) serialize(electronsVector,target,transform,direction);
		}

		void serialize(const std::vector<size_t>& electronsVector,const TargettingType& target,
			      const TransformType& transform,size_t direction)
		{
			FermionSign fs(electronsVector);
			DmrgSerializerType ds(fs,pSprime_,pEprime_,pSE_,target.gs(),transform,direction);
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

		std::string getDirection(size_t dir) const
		{
			if (dir==INFINITE) return  "INFINITE";
			if (dir==EXPAND_ENVIRON) return "EXPAND_ENVIRON";
			return "EXPAND_SYSTEM";
		}

		//! add block X to basis pS and get basis pSprime
		MyBasisWithOperators grow(MyBasisWithOperators &pSprime,const MyBasisWithOperators &pS,BlockType const &X,int dir)
		{
			SparseMatrixType hmatrix;
			BasisDataType q;
			std::vector<OperatorType> creationMatrix;
			model_.setNaturalBasis(creationMatrix,hmatrix,q,X);
			MyBasisWithOperators Xbasis("Xbasis");
			
			Xbasis.setVarious(X,hmatrix,q,creationMatrix);
			pSprime.setToProduct(pS,Xbasis,dir);
			
			SparseMatrixType matrix=pSprime.hamiltonian();

			if (dir==GROW_RIGHT) model_.addHamiltonianConnection(matrix,pSprime,pS,Xbasis,model_.orbitals());
			else		     model_.addHamiltonianConnection(matrix,pSprime,Xbasis,pS,model_.orbitals());

			pSprime.setHamiltonian(matrix);
			
			return Xbasis;
		}

		void updateQuantumSector(size_t sites)
		{
			std::vector<size_t> targetQuantumNumbers(parameters_.targetQuantumNumbers.size());
			for (size_t ii=0;ii<targetQuantumNumbers.size();ii++) 
				targetQuantumNumbers[ii]=int(parameters_.targetQuantumNumbers[ii]*sites);
			if (MyBasis::useSu2Symmetry()) {
				size_t ne = targetQuantumNumbers[0]+targetQuantumNumbers[1];
				if (ne%2==0) {
					if (targetQuantumNumbers[2]%2!=0) targetQuantumNumbers[2]++;
				} else {
					if (targetQuantumNumbers[2]%2==0) targetQuantumNumbers[2]++;
				}
				std::ostringstream msg;
				msg<<"Updating targets to "<<targetQuantumNumbers[0]<<" "<<	targetQuantumNumbers[1];
				msg<<" "<<targetQuantumNumbers[2];
				progress_.printline(msg,std::cout);
			}
			quantumSector_=MyBasis::pseudoQuantumNumber(targetQuantumNumbers);
		}

		//! Truncate basis 
		void changeAndTruncateBasis(MyBasisWithOperators &rSprime,const TargettingType& target,
			MyBasisWithOperators const &pBasis,MyBasisWithOperators const &pBasisSummed,TransformType& ftransform,
   			size_t keptStates,size_t direction)
		{
			DensityMatrixType dmS(target,pBasis,pBasisSummed,pSE_,direction);
			dmS.check(direction);
			
			if (verbose_ && concurrency_.root()) std::cerr<<"Trying to diagonalize density-matrix with size="<<dmS.rank()<<"\n";
			std::vector<RealType> eigs;
			dmS.diag(eigs,'V',concurrency_);
			dmS.check2(direction);
			
			if (verbose_ && concurrency_.root()) std::cerr<<"Done with density-matrix diag.\n";
			
			//! transform basis: dmS^\dagger * operator matrix * dms
			rSprime = pBasis;
			if (verbose_ && concurrency_.root()) std::cerr<<"About to changeBasis...\n";
			
			RealType error = rSprime.changeBasis(ftransform,dmS(),eigs,keptStates,parameters_,concurrency_);
			
			waveFunctionTransformation_.push(ftransform,direction,rSprime,pBasisSummed,pSE_); //,target.m());

			std::ostringstream msg;
			msg<<"new size of basis="<<rSprime.size();
			progress_.printline(msg,std::cout);
			
			std::ostringstream msg2;
			msg2<<"#Error="<<error;
			io_.printline(msg2);
		}
	}; //class DmrgSolver
} // namespace Dmrg

/*@}*/
#endif
