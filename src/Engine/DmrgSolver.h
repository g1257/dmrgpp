// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

namespace Dmrg {

	//!  A class to represent a generic solver for the Dmrg method
	template<
		template<typename,typename> class InternalProductTemplate,
		template<typename,typename,typename,typename> class DensityMatrixTemplate,
		class ModelType,
		class ConcurrencyType,
		class IoType,
		template<typename> class WaveFunctionTransformationTemplate,
  		template<typename> class StackTemplate,
    		template<template<typename,typename> class,template<typename,typename> class,
      			typename,typename,typename,typename,template<typename> class> class TargettingTemplate,
	 	template<typename> class VectorWithOffsetTemplate>
	class DmrgSolver {
			
		typedef typename ModelType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		
		std::string SYSTEM_STACK_STRING;
		std::string ENVIRON_STACK_STRING;
	public:
		typedef typename  OperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename ModelType::MyBasis MyBasis;
		typedef typename MyBasis::RealType RealType;
		typedef typename MyBasis::BlockType BlockType;
		typedef typename ModelType::MyBasisWithOperators MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef StackTemplate<MyBasisWithOperators> StackType;
		typedef WaveFunctionTransformationTemplate<MyBasisWithOperators> WaveFunctionTransformationType;
		typedef TargettingTemplate<LanczosSolver,InternalProductTemplate,WaveFunctionTransformationType,
  				ModelType,ConcurrencyType,IoType,VectorWithOffsetTemplate> TargettingType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
		typedef typename TargettingType::TargettingStructureType TargettingStructureType;
		typedef Diagonalization<ParametersDmrgSolver<RealType>,TargettingType,InternalProductTemplate> DiagonalizationType;
		
		enum {SAVE_TO_DISK=1,DO_NOT_SAVE=0};
		enum {EXPAND_ENVIRON=WaveFunctionTransformationType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransformationType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransformationType::INFINITE};
		
		DmrgSolver(
				ParametersDmrgSolver<RealType> const &parameters,
				ModelType const &model,
				ConcurrencyType &concurrency,
			  	TargettingStructureType& targetStruct) :
				
				SYSTEM_STACK_STRING("SystemStack"),
				ENVIRON_STACK_STRING("EnvironStack"),
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
				quantumSector_(-1),
				stepCurrent_(0),
				systemStack_(SYSTEM_STACK_STRING+parameters_.filename),
				envStack_(ENVIRON_STACK_STRING+parameters_.filename),
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

		void main()
		{
			BlockType S,E;
			std::vector<BlockType> X,Y;
			BasisDataType q;
			SparseMatrixType hmatrix;
			
			std::ostringstream msg;
			msg<<"Turning the engine on";
			progress_.printline(msg,std::cout);
				
			
			io_.print(model_);

			model_.setBlocksOfSites(S,X,Y,E); //! split sites into system, environment, X and Y
			for (size_t i=0;i<X.size();i++) 
				sitesIndices_.push_back(X[i]);
			for (size_t i=0;i<Y.size();i++) sitesIndices_.push_back(Y[Y.size()-i-1]);

			std::vector<OperatorType> creationMatrix;
			model_.setNaturalBasis(creationMatrix,hmatrix,q,S);
			MyBasisWithOperators pS("pS",S,hmatrix,q);
			printOneSiteElectrons(pS);
			pS.setOperators(creationMatrix);
			waveFunctionTransformation_.init(hmatrix.rank());
			if (parameters_.options.find("nowft")!=std::string::npos) waveFunctionTransformation_.disable();

			model_.setNaturalBasis(creationMatrix,hmatrix,q,E);
			MyBasisWithOperators pE("pE",E,hmatrix,q);
			pE.setOperators(creationMatrix);
			
			TargettingType psi(pSprime_,pEprime_,pSE_,model_,targetStruct_,waveFunctionTransformation_);
			
			if (parameters_.options.find("checkpoint")!=std::string::npos)
				checkpointLoad(pS,pE,parameters_.checkpoint.index);
			else
				infiniteDmrgLoop(S,X,Y,E,pS,pE,psi);

			stepCurrent_=X.size();

			if (parameters_.options.find("nofiniteloops")!=std::string::npos) return;

			useReflection_=false; // disable reflection symmetry for finite loop if it was enabled:

			finiteDmrgLoops(S,E,pS,pE,X.size(),psi);
			
			std::ostringstream msg2;
			msg2<<"Turning off the engine.";
			progress_.printline(msg2,std::cout);
		}

	private:
		ParametersDmrgSolver<RealType> parameters_;
		const ModelType& model_;
		ConcurrencyType& concurrency_;
		const TargettingStructureType& targetStruct_;
		bool verbose_,useReflection_;
		MyBasisWithOperators pSprime_,pEprime_;
		MyBasis pSE_;
		typename IoType::Out io_;
		typename IoType::In ioIn_;
		ProgressIndicator progress_;
		int quantumSector_;
		int stepCurrent_;
		StackType systemStack_,envStack_;
		WaveFunctionTransformationType waveFunctionTransformation_;
		std::vector<BlockType> sitesIndices_;
		DiagonalizationType diagonalization_;
		
		void infiniteDmrgLoop(BlockType const &S,std::vector<BlockType> const &X,std::vector<BlockType> const &Y,BlockType const &E,
				      MyBasisWithOperators &pS,MyBasisWithOperators &pE,TargettingType& psi)
		{
			int ns,ne;

			// empty the stacks
			systemStack_.empty();
			envStack_.empty();
			
			// infinite dmrg loop
			systemStack_.push(pS);
			envStack_.push(pE);
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
				
				changeAndTruncateBasis(pS,psi,pSprime_,pEprime_,pSE_,parameters_.keptStatesInfinite,EXPAND_SYSTEM);
				changeAndTruncateBasis(pE,psi,pEprime_,pSprime_,pSE_,parameters_.keptStatesInfinite,EXPAND_ENVIRON);
				
				systemStack_.push(pS); 
				envStack_.push(pE);
			}
			progress_.print("Infinite dmrg loop has been done!\n",std::cout);
		}
		
		void finiteDmrgLoops(
					BlockType const &S,
     					BlockType const &E,
					MyBasisWithOperators &pS,
     					MyBasisWithOperators &pE,
	  				int l,
       					TargettingType& psi)
		{
			
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
			int saveOption = parameters_.finiteLoop[loopIndex].saveOption;
			RealType gsEnergy=0;
			
			size_t direction=EXPAND_SYSTEM;
			if (stepLength<0) direction=EXPAND_ENVIRON;
			//std::cerr<<"PUSHING DIRECTION="<<getDirection(direction)<<"\n";
			int resetCounter = WaveFunctionTransformationType::RESET_COUNTER;
			if (prevDirection ==  direction)
				resetCounter = WaveFunctionTransformationType::DO_NOT_RESET_COUNTER;
			prevDirection = direction;

			waveFunctionTransformation_.setStage(direction,resetCounter); 

			int stepFinal = stepCurrent_+stepLength;
			
			while(true) {
				
				std::ostringstream msg;
				if (size_t(stepCurrent_)>=sitesIndices_.size()) throw std::runtime_error("stepCurrent_ too large!\n");
				if (direction==EXPAND_SYSTEM) {
					grow(pSprime_,pS,sitesIndices_[stepCurrent_],GROW_RIGHT);             //grow system
					shrink(pEprime_,envStack_); //shrink env
				} else {
					grow(pEprime_,pE,sitesIndices_[stepCurrent_],GROW_LEFT);   // grow env.
					shrink(pSprime_,systemStack_); // shrink system
				}
				
				msg<<"finite (dir="<<direction<<"): sys-env: "<<pSprime_.size()<<"x"<<pEprime_.size();
				msg<<" and block="<<pSprime_.block().size()<<"+"<<pEprime_.block().size();
				if (verbose_) {
					msg<<" stackS="<<systemStack_.size()<<" stackE="<<envStack_.size()<< " step="<<stepCurrent_;
					msg<<" loopIndex="<<loopIndex<<" length="<<stepLength<<" StepFinal="<<stepFinal;
				}
				progress_.printline(msg,std::cout);
				
				updateQuantumSector(pSprime_.block().size()+pEprime_.block().size());
				
				pSE_.setToProduct(pSprime_,pEprime_,quantumSector_);
				//if (target.gs().size()==0) throw std:runtime_error("DmrgSolver:: target.size==0 before\n");
				bool needsPrinting = (saveOption==SAVE_TO_DISK);
				gsEnergy =diagonalization_(target,direction,sitesIndices_[stepCurrent_],loopIndex,needsPrinting);
				//if (target.gs().size()==0) throw std:runtime_error("DmrgSolver:: target.size==0 after\n");
				
				if (saveOption==SAVE_TO_DISK) {
					FermionSign fs(pS.electronsVector());
					saveToDiskForObserver(fs,pSprime_,pEprime_,pSE_,target);
				}
				
				progress_.print("Truncating (env) basis now...\n",std::cout);
				if (direction==EXPAND_SYSTEM) {
					changeAndTruncateBasis(pS,target,pSprime_,pEprime_,pSE_,keptStates,direction,saveOption);
					systemStack_.push(pS);
				} else {
					changeAndTruncateBasis(pE,target,pEprime_,pSprime_,pSE_,keptStates,direction,saveOption);
					envStack_.push(pE);
				}
				
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
			checkpointSave(pS,pE,loopIndex);
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
			MyBasisWithOperators Xbasis("Xbasis",X,hmatrix,q);
			
			Xbasis.setOperators(creationMatrix);
			pSprime.setToProduct(pS,Xbasis,dir);
			
			SparseMatrixType matrix=pSprime.hamiltonian();

			if (dir==GROW_RIGHT) model_.addHamiltonianConnection(matrix,pSprime,pS,Xbasis,model_.dof(),model_.orbitals());
			else		     model_.addHamiltonianConnection(matrix,pSprime,Xbasis,pS,model_.dof(),model_.orbitals());
			
			pSprime.setHamiltonian(matrix);
			
			return Xbasis;
		}

		//! shrink pSprime (we don't really shrink, we just undo the growth)
		void shrink(MyBasisWithOperators &pSprime,StackType& thisStack)
		{
			thisStack.pop();
			pSprime=thisStack.top();
		}

		void updateQuantumSector(size_t sites)
		{
			if (parameters_.options.find("hasQuantumNumbers")!=std::string::npos) {
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
			} else quantumSector_= -1;
		}

		//! Truncate basis 
		template<typename TargetType>
		void changeAndTruncateBasis(MyBasisWithOperators &rSprime,const TargetType& target,
			MyBasisWithOperators const &pBasis,MyBasisWithOperators const &pBasisSummed,MyBasis const &pSE,
   			size_t keptStates,size_t direction,int saveOption=DO_NOT_SAVE)
		{
			typedef DensityMatrixTemplate<RealType,MyBasis,MyBasisWithOperators,TargetType> DensityMatrixType;
			DensityMatrixType dmS(target,pBasis,pBasisSummed,pSE,direction);
			dmS.check(direction);
			
			if (verbose_ && concurrency_.root()) std::cerr<<"Trying to diagonalize density-matrix with size="<<dmS.rank()<<"\n";
			std::vector<RealType> eigs;
			dmS.diag(eigs,'V',concurrency_);
			dmS.check2(direction);
			
			if (verbose_ && concurrency_.root()) std::cerr<<"Done with density-matrix diag.\n";
			
			//! transform basis: dmS^\dagger * operator matrix * dms
			rSprime = pBasis;
			typename DensityMatrixType::BuildingBlockType ftransform;
			if (verbose_ && concurrency_.root()) std::cerr<<"About to changeBasis...\n";
			
			RealType error = rSprime.changeBasis(ftransform,dmS(),eigs,keptStates,parameters_,concurrency_);
			if (saveOption==SAVE_TO_DISK) saveToDisk(ftransform,rSprime,direction);
			waveFunctionTransformation_.push(ftransform,direction,rSprime,pBasisSummed,pSE); //,target.m());

			std::ostringstream msg;
			msg<<"new size of basis="<<rSprime.size();
			progress_.printline(msg,std::cout);
			
			std::ostringstream msg2;
			msg2<<"#Error="<<error;
			io_.printline(msg2);
		}

		//! checkpoint save
		void checkpointSave(MyBasisWithOperators &pS,MyBasisWithOperators &pE,size_t loop) 
		{
			pS.save(io_,"#CHKPOINTSYSTEM");
			pE.save(io_,"#CHKPOINTENVIRON");
			
		}

		//! checkpoint load
		void checkpointLoad(MyBasisWithOperators &pS,MyBasisWithOperators &pE,size_t loop)
		{
			typename IoType::In ioTmp(parameters_.checkpoint.filename);
			
			pS.load(ioTmp,"#CHKPOINTSYSTEM",loop);
			pE.load(ioTmp,"#CHKPOINTENVIRON");
			
			//load also the stacks here!!!
			std::string s = appendWithDir(ENVIRON_STACK_STRING,parameters_.checkpoint.filename);
			envStack_.load(s);
			s = appendWithDir(SYSTEM_STACK_STRING,parameters_.checkpoint.filename);
			systemStack_.load(s);
		}

		// Save to disk everything needed to compute any observable (OBSOLETE!!)
		void saveToDiskForObserver(
			const FermionSign& fs,
			MyBasis const &pS,
			MyBasis const &pE,
			MyBasis const &pSE,
			TargettingType const  &target) 
		{
			fs.save(io_);
			pS.save(io_);
			pE.save(io_);
			pSE.save(io_);

			// save wavefunction
			std::string label = "#WAVEFUNCTION_sites=";
			for (size_t i=0;i<pSE.block().size();i++) {
				label += utils::ttos(pSE.block()[i])+",";
			}
			//SparseVector<typename TargettingType::TargetVectorType::value_type> psiSparse(target.gs());
			target.gs().save(io_,label);

		}

		// Save to disk transform
		template<typename SomeMatrixType>
		void saveToDisk(const SomeMatrixType& ftransform,
			MyBasisWithOperators const &pS,size_t direction)
		{
			std::string label = "#TRANSFORM_sites=";
			for (size_t i=0;i<pS.block().size();i++) {
				label += utils::ttos(pS.block()[i])+",";
			}
			io_.printMatrix(ftransform,label);
			std::string s = "#DIRECTION="+utils::ttos(direction);
			io_.printline(s);
		}

		//! Move elsewhere
		//! returns s1+s2 if s2 has no '/', 
		//! if s2 = s2a + '/' + s2b return s2a + '/' + s1 + s2b
		std::string appendWithDir(const std::string& s1,const std::string& s2) const
		{
			size_t x = s2.find("/");
			if (x==std::string::npos) return s1 + s2;
			std::string suf = s2.substr(x+1,s2.length());
			std::string dir = s2.substr(0,s2.length()-suf.length());
			//throw std::runtime_error("testing\n");
			return dir + s1 + suf;
					
		}
		
		void printOneSiteElectrons(const MyBasis& b)
		{
			std::string s = "#ONE_SITE_ELECTRONS\n";
			if (b.block().size()!=1) throw std::runtime_error("printOneSiteElectrons failed\n");
			b.save(io_);
		}
		
	}; //class DmrgSolver
} // namespace Dmrg

/*@}*/
#endif
