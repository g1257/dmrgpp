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

/*! \file Diagonalization.h
 *
 *  FIXME needs doc
 */
#ifndef DIAGONALIZATION_HEADER_H
#define DIAGONALIZATION_HEADER_H

#include "Utils.h"
#include "ProgressIndicator.h"
#include "VectorWithOffset.h" // includes the std::norm functions
#include "VectorWithOffsets.h" // includes the std::norm functions

namespace Dmrg {
	
	template<
		typename ParametersType,
  		typename TargettingType,
    		template<typename,typename> class InternalProductTemplate
    		>
	class Diagonalization {
 public:
	 
	 	typedef typename TargettingType::WaveFunctionTransformationType WaveFunctionTransformationType;
		typedef typename TargettingType::ModelType ModelType;
		typedef typename TargettingType::ConcurrencyType ConcurrencyType;
		typedef typename TargettingType::IoType IoType;
		typedef typename TargettingType::BasisType BasisType;
		typedef typename TargettingType::BasisWithOperatorsType BasisWithOperatorsType;
		typedef typename TargettingType::BlockType BlockType;
		typedef typename TargettingType::TargetVectorType TargetVectorType;
		typedef typename TargettingType::RealType RealType;
		typedef typename IoType::Out IoOutType;
		typedef typename ModelType::OperatorsType OperatorsType;
		typedef typename  OperatorsType::SparseMatrixType SparseMatrixType;
		
		Diagonalization(const ParametersType& parameters,
				const ModelType& model,
    				ConcurrencyType& concurrency,
				const bool& verbose,
    				const bool& useReflection,
				IoOutType& io,
    				const int& quantumSector,
    				WaveFunctionTransformationType& waveFunctionTransformation) 
			:
			parameters_(parameters),
			model_(model),
			concurrency_(concurrency),
			verbose_(verbose),
			useReflection_(useReflection),
			io_(io),
			progress_("Diagonalization",0),
			quantumSector_(quantumSector),
			waveFunctionTransformation_(waveFunctionTransformation)
		{}
		
		RealType operator()(TargettingType& target,size_t direction,const BlockType& block,size_t loopIndex=0,
				   bool needsPrinting = false)
		{
			const BasisType& pSE= target.basisSE();
			const BasisWithOperatorsType& pSprime= target.basisS();
			const BasisWithOperatorsType& pEprime= target.basisE();
			
			
			std::ostringstream msg;
			msg<<"Setting up ham. pse.size="<<pSE.size();
			progress_.printline(msg,std::cout);
		
			TargetVectorType tmpVec;
			std::vector<TargetVectorType> vecSaved;
			std::vector<RealType> energySaved;
			size_t i,j,bs,counter;
			RealType gsEnergy;
			size_t total = pSE.partition()-1;

			energySaved.resize(total);
			vecSaved.resize(total);
			std::vector<size_t> weights(total);

			counter=0;
			for (i=0;i<total;i++) {
				bs = pSE.partition(i+1)-pSE.partition(i);
				if (verbose_) {
					j = pSE.qn(pSE.partition(i));
					std::vector<size_t> qns = BasisType::decodeQuantumNumber(j);
					//std::cerr<<"partition "<<i<<" of size="<<bs<<" has qns=";
					for (size_t k=0;k<qns.size();k++) std::cerr<<qns[k]<<" ";
					//std::cerr<<"\n";
				}

				size_t tmp = pSE.electrons(counter);
				if (quantumSector_<0 && model_.density()*pSE.block().size()!=tmp) weights[i]=0;
				else weights[i]=bs;
				
				// Do only one sector unless we're doing search 
				// Note if we're doing search then quantumSector_ will be negative
				if (quantumSector_>=0 &&
						pSE.pseudoEffectiveNumber(pSE.partition(i))!=size_t(quantumSector_) && 
								weights[i]>0) {
					weights[i]=0;
				}
				
				counter+=bs;
				vecSaved[i].resize(weights[i]);
			}

			// legacy part, need to update code to remove support for search sectors
			if (parameters_.options.find("hasQuantumNumbers")==std::string::npos){
				throw std::runtime_error("DmrgSolver:: You must have quantum number in the input file!\n");
			}
			
			typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
			VectorWithOffsetType initialVector(weights,pSE);
			
			waveFunctionTransformation_.triggerOn(pSprime,pEprime,pSE);
			target.initialGuess(initialVector);
			
			concurrency_.loopCreate(total,weights);

			while (concurrency_.loop(i)) {
				if (weights[i]==0) continue;
				std::ostringstream msg;
				msg<<"About to diag. sector with quantum numbs. ";
				j = pSE.qn(pSE.partition(i));
				std::vector<size_t> qns = BasisType::decodeQuantumNumber(j);
				for (size_t k=0;k<qns.size();k++) msg<<qns[k]<<" ";
				msg<<" pseudo="<<pSE.pseudoEffectiveNumber(pSE.partition(i));
				msg<<" quantumSector="<<quantumSector_;
				
				if (verbose_) {
					msg<<" diagonaliseOneBlock, i="<<i<<" and proc="<<concurrency_.rank()<<" and weight="<<weights[i];
				}
				progress_.printline(msg,std::cout);
				TargetVectorType initialVectorBySector(weights[i]);
				initialVector.extract(initialVectorBySector,i);
				diagonaliseOneBlock(i,tmpVec,gsEnergy,pSprime,pEprime,pSE,initialVectorBySector);
				vecSaved[i] = tmpVec;
				energySaved[i]=gsEnergy;
			}
			
			concurrency_.gather(energySaved);
			concurrency_.gather(vecSaved);
			
			concurrency_.broadcast(energySaved);
			concurrency_.broadcast(vecSaved);
				
			// calc gs energy
			if (verbose_ && concurrency_.root()) std::cerr<<"About to calc gs energy\n";
			gsEnergy=1e6;
			for (i=0;i<total;i++) {
				if (weights[i]==0) continue;
				if (energySaved[i]<gsEnergy) gsEnergy=energySaved[i];
			}
			
			if (verbose_ && concurrency_.root()) std::cerr<<"About to calc gs vector\n";
			//target.reset();
			counter=0;
			for (i=0;i<pSE.partition()-1;i++) {
				if (weights[i]==0) continue;
				if (quantumSector_<0) // && fabs(energySaved[i]-gsEnergy)>1e-6) 
					throw std::runtime_error("DmrgSolver: quantumSector_ must be greater than 0\n");
				
				j = pSE.qn(pSE.partition(i));
				std::vector<size_t> qns = BasisType::decodeQuantumNumber(j);
				msg<<"Found target in partition "<<i<<" of size="<<vecSaved[i].size();
				msg<<" with qns=";
				for (size_t k=0;k<qns.size();k++) msg<<qns[k]<<" ";
				progress_.printline(msg,std::cout);
				counter++;
			}
			
			target.setGs(vecSaved,pSE);

			if (concurrency_.root()) {
				std::ostringstream msg;
				msg.precision(8);
				msg<<"#Energy="<<gsEnergy;
				if (counter>1) msg<<" attention: found "<<counter<<" matrix blocks";
				io_.printline(msg);
			}
			
			// time step targetting: 
			target.evolve(gsEnergy,direction,block,loopIndex,needsPrinting);
			waveFunctionTransformation_.triggerOff(pSprime,pEprime,pSE); //,m);
			return gsEnergy;
		}
 private:
		//! Diagonalise the i-th block of the matrix, return its eigenvectors in tmpVec and its eigenvalues in energyTmp
		template<typename SomeVectorType>
		void diagonaliseOneBlock(
					int i,
					SomeVectorType &tmpVec,
     					double &energyTmp,
					BasisWithOperatorsType const &pSprime,
					BasisWithOperatorsType const &pEprime,
     					BasisType const &pSE,
					const SomeVectorType& initialVector)
		{
			RealType eps=ProgramGlobals::LanczosTolerance;
			int iter=ProgramGlobals::LanczosSteps;
			std::vector<RealType> tmpVec1,tmpVec2;
			//srand48(7123443);
			
			typename ModelType::ModelHelperType modelHelper(i,pSE,pSprime,pEprime,model_.orbitals(),useReflection_);

			if (parameters_.options.find("debugmatrix")!=std::string::npos) {
				SparseMatrixType fullm;
				
				model_.fullHamiltonian(fullm,modelHelper);
				
				
				if (!isHermitian(fullm)) throw std::runtime_error("Not hermitian matrix block\n");
				
				psimag::Matrix<typename SparseMatrixType::value_type> fullm2;
				crsMatrixToFullMatrix(fullm2,fullm);
				if (isZero(fullm2)) std::cerr<<"Matrix is zero\n";
				std::cerr<<fullm2;
				std::vector<RealType> eigs(fullm2.n_row());
				utils::diag(fullm2,eigs,'V');
				std::cerr<<"eigs[0]="<<eigs[0]<<"\n";
				if (parameters_.options.find("test")!=std::string::npos)
					throw std::logic_error("Exiting due to option test in the input file\n");
			}
			std::ostringstream msg;
			msg<<"I will now diagonalize a matrix of size="<<modelHelper.size();
			progress_.printline(msg,std::cout);
			diagonaliseOneBlock(i,tmpVec,energyTmp,modelHelper,initialVector,iter,eps);
		}
		
		template<typename SomeVectorType>
		void diagonaliseOneBlock(
					int i,
     					SomeVectorType &tmpVec,
	  				double &energyTmp,
					typename ModelType::ModelHelperType& modelHelper,
     					const SomeVectorType& initialVector,
					size_t iter,
     					RealType eps,
       					int reflectionSector= -1)
		{
			if (reflectionSector>=0) modelHelper.setReflectionSymmetry(reflectionSector);
			int n = model_.getSize(modelHelper);
			if (verbose_) std::cerr<<"Lanczos: About to do block number="<<i<<" of size="<<n<<"\n";
			typedef InternalProductTemplate<typename SomeVectorType::value_type,ModelType> MyInternalProduct;
			typedef LanczosSolver<MyInternalProduct,SomeVectorType> LanczosSolverType;
			typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);
			size_t mode = LanczosSolverType::WITH_INFO;
			if (parameters_.options.find("lanczosdebug")!=std::string::npos) mode =  LanczosSolverType::DEBUG;
					
			LanczosSolverType lanczosSolver(lanczosHelper,iter,eps,concurrency_.rank(),mode);
			
			tmpVec.resize(lanczosHelper.rank());
			if (lanczosHelper.rank()==0) {
				energyTmp=10000;
				return;
			}
			/*std::ostringstream msg;
			msg<<"Calling computeGroundState...\n";
			progress_.printline(msg,std::cerr);
			*/	
			lanczosSolver.computeGroundState(energyTmp,tmpVec,initialVector);
		}

	 	const ParametersType& parameters_;
		const ModelType& model_;
		ConcurrencyType& concurrency_;
		//const TargettingStructureType& targetStruct_;
		//DiagonalizationType diagonalization_;
		const bool& verbose_;
		const bool& useReflection_;
		IoOutType& io_;
		ProgressIndicator progress_;
		const int& quantumSector_;
		WaveFunctionTransformationType& waveFunctionTransformation_;
	}; // class Diagonalization
} // namespace Dmrg 

/*@}*/
#endif
