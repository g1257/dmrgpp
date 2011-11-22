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
#include "ParametersForSolver.h"

namespace Dmrg {
	
	template<typename ParametersType,
	         typename TargettingType,
	         template<typename,typename> class InternalProductTemplate
    >
	class Diagonalization {

	public:
	 
	 	typedef typename TargettingType::WaveFunctionTransfType WaveFunctionTransfType;
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
		typedef typename ModelType::ModelHelperType ModelHelperType;
		typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;

		Diagonalization(const ParametersType& parameters,
                        const ModelType& model,
                        ConcurrencyType& concurrency,
                        const bool& verbose,
                        const bool& useReflection,
                        IoOutType& io,
                        const size_t& quantumSector,
                       WaveFunctionTransfType& waveFunctionTransformation)
		: parameters_(parameters),
		  model_(model),
		  concurrency_(concurrency),
		  verbose_(verbose),
		  useReflection_(useReflection),
		  io_(io),
		  progress_("Diag.",0),
		  quantumSector_(quantumSector),
		  wft_(waveFunctionTransformation),
		  oldEnergy_(0)
		{}

		RealType operator()(TargettingType& target,
		                    size_t direction,
		                    const BlockType& blockLeft,
		                    const BlockType& blockRight)
		{
			if (direction!=WaveFunctionTransfType::INFINITE) throw std::runtime_error(
				"Diagonalization::operator(): expecting INFINITE direction\n");
			size_t loopIndex = 0;
			size_t nk = 0; // bogus
			RealType gsEnergy = internalMain_(target,direction,loopIndex,false,nk);
			//  targetting: 
			target.evolve(gsEnergy,direction,blockLeft,blockRight,loopIndex);
			wft_.triggerOff(target.leftRightSuper()); //,m);
			return gsEnergy;
		}

		RealType operator()(TargettingType& target,
		                    size_t direction,
		                    const BlockType& block,
		                    size_t loopIndex,
		                    bool needsPrinting)
		{
			assert(direction!=WaveFunctionTransfType::INFINITE);
			assert(block.size()==1);

			size_t nk = model_.hilbertSize(block[0]);
			RealType gsEnergy = internalMain_(target,direction,loopIndex,false,nk);
			//  targetting: 
			target.evolve(gsEnergy,direction,block,block,loopIndex);
			wft_.triggerOff(target.leftRightSuper()); //,m);
			return gsEnergy;
		}

	private:

		RealType internalMain_(TargettingType& target,
		                       size_t direction,
		                       size_t loopIndex,
		                       bool needsPrinting,
		                       size_t nk)

		{
			const LeftRightSuperType& lrs= target.leftRightSuper();
			wft_.triggerOn(lrs);

			RealType gsEnergy = 0;
			if (!target.includeGroundStage()) return gsEnergy;

			bool onlyWft = ((parameters_.finiteLoop[loopIndex].
					saveOption & 2)>0) ? true : false;
			
			std::ostringstream msg;
			msg<<"Setting up Hamiltonian basis of size="<<lrs.super().size();
			progress_.printline(msg,std::cout);
		
			TargetVectorType tmpVec;
			std::vector<TargetVectorType> vecSaved;
			std::vector<RealType> energySaved;
			
			size_t total = lrs.super().partition()-1;

			energySaved.resize(total);
			vecSaved.resize(total);
			std::vector<size_t> weights(total);

			size_t counter=0;
			for (size_t i=0;i<total;i++) {
				size_t bs = lrs.super().partition(i+1)-lrs.super().partition(i);
				if (verbose_) {
					size_t j = lrs.super().qn(lrs.super().partition(i));
					std::vector<size_t> qns = BasisType::decodeQuantumNumber(j);
					//std::cerr<<"partition "<<i<<" of size="<<bs<<" has qns=";
					for (size_t k=0;k<qns.size();k++) std::cerr<<qns[k]<<" ";
					//std::cerr<<"\n";
				}

				weights[i]=bs;

				// Do only one sector unless doing su(2) with j>0, then do all m's
				if (lrs.super().pseudoEffectiveNumber(
						lrs.super().partition(i))!=quantumSector_ )
					weights[i]=0;
				
				counter+=bs;
				vecSaved[i].resize(weights[i]);
			}

			typedef typename TargettingType::VectorWithOffsetType
					VectorWithOffsetType;
			VectorWithOffsetType initialVector(weights,lrs.super());
			
			target.initialGuess(initialVector,nk);
			

			for (size_t i=0;i<total;i++) {
				if (weights[i]==0) continue;
				std::ostringstream msg;
				msg<<"About to diag. sector with quantum numbs. ";
				size_t j = lrs.super().qn(lrs.super().partition(i));
				std::vector<size_t> qns = BasisType::decodeQuantumNumber(j);
				for (size_t k=0;k<qns.size();k++) msg<<qns[k]<<" ";
				msg<<" pseudo="<<lrs.super().pseudoEffectiveNumber(
						lrs.super().partition(i));
				msg<<" quantumSector="<<quantumSector_;
				
				if (verbose_ && concurrency_.root()) {
					msg<<" diagonaliseOneBlock, i="<<i;
					msg<<" and weight="<<weights[i];
				}
				progress_.printline(msg,std::cout);
				TargetVectorType initialVectorBySector(weights[i]);
				initialVector.extract(initialVectorBySector,i);
				if (onlyWft) {
					vecSaved[i]=initialVectorBySector;
					gsEnergy = oldEnergy_;
				} else {
					diagonaliseOneBlock(i,tmpVec,gsEnergy,lrs,initialVectorBySector);
					vecSaved[i] = tmpVec;
				}
				energySaved[i]=gsEnergy;
			}

			// calc gs energy
			if (verbose_ && concurrency_.root()) std::cerr<<"About to calc gs energy\n";
			gsEnergy=1e6;
			for (size_t i=0;i<total;i++) {
				if (weights[i]==0) continue;
				if (energySaved[i]<gsEnergy) gsEnergy=energySaved[i];
			}

			if (verbose_ && concurrency_.root()) std::cerr<<"About to calc gs vector\n";
			//target.reset();
			counter=0;
			for (size_t i=0;i<lrs.super().partition()-1;i++) {
				if (weights[i]==0) continue;

				size_t j = lrs.super().qn(lrs.super().partition(i));
				std::vector<size_t> qns = BasisType::decodeQuantumNumber(j);
				msg<<"Found targetted symmetry sector in partition "<<i;
				msg<<" of size="<<vecSaved[i].size();
				msg<<" with qns=";
				for (size_t k=0;k<qns.size();k++) msg<<qns[k]<<" ";
				progress_.printline(msg,std::cout);
				counter++;
			}

			target.setGs(vecSaved,lrs.super());

			if (concurrency_.root()) {
				std::ostringstream msg;
				msg.precision(8);
				msg<<"#Energy="<<gsEnergy;
				if (counter>1) msg<<" attention: found "<<counter<<" matrix blocks";
				io_.printline(msg);
				oldEnergy_=gsEnergy;
			}
			return gsEnergy;
		}

		//! Diagonalise the i-th block of the matrix, return its eigenvectors in tmpVec and its eigenvalues in energyTmp
		template<typename SomeVectorType>
		void diagonaliseOneBlock(
					int i,
					SomeVectorType &tmpVec,
     				double &energyTmp,
					const LeftRightSuperType& lrs,
					const SomeVectorType& initialVector)
		{
			RealType eps=ProgramGlobals::LanczosTolerance;
			int iter=ProgramGlobals::LanczosSteps;
			std::vector<RealType> tmpVec1,tmpVec2;
			//srand48(7123443);

			typename ModelType::ModelHelperType modelHelper(i,lrs,useReflection_);

			if (parameters_.options.find("debugmatrix")!=std::string::npos) {
				SparseMatrixType fullm;

				model_.fullHamiltonian(fullm,modelHelper);

				if (!isHermitian(fullm))
					throw std::runtime_error("Not hermitian matrix block\n");

				PsimagLite::Matrix<typename SparseMatrixType::value_type> fullm2;
				crsMatrixToFullMatrix(fullm2,fullm);
				if (PsimagLite::isZero(fullm2)) std::cerr<<"Matrix is zero\n";
				printNonZero(fullm2,std::cerr);
				std::vector<RealType> eigs(fullm2.n_row());
				PsimagLite::diag(fullm2,eigs,'V');
				std::cerr<<"eigs[0]="<<eigs[0]<<"\n";
				if (parameters_.options.find("test")!=std::string::npos)
					throw std::logic_error
					         ("Exiting due to option test in the input file\n");
			}
			std::ostringstream msg;
			msg<<"I will now diagonalize a matrix of size="<<modelHelper.size();
			progress_.printline(msg,std::cout);
			diagonaliseOneBlock(i,tmpVec,energyTmp,modelHelper,
					initialVector,iter,eps);
		}
		
		template<typename SomeVectorType>
		void diagonaliseOneBlock(int i,
     		SomeVectorType &tmpVec,
	  		double &energyTmp,
			typename ModelType::ModelHelperType& modelHelper,
     		const SomeVectorType& initialVector,
			size_t iter,
     		RealType eps,
       		int reflectionSector= -1)
		{
			if (reflectionSector>=0) modelHelper.setReflectionSymmetry(reflectionSector);
			int n = modelHelper.size();
			if (verbose_) std::cerr<<"Lanczos: About to do block number="<<i<<" of size="<<n<<"\n";

			typedef InternalProductTemplate<typename SomeVectorType::value_type,ModelType> MyInternalProduct;
			typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
			typedef PsimagLite::LanczosSolver<ParametersForSolverType,MyInternalProduct,SomeVectorType> LanczosSolverType;
			typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);
			
			ParametersForSolverType params;
			params.steps = iter;
			params.tolerance = eps;
			params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;
			params.options= parameters_.options;

			LanczosSolverType lanczosSolver(lanczosHelper,params);

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
		const bool& verbose_;
		const bool& useReflection_;
		IoOutType& io_;
		PsimagLite::ProgressIndicator progress_;
		const size_t& quantumSector_; // this needs to be a reference since DmrgSolver will change it
		WaveFunctionTransfType& wft_;
		double oldEnergy_;
	}; // class Diagonalization
} // namespace Dmrg 

/*@}*/
#endif
