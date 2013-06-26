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

/*! \file CorrectionTargetting.h
 *
 *  corrects the finite-size algorithm
 *  following PRB 72, 180403(R) (2005)
 *
 */
 
#ifndef CORRECTION_TARGETTING_H
#define CORRECTION_TARGETTING_H
#include <iostream>
#include "String.h"
#include "CorrectionParams.h"
#include "CommonTargetting.h"
#include <stdexcept>

namespace Dmrg {
	
	template<
			template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
     			template<typename,typename> class WaveFunctionTransfTemplate,
     			typename ModelType_,
   			typename IoType_,
      			template<typename> class VectorWithOffsetTemplate>
	class CorrectionTargetting  {
		public:
			typedef ModelType_ ModelType;
			typedef IoType_ IoType;
			
			typedef PsimagLite::IoSimple::In IoInputType;

			typedef typename ModelType::RealType RealType;
			typedef InternalProductTemplate<RealType,ModelType> InternalProductType;
			typedef typename PsimagLite::Vector<RealType>::Type VectorType;
			typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
			typedef LanczosSolverTemplate<ParametersForSolverType,InternalProductType,VectorType> LanczosSolverType;
			typedef typename ModelType::ModelHelperType ModelHelperType;
			typedef typename ModelHelperType::LeftRightSuperType
				LeftRightSuperType;
			typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
			typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
			//typedef psimag::Matrix<RealType> MatrixType;
			//typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
			typedef VectorType TargetVectorType;
			typedef CorrectionParams<ModelType> TargettingParamsType;
			typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
			typedef CommonTargetting<ModelType,TargettingParamsType,WaveFunctionTransfType,VectorWithOffsetType,LanczosSolverType>
			CommonTargettingType;
			typedef typename CommonTargettingType::ApplyOperatorType ApplyOperatorType;
			
			enum {DISABLED,ENABLED};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};

			CorrectionTargetting(const LeftRightSuperType& lrs,
 			                     const ModelType& model,
			                     const TargettingParamsType& correctionStruct,
     		                     const WaveFunctionTransfType& wft, // wft is ignored here
			                     const SizeType& quantumSector) // quantumSector ignored here
			: lrs_(lrs),
			  model_(model),
			  correctionStruct_(correctionStruct),
			  waveFunctionTransformation_(wft),
			  progress_("CorrectionTargetting"),
			  stage_(DISABLED),
			  targetVectors_(1),
			  applyOpLocal_(lrs),
			  commonTargetting_(lrs,model,correctionStruct)
			{
			}
			
			const ModelType& model() const { return model_; }

			RealType normSquared(SizeType i) const
			{
				return std::real(multiply(targetVectors_[i],targetVectors_[i]));
			}
			
			RealType weight(SizeType i) const
			{
				if (stage_ == DISABLED) throw PsimagLite::RuntimeError(
					"CorrectionTargetting: What are you doing here?\n");
				return correctionStruct_.correctionA;
			}
			
			RealType gsWeight() const
			{
				return 1;
			}
			
			template<typename SomeBasisType>
			void setGs(const typename PsimagLite::Vector<VectorType>::Type& v,//const typename PsimagLite::Vector<SizeType>::Type& weights,
				   const SomeBasisType& someBasis)
			{
				psi_.set(v,someBasis);
			}
			
			const VectorWithOffsetType& gs() const 
			{
				return psi_;
			}
					
			bool includeGroundStage() const {return true; }
			
			SizeType size() const 
			{
				if (stage_==DISABLED) return 0;
				return targetVectors_.size();
			}
			
			const VectorWithOffsetType& operator()(SizeType i) const
			{
				return targetVectors_[i];
			}

			RealType time() const { return 0; }
			
			void evolve(RealType Eg,
			            SizeType direction,
			            const BlockType& block1,
			            const BlockType& block2,
			            SizeType loopNumber)
			{
				cocoon(block1,direction);

				if (direction==INFINITE) return;
				stage_ = ENABLED;

				commonTargetting_.computeCorrection(targetVectors_[0],direction,block1,psi_);
			}
			
			const LeftRightSuperType& leftRightSuper() const
			{
				return lrs_;
			}

			void initialGuess(VectorWithOffsetType& initialVector,
			                  const typename PsimagLite::Vector<SizeType>::Type& block) const
			{
				RealType eps = 1e-6;
				if (psi_.size()>0 && std::norm(psi_)<eps)
					throw PsimagLite::RuntimeError("psi's norm is zero\n");
				typename PsimagLite::Vector<SizeType>::Type nk;
				commonTargetting_.setNk(nk,block);
				waveFunctionTransformation_.setInitialVector(initialVector,psi_,lrs_,nk);
			}
			
			template<typename IoOutputType>
			void save(const typename PsimagLite::Vector<SizeType>::Type& block,IoOutputType& io) const
			{
				PsimagLite::OstringStream msg;
				msg<<"Saving state...";
				progress_.printline(msg,std::cout);

				if (block.size()!=1) throw PsimagLite::RuntimeError(
						"GST only supports blocks of size 1\n");
				//io.print("#TCENTRALSITE=",block[0]);
				PsimagLite::String s = "#TCENTRALSITE=" + ttos(block[0]);
				io.printline(s);
				psi_.save(io,"PSI");
			}

			void load(const PsimagLite::String& f)
			{
				IoInputType io(f);
				int site=0;
				io.readline(site,"#TCENTRALSITE=",IoType::In::LAST_INSTANCE);
				if (site<0) throw PsimagLite::RuntimeError(
						"GST::load(...): site cannot be negative\n");
				psi_.load(io,"PSI");
			}

			void print(std::ostream& os) const
			{
				os<<"CorrectionWeightGroundState=1\n";
				os<<"CorrectionWeightCorrection="<<correctionStruct_.correctionA<<"\n";
			}

			void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
			{
				// nothing to do here
			}

		private:

			void cocoon(const BlockType& block1,SizeType direction) const
			{
				if (model_.params().insitu=="") return;

				if (BasisType::useSu2Symmetry()) {
					commonTargetting_.noCocoon("not when SU(2) symmetry is in use");
					return;
				}

				try {
					assert(block1.size()>0);
					commonTargetting_.cocoon(direction,block1[0],psi_);
				} catch (std::exception& e) {
					commonTargetting_.noCocoon("unsupported by the model");
				}
			}

			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& correctionStruct_;
			const WaveFunctionTransfType& waveFunctionTransformation_;
			PsimagLite::ProgressIndicator progress_;
			SizeType stage_;
			VectorWithOffsetType psi_;
			typename PsimagLite::Vector<VectorWithOffsetType>::Type targetVectors_;
			ApplyOperatorType applyOpLocal_;
			CommonTargettingType commonTargetting_;
	};     //class CorrectionTargetting

	template<
	template<typename,typename,typename> class LanczosSolverTemplate,
	template<typename,typename> class InternalProductTemplate,
	template<typename,typename> class WaveFunctionTransfTemplate,
	typename ModelType_,
	typename IoType_,
	template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const CorrectionTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		tst.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif
