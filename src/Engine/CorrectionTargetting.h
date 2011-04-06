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

/*! \file CorrectionTargetting.h
 *
 *  corrects the finite-size algorithm
 *  following PRB 72, 180403(R) (2005)
 *
 */
 
#ifndef CORRECTION_TARGETTING_H
#define CORRECTION_TARGETTING_H
#include <iostream>
#include <string>
#include "CorrectionParams.h"
#include "ApplyOperatorLocal.h"
#include <stdexcept>

namespace Dmrg {
	
	template<
			template<typename,typename,typename> class LanczosSolverTemplate,
  			template<typename,typename> class InternalProductTemplate,
     			template<typename,typename> class WaveFunctionTransfTemplate,
     			typename ModelType_,
			typename ConcurrencyType_,
   			typename IoType_,
      			template<typename> class VectorWithOffsetTemplate>
	class CorrectionTargetting  {
		public:
			typedef ModelType_ ModelType;
			typedef ConcurrencyType_ ConcurrencyType;
			typedef IoType_ IoType;
			
			typedef PsimagLite::IoSimple::In IoInputType;

			typedef typename ModelType::RealType RealType;
			typedef InternalProductTemplate<RealType,ModelType> InternalProductType;
			typedef std::vector<RealType> VectorType;
			typedef LanczosSolverTemplate<RealType,InternalProductType,VectorType> LanczosSolverType;
			typedef typename ModelType::ModelHelperType ModelHelperType;
			typedef typename ModelHelperType::LeftRightSuperType
				LeftRightSuperType;
			typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
			typedef typename BasisWithOperatorsType::BasisDataType BasisDataType;
			typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
			//typedef psimag::Matrix<RealType> MatrixType;
			//typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
			typedef VectorType TargetVectorType;
			typedef CorrectionParams<ModelType> TargettingParamsType;
			typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
			typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
			
			enum {DISABLED,ENABLED};
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};
			
			CorrectionTargetting(
	  				const LeftRightSuperType& lrs,
 					const ModelType& model,
					const TargettingParamsType& correctionStruct,
     				const WaveFunctionTransfType& wft) // wft is ignored here
				: lrs_(lrs),
				  model_(model),
				  correctionStruct_(correctionStruct),
				  waveFunctionTransformation_(wft),
				  progress_("CorrectionTargetting",0),
				  stage_(DISABLED),
				  targetVectors_(1),
				  applyOpLocal_(lrs)
			{
			}
			
			RealType normSquared(size_t i) const
			{
				return std::real(multiply(targetVectors_[i],targetVectors_[i]));
			}
			
			RealType weight(size_t i) const
			{
				if (stage_ == DISABLED) throw std::runtime_error(
					"CorrectionTargetting: What are you doing here?\n");
				return correctionStruct_.a;;
			}
			
			RealType gsWeight() const
			{
				return 1;
			}
			
			template<typename SomeBasisType>
			void setGs(const std::vector<VectorType>& v,//const std::vector<size_t>& weights,
				   const SomeBasisType& someBasis)
			{
				psi_.set(v,someBasis);
			}
			
			const VectorWithOffsetType& gs() const 
			{
				return psi_;
			}
					
			bool includeGroundStage() const {return true; }
			
			size_t size() const 
			{
				if (stage_==DISABLED) return 0;
				return targetVectors_.size();
			}
			
			const VectorWithOffsetType& operator()(size_t i) const
			{
				return targetVectors_[i];
			}
			
			void evolve(RealType Eg,size_t direction,const BlockType& block,size_t loopNumber)
			{
				if (direction==INFINITE) return;
				stage_ = ENABLED;
				// block.size()==1 or throw
//				size_t i = block[0]; // center site
//				FermionSign fs(lrs_.left(),correctionStruct_.electrons);
//				const BasisWithOperatorsType *basis = (direction==EXPAND_SYSTEM) ?
//						&(lrs_.left()) : &(lrs_.right());

				// operators in the one-site basis:
				std::vector<OperatorType> creationMatrix;
				SparseMatrixType hmatrix;
				BasisDataType q;

				model_.setNaturalBasis(creationMatrix,hmatrix,q,block);
				std::vector<size_t> electronsOneSite(q.electronsUp.size());
				for (size_t i=0;i<electronsOneSite.size();i++)
					electronsOneSite[i] = q.electronsUp[i] + q.electronsDown[i];

				FermionSign fs(lrs_.left(),electronsOneSite);
				for (size_t j=0;j<creationMatrix.size();j++) {
					VectorWithOffsetType phiTemp;
					applyOpLocal_(phiTemp,psi_,creationMatrix[j],
							fs,direction);
					if (j==0) targetVectors_[0] = phiTemp;
					else targetVectors_[0] += phiTemp;
				}
			}
			
			const LeftRightSuperType& leftRightSuper() const
			{
				return lrs_;
			}

			void initialGuess(VectorWithOffsetType& initialVector) const
			{
				RealType eps = 1e-6;
				if (psi_.size()>0 && std::norm(psi_)<eps) throw std::runtime_error("psi's norm is zero\n");
				waveFunctionTransformation_.setInitialVector(initialVector,psi_,lrs_);
			}
			
			template<typename IoOutputType>
			void save(const std::vector<size_t>& block,IoOutputType& io) const
			{
				std::ostringstream msg;
				msg<<"Saving state...";
				progress_.printline(msg,std::cout);

				if (block.size()!=1) throw std::runtime_error(
						"GST only supports blocks of size 1\n");
				std::string s = "#TCENTRALSITE=" + ttos(block[0]);
				io.printline(s);

				psi_.save(io,"PSI");
			}

			void load(const std::string& f)
			{
				IoInputType io(f);
				int site=0;
				io.readline(site,"#TCENTRALSITE=",IoType::In::LAST_INSTANCE);
				if (site<0) throw std::runtime_error(
						"GST::load(...): site cannot be negative\n");
				psi_.load(io,"PSI");
			}

			void print(std::ostream& os) const
			{
				os<<"CorrectionWeightGroundState=1\n";
				os<<"CorrectionWeightCorrection="<<correctionStruct_.a<<"\n";
			}

		private:
			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const TargettingParamsType& correctionStruct_;
			const WaveFunctionTransfType& waveFunctionTransformation_;
			PsimagLite::ProgressIndicator progress_;
			size_t stage_;
			VectorWithOffsetType psi_;
			std::vector<VectorWithOffsetType> targetVectors_;
			ApplyOperatorType applyOpLocal_;
	};     //class CorrectionTargetting

	template<
	template<typename,typename,typename> class LanczosSolverTemplate,
	template<typename,typename> class InternalProductTemplate,
	template<typename,typename> class WaveFunctionTransfTemplate,
	typename ModelType_,
	typename ConcurrencyType_,
	typename IoType_,
	template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const CorrectionTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,ConcurrencyType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		tst.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif
