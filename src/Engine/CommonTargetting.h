
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

/*! \file CommonTargetting.h
 *
 * Functionality used by many targetting classes
 *
 */

#ifndef COMMON_TARGETTING_H
#define COMMON_TARGETTING_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "DynamicSerializer.h"
#include "DynamicDmrgParams.h"
#include "VectorWithOffsets.h"
#include "ContinuedFraction.h"
#include <cassert>

namespace Dmrg {

	template<typename ModelType,
	         typename TargettingParamsType,
	         typename WaveFunctionTransfType,
	         typename VectorWithOffsetType,
			 typename LanczosSolverType>
	class CommonTargetting  {

	public:

		typedef typename ModelType::RealType RealType;
		typedef typename ModelType::ModelHelperType ModelHelperType;
		typedef typename ModelHelperType::LeftRightSuperType
		                                  LeftRightSuperType;
		typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
		typedef typename LanczosSolverType::PostProcType PostProcType;
		typedef typename LanczosSolverType::DenseMatrixType DenseMatrixType;
		typedef typename VectorWithOffsetType::VectorType VectorType;
		typedef DynamicSerializer<RealType,VectorWithOffsetType,PostProcType>
		        DynamicSerializerType;

		enum {DISABLED,OPERATOR,CONVERGING};
		enum {
			EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE
		};

		static const size_t parallelRank_ = 0; // DYNT needs to support concurrency FIXME

		CommonTargetting(const LeftRightSuperType& lrs,
		                 const ModelType& model,
		                 const TargettingParamsType& tstStruct,
		                 std::vector<VectorWithOffsetType>& targetVectors) 
		: lrs_(lrs),
		  model_(model),
		  tstStruct_(tstStruct),
		  targetVectors_(targetVectors),
		  progress_("CommonTargetting",0)
//		  applyOpLocal_(lrs),
//		  gsWeight_(1.0)
		{
		}

		RealType normSquared(size_t i) const
		{
			// call to mult will conjugate one of the vector
			return std::real(multiply(targetVectors_[i],targetVectors_[i]));
		}
		
		template<typename IoOutputType>
		void save(const std::vector<size_t>& block,
		          IoOutputType& io,
		          const PostProcType& cf) const
		{
			DynamicSerializerType dynS(cf,block[0],targetVectors_);
			dynS.save(io);
		}

		template<typename IoInputType>
		void load(IoInputType& io)
		{
			DynamicSerializerType dynS(io,IoInputType::In::LAST_INSTANCE);
			for (size_t i=0;i<targetVectors_.size();i++)
				targetVectors_[i] = dynS.vector(i);
			
		}

		void checkOrder(size_t i,const std::vector<size_t>& stage) const
		{
			if (i==0) return;
			for (size_t j=0;j<i;j++) {
				if (stage[j] == DISABLED) {
					std::string s ="TST:: Seeing dynamic site "+ttos(tstStruct_.sites[i]);
					s =s + " before having seen";
					s = s + " site "+ttos(j);
					s = s +". Please order your dynamic sites in order of appearance.\n";
					throw std::runtime_error(s);
				}
			}
		}

		bool allStages(size_t x,const std::vector<size_t>& stage) const
		{
			for (size_t i=0;i<stage.size();i++)
				if (stage[i]!=x) return false;
			return true;
		}

		bool noStageIs(size_t x,const std::vector<size_t>& stage) const
		{
			for (size_t i=0;i<stage.size();i++)
				if (stage[i]==x) return false;
			return true;
		}

		std::string getStage(size_t i,const std::vector<size_t>& stage) const
		{
			switch (stage[i]) {
			case DISABLED:
				return "Disabled";
				break;
			case OPERATOR:
				return "Applying operator for the first time";
				break;
			case CONVERGING:
				return "Converging DDMRG";
				break;
			}
			return "undefined";
		}

		void initialGuess(VectorWithOffsetType& v,
		                  const WaveFunctionTransfType& wft,
		                  const VectorWithOffsetType& psi,
		                  const std::vector<size_t>& stage,
		                  const std::vector<RealType>& weights,
		                  size_t nk) const
		{
			wft.setInitialVector(v,psi,lrs_,nk);
			if (!allStages(CONVERGING,stage)) return;
			std::vector<VectorWithOffsetType> vv(targetVectors_.size());
			for (size_t i=0;i<targetVectors_.size();i++) {
				wft.setInitialVector(vv[i],targetVectors_[i],lrs_,nk);
				if (std::norm(vv[i])<1e-6) continue;
				VectorWithOffsetType w= weights[i]*vv[i];
				v += w;
			}
		}
// 		void guessPhiSectors(VectorWithOffsetType& phi,
// 		                     size_t i,
// 		                     size_t systemOrEnviron)
// 		{
// 			FermionSign fs(lrs_.left(),tstStruct_.electrons);
// 			if (allStages(CONVERGING)) {
// 				VectorWithOffsetType tmpVector = psi_;
// 				for (size_t j=0;j<tstStruct_.aOperators.size();j++) {
// 					applyOpLocal_(phi,tmpVector,tstStruct_.aOperators[j],fs,
// 					              systemOrEnviron);
// 					tmpVector = phi;
// 				}
// 				return;
// 			}
// 			applyOpLocal_(phi,psi_,tstStruct_.aOperators[i],fs,
// 			              systemOrEnviron);
// 		}

	private:

// 		void zeroOutVectors()
// 		{
// 			for (size_t i=0;i<targetVectors_.size();i++)
// 				targetVectors_[i].resize(lrs_.super().size());
// 		}

// 		std::vector<size_t> stage_;
// 		VectorWithOffsetType psi_;
		const LeftRightSuperType& lrs_;
		const ModelType& model_;
		const TargettingParamsType& tstStruct_;
		std::vector<VectorWithOffsetType>& targetVectors_;
// 		const WaveFunctionTransfType& waveFunctionTransformation_;
		PsimagLite::ProgressIndicator progress_;
// 		ApplyOperatorType applyOpLocal_;
// 		RealType gsWeight_;
// 		std::vector<VectorWithOffsetType> targetVectors_;
// 		std::vector<RealType> weight_;
	}; // class CommonTargetting

	template<typename ModelType,
	         typename TargettingParamsType,
	         typename WaveFunctionTransfType,
	         typename VectorWithOffsetType,
	         typename LanczosSolverType
	>
	std::ostream& operator<<(std::ostream& os,
	                         const CommonTargetting<ModelType,TargettingParamsType,WaveFunctionTransfType,VectorWithOffsetType,LanczosSolverType>& tst)
	{
		os<<"DT=NothingToSeeHereYet\n";
		return os;
	}

} // namespace
/*@}*/
#endif // COMMON_TARGETTING_H

