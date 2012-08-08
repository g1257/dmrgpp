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

/*! \file GroundStateTargetting.h
 *
 *  targets the ground state
 *
 */

#ifndef GS_TARGETTING_H
#define GS_TARGETTING_H
#include <iostream>
#include <string>
#include "GroundStateParams.h"
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
	class GroundStateTargetting {
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
			typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
			//typedef psimag::Matrix<RealType> MatrixType;
			//typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef typename BasisType::BlockType BlockType;
			typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
			typedef VectorType TargetVectorType;
			typedef GroundStateParams<ModelType> TargettingParamsType;
			typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType,TargetVectorType> ApplyOperatorType;
			typedef WaveFunctionTransfTemplate<LeftRightSuperType,VectorWithOffsetType> WaveFunctionTransfType;
			
			enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
			EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
			INFINITE=WaveFunctionTransfType::INFINITE};
			
			GroundStateTargetting(const LeftRightSuperType& lrs,
			                      const ModelType& model,
			                      const TargettingParamsType& t,
			                      const WaveFunctionTransfType& wft,  // wft is ignored here
			                      const size_t& quantumSector) // quantumSector is ignored here
			: lrs_(lrs),
			  model_(model),
			  waveFunctionTransformation_(wft),
			  progress_("GroundStateTargetting",0),
			  applyOpLocal_(lrs)
			{
			}

			RealType normSquared(size_t i) const
			{
				throw std::runtime_error("GST: What are you doing here?\n");
			}

			RealType weight(size_t i) const
			{
				throw std::runtime_error("GST: What are you doing here?\n");
				return 0;
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
				return 0;
			}

			const VectorWithOffsetType& operator()(size_t i) const
			{
				throw std::runtime_error("GroundStateTargetting::operator()(...)\n");
			}

			void evolve(RealType Eg,
			            size_t direction,
			            const BlockType& block1,
			            const BlockType& block2,
			            size_t loopNumber)
			{
				if (!BasisType::useSu2Symmetry())
					cocoon(direction,block1);
				else
					noCocoon();
			}

			const LeftRightSuperType& leftRightSuper() const
			{
				return lrs_;
			}

			void initialGuess(VectorWithOffsetType& initialVector,size_t nk) const
			{
				//RealType eps = 1e-6;
				//if (psi_.size()>0 && std::norm(psi_)<eps) throw std::runtime_error("psi's norm is zero\n");
				waveFunctionTransformation_.setInitialVector(initialVector,psi_,lrs_,nk);
			}

			template<typename IoOutputType>
			void save(const std::vector<size_t>& block,IoOutputType& io) const
			{
				std::ostringstream msg;
				msg<<"Saving state...";
				progress_.printline(msg,std::cout);

				if (block.size()!=1) throw std::runtime_error(
						"GST only supports blocks of size 1\n");
				io.print("#TCENTRALSITE=",block[0]);

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

			RealType time() const { return 0; }

			void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
			{
				// nothing to do here
			}

		private:

			void noCocoon() const
			{
				std::cout<<"-------------&*&*&* In-situ measurements start\n";
				std::cout<<"----- NO COCOON WHEN SU(2) SYMMETRY IS IN USE (sorry)\n";
				std::cout<<"-------------&*&*&* In-situ measurements end\n";
			}

			// in situ computation:
			void cocoon(size_t direction,const BlockType& block) const
			{
				size_t site = block[0];
				PsimagLite::CrsMatrix<RealType> tmpC(model_.naturalOperator("nup",0,0));
				int fermionSign1 = 1;
				const std::pair<size_t,size_t> jm1(0,0);
				RealType angularFactor1 = 1.0;
				typename OperatorType::Su2RelatedType su2Related1;
				OperatorType nup(tmpC,fermionSign1,jm1,angularFactor1,su2Related1);

				nup.data = tmpC;
				nup.fermionSign = 1;

				std::cout<<"-------------&*&*&* In-situ measurements start\n";
				test(psi_,psi_,direction,"<PSI|nup|PSI>",site,nup);
//				std::string s = "<P0|nup|P0>";
//				test(targetVectors_[0],targetVectors_[0],direction,s,site,nup);

				PsimagLite::CrsMatrix<RealType> tmpC2(model_.naturalOperator("ndown",0,0));
				OperatorType ndown(tmpC2,fermionSign1,jm1,angularFactor1,su2Related1);
				test(psi_,psi_,direction,"<PSI|ndown|PSI>",site,ndown);
//				s = "<P0|ndown|P0>";
//				test(psi_,targetVectors_[0],direction,s,site,ndown);

				PsimagLite::CrsMatrix<RealType> tmpC3 = tmpC * tmpC2;
				OperatorType doubleOcc(tmpC3,fermionSign1,jm1,angularFactor1,su2Related1);
				test(psi_,psi_,direction,"<PSI|doubleOcc|PSI>",site,doubleOcc);

				std::cout<<"-------------&*&*&* In-situ measurements end\n";
			}

			void test(const VectorWithOffsetType& src1,
					  const VectorWithOffsetType& src2,
					  size_t systemOrEnviron,
					  const std::string& label,
					  size_t site,
					  const OperatorType& A) const
			{
				std::vector<size_t> electrons;
				model_.findElectronsOfOneSite(electrons,site);
				FermionSign fs(lrs_.left(),electrons);
				VectorWithOffsetType dest;
				applyOpLocal_(dest,src1,A,fs,systemOrEnviron);

				RealType sum = 0;
				for (size_t ii=0;ii<dest.sectors();ii++) {
					size_t i = dest.sector(ii);
					size_t offset1 = dest.offset(i);
					for (size_t jj=0;jj<src2.sectors();jj++) {
						size_t j = src2.sector(jj);
						size_t offset2 = src2.offset(j);
						if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
						for (size_t k=0;k<dest.effectiveSize(i);k++)
							sum+= dest[k+offset1] * std::conj(src2[k+offset2]);
					}
				}
				std::cout<<site<<" "<<sum<<" "<<" 0";
				std::cout<<" "<<label<<" "<<(src1*src2)<<"\n";
			}

			const LeftRightSuperType& lrs_;
			const ModelType& model_;
			const WaveFunctionTransfType& waveFunctionTransformation_;
			VectorWithOffsetType psi_;
			PsimagLite::ProgressIndicator progress_;
			ApplyOperatorType applyOpLocal_;

	};     //class GroundStateTargetting

	template<
	template<typename,typename,typename> class LanczosSolverTemplate,
	template<typename,typename> class InternalProductTemplate,
	template<typename,typename> class WaveFunctionTransfTemplate,
	typename ModelType_,
	typename ConcurrencyType_,
	typename IoType_,
	template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const GroundStateTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,ConcurrencyType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		os<<"GSTWeightGroundState=1\n";
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif
