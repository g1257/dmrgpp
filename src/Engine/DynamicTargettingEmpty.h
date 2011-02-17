
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

#ifndef DYNAMICTARGETTING_H
#define DYNAMICTARGETTING_H

#include "ProgressIndicator.h"
#include "GroundStateParams.h"

namespace Dmrg {
	
	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
		template<typename,typename> class InternalProductTemplate,
		template<typename,typename> class WaveFunctionTransfTemplate,
		typename ModelType_,
		typename ConcurrencyType_,
		typename IoType_,
		template<typename> class VectorWithOffsetTemplate>
	class DynamicTargetting  {
	public:
		typedef ModelType_ ModelType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef IoType_ IoType;
		typedef typename ModelType::RealType RealType;
		typedef typename ModelType::ModelHelperType ModelHelperType;
		typedef typename ModelHelperType::LeftRightSuperType
						LeftRightSuperType;
		typedef typename LeftRightSuperType::BasisWithOperatorsType
				BasisWithOperatorsType;
		typedef GroundStateParams<ModelType> TargettingParamsType; //3a-01
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename BasisType::BlockType BlockType;
		typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
		typedef typename VectorWithOffsetType::VectorType VectorType;
		typedef VectorType TargetVectorType;
		typedef WaveFunctionTransfTemplate<LeftRightSuperType,
				VectorWithOffsetType> WaveFunctionTransfType;

		enum {	EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
						EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
						INFINITE=WaveFunctionTransfType::INFINITE};
		
		DynamicTargetting(
				const LeftRightSuperType& lrs,
				const ModelType& model,
				const TargettingParamsType& tstStruct,
				const WaveFunctionTransfType& wft)
		:	lrs_(lrs)
		{
			bogus();
		}

		RealType weight(size_t i) const
		{
			return 1.0;
		}

		RealType gsWeight() const
		{
			return 1.0;  // bogus!!
		}

		RealType normSquared(size_t i) const
		{
			return 0.0; // bogus
		}

		void load(const std::string& f) { }

		template<typename IoOutputType>
		void save(const std::vector<size_t>& block,IoOutputType& io) const
		{}

		template<typename SomeBasisType>
		void setGs(const std::vector<TargetVectorType>& v,
				const SomeBasisType& someBasis)
		{
			bogus();
		}

		const VectorType& operator[](size_t i) const { return psi_[i]; } // bogus!!

		VectorType& operator[](size_t i) { return psi_[i]; } // bogus!!

		const VectorWithOffsetType& gs() const { return psi_; } // bogus!!

		bool includeGroundStage() const {return true; }  // bogus!!

		size_t size() const
		{
			return 0;  // bogus!!
		}

		const VectorWithOffsetType& operator()(size_t i) const
		{
			return psi_; // bogus!!
		}

		void evolve(RealType Eg,size_t direction,const BlockType& block,
				size_t loopNumber)
		{
			bogus();
		}

		void initialGuess(VectorWithOffsetType& v) const
		{
			bogus();
		}

		const LeftRightSuperType& leftRightSuper() const
		{
			return lrs_;
		}
	private:
		
		void bogus() const
		{
			throw std::runtime_error("This class is bogus!! You need to use"
								" DynamicTargetting.h instead (which has GSL dependencies)!\n");
		}
		
		VectorWithOffsetType psi_;
		const LeftRightSuperType& lrs_;
	}; // class DynamicTargettingEmpty

	template<
		template<typename,typename,typename> class LanczosSolverTemplate,
			template<typename,typename> class InternalProductTemplate,
 		template<typename,typename> class WaveFunctionTransfTemplate,
			typename ModelType_,
 		typename ConcurrencyType_,
			typename IoType_,
   			template<typename> class VectorWithOffsetTemplate>
	std::ostream& operator<<(std::ostream& os,
			const DynamicTargetting<LanczosSolverTemplate,
			InternalProductTemplate,
			WaveFunctionTransfTemplate,ModelType_,ConcurrencyType_,IoType_,
			VectorWithOffsetTemplate>& tst)
	{
		os<<"DTEmpty=NothingToSeeHere\n";
		return os;
	}
} // namespace
#endif // DYNAMICTARGETTING_H

