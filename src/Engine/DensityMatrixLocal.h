
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

#ifndef DENSITY_MATRIX_LOCAL_H
#define DENSITY_MATRIX_LOCAL_H
#include "ProgressIndicator.h"
#include "TypeToString.h"
#include "BlockMatrix.h"
#include "DensityMatrixBase.h"
#include "ParallelDensityMatrix.h"
#include "NoPthreads.h"

namespace Dmrg {
	//!
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
		typename TargettingType
		>

	class DensityMatrixLocal : public DensityMatrixBase<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType> {
		typedef typename DmrgBasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename TargettingType::VectorWithOffsetType TargetVectorType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;

		typedef typename DmrgBasisType::FactorsType FactorsType;
		typedef PsimagLite::ProgressIndicator ProgressIndicatorType;

		enum {EXPAND_SYSTEM = TargettingType::EXPAND_SYSTEM };

		static int const parallelRank_ = 0; // FIXME <-- ADD PARALLELISM

	public:

		typedef BlockMatrix<DensityMatrixElementType,PsimagLite::Matrix<DensityMatrixElementType> > BlockMatrixType;
		typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;

		DensityMatrixLocal(
			const TargettingType& target,
			const DmrgBasisWithOperatorsType& pBasis,
			const DmrgBasisWithOperatorsType& pBasisSummed,
			const DmrgBasisType& pSE,
			size_t direction,bool debug=false,bool verbose=false) 
		:
			progress_("DensityMatrixLocal",parallelRank_),
			data_(pBasis.size(),
			pBasis.partition()-1),
			debug_(debug),verbose_(verbose)
		{
		}


		virtual BlockMatrixType& operator()()
		{
			return data_;
		}

		virtual size_t rank() { return data_.rank(); }

		virtual void check(int direction)
		{
		}

		virtual void check2(int direction)
		{
		}

		template<typename ConcurrencyType>
		void diag(std::vector<RealType>& eigs,char jobz,ConcurrencyType& concurrency)
		{
			diagonalise<DensityMatrixElementType,RealType,ConcurrencyType>(data_,eigs,jobz,concurrency);
		}

		virtual void init(
				const TargettingType& target,
				DmrgBasisWithOperatorsType const &pBasis,
				const DmrgBasisWithOperatorsType& pBasisSummed,
				DmrgBasisType const &pSE,
				int direction)
		{
			typedef ParallelDensityMatrix<RealType,BlockMatrixType,DmrgBasisWithOperatorsType,TargettingType> ParallelDensityMatrixType;
			{
				std::ostringstream msg;
				msg<<"Init partition for all targets";
				progress_.printline(msg,std::cout);
			}
			//loop over all partitions:
			for (size_t m=0;m<pBasis.partition()-1;m++) {
				// size of this partition
				size_t bs = pBasis.partition(m+1)-pBasis.partition(m);
				
				// density matrix block for this partition:
				BuildingBlockType matrixBlock(bs,bs);
				
				// weight of the ground state:
				RealType w = target.gsWeight();
				
				// if we are to target the ground state do it now:
				ParallelDensityMatrixType helperDm(target,pBasis,pBasisSummed,pSE,direction,m,matrixBlock);

				if (target.includeGroundStage())
					helperDm.initPartition(matrixBlock,pBasis,m,target.gs(),
							pBasisSummed,pSE,direction,w);
				
				// target all other states if any:
				if (target.size()>0) {

					PsimagLite::NoPthreads<ParallelDensityMatrixType> threadedDm;
					PsimagLite::NoPthreads<ParallelDensityMatrixType>::setThreads(target.model().params().nthreads);
					if (threadedDm.name()=="pthreads") {
						std::ostringstream msg;
						msg<<"Threading with "<<threadedDm.threads();
						progress_.printline(msg,std::cout);
					}

					threadedDm.loopCreate(target.size(),helperDm,target.model().concurrency());

//					for (size_t i=0;i<target.size();i++) {
//						w = target.weight(i)/target.normSquared(i);
//						initPartition(matrixBlock,pBasis,m,target(i),
//						              pBasisSummed,pSE,direction,w);
//					}
				}
				// set this matrix block into data_
				data_.setBlock(m,pBasis.partition(m),matrixBlock);
			}
			{
				std::ostringstream msg;
				msg<<"Done with init partition";
				progress_.printline(msg,std::cout);
			}
		}

		template<
			typename RealType_,
			typename DmrgBasisType_,
			typename DmrgBasisWithOperatorsType_,
   			typename TargettingType_
			> 
		friend std::ostream& operator<<(std::ostream& os,
				const DensityMatrixLocal<RealType_,
    					DmrgBasisType_,DmrgBasisWithOperatorsType_,TargettingType_>& dm);

	private:

		ProgressIndicatorType progress_;
		BlockMatrixType data_;
		bool debug_,verbose_;

	}; // class DensityMatrixLocal

	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
  		typename TargettingType
		> 
	std::ostream& operator<<(std::ostream& os,
				const DensityMatrixLocal<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType>& dm)
	{
		for (size_t m=0;m<dm.data_.blocks();m++) {
			size_t ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<" ne="<<ne<<"\n"; 
			os<<dm.data_(m)<<"\n";
		}
		return os;
	}
} // namespace Dmrg

#endif
