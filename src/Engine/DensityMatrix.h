
#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include "Utils.h"

#include "BlockMatrix.h"

#include "DensityMatrixLocal.h"
#include "DensityMatrixSu2.h"

namespace Dmrg {

	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
		typename TargettingType
		>
	class DensityMatrix {

		enum {EXPAND_SYSTEM = TargettingType::EXPAND_SYSTEM };

		typedef typename DmrgBasisWithOperatorsType::SparseMatrixType
			SparseMatrixType;
		typedef typename TargettingType::TargetVectorType::value_type
			DensityMatrixElementType;
		typedef BlockMatrix<DensityMatrixElementType,
			psimag::Matrix<DensityMatrixElementType> > BlockMatrixType;
		typedef typename DmrgBasisType::FactorsType FactorsType;
		typedef DensityMatrixLocal<RealType,DmrgBasisType,
			DmrgBasisWithOperatorsType, TargettingType>
			DensityMatrixLocalType;
		typedef DensityMatrixSu2<RealType,DmrgBasisType,
			DmrgBasisWithOperatorsType,TargettingType>
			DensityMatrixSu2Type;
		typedef DensityMatrixBase<RealType,DmrgBasisType,
			DmrgBasisWithOperatorsType,TargettingType>
			DensityMatrixBaseType;

	public:
		typedef typename BlockMatrixType::BuildingBlockType
			BuildingBlockType;

		DensityMatrix(
			const TargettingType& target,
			const DmrgBasisWithOperatorsType& pBasis,
			const DmrgBasisWithOperatorsType& pBasisSummed,
			const DmrgBasisType& pSE,
			size_t direction,
			bool debug=false,
			bool verbose=false)
			

			: densityMatrixLocal_(target,pBasis,pBasisSummed,pSE,
					direction,debug,verbose),
				densityMatrixSu2_(target,pBasis,pBasisSummed,pSE,
					direction,debug,verbose)
		{

			if (DmrgBasisType::useSu2Symmetry()) {
				densityMatrixImpl_ = &densityMatrixSu2_;
			} else {
				densityMatrixImpl_ = &densityMatrixLocal_;
			}

			densityMatrixImpl_->init(target,pBasis,pBasisSummed,pSE,direction);
		}


		BlockMatrixType& operator()()
		{
			return densityMatrixImpl_->operator()();
		}


		size_t rank() { return densityMatrixImpl_->rank(); }
		

		void check(int direction)
		{
			return densityMatrixImpl_->check(direction);
		}
		

		void check2(int direction)
		{
			densityMatrixImpl_->check2(direction);
		}
		

		template<typename ConcurrencyType>
		void diag(std::vector<RealType>& eigs,char jobz,
				ConcurrencyType& concurrency)
		{
			if (!DmrgBasisType::useSu2Symmetry()) {
				densityMatrixLocal_.diag(eigs,jobz,concurrency);
			} else {
				densityMatrixSu2_.diag(eigs,jobz,concurrency);
			}
			
		}


		template<
			typename RealType_,
			typename DmrgBasisType_,
			typename DmrgBasisWithOperatorsType_,
   			typename TargettingType_
			> 
		friend std::ostream& operator<<(std::ostream& os,
			const DensityMatrix<RealType_,DmrgBasisType_,
				DmrgBasisWithOperatorsType_,TargettingType_>&
						dm);


	private:
		DensityMatrixLocalType densityMatrixLocal_;
		DensityMatrixSu2Type densityMatrixSu2_;
		DensityMatrixBaseType* densityMatrixImpl_;

		
	}; // class DensityMatrix

	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
  		typename TargettingType
		> 
	std::ostream& operator<<(std::ostream& os,
				const DensityMatrix<RealType,DmrgBasisType,
				DmrgBasisWithOperatorsType,TargettingType>& dm)
	{
		os<<(*dm.densityMatrixImpl_);
		return os;
	}
} // namespace Dmrg

#endif

