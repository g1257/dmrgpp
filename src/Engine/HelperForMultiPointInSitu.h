#ifndef HELPERFORMULTIPOINTINSITU_H
#define HELPERFORMULTIPOINTINSITU_H
#include "FermionSign.h"
#include "Matrix.h"
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"
#include "DmrgSerializer.h"

namespace Dmrg {

template<typename CheckpointType>
class HelperForMultiPointInSitu {

public:

	typedef typename CheckpointType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename WaveFunctionTransfType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef FermionSign FermionSignType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename PsimagLite::Vector<LeftRightSuperType*>::Type VectorLeftRightSuperType;
	typedef DmrgSerializer<LeftRightSuperType, VectorWithOffsetType> DmrgSerializerType;
	typedef typename DmrgSerializerType::BlockDiagonalMatrixType BlockDiagonalMatrixType;

	class BogusInput {

	public:

		BogusInput(SizeType numberOfSites,
		           const CheckpointType& checkPoint,
		           const WaveFunctionTransfType& wft,
		           ProgramGlobals::DirectionEnum dir)
		    : numberOfSites_(numberOfSites), checkPoint_(checkPoint), wft_(wft), dir_(dir) {}

		~BogusInput()
		{
			const SizeType n = garbage_.size();
			for (SizeType i = 0; i < n; ++i) {
				delete garbage_[i];
				garbage_[i] = nullptr;
			}
		}

		SizeType numberOfSites() const { return numberOfSites_; }

		const LeftRightSuperType& hookForMultiInSituLrs(SizeType ind)
		{
			std::pair<BasisWithOperatorsType, BasisWithOperatorsType> pair =
			        checkPoint_.hookForMultiInSituLrs(ind);

			BasisType super("superForMultiPointInSitu");
			super.setToProduct(pair.first, pair.second);
			LeftRightSuperType* lrsPtr = new LeftRightSuperType(pair.first, pair.second, super);
			garbage_.push_back(lrsPtr);
			return *lrsPtr;
		}

		const BlockDiagonalMatrixType& getTransform(SizeType ind,
		                                            ProgramGlobals::DirectionEnum dir) const
		{
			return wft_.multiPointGetTransform(ind, dir);
		}

		// the direction remains the same for all the stack depth
		ProgramGlobals::DirectionEnum direction() const { return dir_; }

	private:

		SizeType numberOfSites_;
		const CheckpointType& checkPoint_;
		const WaveFunctionTransfType& wft_;
		ProgramGlobals::DirectionEnum dir_;
		VectorLeftRightSuperType garbage_;
	};

	typedef BogusInput IoInputType;

	HelperForMultiPointInSitu(BogusInput& io,
	                          SizeType start,
	                          SizeType nf,
	                          SizeType trail,
	                          bool withLegacyBugs,
	                          bool readOnDemand)
	    : io_(io), ind_(1 + numberOfSites())
	{}

	const LeftRightSuperType& leftRightSuper(SizeType ind) const
	{
		return io_.hookForMultiInSituLrs(ind);
	}

	SizeType numberOfSites() const { return io_.numberOfSites(); }

	const VectorWithOffsetType& getVectorFromBracketId(const PsimagLite::GetBraOrKet& braOrKet,
	                                                   SizeType index) const
	{
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::getVectorFromBracketId() " +
		                               PsimagLite::String("unimplemented\n"));
	}

	SizeType cols(SizeType ind) const
	{
		if (ind != ind_)
			computeAndSaveTransform(ind);
		return transform_.cols();
	}

	SizeType rows(SizeType ind) const
	{
		if (ind != ind_)
			computeAndSaveTransform(ind);
		return transform_.rows();
	}

	// transform O2 by the transformation in location ind, and put the result in ret
	void transform(SparseMatrixType& ret,
	               const SparseMatrixType& O2,
	               SizeType ind) const
	{
		if (ind != ind_)
			computeAndSaveTransform(ind);
		DmrgSerializerType::transform(ret, O2, transform_);
	}

	ProgramGlobals::DirectionEnum direction(SizeType) const
	{
		return io_.direction();
	}

	short int signsOneSite(SizeType) const
	{
		// given site in the one-site basis return the signs on site
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::signsOneSite() unimplemented\n");
	}

	const FermionSignType& fermionicSignLeft(SizeType) const
	{
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::fermionicSignLeft() unimplemented\n");
	}

private:

	void computeAndSaveTransform(SizeType ind) const
	{
		const ProgramGlobals::DirectionEnum dir = direction(ind);
		transform_ = io_.getTransform(ind, dir);
		ind_ = ind;
	}

	BogusInput& io_;
	mutable BlockDiagonalMatrixType transform_;
	mutable SizeType ind_;
};
}
#endif // HELPERFORMULTIPOINTINSITU_H
