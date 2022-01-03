#ifndef HELPERFORMULTIPOINTINSITU_H
#define HELPERFORMULTIPOINTINSITU_H
#include "FermionSign.h"
#include "Matrix.h"
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"

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

	class BogusInput {

	public:

		BogusInput(SizeType numberOfSites,
		           const CheckpointType& checkPoint,
		           const WaveFunctionTransfType& wft)
		    : numberOfSites_(numberOfSites), checkPoint_(checkPoint), wft_(wft) {}

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

	private:

		SizeType numberOfSites_;
		const CheckpointType& checkPoint_;
		const WaveFunctionTransfType& wft_;
		VectorLeftRightSuperType garbage_;
	};

	typedef BogusInput IoInputType;

	HelperForMultiPointInSitu(BogusInput& io,
	                          SizeType start,
	                          SizeType nf,
	                          SizeType trail,
	                          bool withLegacyBugs,
	                          bool readOnDemand)
	    : io_(io)
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

	SizeType cols(SizeType) const
	{
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::cols() unimplemented\n");
		//return transform_.cols();
	}

	SizeType rows() const
	{
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::rows() unimplemented\n");
		//return transform_.rows();
	}

	void transform(SparseMatrixType& ret,
	               const SparseMatrixType& O2,
	               SizeType ind) const
	{
		// transform O2 by the transformation in location ind, and put the result in ret
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::transform() unimplemented\n");
	}

	ProgramGlobals::DirectionEnum direction(SizeType) const
	{
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::direction() unimplemented\n");
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

	BogusInput& io_;
};
}
#endif // HELPERFORMULTIPOINTINSITU_H
