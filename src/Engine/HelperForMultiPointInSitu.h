#ifndef HELPERFORMULTIPOINTINSITU_H
#define HELPERFORMULTIPOINTINSITU_H
#include "FermionSign.h"
#include "Matrix.h"
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename VectorWithOffsetType_, typename LeftRightSuperType>
class HelperForMultiPointInSitu {

public:

	class BogusInput {

	public:

		BogusInput(SizeType numberOfSites) : numberOfSites_(numberOfSites) {}

		SizeType numberOfSites() const { return numberOfSites_; }

	private:

		SizeType numberOfSites_;

	};

	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef BogusInput IoInputType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef FermionSign FermionSignType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;

	HelperForMultiPointInSitu(BogusInput& io,
	                          SizeType start,
	                          SizeType nf,
	                          SizeType trail,
	                          bool withLegacyBugs,
	                          bool readOnDemand)
	    : io_(io)
	{}

	const LeftRightSuperType& leftRightSuper(SizeType) const
	{
		throw PsimagLite::RuntimeError("HelperForMultiPointInSitu::lrs() unimplemented\n");
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
