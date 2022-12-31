#ifndef OPSFORLINK_HH
#define OPSFORLINK_HH
#include "Vector.h"
#include "Link.h"
#include "OperatorsCached.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class OpsForLink {

public:

	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename LeftRightSuperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef Link<ComplexOrRealType> LinkType;
	typedef typename PsimagLite::Vector<LinkType>::Type VectorLinkType;
	typedef OperatorsCached<LeftRightSuperType> OperatorsCachedType;

	OpsForLink(const OperatorsCachedType& operatorsCached,
	           const VectorLinkType& lps)
	    : lps_(lps),
	      operatorsCached_(operatorsCached),
	      link2_(nullptr),
	      A_(nullptr),
	      B_(nullptr)
	{}

	void setPointer(SizeType xx)
	{
		assert(xx < lps_.size());
		link2_ = &lps_[xx];

		assert(link2_->type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON ||
		       link2_->type == ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM);

		const ProgramGlobals::SysOrEnvEnum sysOrEnv =
		        (link2_->type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::SYSTEM : ProgramGlobals::SysOrEnvEnum::ENVIRON;
		const ProgramGlobals::SysOrEnvEnum envOrSys =
		        (link2_->type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::ENVIRON : ProgramGlobals::SysOrEnvEnum::SYSTEM;

		A_ = &operatorsCached_.getOpStorage(link2_->pairMetaOps.first,
		                                    sysOrEnv);
		B_ = &operatorsCached_.getOpStorage(link2_->pairMetaOps.second,
		                                    envOrSys);

		assert(A_);
		assert(B_);

		if (A_->invalid() || B_->invalid()) return;

		assert(isNonZeroMatrix(*A_));
		assert(isNonZeroMatrix(*B_));

		A_->checkValidity();
		B_->checkValidity();
	}

	bool invalid() const { return A_->invalid() || B_->invalid(); }

	const OperatorStorageType& A() const
	{
		assert(A_);
		return *A_;
	}

	const OperatorStorageType& B() const
	{
		assert(B_);
		return *B_;
	}

	const LinkType& link() const
	{
		assert(link2_);
		return *link2_;
	}

private:

	VectorLinkType lps_;
	OperatorsCachedType operatorsCached_;
	LinkType* link2_;
	OperatorStorageType const* A_;
	OperatorStorageType const* B_;
};
}
#endif // OPSFORLINK_HH
