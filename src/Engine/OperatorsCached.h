#ifndef OPERATORSCACHED_H
#define OPERATORSCACHED_H
#include "Concurrency.h"
#include "MetaOpForConnection.hh"
#include "ProgramGlobals.h"

namespace Dmrg {

template <typename LeftRightSuperType> class OperatorsCached {

public:

	typedef std::pair<SizeType, SizeType>                           PairType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType     BasisWithOperatorsType;
	typedef typename LeftRightSuperType::BasisType                  BasisType;
	typedef typename BasisType::BlockType                           BlockType;
	typedef typename BasisWithOperatorsType::OperatorsType          OperatorsType;
	typedef typename OperatorsType::OperatorType                    OperatorType;
	typedef typename OperatorType::StorageType                      OperatorStorageType;
	typedef PsimagLite::Concurrency                                 ConcurrencyType;
	typedef typename PsimagLite::Vector<OperatorStorageType*>::Type VectorOperatorStorageType;
	typedef typename PsimagLite::Vector<VectorOperatorStorageType>::Type
	                  VectorVectorOperatorStorageType;
	typedef BlockType VectorSizeType;

	OperatorsCached(const LeftRightSuperType& lrs)
	    : lrs_(lrs)
	    , garbage_(ConcurrencyType::codeSectionParams.npthreads)
	    , seen_(ConcurrencyType::codeSectionParams.npthreads)
	{
		ConcurrencyType::mutexInit(&mutex_);
	}

	~OperatorsCached()
	{
		const SizeType n = garbage_.size();
		for (SizeType i = 0; i < n; ++i) {
			const SizeType m = garbage_[i].size();
			for (SizeType j = 0; j < m; ++j) {
				delete garbage_[i][j];
				garbage_[i][j] = 0;
			}
		}

		ConcurrencyType::mutexDestroy(&mutex_);
	}

	void clearThreadSelves() const { threadSelves_.clear(); }

	const OperatorStorageType& getOpStorage(const MetaOpForConnection&         metaOp,
	                                        const ProgramGlobals::SysOrEnvEnum type) const
	{
		return (metaOp.site >= 0) ? getOpStorageLocal(metaOp.index, metaOp.modifier, type)
		                          : getOpStorageNonLocal(metaOp, type);
	}

private:

	const OperatorStorageType& getOpStorageLocal(SizeType                           iifirst,
	                                             char                               modifier,
	                                             const ProgramGlobals::SysOrEnvEnum type) const
	{
		const OperatorStorageType* m = 0;
		if (type == ProgramGlobals::SysOrEnvEnum::SYSTEM) {
			m = &(lrs_.left().localOperator(iifirst).getStorage());
		} else {
			assert(type == ProgramGlobals::SysOrEnvEnum::ENVIRON);
			m = &(lrs_.right().localOperator(iifirst).getStorage());
		}

		m->checkValidity();
		if (modifier == 'N')
			return *m;

		assert(modifier == 'C');
		SizeType typeIndex = (type == ProgramGlobals::SysOrEnvEnum::SYSTEM) ? 0 : 1;
		SizeType packed    = typeIndex + iifirst * 2;
		const ConcurrencyType::PthreadtType threadSelf = ConcurrencyType::threadSelf();
		const SizeType                      threadNum  = threadNumberFromSelf(threadSelf);

		if (garbage_.size() != seen_.size())
			err("reducedOperator: FATAL: internal error\n");

		if (garbage_.size() <= threadNum || seen_.size() <= threadNum)
			err("reducedOperator: FATAL: " + ttos(threadNum)
			    + " >= " + ttos(garbage_.size()) + "\n");

		int indexOfSeen = PsimagLite::indexOrMinusOne(seen_[threadNum], packed);
		if (indexOfSeen >= 0) {
			assert(static_cast<SizeType>(indexOfSeen) < garbage_[threadNum].size());
			return *(garbage_[threadNum][indexOfSeen]);
		}

		OperatorStorageType* mc = new OperatorStorageType;
		transposeConjugate(*mc, *m);
		garbage_[threadNum].push_back(mc);
		seen_[threadNum].push_back(packed);
		mc->checkValidity();

		return *mc;
	}

	SizeType threadNumberFromSelf(ConcurrencyType::PthreadtType threadSelf) const
	{
		ConcurrencyType::mutexLock(&mutex_);

		int threadPreNum = PsimagLite::indexOrMinusOne(threadSelves_, threadSelf);
		if (threadPreNum < 0) {
			threadPreNum = threadSelves_.size();
			threadSelves_.push_back(threadSelf);
		}

		ConcurrencyType::mutexUnlock(&mutex_);

		return threadPreNum;
	}

	const OperatorStorageType&
	getOpStorageNonLocal(const MetaOpForConnection&         metaOp,
	                     const ProgramGlobals::SysOrEnvEnum type) const
	{
		// Non local ops cannot have modifier different than 'N'
		assert(metaOp.modifier == 'N');
		return (type == ProgramGlobals::SysOrEnvEnum::SYSTEM)
		    ? lrs_.left().getSuperByIndex(metaOp.index).getStorage()
		    : lrs_.right().getSuperByIndex(metaOp.index).getStorage();
	}

	const LeftRightSuperType&                                       lrs_;
	mutable VectorVectorOperatorStorageType                         garbage_;
	mutable typename PsimagLite::Vector<BlockType>::Type            seen_;
	mutable ConcurrencyType::MutexType                              mutex_;
	mutable PsimagLite::Vector<ConcurrencyType::PthreadtType>::Type threadSelves_;
};
}
#endif // OPERATORSCACHED_H
