#ifndef PERSISTENTSVD_H
#define PERSISTENTSVD_H
#include "Vector.h"

namespace Dmrg {

template<typename VectorMatrixType, typename VectorVectorRealType, typename VectorQnType, bool>
class PersistentSvd {

public:

	typedef typename VectorMatrixType::value_type MatrixType;
	typedef typename VectorVectorRealType::value_type VectorRealType;
	typedef typename VectorQnType::value_type QnType;

	PersistentSvd() : qns_(QnType::zero()) {}

	void resize(SizeType) {}

	MatrixType& vts(SizeType)
	{
		vt_.clear();
		return vt_;
	}

	VectorRealType& s(SizeType)
	{
		return s_;
	}

	QnType& qns(SizeType)
	{
		return qns_;
	}

	const VectorMatrixType& vts() const { return vts_; }

	const VectorVectorRealType& s() const { return sv_; }

	const VectorQnType& qns() const { return qnsv_; }

private:

	MatrixType vt_;
	VectorRealType s_;
	QnType qns_;
	VectorMatrixType vts_;
	VectorVectorRealType sv_;
	VectorQnType qnsv_;
};

// needef for WFT
template<typename VectorMatrixType, typename VectorVectorRealType, typename VectorQnType>
class PersistentSvd<VectorMatrixType, VectorVectorRealType, VectorQnType, true> {

public:

	typedef typename VectorMatrixType::value_type MatrixType;
	typedef typename VectorVectorRealType::value_type VectorRealType;
	typedef typename VectorQnType::value_type QnType;

	void resize(SizeType n)
	{
		vts_.resize(n);
		s_.resize(n);
		qns_.resize(n, QnType::zero());
	}

	MatrixType& vts(SizeType igroup)
	{
		assert(igroup < vts_.size());
		return vts_[igroup];
	}

	VectorRealType& s(SizeType igroup)
	{
		assert(igroup < s_.size());
		return s_[igroup];
	}

	QnType& qns(SizeType igroup)
	{
		assert(igroup < qns_.size());
		return qns_[igroup];
	}

	const VectorMatrixType& vts() const { return vts_; }

	const VectorVectorRealType& s() const { return s_; }

	const VectorQnType& qns() const { return qns_; }

private:

	VectorMatrixType vts_; // needed for WFT
	typename PsimagLite::Vector<VectorRealType>::Type s_; // needed for WFT
	VectorQnType qns_; // needed for WFT
};

}
#endif // PERSISTENTSVD_H
