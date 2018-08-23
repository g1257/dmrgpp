#ifndef PERSISTENTSVD_H
#define PERSISTENTSVD_H
#include "Vector.h"

namespace Dmrg {

// needef for WFT
template<typename VectorMatrixType, typename VectorVectorRealType, typename VectorQnType>
class PersistentSvd {

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

	void clear()
	{
		vts_.clear();
		s_.clear();
		qns_.clear();
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
