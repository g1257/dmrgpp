#ifndef OFFSETVECTOR_HPP
#define OFFSETVECTOR_HPP
#include "Qn.h"
#include "VectorWithOffsets.h"

namespace Dmrg {

template <typename ComplexOrRealType> class OffsetVectorBase {
public:

	virtual ~OffsetVectorBase() { }

	virtual const ComplexOrRealType& slowAccess(SizeType) const = 0;
};

template <typename ComplexOrRealType>
class OffsetVectorAny : public OffsetVectorBase<ComplexOrRealType> {
public:

	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType, Qn>;

	OffsetVectorAny(const VectorWithOffsetsType& v)
	    : v_(v)
	{
		if (v_.sectors() == 1)
			std::cerr
			    << "Use OffsetVectorOne instead of OffsetVectorAny for performance\n";
	}

	const ComplexOrRealType& slowAccess(SizeType ind) const final { return v_.slowAccess(ind); }

private:

	const VectorWithOffsetsType& v_;
};

template <typename ComplexOrRealType>
class OffsetVectorOne : public OffsetVectorBase<ComplexOrRealType> {
public:

	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType, Qn>;

	OffsetVectorOne(const VectorWithOffsetsType& v)
	    : v_(v)
	    , sector_(0)
	    , offset_(0)
	{
		if (v_.sectors() != 1)
			err("OffsetVectorOne can only be used when sectors==1\n");
		sector_ = v_.sector(0);
		offset_ = v_.offset(sector_);
	}

	const ComplexOrRealType& slowAccess(SizeType ind) const final
	{
		return v_.fastAccess(sector_, ind - offset_);
	}

private:

	const VectorWithOffsetsType& v_;
	SizeType                     sector_;
	SizeType                     offset_;
};

template <typename ComplexOrRealType> class OffsetVector {
public:

	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType, Qn>;

	const OffsetVectorBase<ComplexOrRealType>& makeOffsetVector(const VectorWithOffsetsType& v)
	{
		OffsetVectorBase<ComplexOrRealType>* base_ptr = nullptr;
		if (v.sectors() == 1) {
			base_ptr = new OffsetVectorOne<ComplexOrRealType>(v);
		} else {
			base_ptr = new OffsetVectorAny<ComplexOrRealType>(v);
		}

		base_ptrs_.push_back(base_ptr); // for garbage collection

		return *base_ptr;
	}

	~OffsetVector()
	{
		for (SizeType i = 0; i < base_ptrs_.size(); ++i) {
			delete base_ptrs_[i];
			base_ptrs_[i] = nullptr;
		}
	}

private:

	std::vector<std::unique_ptr<OffsetVectorBase<ComplexOrRealType>>> base_ptrs_;
};

}
#endif // OFFSETVECTOR_HPP
