#ifndef OFFSETVECTOR_HPP
#define OFFSETVECTOR_HPP
#include "VectorWithOffsets.h"
#include "Qn.h"

namespace Dmrg {

template<typename ComplexOrRealType>
class OffsetVectorBase {
public:
	virtual ~OffsetVectorBase() { }

	virtual const ComplexOrRealType& slowAccess(SizeType) const = 0;
};

template<typename ComplexOrRealType>
class OffsetVectorAny : public OffsetVectorBase<ComplexOrRealType> {
public:

	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType, Qn>;

	OffsetVectorAny(const VectorWithOffsetsType& v)
	    : v_(v)
	{
		if (v_.sectors() == 1)
			std::cerr<<"Use OffsetVectorOne instead of OffsetVectorAny for performance\n";
	}

	const ComplexOrRealType& slowAccess(SizeType ind) const final
	{
		return v_.slowAccess(ind);
	}

private:
	const VectorWithOffsetsType& v_;

};

template<typename ComplexOrRealType>
class OffsetVectorOne :  public OffsetVectorBase<ComplexOrRealType> {
public:

	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType, Qn>;

	OffsetVectorOne(const VectorWithOffsetsType& v)
	    : v_(v),sector_(0),offset_(0)
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
	SizeType sector_;
	SizeType offset_;

};

template<typename ComplexOrRealType>
class OffsetVector {
public:
	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType, Qn>;

        const OffsetVectorBase<ComplexOrRealType>& makeOffsetVector(const VectorWithOffsetsType& v)
	{
		if (v.sectors() == 1) {
			base_ptr_ = new OffsetVectorOne<ComplexOrRealType>(v);
		} else {
			base_ptr_ = new OffsetVectorAny<ComplexOrRealType>(v);
		}

		return *base_ptr_;
	}

	~OffsetVector()
	{
		delete base_ptr_;
		base_ptr_ = 0;
	}

private:
	OffsetVectorBase<ComplexOrRealType>* base_ptr_ = nullptr;
};

}
#endif // OFFSETVECTOR_HPP
