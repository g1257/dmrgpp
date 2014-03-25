#ifndef GEOMETRY_EX_H
#define GEOMETRY_EX_H

#include "Vector.h"
#include "String.h"

namespace PsimagLite {

template<typename GeometryTermType>
class GeometryExBase {

	typedef typename GeometryTermType::RealType RealType;

public:

	typedef typename Vector<RealType>::Type VectorRealType;

	SizeType dimension(SizeType i = 0) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: dimension\n");
	}
	
	void index2Rvector(SizeType r,VectorRealType& rvector) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: index2Rvector\n");
	}
	
	void index2Kvector(SizeType k,VectorRealType& kvector) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: index2Kvector\n");
	}
	
	//! Number of symmetry operations for this K Geometry.
	SizeType nGroupK() const
	{
		throw PsimagLite::RuntimeError("GeometryEx: nGroupK\n");
	}
	
	SizeType ickequ(SizeType j,SizeType op) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: ickequ\n");
	}
};

#ifndef USE_MS_GEOMETRY_H
template<typename GeometryTermType>
class GeometryEx  : public GeometryExBase<GeometryTermType> {

	typedef typename Vector<GeometryTermType*>::Type VectorGeometryTermType;
	
public:
		
	GeometryEx(SizeType linSize,VectorGeometryTermType& terms)
	{}
};
#else

template<typename GeometryTermType>
class GeometryEx  : public GeometryExBase<GeometryTermType> {

	typedef typename GeometryTermType::RealType RealType;
	typedef typename Vector<GeometryTermType*>::Type VectorGeometryTermType;
	
public:

	typedef typename Vector<RealType>::Type VectorRealType;
	
	GeometryEx(SizeType linSize,VectorGeometryTermType& terms)
	: linSize_(linSize), 
	  terms_(terms)
	  nGroupR_(icrequ_data.n_row()),
	  nGroupK_(ickequ_data.n_row())
	{}

	SizeType dimension(SizeType i = 0) const
	{
		assert(terms_.size() > i);
		PsimagLite::String s = terms_[i]->label();
		
		if (s == "chain") return 1;
		
		if (s=="ladder" ||  s=="ladderbath") return 2;
		
		throw PsimagLite::RuntimeError("GeometryEx: dimension\n");
	}
	
	void index2Rvector(SizeType r,VectorRealType& rvector) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: index2Rvector\n");
	}
	
	void index2Kvector(SizeType k,VectorRealType& kvector) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: index2Kvector\n");
	}
	
	//! Number of symmetry operations for this K Geometry.
	SizeType nGroupK() const
	{
		return nGroupK_; 
	}
	
	SizeType ickequ(SizeType j,SizeType op) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: ickequ\n");
	}

private:

	SizeType linSize_;
	VectorGeometryTermType& terms_;
};
#endif

} // namespace PsimagLite

#endif // GEOMETRY_EX_H