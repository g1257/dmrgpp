#ifndef GEOMETRY_EX_H
#define GEOMETRY_EX_H

#include "Vector.h"
#include "String.h"

namespace PsimagLite {

template<typename GeometryTermType>
class GeometryEx {

	typedef typename GeometryTermType::RealType RealType;
	typedef typename Vector<GeometryTermType*>::Type VectorGeometryTermType;
	
public:

	typedef typename Vector<RealType>::Type VectorRealType;
	
	GeometryEx(SizeType linSize,VectorGeometryTermType& terms)
	: linSize_(linSize), terms_(terms)
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
	
	SizeType nGroupK() const
	{
		throw PsimagLite::RuntimeError("GeometryEx: nGroupK\n");
	}
	
	SizeType ickequ(SizeType j,SizeType op) const
	{
		throw PsimagLite::RuntimeError("GeometryEx: ickequ\n");
	}

private:

	SizeType linSize_;
	VectorGeometryTermType& terms_;
};

} // namespace PsimagLite

#endif // GEOMETRY_EX_H