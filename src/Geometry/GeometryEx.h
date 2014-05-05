#ifndef GEOMETRY_EX_H
#define GEOMETRY_EX_H

#include "Vector.h"
#include "String.h"

#ifndef USE_MS_GEOMETRY

namespace PsimagLite {

template<typename RealType,typename InputType>
class GeometryEx {

public:

	typedef typename Vector<RealType>::Type VectorRealType;

	GeometryEx(InputType& io, SizeType meshPoints)
	{
		meshLength_ = sqrt(meshPoints);
		meshStep_ = static_cast<RealType>(2*M_PI/meshLength_);
		PsimagLite::String str;
		io.readline(str,"GeometryKind=");
		if (str != "star")
			throw RuntimeError("GeometryEx: expecting GeometryKind=star\n");
	}

	SizeType dimension(SizeType i = 0) const
	{
		return 2;
	}

	void index2Rvector(SizeType r,VectorRealType& rvector) const
	{
		rvector.resize(2);
		rvector[0] = 0;
		rvector[1] = 0;
	}

	void index2Kvector(SizeType k,VectorRealType& kvector) const
	{
		kvector.resize(2);
		kvector[0] = 0;
		kvector[1] = 0;
	}

	//! Number of symmetry operations for this K Geometry.
	SizeType nGroupK() const
	{
		return 1;	
	}

	SizeType ickequ(SizeType j,SizeType op) const
	{
		return 0;	
	}

	void getMeshVector(VectorRealType& kvector,SizeType k) const
	{
		div_t q = div(k,meshLength_);
		kvector[0] = -M_PI + q.quot*meshStep_;
		kvector[1] = -M_PI + q.rem*meshStep_;	
	}

	SizeType sizeOfMesh() const
	{
		return meshLength_*meshLength_;
	}

private:

	SizeType meshLength_;
	RealType meshStep_;
};

}

#else

#include "msGeometry.h"

namespace PsimagLite {

template<typename RealType,typename InputType>
class GeometryEx  {

public:

	typedef typename Vector<RealType>::Type VectorRealType;

	GeometryEx(InputType& io,SizeType meshPoints)
	{
		SizeType ly = 2;
		SizeType lx = 2;
		std::cerr<<"WARNING: GeometryEx(): lattice of 2x2 hard wired\n";
		msGeometryInit2D(lx,0,0,ly,meshPoints);
	}

	~GeometryEx()
	{
		msGeometryEnd();
	}

	SizeType dimension(SizeType i = 0) const
	{
		return msGeometryDim();
	}

	void index2Rvector(SizeType r,VectorRealType& rvector) const
	{
		return msGeometryIndex2Rvector(r,rvector);
	}

	void index2Kvector(SizeType k,VectorRealType& kvector) const
	{
		return msGeometryIndex2Kvector(k,kvector);
	}

	//! Number of symmetry operations for this K Geometry.
	SizeType nGroupK() const
	{
		return msGeometryNgroupK();
	}

	SizeType ickequ(SizeType j,SizeType op) const
	{
		return msGeometryIckequ(j,op);
	}

	void getMeshVector(VectorRealType& kvector,SizeType k) const
	{
		return msGeometryGetMeshVector(kvector,k);
	}

	SizeType sizeOfMesh() const
	{
		return msGeometrySizeOfMesh();
	}
};

} // namespace PsimagLite
#endif

#endif // GEOMETRY_EX_H

