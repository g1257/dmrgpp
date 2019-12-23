#ifndef SUPERGEOMETRY_H
#define SUPERGEOMETRY_H
#include "Geometry/Geometry.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ComplexOrRealType_,typename InputType_, typename ProgramGlobalsType>
class SuperGeometry {

	typedef PsimagLite::Geometry<ComplexOrRealType_, InputType_, ProgramGlobalsType> GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef typename GeometryType::VectorSizeType VectorSizeType;
	typedef typename PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

public:

	typedef typename GeometryType::RealType RealType;

	SuperGeometry(InputType_& io)
	    : geometry_(io)
	{}

	// const GeometryType& geometry() const { return geometry_; }

	void split(SizeType sitesPerBlock,
	           VectorSizeType& S,
	           VectorVectorSizeType& X,
	           VectorVectorSizeType& Y,
	           VectorSizeType& E,
	           bool allInSystem = false) const
	{
		geometry_.split(sitesPerBlock, S, X, Y, E, allInSystem);
	}

	SizeType numberOfSites() const { return geometry_.numberOfSites(); }

	SizeType terms() const { return geometry_.terms(); }

	void write(PsimagLite::String label, PsimagLite::IoNgSerializer& ioSerializer) const
	{
		geometry_.write(label, ioSerializer);
	}

//	SizeType maxConnections(SizeType termId) const
//	{
//		return geometry_.maxConnections(termId);
//	}

	SizeType maxConnections() const { return geometry_.maxConnections(); }

	PsimagLite::String label(SizeType i) const { return geometry_.label(i); }

	template<typename T>
	typename PsimagLite::EnableIf<PsimagLite::IsComplexNumber<T>::True ||
	Loki::TypeTraits<T>::isStdFloat,
	T>::Type vModifier(SizeType term, T value, RealType time) const
	{
		return geometry_.vModifier(term, value, time);
	}

	ComplexOrRealType operator()(SizeType smax,
	                             SizeType emin,
	                             const VectorSizeType& hItems,
	                             const VectorSizeType& edofs,
	                             SizeType term) const
	{
		checkVectorHasTwoEntries(hItems);
		checkVectorHasTwoEntries(edofs);
		return geometry_(smax, emin, hItems[0], edofs[0], hItems[1], edofs[1], term);
	}

	bool connected(SizeType smax,SizeType emin, const VectorSizeType& hItems) const
	{
		checkVectorHasTwoEntries(hItems);

		return geometry_.connected(smax, emin, hItems[0], hItems[1]);
	}

	typename ProgramGlobals::ConnectionEnum connectionKind(SizeType smax,
	                                                       const VectorSizeType& hItems) const
	{
		checkVectorHasTwoEntries(hItems);

		return geometry_.connectionKind(smax, hItems[0], hItems[1]);
	}

	SizeType overSize(SizeType blockSize) const
	{
		return blockSize*(blockSize/2 + 1); // + superc_.size();
	}

//	SizeType addSuperConnections(VectorVectorSizeType& data,
//	                             SizeType smax,
//	                             SizeType emin,
//	                             const VectorSizeType& block,
//	                             SizeType counter) const
//	{
//		SizeType c = superc_.addSuperConnections(data, smax, emin, block);
//		return c + counter;
//	}

	friend std::ostream& operator<<(std::ostream& os, const SuperGeometry& supergeometry)
	{
		os<<supergeometry.geometry_;
		return os;
	}

private:

	static void checkVectorHasTwoEntries(const VectorSizeType& hItems)
	{
		if (hItems.size() != 2)
			err("SuperGeometry unimplemented\n");
	}

	const GeometryType& geometry_;
};

}
#endif // SUPERGEOMETRY_H
