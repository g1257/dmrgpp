#ifndef SUPERGEOMETRY_H
#define SUPERGEOMETRY_H
#include "Geometry/Geometry.h"
#include "ProgramGlobals.h"
#include "Geometry/GeometryDca.h"

namespace Dmrg {

template<typename ComplexOrRealType_,typename InputType_, typename ProgramGlobalsType>
class SuperGeometry {

	typedef PsimagLite::Geometry<ComplexOrRealType_, InputType_, ProgramGlobalsType> GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef typename GeometryType::VectorSizeType VectorSizeType;
	typedef typename PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

public:

	typedef typename GeometryType::RealType RealType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	SuperGeometry(InputType_& io)
	    : geometry_(io), dcaPtr_(0), hollowOutRadius_(0)
	{
		// add super terms as needed
		const SizeType n = geometry_.terms();
		for (SizeType i = 0; i < n; ++i) {
			if (geometry_.directions(i) > 0) continue;
			// super term found
			// it's gotta be "SuperPlaquette" for now, (only one option, sorry!)
			if (geometry_.options(i) == "SuperPlaquette") {
				if (std::find(superStrings_.begin(),
				              superStrings_.end(),
				              "SuperPlaquette") != superStrings_.end())
					continue;
				superStrings_.push_back("SuperPlaquette");
				hollowOutRadius_ = 4;
			}
		}
	}

	~SuperGeometry()
	{
		delete dcaPtr_;
		dcaPtr_ = 0;
	}

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

	SizeType hollowOutRadius(SizeType maxLeft) const
	{
		return std::max(maxLeft*geometry_.maxConnections(), hollowOutRadius_);
	}

	SizeType orbitals(SizeType term, SizeType site) const
	{
		return geometry_.orbitals(term, site);
	}

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

	const GeometryDcaType& createDcaObject(SizeType orbitals) const
	{
		if (!dcaPtr_)
			dcaPtr_ = new GeometryDcaType(geometry_, orbitals);
		return *dcaPtr_;
	}

	SizeType overSize(SizeType blockSize) const
	{
		return blockSize*(blockSize/2 + 1); // + superc_.size();
	}

	SizeType addSuperConnections(VectorVectorSizeType& data,
	                             SizeType smax,
	                             SizeType emin,
	                             const VectorSizeType& block,
	                             SizeType counter) const
	{
		SizeType c = addSuperConnections_(data, smax, emin, block);
		return c + counter;
	}

	friend std::ostream& operator<<(std::ostream& os, const SuperGeometry& supergeometry)
	{
		os<<supergeometry.geometry_;
		return os;
	}

private:

	SizeType addSuperConnections_(VectorVectorSizeType& data,
	                              SizeType smax,
	                              SizeType emin,
	                              const VectorSizeType& block)
	{
		return (smax + 1 == emin) ? addSuperConnectionsFinite_(data, smax, emin, block) : 0;
		//		                          : addSuperConnectionsInfinite_(data, smax, emin, block);
	}

	SizeType addSuperConnectionsInfinite_(VectorVectorSizeType& data,
	                                      SizeType smax,
	                                      SizeType emin,
	                                      const VectorSizeType& block)
	{
		if (superStrings_.size() == 0) return 0;


		assert(superStrings_.size() == 1);

		assert(superStrings_[0] == "SuperPlaquette");

		const SizeType linSize = geometry_.numberOfSites();

		// smax - 1, smax, emin, emin + 1
		if (smax > 0 && emin + 1 < linSize)
			data.push_back(VectorSizeType{smax - 1, smax, emin, emin + 1});

		// smax, emin, emin + 1, emin + 2
		if (emin + 2 < linSize)
			data.push_back(VectorSizeType{smax, emin, emin + 1, emin + 2});

		// smax - 2, smax -1, smax, emin
		if (smax > 1)
			data.push_back(VectorSizeType{smax - 2, smax - 1, smax, emin});
	}

	static void checkVectorHasTwoEntries(const VectorSizeType& hItems)
	{
		if (hItems.size() != 2)
			err("SuperGeometry unimplemented\n");
	}

	const GeometryType geometry_;
	mutable GeometryDcaType* dcaPtr_;
	SizeType hollowOutRadius_;
	VectorStringType superStrings_;
};

}
#endif // SUPERGEOMETRY_H
