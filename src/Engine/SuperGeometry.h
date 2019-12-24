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

	class SuperPlaquette {

	public:

		SuperPlaquette(InputType_& io) : value_(0)
		{
			io.readline(value_, "SuperPlaquetteValue=");
		}

		bool connected(SizeType emax, SizeType emin, const VectorSizeType& hItems) const
		{
			assert(hItems.size() == 4);

			return (hItems[1] == hItems[0] + 1 &&
			        hItems[2] == hItems[1] + 1 &&
			        hItems[3] == hItems[2] + 1);

		}

		void addSuperConnections(VectorVectorSizeType& data,
		                         SizeType smax,
		                         SizeType emin,
		                         SizeType linSize) const
		{
			return (smax + 1 == emin) ? addSuperConnectionsFinite_(data, smax, emin, linSize)
			                          : addSuperConnectionsInfinite_(data, smax, emin, linSize);
		}

		SizeType holloutRadius() const { return 4; }

		ComplexOrRealType operator()(SizeType smax,
		                             SizeType emin,
		                             const VectorSizeType& hItems,
		                             const VectorSizeType& edofs) const
		{
			return value_;
		}

	private:

		void addSuperConnectionsInfinite_(VectorVectorSizeType& data,
		                                  SizeType smax,
		                                  SizeType emin,
		                                  SizeType linSize) const
		{
			// FIXME: Add here site substitutions for when the lattice is not fully built
			return;
		}

		void addSuperConnectionsFinite_(VectorVectorSizeType& data,
		                                SizeType smax,
		                                SizeType emin,
		                                SizeType linSize) const
		{
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

		typename GeometryType::RealType value_;
	};

public:

	typedef typename GeometryType::RealType RealType;
	typedef PsimagLite::GeometryDca<RealType,GeometryType> GeometryDcaType;
	typedef typename PsimagLite::Vector<SuperPlaquette*>::Type VectorSuperPlaquetteType;

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
				auto ptr = new SuperPlaquette(io);
				superStrings_.push_back(ptr);
				hollowOutRadius_ = std::max(hollowOutRadius_, ptr->holloutRadius());
			} else {
				err("SuperGeometry " + geometry_.options(i) + " unsupported\n");
			}
		}
	}

	~SuperGeometry()
	{
		delete dcaPtr_;
		dcaPtr_ = nullptr;
		const SizeType n = superStrings_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete superStrings_[i];
			superStrings_[i] = nullptr;
		}
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
		const SizeType n = hItems.size();
		assert(n == edofs.size());
		if (n == 2)
			return geometry_(smax, emin, hItems[0], edofs[0], hItems[1], edofs[1], term);

		assert(superStrings_.size() == 1);
		return superStrings_[0]->operator()(smax, emin, hItems, edofs);

	}

	bool connected(SizeType smax,SizeType emin, const VectorSizeType& hItems) const
	{
		if (hItems.size() == 2)
			return geometry_.connected(smax, emin, hItems[0], hItems[1]);

		const SizeType n = superStrings_.size();
		for (SizeType i = 0; i < n; ++i)
			if (superStrings_[i]->connected(smax, emin, hItems))
				return true;

		return false;
	}

	typename ProgramGlobals::ConnectionEnum connectionKind(SizeType smax,
	                                                       const VectorSizeType& hItems) const
	{
		if (hItems.size() == 2)
			return geometry_.connectionKind(smax, hItems[0], hItems[1]);

		const SizeType n = hItems.size();
		SizeType flag = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (hItems[i] <= smax) flag |= 1;
			else flag |= 2;
		}

		return (flag == 3) ? ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON
		                   : ProgramGlobals::ConnectionEnum::SYSTEM_SYSTEM;
	}

	const GeometryDcaType& createDcaObject(SizeType orbitals) const
	{
		if (!dcaPtr_)
			dcaPtr_ = new GeometryDcaType(geometry_, orbitals);
		return *dcaPtr_;
	}

	void addSuperConnections(VectorVectorSizeType& data,
	                         SizeType smax,
	                         SizeType emin) const
	{
		const SizeType n = superStrings_.size();
		for (SizeType i = 0; i < n; ++i)
			superStrings_[i]->addSuperConnections(data, smax, emin, geometry_.numberOfSites());
	}

	friend std::ostream& operator<<(std::ostream& os, const SuperGeometry& supergeometry)
	{
		os<<supergeometry.geometry_;
		return os;
	}

private:

	const GeometryType geometry_;
	mutable GeometryDcaType* dcaPtr_;
	SizeType hollowOutRadius_;
	VectorSuperPlaquetteType superStrings_;
};

}
#endif // SUPERGEOMETRY_H
