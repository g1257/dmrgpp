/*
 */

#ifndef BASISRASHBASOC_H
#define BASISRASHBASOC_H
#include "../../Engine/BasisBase.h"
#include "../HubbardOneOrbital/BasisOneSpin.h"

namespace LanczosPlusPlus {

template <typename GeometryType> class BasisRashbaSOC : public BasisBase<GeometryType> {

	enum
	{
		SPIN_UP   = LanczosGlobals::SPIN_UP,
		SPIN_DOWN = LanczosGlobals::SPIN_DOWN
	};

public:

	typedef LanczosGlobals::PairIntType                     PairIntType;
	typedef BasisOneSpin                                    BasisType;
	typedef BasisBase<GeometryType>                         BaseType;
	typedef typename BaseType::WordType                     WordType;
	typedef std::pair<WordType, WordType>                   PairWordType;
	typedef typename PsimagLite::Vector<PairWordType>::Type VectorPairWordType;
	typedef typename BaseType::VectorWordType               VectorWordType;
	typedef typename BaseType::LabeledOperatorType          LabeledOperatorType;

	BasisRashbaSOC(const GeometryType& geometry, SizeType ne)
	    : ne_(ne)
	{
		const SizeType nsite = geometry.numberOfSites();
		LanczosGlobals::doCombinatorial(2 * nsite + 2);
		const SizeType hilbert = LanczosGlobals::combinatorial(2 * nsite, ne);
		data_.resize(hilbert);
		SizeType k = 0;
		for (SizeType ndown = 0; ndown <= ne; ++ndown) {

			const SizeType nup = ne - ndown;
			BasisOneSpin   basisUp(nsite, nup);
			BasisOneSpin   basisDown(nsite, ndown);
			const SizeType sizeUp   = basisUp.size();
			const SizeType sizeDown = basisDown.size();
			for (SizeType i = 0; i < sizeUp; ++i) {
				for (SizeType j = 0; j < sizeDown; ++j) {
					assert(k < data_.size());
					data_[k++] = PairWordType(basisUp[i], basisDown[j]);
				}
			}
		}

		assert(k == hilbert);
	}

	PairIntType parts() const
	{
		return PairIntType(0, 0); // dummy, we don't care
	}

	static const WordType& bitmask(SizeType i) { return BasisType::bitmask(i); }

	SizeType size() const { return data_.size(); }

	//! Spin up and spin down
	SizeType dofs() const { throw PsimagLite::RuntimeError("dofs() unimplemented\n"); }

	SizeType hilbertOneSite(SizeType) const { return 4; }

	SizeType perfectIndex(const typename PsimagLite::Vector<WordType>::Type& kets) const
	{
		throw PsimagLite::RuntimeError("dofs() perfectIndex\n");
	}

	SizeType perfectIndex(WordType ket1, WordType ket2) const
	{
		const PairWordType p(ket1, ket2);
		auto               it = std::find(data_.begin(), data_.end(), p);
		if (it != data_.end())
			return it - data_.begin();
		throw PsimagLite::RuntimeError("perfectIndex(): Index not found for state\n");
	}

	SizeType perfectIndex(WordType newKet, SizeType ispace, SizeType spinOfNew) const
	{
		throw PsimagLite::RuntimeError("dofs() perfectIndex\n");
	}

	WordType operator()(SizeType i, SizeType spin) const
	{
		assert(i < data_.size());
		return (spin == LanczosGlobals::SPIN_UP) ? data_[i].first : data_[i].second;
	}

	PairIntType getBraIndex(WordType ket1,
	                        WordType ket2,
	                        const LabeledOperatorType&,
	                        SizeType site,
	                        SizeType spin,
	                        SizeType orb) const
	{
		throw PsimagLite::RuntimeError("getBraIndex()\n");
	}

	int doSign(WordType ket1,
	           WordType ket2,
	           SizeType i,
	           SizeType orb1,
	           SizeType j,
	           SizeType orb2,
	           SizeType spin) const
	{
		assert(i <= j);
		return (spin == SPIN_UP) ? BasisOneSpin::doSign(ket1, i, j)
		                         : BasisOneSpin::doSign(ket2, i, j);
	}

	int doSignGf(WordType a, WordType b, SizeType ind, SizeType spin, SizeType orb) const
	{
		throw PsimagLite::RuntimeError("doSignGf()\n");
	}

	int doSignSpSm(WordType a, WordType b, SizeType ind, SizeType spin, SizeType orb) const
	{
		return 1;
	}

	SizeType isThereAnElectronAt(WordType ket1,
	                             WordType ket2,
	                             SizeType site,
	                             SizeType spin,
	                             SizeType) const
	{
		return (spin == LanczosGlobals::SPIN_UP)
		    ? BasisOneSpin::isThereAnElectronAt(ket1, site)
		    : BasisOneSpin::isThereAnElectronAt(ket2, site);
	}

	SizeType orbsPerSite(SizeType i) const { return 1; }

	SizeType orbs() const { return 1; }

	SizeType getN(WordType ket1, WordType ket2, SizeType site, SizeType spin, SizeType) const
	{
		return (spin == LanczosGlobals::SPIN_UP) ? BasisOneSpin::getN(ket1, site)
		                                         : BasisOneSpin::getN(ket2, site);
	}

	bool
	getBra(WordType&, WordType, WordType, const LabeledOperatorType&, SizeType, SizeType) const
	{
		throw PsimagLite::RuntimeError("getBra()\n");
	}

	void print(std::ostream& os, typename BaseType::PrintEnum binaryOrDecimal) const
	{
		const SizeType n = data_.size();
		for (SizeType i = 0; i < n; ++i)
			os << data_[i].first << " " << data_[i].second << "\n";
	}

private:

	SizeType           ne_;
	VectorPairWordType data_;

}; // class BasisRashbaSOC

} // namespace LanczosPlusPlus
#endif
