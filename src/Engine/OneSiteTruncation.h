#ifndef ONESITETRUNCATION_H
#define ONESITETRUNCATION_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "PackIndices.h"

namespace Dmrg {

template<typename ModelType, typename VectorWithOffsetType>
class OneSiteTruncation {

public:

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename ModelType::MatrixType MatrixType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;

	OneSiteTruncation(const LeftRightSuperType& lrs, const ModelType& model)
	    : lrs_(lrs), model_(model) {}

	void update(bool OneSiteTruncActive,
	            const VectorWithOffsetType& psi,
	            ProgramGlobals::DirectionEnum dir)
	{
		if (!OneSiteTruncActive) return;

		// compute U ...
		MatrixType U;
		computeU(U, psi, dir);

		// ... and send it to model
		model_.oneSiteTruncationUpdate(U);
	}

private:

	void computeU(MatrixType& U,
	              const VectorWithOffsetType& psi,
	              ProgramGlobals::DirectionEnum dir)
	{
		const SizeType oneSiteTotal = findHilbertSize(dir);
		const SizeType sectors = psi.sectors();
		const SizeType nsysOrEnv = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		        ? lrs_.left().size() : lrs_.right().size();
		PackIndicesType pack(nsysOrEnv);
		PackIndicesType packSuper(lrs_.super().size());
		for (SizeType i = 0; i < sectors; ++i) {
			const SizeType sector = psi.sector(i);
			const SizeType total = psi.effectiveSize(sector);
			const SizeType offset = psi.offset(sector);
			for (SizeType s = 0; s < total; ++s) {
				SizeType x = 0;
				SizeType y = 0;
				packSuper.unpack(x, y, lrs_.super().permutation(s + offset));
				const ComplexOrRealType value = psi.fastAccess(sector, s);
				if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
					SizeType x0 = 0;
					SizeType x1 = 0;
					pack.unpack(x0, x1, lrs_.left().permutation(x));
					for (SizeType x1prime = 0; x1prime < oneSiteTotal; ++x1prime) {
						const SizeType xprime = pack.pack(x0,
						                                  x1prime,
						                                  lrs_.left().permutationInverse());
						const SizeType sWithOffsetprime =
						        packSuper.pack(xprime, y, lrs_.super().permutationInverse());
						if (sWithOffsetprime < offset) continue; // respect all symmetries
						const SizeType sprime = sWithOffsetprime - offset;
						if (sprime >= total) continue; // respect all symmetries
						U(x1, x1prime) += PsimagLite::conj(psi.fastAccess(sector, sprime))*value;
					}
				} else {
					SizeType y0 = 0;
					SizeType y1 = 0;
					pack.unpack(y0, y1, lrs_.right().permutation(y));
					for (SizeType y0prime = 0; y0prime < oneSiteTotal; ++y0prime) {
						const SizeType yprime = pack.pack(y0prime,
						                                  y1,
						                                  lrs_.right().permutationInverse());
						const SizeType sWithOffsetprime =
						        packSuper.pack(x, yprime, lrs_.super().permutationInverse());
						if (sWithOffsetprime < offset) continue; // respect all symmetries
						const SizeType sprime = sWithOffsetprime - offset;
						if (sprime >= total) continue; // respect all symmetries
						U(y0, y0prime) += PsimagLite::conj(psi.fastAccess(sector, sprime))*value;
					}
				}
			}
		}
	}

	SizeType findHilbertSize(ProgramGlobals::DirectionEnum dir) const
	{
		const BasisWithOperatorsType& basis = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		        ? lrs_.left() : lrs_.right();

		const SizeType ind = basis.numberOfLocalOperators();
		assert(ind > 0);

		return basis.localOperator(ind - 1).getStorage().rows();
	}

	const LeftRightSuperType& lrs_;
	const ModelType& model_;
};
}
#endif // ONESITETRUNCATION_H
