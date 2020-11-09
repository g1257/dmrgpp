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
	typedef typename ModelType::VectorRealType VectorRealType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;

	OneSiteTruncation(const LeftRightSuperType& lrs, const ModelType& model)
	    : lrs_(lrs), model_(model) {}

	void update(SizeType oneSiteTruncSize,
	            const VectorWithOffsetType& psi,
	            ProgramGlobals::DirectionEnum dir)
	{
		if (oneSiteTruncSize == 0) return;

		// compute U ...
		MatrixType U(oneSiteTruncSize, oneSiteTruncSize);
		computeU(U, oneSiteTruncSize, psi, dir);

		// ... and send it to model
		model_.oneSiteTruncationUpdate(U);
	}

private:

	void computeU(MatrixType& U,
	              SizeType oneSiteTruncSize,
	              const VectorWithOffsetType& psi,
	              ProgramGlobals::DirectionEnum dir)
	{
		const SizeType sectors = psi.sectors();
		const SizeType nsysOrEnv = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		        ? lrs_.left().size()/oneSiteTruncSize : oneSiteTruncSize;
		PackIndicesType pack(nsysOrEnv);
		PackIndicesType packSuper(lrs_.left().size());
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
					for (SizeType x1prime = 0; x1prime < oneSiteTruncSize; ++x1prime) {
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
					for (SizeType y0prime = 0; y0prime < oneSiteTruncSize; ++y0prime) {
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

		VectorRealType eigs(U.rows());
		diag(U, eigs, 'V');
		truncateU(U, eigs);
	}

	void truncateU(MatrixType& U, const VectorRealType& eigs) const
	{
		std::cerr<<"truncateU: unimplemented, sorry, todo, fixme\n";
	}

	const LeftRightSuperType& lrs_;
	const ModelType& model_;
};
}
#endif // ONESITETRUNCATION_H
