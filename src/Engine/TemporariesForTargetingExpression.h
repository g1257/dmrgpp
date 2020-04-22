#ifndef TEMPORARIESFORTARGETINGEXPRESSION_H
#define TEMPORARIESFORTARGETINGEXPRESSION_H

namespace Dmrg {

template<typename AuxForTargetingExpression>
class TemporariesForTargetingExpression {

public:

	typedef AuxForTargetingExpression AuxiliaryType;
	typedef typename AuxiliaryType::VectorWithOffsetType VectorWithOffsetType;

	TemporariesForTargetingExpression(const AuxiliaryType& aux) : aux_(aux)
	{}

	const VectorWithOffsetType& getCurrentVector(PsimagLite::String braOrKet) const
	{
		PsimagLite::GetBraOrKet getBraOrKet(braOrKet);

		if (getBraOrKet.isPvector()) {
			const SizeType pIndex = getBraOrKet.pIndex();
			if (pIndex >= aux_.pvectors.size())
				err("getVector: out of range for " + braOrKet + "\n");
			return aux_.pvectors[pIndex];
		}

		const SizeType sectorIndex = getBraOrKet.sectorIndex();
		return *(aux_.psi[sectorIndex][getBraOrKet.levelIndex()]);
	}

	const AuxiliaryType& aux_;
};
}
#endif // TEMPORARIESFORTARGETINGEXPRESSION_H
