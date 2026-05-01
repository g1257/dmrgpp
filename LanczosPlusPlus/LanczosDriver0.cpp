#include "LanczosDriver1.h"

typedef PsimagLite::Geometry<RealType, InputNgType::Readable, LanczosPlusPlus::LanczosGlobals>
    Geometry0Type;

typedef LanczosPlusPlus::ModelSelector<RealType, Geometry0Type, InputNgType::Readable>
                                          ModelSelector0Type;
typedef ModelSelector0Type::ModelBaseType ModelBase0Type;

typedef ModelBase0Type::BasisBaseType BasisBase0Type;

typedef LanczosPlusPlus::TranslationSymmetry<Geometry0Type, BasisBase0Type> Symmetry0Type;

template void
mainLoop2<ModelBase0Type, Symmetry0Type>(const ModelBase0Type&            model,
                                         InputNgType::Readable&           io,
                                         LanczosPlusPlus::LanczosOptions& lanczosOptions);

namespace LanczosPlusPlus {

bool operator==(const LanczosPlusPlus::ReflectionItem& item1,
                const LanczosPlusPlus::ReflectionItem& item2)
{
	if (item1.type != item2.type)
		return false;

	if (item1.i == item2.j && item1.j == item2.i)
		return true;

	return (item1.i == item2.i && item1.j == item2.j);
}
}

std::ostream& operator<<(std::ostream& os, const LanczosPlusPlus::BasisOneSpinFeAs& b)
{
	for (SizeType i = 0; i < b.size(); i++)
		os << i << " " << b[i] << "\n";
	return os;
}

std::ostream& operator<<(std::ostream& os, const LanczosPlusPlus::BasisOneSpinImmm& b)
{
	for (SizeType i = 0; i < b.size(); i++)
		os << i << " " << b[i] << "\n";
	return os;
}
