#ifndef ONELINK_HH
#define ONELINK_HH
#include "ProgramGlobals.h"
#include <functional>

namespace Dmrg
{

template <typename ComplexOrRealType>
class OneLink
{

public:

	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using OldLambdaType = std::function<void(ComplexOrRealType&)>;
	using LambdaType = std::function<void(ComplexOrRealType&, RealType, SizeType)>;
	using VectorSizeType = std::vector<SizeType>;

	OneLink(VectorSizeType indices_,
	    VectorSizeType orbs_,
	    ProgramGlobals::FermionOrBosonEnum fermionOrBoson_,
	    PsimagLite::String mods_,
	    SizeType angularMomentum_,
	    RealType angularFactor_,
	    SizeType category_,
	    LambdaType vModifier_)
	    : indices(indices_)
	    , orbs(orbs_)
	    , fermionOrBoson(fermionOrBoson_)
	    , mods(mods_)
	    , angularMomentum(angularMomentum_)
	    , angularFactor(angularFactor_)
	    , category(category_)
	    , modifier(vModifier_)
	{
	}

	OneLink(VectorSizeType indices_,
	    VectorSizeType orbs_,
	    ProgramGlobals::FermionOrBosonEnum fermionOrBoson_,
	    PsimagLite::String mods_,
	    SizeType angularMomentum_,
	    RealType angularFactor_,
	    SizeType category_,
	    OldLambdaType vModifier_)
	    : indices(indices_)
	    , orbs(orbs_)
	    , fermionOrBoson(fermionOrBoson_)
	    , mods(mods_)
	    , angularMomentum(angularMomentum_)
	    , angularFactor(angularFactor_)
	    , category(category_)
	{
		modifier = [vModifier_](ComplexOrRealType& value, RealType, SizeType) { vModifier_(value); };
	}

	VectorSizeType indices;
	VectorSizeType orbs;
	ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
	PsimagLite::String mods;
	SizeType angularMomentum;
	RealType angularFactor;
	SizeType category;
	LambdaType modifier;
}; // OneLink
}
#endif // ONELINK_HH
