#ifndef EFFECTIVEQUANTUMNUMBER_H
#define EFFECTIVEQUANTUMNUMBER_H
#include "Vector.h"

namespace Dmrg {

class EffectiveQuantumNumber {

public:

	EffectiveQuantumNumber(PsimagLite::String str)
	    : str_(str) {}

	const PsimagLite::String& name() const { return str_;}

private:

	PsimagLite::String str_;
};
}
#endif // EFFECTIVEQUANTUMNUMBER_H
