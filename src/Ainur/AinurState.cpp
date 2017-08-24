#include "AinurState.h"

namespace PsimagLite {

void AinurState::assign(String k, String v)
{
	int x = storageIndexByName(k);
	if (x < 0)
		err(errLabel(ERR_PARSE_UNDECLARED, k));

	assert(static_cast<SizeType>(x) < values_.size());
	values_[x] = v;
}
} // namespace PsimagLite
