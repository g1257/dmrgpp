#include "Concurrency.h"

namespace PsimagLite {

SizeType Concurrency::mode = 0;
MpiDisabled Concurrency::mpiDisabled_;
CodeSectionParams Concurrency::codeSectionParams(1);
} // namespace PsimagLite

