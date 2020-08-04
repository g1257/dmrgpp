#include "Concurrency.h"

namespace PsimagLite {

SizeType Concurrency::mode = 0;
LabelDisabled Concurrency::mpiDisabled_;
CodeSectionParams Concurrency::codeSectionParams(1);
} // namespace PsimagLite

