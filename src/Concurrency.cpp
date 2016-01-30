#include "Concurrency.h"

namespace PsimagLite {

SizeType Concurrency::mode = 0;
SizeType Concurrency::npthreads = 1;
MpiDisabled Concurrency::mpiDisabled_;

} // namespace PsimagLite

