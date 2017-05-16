#include "Concurrency.h"

namespace PsimagLite {

SizeType Concurrency::mode = 0;
SizeType Concurrency::npthreads = 1;
MpiDisabled Concurrency::mpiDisabled_;
bool Concurrency::setAffinitiesDefault = false;
} // namespace PsimagLite

