#ifndef IOSELECTOR_H
#define IOSELECTOR_H

#include "IoSerializerStub.h"

#ifdef USE_IO_SIMPLE
#include "IoSimple.h"
#else
#include "IoNg.h"
#endif

namespace PsimagLite {

#ifdef USE_IO_SIMPLE
typedef PsimagLite::IoSimple IoSelector;
#else
typedef PsimagLite::IoNg IoSelector;
#endif
} // namespace PsimagLite

#endif // IOSELECTOR_H
