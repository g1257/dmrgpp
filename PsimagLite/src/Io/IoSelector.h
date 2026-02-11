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
using IoSelector = PsimagLite::IoSimple;
#else
using IoSelector = PsimagLite::IoNg;
#endif
} // namespace PsimagLite

#endif // IOSELECTOR_H
