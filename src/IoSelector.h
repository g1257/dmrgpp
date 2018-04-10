#ifndef IOSELECTOR_H
#define IOSELECTOR_H

#include "IoSerializerStub.h"

#ifndef USE_IO_NG
#include "IoSimple.h"
#else
#include "IoNg.h"
#endif

namespace PsimagLite {

#ifndef USE_IO_NG
typedef PsimagLite::IoSimple IoSelector;
#else
typedef PsimagLite::IoNg IoSelector;
#endif
} // namespace PsimagLite

#endif // IOSELECTOR_H
