#ifndef IOSERIALIZERSTUB_H
#define IOSERIALIZERSTUB_H

#ifndef USE_IO_SIMPLE

#include "IoNgSerializer.h"
namespace PsimagLite {

typedef IoNgSerializer IoSerializer;

}

#else

#include "IoSerializerEmpty.h"

namespace  PsimagLite {

typedef IoSerializerEmpty IoSerializer;

}

#endif
#endif // IOSERIALIZERSTUB_H
