#ifndef IOSERIALIZERSTUB_H
#define IOSERIALIZERSTUB_H

#ifndef USE_IO_SIMPLE

#include "IoNgSerializer.h"
namespace PsimagLite {

using IoSerializer = IoNgSerializer;

}

#else

#include "IoSerializerEmpty.h"

namespace PsimagLite {

using IoSerializer = IoSerializerEmpty;

}

#endif
#endif // IOSERIALIZERSTUB_H
