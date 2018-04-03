#ifndef IOSERIALIZERSTUB_H
#define IOSERIALIZERSTUB_H

#ifdef USE_IO_NG
#include "IoNgSerializer.h"
namespace PsimagLite {
typedef IoNgSerializer IoSerializer;
}
#else
namespace PsimagLite {
class IoSerializer {

};
}
#endif
#endif // IOSERIALIZERSTUB_H
