#ifndef AINURDOUBLEORFLOAT_H
#define AINURDOUBLEORFLOAT_H

namespace PsimagLite {

#ifdef USE_FLOAT
#define DoubleOrFloatType float
#define DoubleOrFloatUnderscore float_
#else
#define DoubleOrFloatType double
#define DoubleOrFloatUnderscore double_
#endif

}
#endif // AINURDOUBLEORFLOAT_H
