#ifdef USE_IO_NG

#include "TypeToH5.h"

namespace PsimagLite {

template<>
 H5::PredType typeToH5<double>() { return H5::PredType::NATIVE_DOUBLE; }

template<>
 H5::PredType typeToH5<unsigned int>() { return H5::PredType::NATIVE_UINT8; }

template<>
 H5::PredType typeToH5<int>() { return H5::PredType::NATIVE_INT8; }

template<>
 H5::PredType typeToH5<unsigned long>() { return H5::PredType::NATIVE_UINT16; }

template<>
 H5::PredType typeToH5<long>() { return H5::PredType::NATIVE_INT16; }

template<>
 H5::PredType typeToH5<bool>() { return H5::PredType::NATIVE_HBOOL; }

template<>
H5::PredType typeToH5<char>() { return H5::PredType::NATIVE_CHAR; }
} // namespace PsimagLite
#endif

