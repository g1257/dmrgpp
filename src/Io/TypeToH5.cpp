#ifdef USE_IO_NG

#include "TypeToH5.h"

namespace PsimagLite {

template<>
H5::PredType typeToH5<char>() { return H5::PredType::NATIVE_CHAR; }

template<>
H5::PredType typeToH5<bool>() { return H5::PredType::NATIVE_HBOOL; }

template<>
H5::PredType typeToH5<int>() { return H5::PredType::NATIVE_INT; }

template<>
H5::PredType typeToH5<short>() { return H5::PredType::NATIVE_SHORT; }

template<>
H5::PredType typeToH5<long>() { return H5::PredType::NATIVE_LONG; }

template<>
H5::PredType typeToH5<unsigned int>() { return H5::PredType::NATIVE_UINT; }

template<>
H5::PredType typeToH5<unsigned short>() { return H5::PredType::NATIVE_USHORT; }

template<>
H5::PredType typeToH5<unsigned long>() { return H5::PredType::NATIVE_ULONG; }

template<>
H5::PredType typeToH5<float>() { return H5::PredType::NATIVE_FLOAT; }

template<>
H5::PredType typeToH5<double>() { return H5::PredType::NATIVE_DOUBLE; }




} // namespace PsimagLite
#endif

