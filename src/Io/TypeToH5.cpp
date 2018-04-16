#ifdef USE_IO_NG

#include "TypeToH5.h"

namespace PsimagLite {

template<>
const H5::PredType TypeToH5<double>::type = H5::PredType::NATIVE_DOUBLE;

template<>
const H5T_class_t TypeToH5<double>::super = H5T_FLOAT;

template<>
const H5::PredType TypeToH5<unsigned int>::type = H5::PredType::NATIVE_UINT8;

template<>
const H5::PredType TypeToH5<int>::type = H5::PredType::NATIVE_INT8;

template<>
const H5T_class_t TypeToH5<unsigned int>::super = H5T_INTEGER;

template<>
const H5::PredType TypeToH5<bool>::type = H5::PredType::NATIVE_HBOOL;

template<>
const H5::PredType TypeToH5<char>::type = H5::PredType::NATIVE_CHAR;
} // namespace PsimagLite
#endif

