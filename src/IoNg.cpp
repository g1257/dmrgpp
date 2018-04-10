#ifdef USE_IO_NG

#include "IoNg.h"

namespace PsimagLite {

template<>
const H5::PredType ToH5<double>::type = H5::PredType::NATIVE_DOUBLE;

template<>
const H5T_class_t ToH5<double>::super = H5T_FLOAT;

template<>
const H5::PredType ToH5<unsigned int>::type = H5::PredType::NATIVE_UINT8;

template<>
const H5T_class_t ToH5<unsigned int>::super = H5T_INTEGER;

template<>
const H5::PredType ToH5<bool>::type = H5::PredType::NATIVE_HBOOL;

} // namespace PsimagLite
#endif

