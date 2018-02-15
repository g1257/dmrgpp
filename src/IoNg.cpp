#include "IoNg.h"

namespace PsimagLite {

const H5::PredType ToH5<double>::type = H5::PredType::NATIVE_DOUBLE;
const H5T_class_t ToH5<double>::super = H5T_FLOAT;

} // namespace PsimagLite
