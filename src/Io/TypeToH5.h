#ifndef PSI_TYPETOH5_H
#define PSI_TYPETOH5_H
#include "../AllocatorCpu.h"
#include <H5Cpp.h>

namespace PsimagLite {

template <typename T>
typename EnableIf<!std::is_enum<T>::value, H5::PredType>::Type typeToH5();

template <typename T>
typename EnableIf<std::is_enum<T>::value, H5::PredType>::Type typeToH5()
{
	return H5::PredType::NATIVE_UINT8;
}
} // namespace PsimagLite
#endif // PSI_TYPETOH5_H
