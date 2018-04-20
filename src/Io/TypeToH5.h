#ifndef PSI_TYPETOH5_H
#define PSI_TYPETOH5_H
#include "H5Cpp.h"
#include "../AllocatorCpu.h"

namespace PsimagLite {

template<typename T>
typename EnableIf<!IsEnum<T>::True, H5::PredType>::Type typeToH5();

template<typename T>
typename EnableIf<IsEnum<T>::True, H5::PredType>::Type typeToH5()
{
	return H5::PredType::NATIVE_UINT8;
}
} // PsimagLite namespace
#endif // PSI_TYPETOH5_H
