#ifndef PSI_TYPETOH5_H
#define PSI_TYPETOH5_H
#include "H5Cpp.h"

namespace PsimagLite {

template<typename T>
struct TypeToH5 {
	static const H5::PredType type;
	static const H5T_class_t super;
};
} // PsimagLite namespace
#endif // PSI_TYPETOH5_H
