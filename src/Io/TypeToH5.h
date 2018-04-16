#ifndef PSI_TYPETOH5_H
#define PSI_TYPETOH5_H
namespace PsimagLite {

template<typename T>
struct ToH5 {
	static const H5::PredType type;
	static const H5T_class_t super;
};
} // PsimagLite namespace
#endif // PSI_TYPETOH5_H
