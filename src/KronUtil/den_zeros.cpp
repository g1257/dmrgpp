#include "util.h"

template<typename ComplexOrRealType>
void den_zeros(const int,
               const int,
               PsimagLite::Matrix<ComplexOrRealType>& a_)
{
	/*
 * ----------------------
 * set array to all zeros
 * ----------------------
 */
	a_.setTo(0.0);
}
#undef A
