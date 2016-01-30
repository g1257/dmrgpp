#ifndef PSI_PSIMAGLITE_H
#define PSI_PSIMAGLITE_H

#include <iostream>
#include <utility>
#include "Vector.h"

namespace PsimagLite {

std::ostream& operator<<(std::ostream&,const std::pair<SizeType,SizeType>&);

std::istream& operator>>(std::istream&,std::pair<SizeType,SizeType>&);

} // namespace PsimagLite

#endif // PSI_PSIMAGLITE_H

