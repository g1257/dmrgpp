#ifndef INDEXOFITEM_H
#define INDEXOFITEM_H
#include "Vector.h"

namespace PsimagLite {

template<typename T>
int indexOfItemOrMinusOne(const typename Vector<T>::Type& v, const T& x)
{
	SizeType n = v.size();
	for (SizeType i = 0; i < n; ++i)
		if (v[i] == x) return i;

	return -1;
}

template<typename T>
SizeType indexOfItem(const typename Vector<T>::Type& v, const T& x)
{
	int y = indexOfItemOrMinusOne(v, x);
	if (y >= 0) return y;

	throw RuntimeError("indexOfItem(): item not found " + typeToString(x) + "\n");
}

}
#endif // INDEXOFITEM_H
