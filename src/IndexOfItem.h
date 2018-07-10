#ifndef INDEXOFITEM_H
#define INDEXOFITEM_H
#include "Vector.h"

namespace PsimagLite {

static SizeType indexOfItemOrMinusOne(const Vector<SizeType>::Type& v, SizeType x)
{
	SizeType n = v.size();
	for (SizeType i = 0; i < n; ++i)
		if (v[i] == x) return i;

	return -1;
}

static SizeType indexOfItem(const Vector<SizeType>::Type& v, SizeType x)
{
	int y = indexOfItemOrMinusOne(v, x);
	if (y >= 0) return y;

	throw PsimagLite::RuntimeError("indexOfItem(): item not found " + ttos(x) + "\n");
}
}
#endif // INDEXOFITEM_H
