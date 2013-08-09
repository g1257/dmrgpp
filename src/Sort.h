// sort algorithm example
#ifndef SORT_H_H
#define SORT_H_H
#include <iostream>
#include <algorithm>
#include "Vector.h"

namespace PsimagLite {
template<typename ContainerType>
class Sort {
public:

	typedef typename ContainerType::value_type FieldType;
	typedef std::pair<FieldType,SizeType> PairType;
	class Compare {

	public:

		Compare(const typename Vector<PairType>::Type& x) : x_(x)
		{}

		bool operator()(const PairType& x1,const PairType& x2)
		{
			if (x1.first<x2.first) return true;
			return false;
		}

	private:

		const typename Vector<PairType>::Type& x_;
	};

	template<typename A>
	void sort(ContainerType& x,
	          typename std::vector<SizeType,A>& iperm,
	          SizeType smallSize=0)
	{
		SizeType n = x.size();
		if (n==0) return;
		// FIXME: DON'T USE smallSize, just say n=iperm.size()
		if (smallSize!=0) n = smallSize;
		const ContainerType& xread = x;
		PairType onep(xread[0],0);
		typename Vector<PairType>::Type p(n,onep);
		for (SizeType i=0;i<n;i++) {
			p[i].first = xread[i];
			p[i].second = i;
		}
		std::sort(p.begin(),p.end(),Compare(p));
		for (SizeType i=0;i<n;i++) {
			x[i] = p[i].first;
			iperm[i] = p[i].second;
		}
	}

}; //class Sort
} // namespace PsimagLite

#endif // SORT_H_H

