// sort algorithm example
#ifndef SORT_H_H
#define SORT_H_H
#include <iostream>
#include <algorithm>
#include <vector>

template<typename ContainerType>
class Sort {
	public:

	typedef typename ContainerType::value_type FieldType;
	typedef std::pair<FieldType,size_t> PairType;
	class Compare {

		public:
			Compare(const std::vector<PairType>& x) : x_(x)
			{}

			bool operator()(const PairType& x1,const PairType& x2)
			{
				if (x1.first<x2.first) return true;
				return false;
			}
		private:
			const std::vector<PairType>& x_;
	};

	
	void sort(ContainerType& x,std::vector<size_t>& iperm)
	{
		if (x.size()==0) return;
		PairType onep(x[0],0);
		std::vector<PairType> p(x.size(),onep);
		for (size_t i=0;i<x.size();i++) {
			p[i].first = x[i];
			p[i].second = i;
		}
		std::sort(p.begin(),p.end(),Compare(p));
		for (size_t i=0;i<x.size();i++) {
			x[i] = p[i].first;
			iperm[i] = p[i].second;
		}
	}

}; //class Sort

#endif // SORT_H_H
