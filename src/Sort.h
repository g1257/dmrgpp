// sort algorithm example
#ifndef SORT_H_H
#define SORT_H_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>

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
				if (isInteger(x1.first) && isInteger(x2.first))
					return operatorParensInt(x1,x2);
				
				FieldType epsilon = static_cast<FieldType>(1e-99);
				return (x1.first-x2.first<epsilon); //std::numeric_limits<double>::epsilon());
			}
		private:
			const std::vector<PairType>& x_;

			bool operatorParensInt(const PairType& x1,const PairType& x2)
			{
				return (x1.first<x2.first);
			}

			bool isInteger(const double& t)
			{
				return false;
			}

			bool isInteger(const int& t)
			{
				return true;
			}

			bool isInteger(const size_t& t)
			{
				return true;
			} 

	};

	void sort(ContainerType& x,std::vector<size_t>& iperm,size_t smallSize=0)
	{
		size_t n = x.size();
		if (n==0) return;
		// FIXME: DON'T USE smallSize, just say n=iperm.size()
		if (smallSize!=0) n = smallSize;
		PairType onep(x[0],0);
		std::vector<PairType> p(n,onep);
		for (size_t i=0;i<n;i++) {
			p[i].first = x[i];
			p[i].second = i;
		}
		std::sort(p.begin(),p.end(),Compare(p));
		for (size_t i=0;i<n;i++) {
			x[i] = p[i].first;
			iperm[i] = p[i].second;
		}
	}

}; //class Sort

#endif // SORT_H_H
