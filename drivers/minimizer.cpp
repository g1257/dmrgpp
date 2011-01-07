
#include "Minimizer.h"

template<typename T>
std::ostream& operator<<(std::ostream& os,const std::vector<T>& v)
{
	os<<v.size()<<"\n";
	for (size_t i=0;i<v.size();i++)
		os<<v[i]<<" ";
	os<<"\n";
	return os;
}

using namespace PsimagLite;
typedef  double RealType;

class MyFunctionTest
{
public:
	typedef double FieldType;

	template<typename SomeVectorType>
	FieldType operator()(const SomeVectorType &v) const
	{
		return (v[0]-2)*(v[0]-2);
	}
	size_t size() const { return  1; }
};


int main(int argc,char *argv[])
{
	size_t n=1;
	std::vector<RealType> x(n);
	
	// inital guess:
	for (size_t i=0;i<n;i++)
		x[i] = 1.0;
	
	size_t maxIter = 100;
	MyFunctionTest f;
	Minimizer<RealType,MyFunctionTest> min(f,maxIter);

	int iter = min.simplex(x);
	if (iter<0) {
		std::cout<<"No minimum found\n";
		return 1;
	}
	std::cout<<"Minimum found after "<<iter<<" iterations.\n";
	std::cout<<x;
}
