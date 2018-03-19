#include "Vector.h"
#include "Minimizer.h"
//#include "Square.h"

template<typename T>
T square(const T& t1)
{
	return t1*t1;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const typename PsimagLite::Vector<T>::Type& v)
{
	os<<v.size()<<"\n";
	for (SizeType i=0;i<v.size();i++)
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
		return square(v[0]-2)+square(v[1]-3);
	}
	SizeType size() const { return  2; }
};


int main(int argc,char *argv[])
{
	SizeType n=2;
	typename Vector<RealType>::Type x(n);
	
	// inital guess:
	for (SizeType i=0;i<n;i++)
		x[i] = drand48();
	
	SizeType maxIter = 100;
	MyFunctionTest f;
	Minimizer<RealType,MyFunctionTest> min(f,maxIter,true);

	int iter = min.simplex(x,1e-3,1e-5);
	if (iter<0) {
		std::cout<<"No minimum found\n";
		return 1;
	}
	std::cout<<"Minimum found after "<<iter<<" iterations.\n";
	std::cout<<x;
}
