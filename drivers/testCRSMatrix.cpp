//#include "CrsMatrix.h"
#include "CrsMatrix.h"
#include <cstdlib>
#include <fstream>
using namespace PsimagLite;
typedef double RealType;

template<typename T>
std::ostream& operator<<(std::ostream& os,const typename Vector<T>::Type& v)
{
	os<<v.size()<<"\n";
	for (SizeType i=0;i<v.size();i++)
		os<<v[i]<<" ";
	os<<"\n";
	return os;
}

void fillRandomVector(Vector<SizeType>::Type& x, SizeType maxValue)
{
	unsigned int long seed = 7334211;
	srand48(seed);
	for (SizeType i=0;i<x.size();i++)
		x[i] = drand48()*maxValue;
}

void fillRandomVector(Vector<RealType>::Type& x,RealType maxValue)
{
	unsigned int long seed = 7334211;
	srand48(seed);
	for (SizeType i=0;i<x.size();i++)
		x[i] = drand48()*maxValue;
}

template<typename T>
void testMultiply(const CrsMatrix<T>& m,RealType maxValue)
{
	typename Vector<RealType>::Type x(m.rows(),0.0),y(m.rows());
	fillRandomVector(y,maxValue);
	std::cout<<"initial vector:\n";
	std::cout<<y;
	m.matrixVectorProduct(x,y);
	std::cout<<"final vector:\n";
	std::cout<<x;
}

template<typename T>
CrsMatrix<T> createRandomCrs(SizeType rank, SizeType seed, SizeType nonZeros, T maxValue)
{
	typename Vector<SizeType>::Type rows;
	typename Vector<SizeType>::Type cols;
	typename Vector<T>::Type vals;

	srand48(seed);
	for (SizeType i = 0; i < nonZeros; ++i) {
		// pick a row
		const SizeType row = SizeType(drand48()*rank);
		// and a column
		const SizeType col = SizeType(drand48()*rank);
		// and a value
		const T val = drand48()*maxValue;
		rows.push_back(row);
		cols.push_back(col);
		vals.push_back(val);
	}

	// fill the matrix with this data:
	return CrsMatrix<T>(rank, rows, cols, vals);
}

int main(int argc,char *argv[])
{

	if (argc==3) {
		SizeType rank = std::atoi(argv[1]);
		unsigned int long seed = 343981;
		RealType ratio = std::atof(argv[2]);
		SizeType nonZeros = SizeType(ratio * rank *rank);
		RealType maxValue = 10.0;
		CrsMatrix<RealType> m = createRandomCrs(rank, seed, nonZeros, maxValue);
		std::cout<<m;

		testMultiply(m,maxValue);
	} else if (argc==2) {
		std::ifstream fin(argv[1]);
		Matrix<RealType> mdense(fin);
		fin.close();
		std::cout<<mdense;

		CrsMatrix<RealType> m(mdense);
		RealType maxValue = 10.0;
		testMultiply(m, maxValue);
		std::cout<<m;
		std::cout<<"----------\n";
		std::cout<<m.toDense();
	} else {
		throw RuntimeError("Wrong number of arguments\n");
	}


}

