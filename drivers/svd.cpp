#include "Matrix.h"
#include "Vector.h"
#include "Svd.h"

int main()
{
	typedef PsimagLite::Matrix<double> MatrixType;

	MatrixType a(4,2);
	a(0,0)=2; a(0,1)=4;
	a(1,0)=1; a(1,1)=3;

	std::cout<<"A\n";
	std::cout<<a;

	PsimagLite::Vector<double>::Type s;
	MatrixType vt;

	PsimagLite::Svd<double> svd;
	svd('A', a, s, vt);

	std::cout<<"U\n";
	std::cout<<a;
	std::cout<<"S\n";
	std::cout<<s;
}

