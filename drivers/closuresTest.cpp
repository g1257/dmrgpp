#include "Matrix.h"
#include "Vector.h"

void testVector()
{
	int n = 4;
	std::vector<double> v1(n,1.0);
	std::vector<double> v2(n,1.1);
	std::vector<double> v3;

	// multiplication by scalar
	v3 <= v2*1.2;
	std::cout<<v3;

	// plus
	v3 <= v1 + v2;
	std::cout<<v3;
	v3 <= v1 + 1.3*v2;
	std::cout<<v3;
	v3 <= 1.3*v2 + v1;
	std::cout<<v3;

//	// minus
//	v3 <= v1 - v2;
//	std::cout<<v3;
//	v3 <= v1 - 0.5*v2;
//	std::cout<<v3;
//	v3 <= 0.5*v2 - v1;
//	std::cout<<v3;
}

void testMatrix()
{
	SizeType n = 4;
	PsimagLite::Matrix<double> m1(n,n);
	for (SizeType i = 0; i < n; ++i) {
		m1(i,i) = i;
		SizeType iPlus1 = i + 1;
		if (iPlus1 >= n) iPlus1 = 0;
		m1(i,iPlus1) = m1(iPlus1,i) = 0.1;
	}

	PsimagLite::Matrix<double> m2;
	// multiplication by scalar
	m2 = m1*1.2;
	std::cout<<m2;

	// plus
	PsimagLite::Matrix<double> m3;
	m3 = m1 + m2;
	std::cout<<m3;
	m3 = m1 + 1.3*m2;
	std::cout<<m3;
	m3 = 1.3*m2 + m1;
	std::cout<<m3;

	// minus tests omitted here

	// matrix * matrix
	m3 = m1 * m2;
	std::cout<<m3;

	// matrix * vector
	std::vector<double> v1(n,3.0);
	std::vector<double> v2;
	v2 <= m3 * v1;
	std::cout<<v2;
	v2 <= v1 * m3;
	std::cout<<v2;
}

int main(int, char**)
{
	testVector();
	testMatrix();
}

