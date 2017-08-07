#include "CrsMatrix.h"
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
	std::cout<<"--------------\n";
	std::cout<<m3;
	std::cout<<m2;
	m3 = 0.5*(m3-m2);
	std::cout<<m3;
	std::cout<<"--------------\n";
	m3 += m2;
	std::cout<<m3;
	m3 -= m2;
	std::cout<<m3;
	std::cout<<"--------------\n";

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

template<typename T>
void createRandomCrs(PsimagLite::CrsMatrix<T>& crs,
                     SizeType seed,
                     SizeType nonZeros,
                     T maxValue)
{
	srand48(seed);
	PsimagLite::Matrix<T> m(crs.rows(),crs.cols());
	for (SizeType i=0;i<nonZeros;i++) {
		// pick a row
		SizeType row = SizeType(drand48()*m.rows());
		// and a column
		SizeType col = SizeType(drand48()*m.cols());
		// and a value
		T val = drand48()*maxValue;
		m(row,col) = val;
	}

	fullMatrixToCrsMatrix(crs,m);
}

template<typename T>
bool checkCrs(const PsimagLite::CrsMatrix<T>& m,
              const PsimagLite::Matrix<T>& fm)
{
	PsimagLite::Matrix<T> mf;
	crsMatrixToFullMatrix(mf,m);
	for (SizeType i = 0; i < fm.n_row(); ++i) {
		for (SizeType j = 0; j < fm.n_col(); ++j) {
			if (fabs(mf(i,j)-fm(i,j)) > 1e-6) {
				std::cout<<mf;
				std::cout<<fm;
				return false;
			}
		}
	}

	return true;
}

void testCrsMatrix()
{
	SizeType n = 10;
	unsigned int long seed = 343981;
	double ratio = 0.5;
	SizeType nonZeros = SizeType(ratio*n*n);
	double maxValue = 10.0;
	PsimagLite::CrsMatrix<double> m1(n,n);
	createRandomCrs(m1,seed,nonZeros,maxValue);
	PsimagLite::Matrix<double> fm1;
	crsMatrixToFullMatrix(fm1,m1);
	//std::cout<<m1;
	//std::cout<<"-------------\n";

	seed += 111;
	PsimagLite::CrsMatrix<double> m2(n,n);
	createRandomCrs(m2,seed,nonZeros,maxValue);
	PsimagLite::Matrix<double> fm2;
	crsMatrixToFullMatrix(fm2,m2);
	//std::cout<<m2;
	std::cout<<"---operator+=scalar*crs----TEST 1/4\n";
	m1 += 1.2*m2;
//	std::cout<<m1;
	fm1 += 1.2*fm2;
	std::cout<<"CHECK PASSES="<<checkCrs(m1,fm1)<<"\n";
	std::cout<<"-------------\n";

	std::cout<<"---ctor(crs*crs)----TEST 2/4\n";
	PsimagLite::CrsMatrix<double> m3 = m1*m2;
	PsimagLite::Matrix<double> fm3 = fm1*fm2;
	std::cout<<"CHECK PASSES="<<checkCrs(m3,fm3)<<"\n";
//	std::cout<<m3;
//	std::cout<<"-------------\n";

	std::cout<<"---ctor(scalar*crs)----TEST 3/4\n";
	PsimagLite::CrsMatrix<double> m4 = 1.3*m2;
	PsimagLite::Matrix<double> fm4;
	fm4 = 1.3*fm2;
	std::cout<<"CHECK PASSES="<<checkCrs(m4,fm4)<<"\n";
//	std::cout<<m4;
//	std::cout<<"-------------\n";

	std::cout<<"---ctor(scalar*crs)----TEST 4/4\n";
	m4 += 1.4*m1*m3;
	PsimagLite::Matrix<double> fm5;
	fm5 = 1.4*fm1*fm3;
	fm4 += fm5;
//	std::cout<<m4;
	std::cout<<"CHECK PASSES="<<checkCrs(m4,fm4)<<"\n";
	std::cout<<"-------------\n";
}

int main(int, char**)
{
	testVector();
	testMatrix();
	testCrsMatrix();
}

