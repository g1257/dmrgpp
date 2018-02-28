#include <iostream>
#include <fstream>
#include <complex>
#include "Vector.h"
#include "Matrix.h"
#include "RungeKutta.h"
#include "IoSimple.h"

// authored by K.A.A

template<typename RealType,typename VectorType,typename MatrixType>
class GammaCiCj {

	typedef typename MatrixType::value_type ComplexOrRealType;

public:

	GammaCiCj(MatrixType& T,
		VectorType& V,
		VectorType& W,
		const RealType& omega)
		: T_(T),
		  V_(V),
		  W_(W),
		  mV_(V.size(), V.size()),
		  mW_(W.size(), W.size()),
		  omega_(omega)
	{
		assert(V_.size() == W_.size());

		for (SizeType i = 0; i < mV_.n_row(); i++) {
			for (SizeType j = 0; j < mV_.n_col(); j++) {
				mV_(i,j) = V_[j] - V_[i];
				mW_(i,j) = W_[j] - W_[i];
			}
		}
	}

	MatrixType operator()(const RealType& t,
			      const MatrixType& y) const
	{
		ComplexOrRealType c(0., -1.);
		MatrixType tmp = (-1.0)*(T_* y);
		//MatrixType yTranspose;
		//transposeConjugate(yTranspose,y);
		tmp += y * T_;
		for (SizeType i=0;i<tmp.n_row();i++)
			for (SizeType j=0;j<tmp.n_col();j++)
				tmp(i,j) += mV_(i,j)*y(i,j) + mW_(i,j)*cos(omega_*t)*y(i,j);
		return c*tmp;
	}

private:

    const MatrixType& T_;
    const VectorType& V_;
    const VectorType& W_;
    MatrixType mV_;
    MatrixType mW_;
    RealType omega_;
}; // class GammaCiCj

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType
#endif
typedef std::complex<RealType> ComplexOrRealType;
//typedef RealType ComplexOrRealType;
typedef PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
typedef GammaCiCj<RealType,VectorType,MatrixType> GammaCiCjType;
typedef PsimagLite::IoSimple::In IoInType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -f file  -i file2";
	std::cerr<<" -b t1 -e te -s ts \n";
}

int main(int argc, char* argv[])
{
	int opt;
	PsimagLite::String file = "";
	RealType wbegin=0.;
	RealType wend=10.;
	RealType wstep=0.01;
	PsimagLite::String file2 = "";

	while ((opt = getopt(argc, argv,"f:i:b:e:s:")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'b':
			wbegin = atof(optarg);
			break;
		case 'e':
			wend = atof(optarg);
			break;
		case 's':
			wstep = atof(optarg);
			break;
		case 'i':
			file2 = optarg;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	if (file=="" || file2=="" || wend<=wbegin || wstep<0) {
		usage(argv[0]);
		return 1;
	}

	IoInType io(file);

	SizeType N;
	io.readline(N,"TotalNumberOfSites=");

	RealType hopping = 1.0;
	MatrixType T(N, N);
	for(SizeType i = 0; i < N-1; i++)
		T(i,i+1) = T(i+1,i) = hopping;

	VectorType V;
	io.read(V,"potentialV");
	if (V.size()>N) V.resize(N);

	VectorType W;
	io.read(W,"PotentialT");
	assert(W.size()==N);

	RealType omega;
	io.readline(omega,"omega=");

	GammaCiCjType f(T, V, W, omega);
	PsimagLite::RungeKutta<RealType,GammaCiCjType,MatrixType> rk(f, wstep);

	MatrixType y0;
	IoInType io2(file2);
	io2.read(y0, "MatrixCiCj");
	for (SizeType i=0;i<y0.n_row();i++) {
		for (SizeType j=0;j<y0.n_col();j++) {
			//if (i==j) y0(i,j) = 1.-y0(i,j);
			//if (i!=j) y0(i,j) = -y0(i,j);
		}
	}

	PsimagLite::Vector<VectorType>::Type result;
	rk.solveEx(result,wbegin,wend, y0);
	for (SizeType i=0;i<result.size();i++) {
		RealType time = wbegin + wstep*i;
		std::cout<<time<<" ";
		for (SizeType j=0;j<result[i].size();j++)
			std::cout<<PsimagLite::real(result[i][j])<<" ";
		std::cout<<"\n";
	}
}
