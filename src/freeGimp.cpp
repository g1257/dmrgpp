#include "Matsubaras.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "ParamsDmftSolver.h"
#include "FunctionOfFrequency.h"

template<typename ComplexType>
void computeGamma(Dmft::FunctionOfFrequency<ComplexType>& gamma,
				  const typename PsimagLite::Vector<typename PsimagLite::Real<ComplexType>::Type>::Type& bath)
{
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	const SizeType n = gamma.totalMatsubaras();
	const SizeType nBath = bath.size()/2;
	for (SizeType i = 0; i < n; ++i) {
		const RealType wn = gamma.omega(i);
		ComplexType sum = 0;
		for (SizeType alpha = 0; alpha < nBath; ++alpha)
			sum += PsimagLite::conj(bath[alpha])*bath[alpha]/
				ComplexType(-bath[alpha + nBath], wn);
		gamma(i) = sum;
	}
}

template<typename ComplexType>
void computeFreeGimp(Dmft::FunctionOfFrequency<ComplexType>& freeGimp,
					 const Dmft::FunctionOfFrequency<ComplexType>& gamma)
{
	typedef typename PsimagLite::Real<ComplexType>::Type RealType;
	const SizeType n = gamma.totalMatsubaras();
	RealType sum = 0;
	for (SizeType i = 0; i < n; ++i) {
		const RealType wn = gamma.omega(i);
		const ComplexType value(0, wn);
		freeGimp(i) = 1.0/(value - gamma(i));
		sum += PsimagLite::imag(freeGimp(i));
	}

	std::cerr<<"Sum of Gimp: "<<sum<<"\n";
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr<<"USAGE: "<<argv[0]<<" filename\n";
		return 1;
	}

	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::InputNg<Dmft::InputCheck> InputNgType;
	typedef Dmft::ParamsDmftSolver<RealType, InputNgType> ParamsDmftSolverType;
	typedef Dmft::FunctionOfFrequency<ComplexType> FunctionOfFrequencyType;

	PsimagLite::String inputfile = argv[1];
	Dmft::InputCheck inputCheck;
    InputNgType::Writeable ioWriteable(inputfile, inputCheck);
    InputNgType::Readable io(ioWriteable);

    ParamsDmftSolverType params(io);
	Dmft::Matsubaras<RealType> matsubaras(params.ficticiousBeta, params.nMatsubaras);

	PsimagLite::Vector<RealType>::Type bath;
	io.read(bath, "Bath");
	FunctionOfFrequencyType gamma(params.ficticiousBeta, params.nMatsubaras);
	computeGamma(gamma, bath);

	FunctionOfFrequencyType freeGimp(params.ficticiousBeta, params.nMatsubaras);
	computeFreeGimp(freeGimp, gamma);

	std::cout<<freeGimp;
}


