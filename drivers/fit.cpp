#include "Vector.h"
#include <fstream>
#include <cstdlib>
#include "Minimizer.h"
#include <cassert>
#include "PsimagLite.h"

template<typename RealType_>
class OracleData {

public:

	typedef RealType_ RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	OracleData(PsimagLite::String file, RealType kf)
	    : file_(file)
	{
		// kx omega real imag
		std::ifstream fin(file_.c_str());

		while (!fin.eof()) {
			RealType tmp = 0.0;
			fin >> tmp;
			RealType kx = tmp;
			fin >> tmp;
			RealType omega = tmp;
			fin >> tmp;
			fin >> tmp;
			RealType value = tmp;
			if (fabs(kx - kf) >= 1e-6)
				continue;
			omegas_.push_back(omega);
			values_.push_back(value);
		}

		std::cerr<<"#Found "<<omegas_.size()<<" omega values for kf="<<kf<<"\n";
		if (omegas_.size() == 0)
			err("No data found in " + file + "\n");
	}

	const RealType& operator()(SizeType i) const
	{
		assert(i < values_.size());
		return values_[i];
	}

	const VectorRealType& omegas() const
	{
		return omegas_;
	}

private:

	PsimagLite::String file_;
	VectorRealType omegas_;
	VectorRealType values_;
};

template<typename RealType>
class FitData {

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	FitData(const VectorRealType& omegas, RealType mu, RealType kf, int ky)
	    : omegas_(omegas),
	      ekf_(dispersion(kf, ky*M_PI) - mu),
	      initDelta_(1),
	      initGamma_(0.1),
	      anorm_(1.0)
	{
		//anorm_ = 1.0/sum();
		std::cerr<<"#FitData ctor(): ekf_= "<<ekf_;
		std::cerr<<" initDelta_= "<<initDelta_<<" initGamma_= "<<initGamma_;
		std::cerr<<" mu= "<<mu<<"\n";
		std::cout<<"anorm="<<anorm_<<"\n";
	}

	RealType operator()(SizeType i, const VectorRealType& v) const
	{
		assert(v.size() == 2);
		const RealType& omega = omegas_[i];
		const RealType& delta = v[0];
		const RealType& gamma = v[1];
		return finternal(omega, delta, gamma);
	}

	RealType df(SizeType i, const VectorRealType& v, SizeType j) const
	{
		assert(v.size() == 2);
		const RealType& omega = omegas_[i];
		const RealType& delta = v[0];
		const RealType& gamma = v[1];
		return (j == 0) ? dfDelta(omega, delta, gamma) : dfGamma(omega, delta, gamma);
	}

	SizeType size() const { return 2; }

	static void init(VectorRealType& x)
	{
		assert(x.size() == 2);
		x[0] = 1;
		x[1] = 0.1;
	}

private:

	RealType finternal(RealType omega, RealType delta, RealType gamma) const
	{
		RealType gaom = gamma*omega;
		RealType num = (omega + ekf_)*gaom*2.0*anorm_/M_PI;
		RealType omega2 = omega*omega;
		RealType phi2 = ekf_*ekf_ + gamma*gamma + delta*delta;
		RealType den = square(omega2 - phi2) + 4*gaom*gaom;
		return num/den;
	}

	RealType dfDelta(RealType omega, RealType delta, RealType gamma) const
	{
		err("fixme\n");
		return 0;
	}

	RealType dfGamma(RealType omega, RealType delta, RealType gamma) const
	{
		err("fixme\n");
		return 0;
	}

	RealType dispersion(RealType kx, RealType ky) const
	{
		return -2*cos(kx) - cos(ky);
	}

	RealType sum() const
	{
		SizeType n = omegas_.size();
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += finternal(omegas_[i], initDelta_, initGamma_);
		return sum;
	}

	static RealType square(RealType x) { return x*x; }

	VectorRealType omegas_;
	RealType ekf_;
	const RealType initDelta_;
	const RealType initGamma_;
	RealType anorm_;
};

template<typename OracleType, typename FitDataType>
class Fitter {

	typedef typename OracleType::RealType RealType;
	typedef typename OracleType::VectorRealType VectorRealType;

	class MyFunctionTest {

	public:

		typedef RealType FieldType;

		MyFunctionTest(const OracleType& od, const FitDataType& fd)
		    : od_(od), fd_(fd)
		{}

		RealType operator()(const VectorRealType &v) const
		{
			RealType sum = 0.0;
			SizeType n = od_.omegas().size();
			for (SizeType i = 0; i < n; ++i) {
				RealType x = fabs(od_(i) - fd_(i, v));
				sum += x*x;
			}

			return sum;
		}

		void df(VectorRealType& result, const VectorRealType &v) const
		{
			assert(result.size() == size());
			for (SizeType j = 0; j < size(); ++j) {
				RealType sum = 0.0;
				SizeType n = od_.omegas().size();
				for (SizeType i = 0; i < n; ++i) {
					// FIXME CHECK SIGN OF DERIVATIVE HERE
					RealType x = fabs(od_(i) - fd_(i, v))*2.0*fd_.df(i,v,j);
					sum += x;
				}

				result[j] = sum;
			}
		}

		SizeType size() const { return  2; }

	private:

		const OracleType& od_;
		const FitDataType& fd_;
	};

public:

	Fitter(const OracleType& od, const FitDataType& fd)
	    : od_(od), fd_(fd), results_(fd.size(), 0)
	{}

	void fit(SizeType maxIter)
	{
		FitDataType::init(results_);

		MyFunctionTest f(od_, fd_);
		PsimagLite::Minimizer<RealType,MyFunctionTest> min(f, maxIter);

		int iter = min.simplex(results_, 1e-5, 1e-7);
		if (iter<0)
			std::cerr<<"No minimum found\n";
	}

	void fit2(SizeType maxIter)
	{
		FitDataType::init(results_);

		MyFunctionTest f(od_, fd_);
		PsimagLite::Minimizer<RealType,MyFunctionTest> min(f, maxIter);

		int iter = min.conjugateGradient(results_, 1e-5, 1e-5, 1e-7);
		if (iter<0)
			std::cerr<<"No minimum found\n";
	}

	void print(std::ostream& os) const
	{
		if (results_.size() != 2)
			err("print results not size 2\n");

		os<<"delta="<<results_[0]<<"\n";
		os<<"gamma="<<results_[1]<<"\n";
	}

private:

	const OracleType& od_;
	const FitDataType& fd_;
	VectorRealType results_;
};

int main(int argc, char** argv)
{
	if (argc != 5) {
		std::cerr<<"USAGE: "<<argv[0]<<" filename.gnuplot mu kf ky\n";
		return 1;
	}

	double mu = atof(argv[2]);
	double kf = atof(argv[3]);
	int ky = atoi(argv[4]);
	OracleData<double> od(argv[1], kf);

	FitData<double> fit(od.omegas(), mu, kf, ky);

	Fitter<OracleData<double>, FitData<double> > fitter(od, fit);

	SizeType maxIter = 1000;
	fitter.fit(maxIter);

	fitter.print(std::cout);
}
