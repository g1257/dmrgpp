#include "Integrator.h"

template<typename RealType_>
class SquareFunction {

	struct Params {
		Params(RealType_ p_) : p(p_)
		{}

		RealType_ p;
	};

public:

	typedef RealType_ RealType;

	SquareFunction(RealType p) : p_(p)
	{}

	static RealType function(RealType x, void* vp)
	{
		Params* p = static_cast<Params*>(vp);
		return x*x*p->p;
	}

	Params& params() { return p_; }

private:

	Params p_;
};

int main(int argc, char** argv)
{
	if (argc != 4) {
		std::cerr<<"USAGE: "<<argv[0]<<" x0 total xstep\n";
		return 1;
	}

	double x0 = atof(argv[1]);
	SizeType total = atoi(argv[2]);
	double xstep = atof(argv[3]);
	typedef SquareFunction<double>  SquareFunctionType;
	SquareFunctionType squareFunction(3.0);
	PsimagLite::Integrator<SquareFunctionType> integrator(squareFunction);
	PsimagLite::Vector<double>::Type pts(2,0);
	for (SizeType i = 0; i < total; ++i) {
		pts[1] = x0 + i*xstep;
		std::cout<<pts[1]<<" "<<integrator(pts)<<"\n";
	}
}

