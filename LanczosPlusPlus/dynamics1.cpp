#include "LabeledOperator.h"
#include "lanczos.cpp"

template <typename EngineType, bool b> class FindVectorTtoGs {

public:

	typedef typename EngineType::ComplexOrRealType ComplexOrRealType;
	typedef typename EngineType::ModelType         ModelType;

	void findVectorTtoGs(std::vector<ComplexOrRealType>& zvector,
	                     const EngineType&               engine,
	                     const ModelType&                model,
	                     SizeType                        mForK)
	{
		throw PsimagLite::RuntimeError("Only works with useComplex\n");
	}
};

template <typename EngineType> class FindVectorTtoGs<EngineType, true> {

public:

	typedef typename EngineType::ComplexOrRealType ComplexOrRealType;
	typedef typename EngineType::ModelType         ModelType;

	void findVectorTtoGs(std::vector<ComplexOrRealType>& zvector,
	                     const EngineType&               engine,
	                     const ModelType&                model,
	                     SizeType                        mForK)
	{
		const SizeType                        sites    = model.geometry().numberOfSites();
		const std::vector<ComplexOrRealType>& gsVector = engine.eigenvector(0);
		const SizeType                        n        = gsVector.size();
		const SizeType                        spin     = 0; // bogus
		const SizeType                        orb      = 0; // bogus
		zvector.resize(n);
		for (SizeType site = 0; site < sites; ++site) {
			RealType          arg = 2.0 * M_PI * mForK / sites;
			ComplexOrRealType factor(cos(arg), sin(arg));
			engine.accModifiedState_(
			    zvector,
			    LabeledOperator::Label::OPERATOR_CDAGGER_A_UP_C_B_UP,
			    model.basis(),
			    gsVector,
			    model.basis(),
			    site,
			    spin,
			    orb,
			    factor);
		}
	}
};

template <typename ModelType,
          typename SpecialSymmetryType,
          template <typename, typename> class InternalProductTemplate>
void mainLoop3(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	typedef typename ModelType::GeometryType         GeometryType;
	typedef typename GeometryType::ComplexOrRealType ComplexOrRealType;
	typedef LanczosPlusPlus::Engine<ModelType, InternalProductTemplate, SpecialSymmetryType>
	                                                             EngineType;
	typedef typename EngineType::TridiagonalMatrixType           TridiagonalMatrixType;
	typedef PsimagLite::ContinuedFraction<TridiagonalMatrixType> ContinuedFractionType;

	EngineType engine(model, io);

	//! get the g.s.:
	RealType Eg = engine.energies(0);
	std::cout.precision(8);
	std::cout << "Energy=" << Eg << "\n";

	std::vector<ComplexOrRealType>                                                    tTogs;
	FindVectorTtoGs<EngineType, PsimagLite::IsComplexNumber<ComplexOrRealType>::True> fv;
	fv.findVectorTtoGs(tTogs, engine, model, lanczosOptions.split);

	SpecialSymmetryType symm(model.basis(), model.geometry(), "");
	InternalProductTemplate<ModelType, SpecialSymmetryType> matrix(model, model.basis(), symm);
	ContinuedFractionType                                   cf;

	if (PsimagLite::norm(tTogs) < 1e-10)
		std::cerr << "spectralFunction: modifVector==0\n";

	static const bool isFermionic = false;
	static const bool isDiagonal  = true;
	const SizeType    type        = 0;
	const SizeType    spin        = 0;
	engine.calcSpectral(cf, isFermionic, tTogs, matrix, type, spin, isDiagonal);
	PsimagLite::IoSimple::Out ioOut(std::cout);
	cf.write(ioOut, "SPECTRAL");
}

template <typename ModelType, typename SpecialSymmetryType>
void mainLoop2(const ModelType&                 model,
               InputNgType::Readable&           io,
               LanczosPlusPlus::LanczosOptions& lanczosOptions)
{
	PsimagLite::String tmp;
	io.readline(tmp, "SolverOptions=");
	bool onthefly = (tmp.find("InternalProductOnTheFly") != PsimagLite::String::npos);

	if (onthefly) {
		mainLoop3<ModelType, SpecialSymmetryType, LanczosPlusPlus::InternalProductOnTheFly>(
		    model, io, lanczosOptions);
	} else {
		mainLoop3<ModelType, SpecialSymmetryType, LanczosPlusPlus::InternalProductStored>(
		    model, io, lanczosOptions);
	}
}
