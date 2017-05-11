#include "ObserveDriver.h"

template<typename ModelHelperType>
bool Dmrg::LinkProductHeisenbergAncillaC<ModelHelperType>::hot_ = false;

template<typename ModelHelperType>
bool Dmrg::LinkProductHubbardAncillaExtended<ModelHelperType>::hot_ = false;

template<typename ModelHelperType>
bool Dmrg::LinkProductTjAncillaC2<ModelHelperType>::hot_ = false;

template<typename ModelHelperType>
bool Dmrg::LinkProdExtendedSuperHubbard1Orb<ModelHelperType>::hasSpinOrbit_ = false;

using namespace Dmrg;

template<typename GeometryType,
         typename ModelHelperType,
         typename VectorWithOffsetType>
void mainLoop(GeometryType& geometry,
              const PsimagLite::String& targetting,
              InputNgType::Readable& io,
              const ParametersDmrgSolverType& params,
              const PsimagLite::String& list)
{
	typedef ModelBase<ModelHelperType,
	        ParametersDmrgSolverType,
	        InputNgType::Readable,
	        GeometryType> ModelBaseType;

	SizeType orbitals = 1.0;
	try {
		io.readline(orbitals,"Orbitals=");
	} catch (std::exception&) {}

	ModelSelector<ModelBaseType> modelSelector(params.model);
	const ModelBaseType& model = modelSelector(params,io,geometry);

	bool moreData = true;
	const PsimagLite::String& datafile = params.filename;
	ArchiveFiles<ParametersDmrgSolverType>::unpackIfNeeded(datafile);
	IoInputType dataIo(datafile);
	bool hasTimeEvolution = (targetting != "GroundStateTargetting");

	while (moreData) {
		try {
			moreData = !observeOneFullSweep<VectorWithOffsetType,ModelBaseType>
			        (dataIo,model,list,hasTimeEvolution,orbitals);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"There's no more data\n";
			break;
		}
	}
}

template<typename GeometryType,
         template<typename> class ModelHelperTemplate,
         typename MySparseMatrix>
void mainLoop1(GeometryType& geometry,
               const PsimagLite::String& targetting,
               InputNgType::Readable& io,
               const ParametersDmrgSolverType& params,
               const PsimagLite::String& list)
{
	typedef Basis<MySparseMatrix, CvectorSizeType> BasisType;
	typedef Operators<BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef typename MySparseMatrix::value_type ComplexOrRealType;

	if (params.options.find("vectorwithoffsets")!=PsimagLite::String::npos) {
		typedef VectorWithOffsets<ComplexOrRealType> VectorWithOffsetType;
		mainLoop<GeometryType,ModelHelperType,VectorWithOffsetType>
		        (geometry,targetting,io,params, list);
	} else {
		typedef VectorWithOffset<ComplexOrRealType> VectorWithOffsetType;
		mainLoop<GeometryType,ModelHelperType,VectorWithOffsetType>
		        (geometry,targetting,io,params, list);
	}
}

template<typename MySparseMatrix>
void mainLoop0(InputNgType::Readable& io,
               ParametersDmrgSolverType& dmrgSolverParams,
               InputCheck& inputCheck,
               const PsimagLite::String& list)
{
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef PsimagLite::Geometry<ComplexOrRealType,
	        InputNgType::Readable,
	        ProgramGlobals> GeometryType;

	GeometryType geometry(io);
	int tmp = 0;
	try {
		io.readline(tmp,"UseSu2Symmetry=");
	} catch (std::exception&) {}

	bool su2 = (tmp > 0);

	PsimagLite::String targetting=inputCheck.getTargeting(dmrgSolverParams.options);

	if (targetting!="GroundStateTargetting" && su2)
		throw PsimagLite::RuntimeError("SU(2) supports only GroundStateTargetting");

	if (su2) {
		mainLoop1<GeometryType,ModelHelperSu2,MySparseMatrix>
		        (geometry,targetting,io,dmrgSolverParams, list);

		return;
	}

	if (targetting=="GroundStateTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,MySparseMatrix>
		        (geometry,targetting,io,dmrgSolverParams, list);
	} else if (targetting=="TimeStepTargetting") {
		mainLoop1<GeometryType,ModelHelperLocal,MySparseMatrixComplex>
		        (geometry,targetting,io,dmrgSolverParams, list);
	} else {
		mainLoop1<GeometryType,ModelHelperLocal,MySparseMatrix>
		        (geometry,targetting,io,dmrgSolverParams, list);
	}
}

void usage(const char* name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [-p precision] [-F fileoption]";
	std::cerr<<" [-V] whatToMeasure\n";
}

int main(int argc,char **argv)
{
	using namespace Dmrg;
	PsimagLite::PsiApp application("observe",&argc,&argv,1);
	PsimagLite::String filename;
	PsimagLite::String filesOption;
	int opt = 0;
	int precision = 6;
	bool versionOnly = false;
	while ((opt = getopt(argc, argv,"f:o:p:F:V")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'o':
			std::cerr<<argv[0]<<": Omit the \"-o\". It's not needed anymore.\n";
			std::cerr<<"\t Write the insitu measurements at the end of the command line\n";
			return 1;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'F':
			filesOption = optarg;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	PsimagLite::String list = (optind < argc) ? argv[optind] : "";

	//sanity checks here
	if (filename=="" || (filesOption != "keep" && filesOption != "")) {
		if (!versionOnly) {
			usage(argv[0]);
			return 1;
		}
	}

	typedef PsimagLite::Concurrency ConcurrencyType;

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;

	//Setup the Geometry
	InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(filename,
	                                   inputCheck,
	                                   "#InputStartsHere",
	                                   "#InputEndsHere");
	InputNgType::Readable io(ioWriteable);

	//! Read the parameters for this run
	//ParametersModelType mp(io);
	ParametersDmrgSolverType dmrgSolverParams(io,false,true);

	ConcurrencyType::npthreads = dmrgSolverParams.nthreads;

	if (dmrgSolverParams.options.find("useComplex") != PsimagLite::String::npos) {
		mainLoop0<MySparseMatrixComplex>(io,dmrgSolverParams,inputCheck, list);
	} else {
		mainLoop0<MySparseMatrixReal>(io,dmrgSolverParams,inputCheck, list);
	}

	if (filesOption != "keep")
		ArchiveFiles<ParametersDmrgSolverType>::staticDelete();
} // main

