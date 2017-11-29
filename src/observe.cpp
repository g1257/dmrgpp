#include "ObserveDriver.h"

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
	bool hasTimeEvolution = (targetting != "GroundStateTargetting" &&
	        targetting != "CorrectionTargetting");

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

void usage(const PsimagLite::String& name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [-p precision] [-F fileoption]";
	std::cerr<<" [-V] whatToMeasure\n";
}

/* PSIDOC ObserveDriver
 The observe driver may be used to measure after the main dmrg run has finished.
The observe driver is used as follows,
\begin{verbatim}
./observe -f input.inp whatToMeasure
\end{verbatim}
  The command line arguments of observe are the following.
\begin{itemize}
\item[-f] [Mandatory, String] Input to use. Files
 referred to by \verb!OutputFile=! are now inputs, and
must be present.
\item[whatToMeasure] {[}Mandatory, String{]} What to measure post process.
This is a comma-separated list of braket specifications.
Braket specifications can be bare or dressed, and are explained elsewhere.
\item[-p] [Optional, Integer] Digits of precision for printing.
\item[-o] {[}Optional, String{]} Extra options for SolverOptions
\item[-F] [Optional, string] TBW
\item[-V] [Optional] Print version and exit
\end{itemize}
  */
int main(int argc,char **argv)
{
	using namespace Dmrg;
	PsimagLite::PsiApp application("observe",&argc,&argv,1);
	PsimagLite::String filename;
	PsimagLite::String filesOption;
	int opt = 0;
	int precision = 6;
	bool versionOnly = false;
	PsimagLite::String sOptions("");
	while ((opt = getopt(argc, argv,"f:p:o:F:V")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'o':
			sOptions += optarg;
			break;
		case 'F':
			filesOption = optarg;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			usage(application.name());
			return 1;
		}
	}

	PsimagLite::String list = (optind < argc) ? argv[optind] : "";

	//sanity checks here
	if (filename=="" || (filesOption != "keep" && filesOption != "")) {
		if (!versionOnly) {
			usage(application.name());
			return 1;
		}
	}

	typedef PsimagLite::Concurrency ConcurrencyType;

	// print license
	if (ConcurrencyType::root()) {
		Provenance provenance;
		std::cout<<provenance;
		std::cerr<<application.name()<<"\x1b[38;5;124m";
		std::cerr<<" [features]"<<PsimagLite::AnsiColor::reset<<"\n";
	}

	if (versionOnly) return 0;

	//Setup the Geometry
	InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(filename,
	                                   inputCheck,
	                                   "#InputStartsHere",
	                                   "#InputEndsHere");
	InputNgType::Readable io(ioWriteable);

	ParametersDmrgSolverType dmrgSolverParams(io,false,true);
	dmrgSolverParams.options += sOptions;

	bool setAffinities = (dmrgSolverParams.options.find("setAffinities")
                                  != PsimagLite::String::npos);
	ConcurrencyType::setOptions(dmrgSolverParams.nthreads, setAffinities);

	if (dmrgSolverParams.options.find("useComplex") != PsimagLite::String::npos) {
		mainLoop0<MySparseMatrixComplex>(io,dmrgSolverParams,inputCheck, list);
	} else {
		mainLoop0<MySparseMatrixReal>(io,dmrgSolverParams,inputCheck, list);
	}

	if (filesOption != "keep")
		ArchiveFiles<ParametersDmrgSolverType>::staticDelete();
} // main

