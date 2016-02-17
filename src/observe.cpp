#include <unistd.h>
#include "Observer.h"
#include "ObservableLibrary.h"
#include "IoSimple.h"
#include "Operators.h"
#include "Geometry/Geometry.h"
#include "CrsMatrix.h"
#include "ModelHelperLocal.h"
#include "ModelHelperSu2.h"
#include "VectorWithOffset.h"
#include "VectorWithOffsets.h"
#include "DmrgSolver.h" // only used for types
#include "TargetingGroundState.h"
#include "TargetingTimeStep.h"
#include "TargetingDynamic.h"
#include "TargetingAdaptiveDynamic.h"
#include "TargetingCorrection.h"
#include "TargetingCorrectionVector.h"
#include "TargetingMetts.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "InputNg.h"
#include "Provenance.h"
#include "InputCheck.h"
#include "ModelSelector.h"
#include "ArchiveFiles.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

using namespace Dmrg;

typedef std::complex<RealType> ComplexType;

typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;
typedef PsimagLite::IoSimple::In IoInputType;
typedef PsimagLite::InputNg<InputCheck> InputNgType;
typedef ParametersDmrgSolver<RealType,InputNgType::Readable> ParametersDmrgSolverType;

template<typename ModelType>
SizeType dofsFromModelName(const ModelType& model)
{
	const PsimagLite::String& modelName = model.params().model;
	SizeType site = 0; // FIXME : account for Hilbert spaces changing with site
	SizeType dofs = SizeType(log(model.hilbertSize(site))/log(2.0));
	std::cerr<<"DOFS= "<<dofs<<" <------------------------------------\n";
	if (modelName.find("FeAsBasedSc")!=PsimagLite::String::npos) return dofs;
	if (modelName.find("FeAsBasedScExtended")!=PsimagLite::String::npos) return dofs;
	if (modelName.find("HubbardOneBand")!=PsimagLite::String::npos) return dofs;

	// max number here, site dependence taken into account elsewhere
	if (modelName.find("Immm")!=PsimagLite::String::npos) return 4;
	return 0;
}

template<typename VectorWithOffsetType,
         typename ModelType>
bool observeOneFullSweep(IoInputType& io,
                         const ModelType& model,
                         const PsimagLite::String& list2,
                         bool hasTimeEvolution)
{
	typedef typename ModelType::GeometryType GeometryType;
	const GeometryType& geometry = model.geometry();
	bool verbose = false;
	PsimagLite::String obsOptions("");
	PsimagLite::String list = list2;

	if (list2.length() > 0 && list2[0] != '<') {
        obsOptions = list2;
		list = "";
	}

	typedef Observer<VectorWithOffsetType,ModelType,IoInputType> ObserverType;
	typedef ObservableLibrary<ObserverType> ObservableLibraryType;

	SizeType n  = geometry.numberOfSites();

	ObservableLibraryType observerLib(io,n,hasTimeEvolution,model,verbose);

	bool ot = false;
	if (obsOptions.find("ot")!=PsimagLite::String::npos ||
	        obsOptions.find("time")!=PsimagLite::String::npos) ot = true;

	if (hasTimeEvolution && ot) {
		observerLib.measureTime("superDensity");
		observerLib.measureTime("nupNdown");
		observerLib.measureTime("nup+ndown");
		if (obsOptions.find("sz")!=PsimagLite::String::npos)
			observerLib.measureTime("sz");
	}

	if (hasTimeEvolution) observerLib.setBrakets("time","time");

	const PsimagLite::String& modelName = model.params().model;
	SizeType rows = n; // could be n/2 if there's enough symmetry

	observerLib.interpret(list,rows,n);

	// Immm supports only onepoint:
	if (modelName=="Immm" && obsOptions!="onepoint") {
		PsimagLite::String str(__FILE__);
		str += " "  + ttos(__LINE__) + "\n";
		str += "Model Immm only supports onepoint\n";
		throw PsimagLite::RuntimeError(str.c_str());
	}

	SizeType numberOfDofs = dofsFromModelName(model);

	if (!hasTimeEvolution && obsOptions.find("onepoint")!=PsimagLite::String::npos) {
		observerLib.measureTheOnePoints(numberOfDofs);
	}

	if (modelName.find("Heisenberg")==PsimagLite::String::npos) {
		if (obsOptions.find("cc")!=PsimagLite::String::npos) {
			observerLib.measure("cc",rows,n);
		}

		if (obsOptions.find("nn")!=PsimagLite::String::npos) {
			observerLib.measure("nn",rows,n);
		}
	}
	if (obsOptions.find("szsz")!=PsimagLite::String::npos) {
		observerLib.measure("szsz",rows,n);
	}

	if (modelName.find("FeAsBasedSc")!=PsimagLite::String::npos ||
	        modelName.find("FeAsBasedScExtended")!=PsimagLite::String::npos ||
	        modelName.find("HubbardOneBand")!=PsimagLite::String::npos) {
		bool dd4 = (obsOptions.find("dd4")!=PsimagLite::String::npos);

		if (obsOptions.find("dd")!=PsimagLite::String::npos && !dd4) { // &&
			//geometry.label(0).find("ladder")!=PsimagLite::String::npos) {
			observerLib.measure("dd",rows,n);
		}

		// FOUR-POINT DELTA-DELTA^DAGGER:
		if (dd4 && geometry.label(0).find("ladder")!=PsimagLite::String::npos) {
			observerLib.measure("dd4",rows,n);
		} // if dd4
	}

	if (modelName.find("HubbardOneBand")!=PsimagLite::String::npos &&
	        obsOptions.find("multi")!=PsimagLite::String::npos) {
		observerLib.measure("multi",rows,n);
	}

	if (obsOptions.find("s+s-")!=PsimagLite::String::npos) {
		observerLib.measure("s+s-",rows,n);
	}

	if (obsOptions.find("ss")!=PsimagLite::String::npos) {
		observerLib.measure("ss",rows,n);
	}

	return observerLib.endOfData();
}

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

	ModelSelector<ModelBaseType> modelSelector(params.model);
	const ModelBaseType& model = modelSelector(params,io,geometry);

	bool moreData = true;
	const PsimagLite::String& datafile = params.filename;
	ArchiveFiles<ParametersDmrgSolverType>::unpackIfNeeded(datafile);
	IoInputType dataIo(datafile);
	bool hasTimeEvolution = (targetting == "TimeStepTargetting" ||
	                         targetting=="MettsTargetting" ||
	                         targetting=="TargetingAncilla");
	while (moreData) {
		try {
			moreData = !observeOneFullSweep<VectorWithOffsetType,ModelBaseType>
			        (dataIo,model,list,hasTimeEvolution);
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
	typedef Basis<MySparseMatrix> BasisType;
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

int main(int argc,char *argv[])
{
	using namespace Dmrg;

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
	ConcurrencyType concurrency(&argc,&argv,1);

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;

	//Setup the Geometry
	InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	//! Read the parameters for this run
	//ParametersModelType mp(io);
	ParametersDmrgSolverType dmrgSolverParams(io,false,true);

	ConcurrencyType::npthreads = dmrgSolverParams.nthreads;

	if (dmrgSolverParams.options.find("useComplex") != PsimagLite::String::npos) {
		std::cerr<<argv[0]<<" EXPERIMENTAL option complex is in use\n";
		mainLoop0<MySparseMatrixComplex>(io,dmrgSolverParams,inputCheck, list);
	} else {
		mainLoop0<MySparseMatrixReal>(io,dmrgSolverParams,inputCheck, list);
	}

	if (filesOption != "keep")
		ArchiveFiles<ParametersDmrgSolverType>::staticDelete();
} // main

