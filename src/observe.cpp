#include "String.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

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
#include "MettsTargetting.h"
#include "BasisWithOperators.h"
#include "LeftRightSuper.h"
#include "InputNg.h"
#include "Provenance.h"
#include "InputCheck.h"
#include "ModelSelector.h"
#include "ObserverInterpreter.h"

using namespace Dmrg;

typedef std::complex<RealType> ComplexType;

typedef  PsimagLite::CrsMatrix<ComplexType> MySparseMatrixComplex;
typedef  PsimagLite::CrsMatrix<RealType> MySparseMatrixReal;
#ifdef USE_COMPLEX
typedef MySparseMatrixComplex MySparseMatrixC;
#else
typedef MySparseMatrixReal MySparseMatrixC;
#endif
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

template<typename GeometryType,
         typename VectorWithOffsetType,
         typename ModelType,
         typename SparseMatrixType>
bool observeOneFullSweep(IoInputType& io,
                         const ModelType& model,
                         const PsimagLite::String& obsOptions,
                         const PsimagLite::String& list,
                         bool hasTimeEvolution)
{
	const GeometryType& geometry = model.geometry();
	bool verbose = false;

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

	if (hasTimeEvolution) observerLib.setBrackets("time","time");

	const PsimagLite::String& modelName = model.params().model;
	SizeType rows = n; // could be n/2 if there's enough symmetry

	ObserverInterpreter<ObservableLibraryType> observerInterpreter(observerLib);

	observerInterpreter(list,rows,n);

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
         template<typename> class ModelHelperTemplate,
         template<typename> class VectorWithOffsetTemplate,
         typename MySparseMatrix>
void mainLoop(GeometryType& geometry,
              const PsimagLite::String& targetting,
              InputNgType::Readable& io,
              const ParametersDmrgSolverType& params,
              const PsimagLite::String& obsOptions,
              const PsimagLite::String& list)
{
	typedef Basis<MySparseMatrix> BasisType;
	typedef Operators<BasisType> OperatorsType;
	typedef BasisWithOperators<OperatorsType> BasisWithOperatorsType;
	typedef LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
	typedef ModelHelperTemplate<LeftRightSuperType> ModelHelperType;
	typedef ModelBase<ModelHelperType,
	                  ParametersDmrgSolverType,
	                  InputNgType::Readable,
	                  GeometryType> ModelBaseType;
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef VectorWithOffsetTemplate<ComplexOrRealType> VectorWithOffsetType;

	ModelSelector<ModelBaseType> modelSelector(params.model);
	const ModelBaseType& model = modelSelector(params,io,geometry);

	bool moreData = true;
	const PsimagLite::String& datafile = params.filename;
	IoInputType dataIo(datafile);
	bool hasTimeEvolution = (targetting == "TimeStepTargetting" || targetting=="MettsTargetting") ? true : false;
	while (moreData) {
		try {
			moreData = !observeOneFullSweep<GeometryType,VectorWithOffsetType,ModelBaseType,
			            MySparseMatrix>
			(dataIo,model,obsOptions,list,hasTimeEvolution);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"There's no more data\n";
			break;
		}
	}
}

template<typename MySparseMatrix>
void mainLoop0(InputNgType::Readable& io,
               ParametersDmrgSolverType& dmrgSolverParams,
               InputCheck& inputCheck,
               const PsimagLite::String& options,
               const PsimagLite::String& list)
{
	typedef typename MySparseMatrix::value_type ComplexOrRealType;
	typedef PsimagLite::Geometry<ComplexOrRealType,
	                             InputNgType::Readable,
	                             ProgramGlobals> GeometryType;

	GeometryType geometry(io);
	bool su2 = (dmrgSolverParams.options.find("useSu2Symmetry")!=PsimagLite::String::npos);

	PsimagLite::String targetting=inputCheck.getTargeting(dmrgSolverParams.options);

	if (targetting!="GroundStateTargetting" && su2)
		throw PsimagLite::RuntimeError("SU(2) supports only GroundStateTargetting");

	if (su2) {
		if (dmrgSolverParams.targetQuantumNumbers[2]>0) {
			mainLoop<GeometryType,ModelHelperSu2,VectorWithOffsets,MySparseMatrix>
			        (geometry,targetting,io,dmrgSolverParams, options, list);
		} else {
			mainLoop<GeometryType,ModelHelperSu2,VectorWithOffset,MySparseMatrix>
			        (geometry,targetting,io,dmrgSolverParams, options, list);
		}
		return;
	}

	if (targetting=="GroundStateTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffset,MySparseMatrix>
		        (geometry,targetting,io,dmrgSolverParams, options, list);
	} else if (targetting=="TimeStepTargetting") {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,MySparseMatrixComplex>
		        (geometry,targetting,io,dmrgSolverParams, options, list);
	} else {
		mainLoop<GeometryType,ModelHelperLocal,VectorWithOffsets,MySparseMatrix>
		        (geometry,targetting,io,dmrgSolverParams, options, list);
	}
}

void usage(const char* name)
{
	std::cerr<<"USAGE is "<<name<<" -f filename [-o options]\n";
}

int main(int argc,char *argv[])
{
	using namespace Dmrg;

	PsimagLite::String filename="";
	PsimagLite::String options("");
	int opt = 0;
	int precision = 6;
	while ((opt = getopt(argc, argv,"f:o:p:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'o':
			options = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	PsimagLite::String list = (optind < argc) ? argv[optind] : "";

	//sanity checks here
	if (filename=="") {
		usage(argv[0]);
		return 1;
	}

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(&argc,&argv,1);

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	//Setup the Geometry
	InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	//! Read the parameters for this run
	//ParametersModelType mp(io);
	ParametersDmrgSolverType dmrgSolverParams(io);

	ConcurrencyType::npthreads = dmrgSolverParams.nthreads;

#ifdef USE_COMPLEX
		std::cerr<<argv[0]<<" EXPERIMENTAL option complex is in use\n";
		mainLoop0<MySparseMatrixC>(io,dmrgSolverParams,inputCheck, options, list);
#else
		mainLoop0<MySparseMatrixReal>(io,dmrgSolverParams,inputCheck, options, list);
#endif
} // main

