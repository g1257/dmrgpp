#include "DmrgRunner.h"
#include "MatrixVectorTypes.hpp"
#include <type_traits>

namespace Dmrg {

template <typename RealType>
DmrgRunner<RealType>::DmrgRunner(const ApplicationType&    app,
                                 const PsimagLite::String& data,
                                 const CmdLineOptions&     cmd_line)
    : application_(app)
    , io_(nullptr)
    , dmrg_solver_params_(nullptr)
{
	// io_ doesn't need ioWritable past its construction, else
	// we would have to store it as member data
	{
		InputNgType::Writeable ioWriteable(inputCheck_, data);
		io_ = new InputNgType::Readable(ioWriteable);
	}

	assert(io_);
	dmrg_solver_params_ = new ParametersDmrgSolverType(*io_, cmd_line.solver_options);

	bool introspect_only = (dmrg_solver_params_->options.isSet("introspect"));

	std::string insitu   = cmd_line.in_situ_measurements;
	SizeType    nthreads = cmd_line.number_of_threads;
	if (introspect_only) {
		insitu   = "";
		nthreads = 1;
	}

	adjustDmrgSolverParams(insitu, cmd_line.precision, nthreads);

	if (!introspect_only) {
		if (endDueToClobber(dmrg_solver_params_->options.isSet("noClobber"),
		                    cmd_line.logfile)) {
			return;
		}

		dealWithConsoleOutput(cmd_line.logfile, cmd_line.unbuffered_output);

		constexpr bool PRINT_HEADER = true;
		application_.base64encode(std::cout, data, PRINT_HEADER);

		application_.printCmdLine(std::cout);

		setCorrectThreading();
	}
}

template <typename RealType> DmrgRunner<RealType>::~DmrgRunner()
{
	delete dmrg_solver_params_;
	dmrg_solver_params_ = nullptr;

	delete io_;
	io_ = nullptr;
}

template <typename RealType> void DmrgRunner<RealType>::doOneRun() const
{
	OptionsForIntrospect op_options;
	this->doOneRun((op_options));
}

template <typename RealType>
void DmrgRunner<RealType>::doOneRun(OptionsForIntrospect& op_options) const
{
	assert(dmrg_solver_params_);
	bool isComplex = (dmrg_solver_params_->options.isSet("useComplex")
	                  || dmrg_solver_params_->options.isSet("TimeStepTargeting"));

	if (isComplex)
		doOneRun2<std::complex<RealType>>(op_options);
	else
		doOneRun2<RealType>(op_options);
}

template <typename RealType>
template <typename ComplexOrRealType>
void DmrgRunner<RealType>::doOneRun2(const OptionsForIntrospect& op_options) const
{
	using SuperGeometryType
	    = SuperGeometry<ComplexOrRealType, InputNgType::Readable, ProgramGlobals>;
	using MySparseMatrixComplex  = PsimagLite::CrsMatrix<ComplexOrRealType>;
	using BasisType              = Basis<MySparseMatrixComplex>;
	using BasisWithOperatorsType = BasisWithOperators<BasisType>;
	using LeftRightSuperType     = LeftRightSuper<BasisWithOperatorsType, BasisType>;
	using ModelHelperType        = ModelHelperLocal<LeftRightSuperType>;
	using ModelBaseType          = ModelBase<ModelHelperType,
	                                         ParametersDmrgSolverType,
	                                         InputNgType::Readable,
	                                         SuperGeometryType>;
	using MatrixVectorModelType  = typename MatrixVectorTypes<ComplexOrRealType>::ModelType;
	static_assert(std::is_same<ModelBaseType, MatrixVectorModelType>::value,
	              "DmrgRunner and MatrixVectorTypes must use the same model type");

	assert(dmrg_solver_params_);
	if (dmrg_solver_params_->options.isSet("MatrixVectorStored")) {
		doOneRun3<MatrixVectorStored<ComplexOrRealType>>(op_options);
	} else if (dmrg_solver_params_->options.isSet("MatrixVectorOnTheFly")) {
		doOneRun3<MatrixVectorOnTheFly<ComplexOrRealType>>(op_options);
	} else {
		doOneRun3<MatrixVectorKron<ComplexOrRealType>>(op_options);
	}
}

template <typename RealType>
template <typename MatrixVectorType>
void DmrgRunner<RealType>::doOneRun3(const OptionsForIntrospect& op_options) const
{
	using ComplexOrRealType = typename MatrixVectorType::ComplexOrRealType;
	using SuperGeometryType
	    = SuperGeometry<ComplexOrRealType, InputNgType::Readable, ProgramGlobals>;
	using VectorWithOffsetsType = VectorWithOffsets<ComplexOrRealType>;

	assert(io_);
	assert(dmrg_solver_params_);

	SuperGeometryType geometry(*io_);
	if (dmrg_solver_params_->options.isSet("printgeometry"))
		std::cout << geometry;

	using SolverType    = PsimagLite::LanczosSolver<MatrixVectorType>;
	using ModelBaseType = typename SolverType::MatrixType::ModelType;

	//! Setup the Model
	ModelSelector<ModelBaseType> modelSelector(dmrg_solver_params_->model);
	ModelBaseType&               model = modelSelector(*dmrg_solver_params_, *io_, geometry);

	Introspect introspect(dmrg_solver_params_->options.isSet("introspect"));

	if (introspect(model, op_options)) {
		return; // <--- EARLY EXIT HERE
	}

	//! Setup the dmrg solver
	using DmrgSolverType = DmrgSolver<SolverType, VectorWithOffsetsType>;
	DmrgSolverType dmrgSolver(model, *io_);

	//! Calculate observables:
	dmrgSolver.main(geometry);

	std::cout.flush();
}

template <typename RealType>
void DmrgRunner<RealType>::dealWithConsoleOutput(const std::string& logfile, bool unbuffered)
{
	if (logfile == "-") {
		return;
	}

	std::ios_base::openmode open_mode
	    = (dmrg_solver_params_->autoRestart) ? std::ofstream::app : std::ofstream::out;

	redirect_output_.doIt(logfile, open_mode, unbuffered);
}

template <typename RealType>
void DmrgRunner<RealType>::adjustDmrgSolverParams(const std::string& insitu,
                                                  SizeType           precision,
                                                  SizeType           threadsInCmd)
{
	assert(dmrg_solver_params_);
	if (threadsInCmd > 0)
		dmrg_solver_params_->nthreads = threadsInCmd;

	if (precision > 0)
		dmrg_solver_params_->precision = precision;

	dmrg_solver_params_->insitu = insitu;

	if (dmrg_solver_params_->options.isSet("minimizeDisk")) {
		dmrg_solver_params_->options += ",noSaveWft,noSaveStacks,noSaveData";
	}

	// The below doesn't change dmrgSolverParams
	if (dmrg_solver_params_->options.isSet("hd5DontPrint"))
		PsimagLite::IoNg::dontPrintDebug();

	if (dmrg_solver_params_->autoRestart) {
		std::cout << "\nAutoRestart possible\n";
	}
}

template <typename RealType> void DmrgRunner<RealType>::setCorrectThreading() const
{
	typedef PsimagLite::Concurrency ConcurrencyType;

	assert(io_);
	SizeType threadsStackSize = 0;
	try {
		io_->readline(threadsStackSize, "ThreadsStackSize=");
	} catch (std::exception&) { }

	constexpr bool setAffinities = false; // no longer supported, I mean, deleted

	assert(dmrg_solver_params_);
	PsimagLite::CodeSectionParams codeSection(dmrg_solver_params_->nthreads,
	                                          dmrg_solver_params_->nthreads2,
	                                          setAffinities,
	                                          threadsStackSize);
	ConcurrencyType::setOptions(codeSection);
}

template <typename RealType>
bool DmrgRunner<RealType>::endDueToClobber(bool is_no_clobber_set, const std::string& logfile)
{
	if (logfile == "-")
		return false;

	RunFinished runFinished(is_no_clobber_set);
	if (runFinished.OK(logfile.c_str())) {
		runFinished.printTermination(logfile.c_str());
		return true;
	}

	return false;
}

template class DmrgRunner<double>;
}
