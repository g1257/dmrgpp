#ifndef DMRGRUNNER_H
#define DMRGRUNNER_H
#include "BasisWithOperators.h"
#include "CrsMatrix.h"
#include "DmrgSolver.h"
#include "InputCheck.h"
#include "LeftRightSuper.h"
#include "MatrixVectorKron/MatrixVectorKron.h"
#include "MatrixVectorOnTheFly.h"
#include "MatrixVectorStored.h"
#include "ModelHelperLocal.h"
#include "ModelSelector.h"
#include "ParametersDmrgSolver.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "Qn.h"
#include "SuperGeometry.h"
#include "VectorWithOffset.h"

namespace Dmrg
{

template <typename ComplexOrRealType>
class DmrgRunner
{

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef Dmrg::ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn>
	    ParametersDmrgSolverType;
	typedef Dmrg::SuperGeometry<ComplexOrRealType,
	    InputNgType::Readable,
	    Dmrg::ProgramGlobals>
	    SuperGeometryType;
	typedef Dmrg::VectorWithOffset<ComplexOrRealType, Dmrg::Qn> VectorWithOffsetType;
	typedef PsimagLite::PsiApp ApplicationType;

	DmrgRunner(RealType precision, const ApplicationType& application)
	    : precision_(precision)
	    , application_(application)
	{
	}

	void doOneRun(PsimagLite::String data,
	    PsimagLite::String insitu,
	    PsimagLite::String logfile) const
	{
		typedef PsimagLite::CrsMatrix<ComplexOrRealType> MySparseMatrixComplex;
		typedef Dmrg::Basis<MySparseMatrixComplex> BasisType;
		typedef Dmrg::BasisWithOperators<BasisType> BasisWithOperatorsType;
		typedef Dmrg::LeftRightSuper<BasisWithOperatorsType, BasisType> LeftRightSuperType;
		typedef Dmrg::ModelHelperLocal<LeftRightSuperType> ModelHelperType;
		typedef Dmrg::ModelBase<ModelHelperType,
		    ParametersDmrgSolverType,
		    InputNgType::Readable,
		    SuperGeometryType>
		    ModelBaseType;

		std::streambuf* globalCoutBuffer = 0;
		std::ofstream globalCoutStream;
		if (logfile != "-") {
			globalCoutBuffer = std::cout.rdbuf(); // save old buf
			globalCoutStream.open(logfile.c_str(), std::ofstream::out);
			if (!globalCoutStream || globalCoutStream.bad() || !globalCoutStream.good()) {
				PsimagLite::String str(application_.name());
				str += ": Could not redirect std::cout to " + logfile + "\n";
				err(str);
			}

			std::cout.rdbuf(globalCoutStream.rdbuf()); // redirect std::cout to file

			printOutputChange(logfile, data);
		}

		Dmrg::InputCheck inputCheck;
		InputNgType::Writeable ioWriteable(inputCheck, data);
		InputNgType::Readable io(ioWriteable);

		ParametersDmrgSolverType dmrgSolverParams(io, "", false);
		if (dmrgSolverParams.options.isSet("hd5DontPrint"))
			PsimagLite::IoNg::dontPrintDebug();

		if (precision_ > 0)
			dmrgSolverParams.precision = precision_;

		dmrgSolverParams.insitu = insitu;

		if (dmrgSolverParams.options.isSet("MatrixVectorStored")) {
			doOneRun2<Dmrg::MatrixVectorStored<ModelBaseType>>(dmrgSolverParams, io);
		} else if (dmrgSolverParams.options.isSet("MatrixVectorOnTheFly")) {
			doOneRun2<Dmrg::MatrixVectorOnTheFly<ModelBaseType>>(dmrgSolverParams, io);
		} else {
			doOneRun2<Dmrg::MatrixVectorKron<ModelBaseType>>(dmrgSolverParams, io);
		}

		if (logfile == "-" || globalCoutBuffer == 0)
			return;
		globalCoutStream.close();
		std::cout.rdbuf(globalCoutBuffer);
	}

	void printOutputChange(PsimagLite::String logfile, PsimagLite::String data) const
	{
		std::cerr << Provenance::logo(application_.name());
		std::cerr << "Standard output sent to ";
		std::cerr << logfile << "\n";
		std::cerr.flush();

		application_.printCmdLine(std::cout);
		application_.base64encode(std::cout, data, true);
	}

	const ApplicationType& application() const { return application_; }

private:

	template <typename MatrixVectorType>
	void doOneRun2(const ParametersDmrgSolverType& dmrgSolverParams,
	    InputNgType::Readable& io) const
	{
		SuperGeometryType geometry(io);
		if (dmrgSolverParams.options.isSet("printgeometry"))
			std::cout << geometry;

		typedef PsimagLite::ParametersForSolver<typename MatrixVectorType::RealType>
		    ParametersForSolverType;
		typedef PsimagLite::LanczosSolver<ParametersForSolverType,
		    MatrixVectorType,
		    typename MatrixVectorType::VectorType>
		    SolverType;
		typedef typename SolverType::MatrixType::ModelType ModelBaseType;

		//! Setup the Model
		Dmrg::ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
		ModelBaseType& model = modelSelector(dmrgSolverParams, io, geometry);

		//! Setup the dmrg solver: (vectorwithoffset.h only):
		typedef Dmrg::DmrgSolver<SolverType, VectorWithOffsetType> DmrgSolverType;
		DmrgSolverType dmrgSolver(model, io);

		//! Calculate observables:
		dmrgSolver.main(geometry);

		std::cout.flush();
	}

	static void echoBase64(std::ostream& os, const PsimagLite::String& str)
	{
		os << "ImpuritySolver::echoBase64: Echo of [[data]] in base64\n";
		PsimagLite::PsiBase64::Encode base64(str);
		os << base64() << "\n";
	}

	RealType precision_;
	const ApplicationType& application_;
};
}
#endif // DMRGRUNNER_H
