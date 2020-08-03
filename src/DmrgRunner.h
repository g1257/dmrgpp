#ifndef DMRGRUNNER_H
#define DMRGRUNNER_H
#include "CrsMatrix.h"
#include "../../dmrgpp/src/Engine/InputCheck.h"
#include "../../dmrgpp/src/Engine/Qn.h"
#include "../../dmrgpp/src/Engine/ProgramGlobals.h"
#include "../../dmrgpp/src/Engine/SuperGeometry.h"
#include "../../dmrgpp/src/Engine/ParametersDmrgSolver.h"
#include "../../dmrgpp/src/Engine/ModelSelector.h"
#include "../../dmrgpp/src/Engine/ModelHelperLocal.h"
#include "../../dmrgpp/src/Engine/MatrixVectorKron/MatrixVectorKron.h"
#include "../../dmrgpp/src/Engine/MatrixVectorOnTheFly.h"
#include "../../dmrgpp/src/Engine/MatrixVectorStored.h"
#include "../../dmrgpp/src/Engine/LeftRightSuper.h"
#include "../../dmrgpp/src/Engine/BasisWithOperators.h"
#include "../../dmrgpp/src/Engine/DmrgSolver.h"
#include "../../dmrgpp/src/Engine/VectorWithOffset.h"

namespace Dmft {

template<typename  ComplexOrRealType>
class DmrgRunner {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef Dmrg::ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn>
	ParametersDmrgSolverType;
	typedef Dmrg::SuperGeometry<ComplexOrRealType,
	InputNgType::Readable,
	Dmrg::ProgramGlobals> SuperGeometryType;
	typedef Dmrg::VectorWithOffset<ComplexOrRealType, Dmrg::Qn> VectorWithOffsetType;

	DmrgRunner(RealType precision) : precision_(precision)
	{}

	void doOneRun(PsimagLite::String data, PsimagLite::String sOptions)
	{
		typedef  PsimagLite::CrsMatrix<std::complex<RealType> > MySparseMatrixComplex;
		typedef Dmrg::Basis<MySparseMatrixComplex> BasisType;
		typedef Dmrg::BasisWithOperators<BasisType> BasisWithOperatorsType;
		typedef Dmrg::LeftRightSuper<BasisWithOperatorsType,BasisType> LeftRightSuperType;
		typedef Dmrg::ModelHelperLocal<LeftRightSuperType> ModelHelperType;
		typedef Dmrg::ModelBase<ModelHelperType,
		        ParametersDmrgSolverType,
		        InputNgType::Readable,
		        SuperGeometryType> ModelBaseType;

		Dmrg::InputCheck inputCheck;
		InputNgType::Writeable ioWriteable(inputCheck, data);
		InputNgType::Readable io(ioWriteable);

		ParametersDmrgSolverType dmrgSolverParams(io, sOptions, false);

		if (precision_ > 0) dmrgSolverParams.precision = precision_;

		if (dmrgSolverParams.options.isSet("MatrixVectorStored")) {
			doOneRun2<Dmrg::MatrixVectorStored<ModelBaseType> >(dmrgSolverParams, io);
		} else if (dmrgSolverParams.options.isSet("MatrixVectorOnTheFly")) {
			doOneRun2<Dmrg::MatrixVectorOnTheFly<ModelBaseType> >(dmrgSolverParams, io);
		} else {
			doOneRun2<Dmrg::MatrixVectorKron<ModelBaseType> >(dmrgSolverParams, io);
		}
	}

	template<typename MatrixVectorType>
	void doOneRun2(const ParametersDmrgSolverType& dmrgSolverParams, InputNgType::Readable& io)
	{
		SuperGeometryType geometry(io);
		if (dmrgSolverParams.options.isSet("printgeometry"))
			std::cout<<geometry;

		typedef PsimagLite::ParametersForSolver<typename MatrixVectorType::RealType>
		        ParametersForSolverType;
		typedef PsimagLite::LanczosSolver<ParametersForSolverType,
		        MatrixVectorType, typename MatrixVectorType::VectorType> SolverType;
		typedef typename SolverType::MatrixType::ModelType ModelBaseType;

		//! Setup the Model
		Dmrg::ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
		const ModelBaseType& model = modelSelector(dmrgSolverParams, io, geometry);

		//! Setup the dmrg solver: (vectorwithoffset.h only):
		typedef Dmrg::DmrgSolver<SolverType, VectorWithOffsetType> DmrgSolverType;
		DmrgSolverType dmrgSolver(model,io);

		//! Calculate observables:
		dmrgSolver.main(geometry);

		std::cout.flush();
	}

private:

	RealType precision_;
};
}
#endif // DMRGRUNNER_H
