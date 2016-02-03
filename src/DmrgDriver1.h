#ifndef DMRG_DMRGDRIVER_1_H
#define DMRG_DMRGDRIVER_1_H

#include "DmrgDriver.h"


template<typename SolverType>
void mainLoop3(SolverType::LanczosMatrixType::ModelType::GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io,
               const OperatorOptions& opOptions,
               PsimagLite::String targeting)
{
	typedef typename MatrixVectorType::ModelType ModelBaseType;

	//! Setup the Model
	Dmrg::ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
	const ModelBaseType& model = modelSelector(dmrgSolverParams,io,geometry);

	if (opOptions.enabled) {
		operatorDriver(model,opOptions);
		return;
	}

	//! Setup the dmrg solver:
	typedef Dmrg::DmrgSolver<SolverType> DmrgSolverType;
	DmrgSolverType dmrgSolver(model,io);

	//! Calculate observables:
	dmrgSolver.main(geometry,targeting);
}

#endif // DMRG_DMRGDRIVER_1_H

