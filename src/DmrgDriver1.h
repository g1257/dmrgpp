#ifndef DMRG_DMRGDRIVER_1_H
#define DMRG_DMRGDRIVER_1_H

#include "DmrgDriver.h"

template <typename SolverType, typename VectorWithOffsetType>
void mainLoop4(typename SolverType::MatrixType::ModelType::SuperGeometryType& geometry,
               const ParametersDmrgSolverType&                                dmrgSolverParams,
               InputNgType::Readable&                                         io,
               const OperatorOptions&                                         opOptions)
{
	using ModelBaseType = typename SolverType::MatrixType::ModelType;

	//! Setup the Model
	Dmrg::ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
	ModelBaseType&                     model = modelSelector(dmrgSolverParams, io, geometry);

	if (opOptions.enabled) {
		operatorDriver(model, opOptions);
		return;
	}

	//! Setup the dmrg solver:
	using DmrgSolverType = Dmrg::DmrgSolver<SolverType, VectorWithOffsetType>;
	DmrgSolverType dmrgSolver(model, io);

	//! Calculate observables:
	dmrgSolver.main(geometry);
}

#endif // DMRG_DMRGDRIVER_1_H
