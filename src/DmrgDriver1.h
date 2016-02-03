#ifndef DMRG_DMRGDRIVER_1_H
#define DMRG_DMRGDRIVER_1_H

#include "DmrgDriver.h"

template<typename GeometryType, typename TargettingType>
void mainLoop3(GeometryType& geometry,
               const ParametersDmrgSolverType& dmrgSolverParams,
               InputNgType::Readable& io,
               const OperatorOptions& opOptions,
               PsimagLite::String targeting)
{
	typedef typename TargettingType::TargetParamsType TargetParamsType;
	typedef typename TargettingType::MatrixVectorType::ModelType ModelBaseType;

	//! Setup the Model
	Dmrg::ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
	const ModelBaseType& model = modelSelector(dmrgSolverParams,io,geometry);

	if (opOptions.enabled) {
		operatorDriver(model,opOptions);
		return;
	}

	//! Read TimeEvolution if applicable:
	TargetParamsType tsp(io,model);

	//! Setup the dmrg solver:
	typedef Dmrg::DmrgSolver<TargettingType> SolverType;
	SolverType dmrgSolver(model,tsp,io);

	//! Calculate observables:
	dmrgSolver.main(geometry,targeting);
}

#endif // DMRG_DMRGDRIVER_1_H

