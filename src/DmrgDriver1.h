#ifndef DMRG_DMRGDRIVER_1_H
#define DMRG_DMRGDRIVER_1_H

#include "DmrgDriver.h"

template <typename RealOrComplexType,
          template <typename> typename MatrixVectorTemplate,
          template <typename> typename LanczosSolverTemplate>
struct ExplicitInstantiationHelper {
	using SparseMatrixType = PsimagLite::CrsMatrix<RealOrComplexType>;
	using GeometryType     = Dmrg::SuperGeometry<RealOrComplexType,
	                                             PsimagLite::InputNg<Dmrg::InputCheck>::Readable,
	                                             Dmrg::ProgramGlobals>;
	using MatrixVectorType = MatrixVectorTemplate<
	    Dmrg::ModelBase<Dmrg::ModelHelperLocal<Dmrg::LeftRightSuper<
	                        Dmrg::BasisWithOperators<Dmrg::Basis<SparseMatrixType>>,
	                        Dmrg::Basis<SparseMatrixType>>>,
	                    ParametersDmrgSolverType,
	                    InputNgType::Readable,
	                    GeometryType>>;

	using SolverType           = LanczosSolverTemplate<MatrixVectorType>;
	using VectorWithOffsetType = Dmrg::VectorWithOffset<RealOrComplexType, Dmrg::Qn>;

	using VectorWithOffsetsType = Dmrg::VectorWithOffsets<RealOrComplexType, Dmrg::Qn>;
};

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
