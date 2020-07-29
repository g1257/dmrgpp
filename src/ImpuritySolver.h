#ifndef IMPURITYSOLVER_H
#define IMPURITYSOLVER_H

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

#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "ParamsDmftSolver.h"

namespace Dmft {

template<typename ParamsDmftSolverType>
class ImpuritySolver {

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef Dmrg::ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn>
	ParametersDmrgSolverType;
	typedef Dmrg::SuperGeometry<ComplexOrRealType,
	        InputNgType::Readable,
	        Dmrg::ProgramGlobals> SuperGeometryType;
	typedef Dmrg::VectorWithOffset<ComplexOrRealType, Dmrg::Qn> VectorWithOffsetType;

	ImpuritySolver(const ParamsDmftSolverType& params)
	    : params_(params)
	{}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{

		PsimagLite::String data;
		InputNgType::Writeable::readFile(data, params_.gsTemplate);
		PsimagLite::String data2 = modifyBathParams(data, bathParams, params_.echoInput);
		PsimagLite::String sOptions = "";

		doOneRun(data2, sOptions);
	}

	ComplexOrRealType gimp(SizeType i)
	{
		throw PsimagLite::RuntimeError("gimp not ready yet\n");
	}

private:

	static PsimagLite::String modifyBathParams(PsimagLite::String data,
	                                           const VectorRealType& bathParams,
	                                           bool echoInput)
	{
		const SizeType nBath = int(bathParams.size() / 2);
		static const PsimagLite::String label = "dir0:Connectors=";
		size_t pos1 = data.find(label);
		if (pos1 == PsimagLite::String::npos)
			err("modifyBathParams(): cannot find " + label + "\n");

		PsimagLite::String connectors = findBathParams(0, nBath, bathParams);
		PsimagLite::String buffer = data.substr(0, pos1) + label + "[" + connectors + "]";

		size_t pos2 = data.find(";", pos1);
		if (pos2 < pos1 || pos2 == PsimagLite::String::npos)
			err("modifyBathParams(): Internal error\n");

		static const PsimagLite::String label2 = "potentialV=";
		size_t pos3 = data.find(label2);
		if (pos3 < pos2 + 1 || pos3 == PsimagLite::String::npos)
			err("modifyBathParams(): cannot find " + label2 + "\n");

		SizeType len = pos3 - pos2 - 1;
		if (bathParams.size() % 2)
			err("BathParams not even\n");

		PsimagLite::String potentialV = findBathParams(nBath, 2*nBath, bathParams);
		buffer += data.substr(pos2, len) + label2 + "[0, " + potentialV +
		        ", 0, " + potentialV + "]";

		size_t pos4 = data.find(";", pos3);
		if (pos4 < pos3 + 1 || pos4 == PsimagLite::String::npos)
			err("modifyBathParams(): internal error(2)\n");

		const SizeType len2 = data.length() - pos4 - 1;
		buffer += data.substr(pos4, len2);
		if (echoInput) echoBase64(std::cout, data);
		else std::cout<<data;
		return buffer;
	}

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

		if (params_.precision > 0) dmrgSolverParams.precision = params_.precision;

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

	static PsimagLite::String findBathParams(SizeType start,
	                                         SizeType end,
	                                         const VectorRealType& bathParams)
	{
		PsimagLite::String buffer = ttos(bathParams[start]);
		for (SizeType i = start + 1; i < end; ++i)
			buffer += "," + ttos(bathParams[i]);

		return buffer;
	}

	static void echoBase64(std::ostream& os, const PsimagLite::String& str)
	{
		os<<"ImpuritySolver::echoBase64: Echo of [[data]] in base64\n";
		PsimagLite::PsiBase64::Encode base64(str);
		os<<base64()<<"\n";
	}

	const ParamsDmftSolverType& params_;
};
}
#endif // IMPURITYSOLVER_H
