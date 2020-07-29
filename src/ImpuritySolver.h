#ifndef IMPURITYSOLVER_H
#define IMPURITYSOLVER_H
#include "Vector.h"
#include "../../dmrgpp/src/Engine/InputCheck.h"
#include "../../dmrgpp/src/Engine/Qn.h"
#include "../../dmrgpp/src/Engine/ProgramGlobals.h"
#include "../../dmrgpp/src/Engine/SuperGeometry.h"
#include "../../dmrgpp/src/Engine/ParametersDmrgSolver.h"
#include "InputNg.h"

namespace Dmft {

template<typename ComplexOrRealType>
class ImpuritySolver {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef Dmrg::ParametersDmrgSolver<RealType, InputNgType::Readable, Dmrg::Qn>
	ParametersDmrgSolverType;
	typedef Dmrg::SuperGeometry<ComplexOrRealType,
	        InputNgType::Readable,
	        Dmrg::ProgramGlobals> SuperGeometryType;

	ImpuritySolver(PsimagLite::String gsTemplate, PsimagLite::String omegaTemplate)
	    : gsTemplate_(gsTemplate), omegaTemplate_(omegaTemplate)
	{}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{

		PsimagLite::String data;
		InputNgType::Writeable::readFile(data, gsTemplate_);
		PsimagLite::String data2 = modifyBathParams(data, bathParams);
		PsimagLite::String sOptions = "";

		doOneRun(data2, sOptions);
	}

	ComplexOrRealType gimp(SizeType i)
	{
		return 0.0;
	}

private:

	static PsimagLite::String modifyBathParams(PsimagLite::String data,
	                                           const VectorRealType& bathParams)
	{
		const SizeType nBath = int(bathParams.size() / 2);
		static const PsimagLite::String label = "dir0:Connectors=";
		size_t pos1 = data.find(label);
		if (pos1 == PsimagLite::String::npos)
			err("modifyBathParams(): cannot find " + label + "\n");

		PsimagLite::String connectors = findBathParams(0, nBath, bathParams);
		PsimagLite::String buffer = data.substr(0, pos1) + label + connectors + ";";

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
		buffer += data.substr(pos2 + 1, len) + label2 + potentialV;

		size_t pos4 = data.find(";", pos3);
		if (pos4 < pos3 + 1 || pos4 == PsimagLite::String::npos)
			err("modifyBathParams(): internal error(2)\n");

		const SizeType len2 = data.length() - pos4 - 1;
		buffer += data.substr(pos4, len2);
		std::cerr<<buffer<<"\n";
		return buffer;
	}


	void doOneRun(PsimagLite::String data, PsimagLite::String sOptions)
	{
		Dmrg::InputCheck inputCheck;
		InputNgType::Writeable ioWriteable(inputCheck, data);
		InputNgType::Readable io(ioWriteable);

		ParametersDmrgSolverType dmrgSolverParams(io, sOptions, false);

		SuperGeometryType geometry(io);
		if (dmrgSolverParams.options.isSet("printgeometry"))
			std::cout<<geometry;

//		//! Setup the Model
//		Dmrg::ModelSelector<ModelBaseType> modelSelector(dmrgSolverParams.model);
//		const ModelBaseType& model = modelSelector(dmrgSolverParams,io,geometry);

//		//! Setup the dmrg solver:
//		typedef Dmrg::DmrgSolver<SolverType, VectorWithOffsetType> DmrgSolverType;
//		DmrgSolverType dmrgSolver(model,io);

//		//! Calculate observables:
//		dmrgSolver.main(geometry);
	}

	static PsimagLite::String findBathParams(SizeType start,
	                                         SizeType end,
	                                         const VectorRealType& bathParams)
	{
		PsimagLite::String buffer = "[" + ttos(bathParams[start]);
		for (SizeType i = start + 1; i < end; ++i)
			buffer += "," + ttos(bathParams[i + start]);

		buffer += "]";
		return buffer;
	}

	PsimagLite::String gsTemplate_;
	PsimagLite::String omegaTemplate_;
};
}
#endif // IMPURITYSOLVER_H
