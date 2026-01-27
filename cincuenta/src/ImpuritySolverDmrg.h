#ifndef IMPURITYSOLVER_DMRG_H
#define IMPURITYSOLVER_DMRG_H

#include "../../dmrgpp/src/Engine/DmrgRunner.h"
#include "../../dmrgpp/src/Engine/ManyOmegas.h"
#include "../../dmrgpp/src/Engine/ProcOmegas.h"
#include "Geometry/Star.h"
#include "ImpuritySolverBase.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Matsubaras.h"
#include "ParamsDmftSolver.h"
#include "PsiBase64.h"
#include "PsimagLite.h"
#include "Vector.h"

namespace Dmft {

template <typename ParamsDmftSolverType>
class ImpuritySolverDmrg : public ImpuritySolverBase<ParamsDmftSolverType> {

	enum class DmrgType
	{
		TYPE_0,
		TYPE_1
	};

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType   ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType>                             ComplexType;
	typedef typename PsimagLite::Vector<RealType>::Type        VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type     VectorComplexType;
	typedef Dmrg::DmrgRunner<RealType>                         DmrgRunnerType;
	typedef typename DmrgRunnerType::InputNgType               InputNgType;
	typedef PsimagLite::PsiApp                                 ApplicationType;
	typedef Matsubaras<RealType>                               MatsubarasType;
	typedef Dmrg::ManyOmegas<RealType, MatsubarasType>         ManyOmegasType;
	typedef Dmrg::ProcOmegas<RealType, MatsubarasType>         ProcOmegasType;

	static const SizeType CENTER = PsimagLite::Star<ComplexOrRealType, int>::CENTER;

	ImpuritySolverDmrg(const ParamsDmftSolverType& params, const ApplicationType& app)
	    : params_(params)
	    , runner_(params_.precision, app)
	{ }

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);

		if (mpiRank == 0) {
			PsimagLite::String data;
			InputNgType::Writeable::readFile(data, params_.gsTemplate);
			PsimagLite::String data2  = addBathParams(data, bathParams);
			PsimagLite::String insitu = "<gs|nup|gs>";

			runner_.doOneRun(data2, insitu, "-");
		}

		PsimagLite::MPI::barrier(PsimagLite::MPI::COMM_WORLD);

		PsimagLite::String data3;
		InputNgType::Writeable::readFile(data3, params_.omegaTemplate);
		PsimagLite::String data4 = addBathParams(data3, bathParams);

		doType(DmrgType::TYPE_0, data4, mpiRank);

		doType(DmrgType::TYPE_1, data4, mpiRank);

		if (mpiRank == 0) {
			scaleGimp();

			std::cerr << "Sum of Gimp=" << density() << "\n";
			writeGimpForDebugOnly();
		}

		PsimagLite::MPI::barrier(PsimagLite::MPI::COMM_WORLD);
	}

	const VectorComplexType& gimp() const { return gimp_; }

private:

	static PsimagLite::String addBathParams(PsimagLite::String    data,
	                                        const VectorRealType& bathParams)
	{
		const SizeType           nBath      = int(bathParams.size() / 2);
		const PsimagLite::String connectors = findBathParams(0, nBath, bathParams);
		const PsimagLite::String label      = "dir0:Connectors=[" + connectors + "];\n";
		const PsimagLite::String potentialV = findBathParams2(nBath, 2 * nBath, bathParams);
		const PsimagLite::String label2
		    = "potentialV=[" + potentialV + "," + potentialV + "];\n";

		return data + label + label2;
	}

	static PsimagLite::String
	findBathParams(SizeType start, SizeType end, const VectorRealType& bathParams)
	{
		PsimagLite::String buffer = ttos(bathParams[start]);
		for (SizeType i = start + 1; i < end; ++i)
			buffer += "," + ttos(bathParams[i]);

		return buffer;
	}

	static PsimagLite::String
	findBathParams2(SizeType start, SizeType end, const VectorRealType& bathParams)
	{
		PsimagLite::String buffer = ttos(bathParams[start]);
		SizeType           j      = start + 1;
		for (SizeType i = start + 1; i < end + 1; ++i) {
			assert(j < bathParams.size());
			PsimagLite::String tmp = ttos(bathParams[j]);
			if (i - start == CENTER)
				tmp = "0";
			else
				++j;
			buffer += "," + tmp;
		}

		return buffer;
	}

	static PsimagLite::String addTypeAndObs(DmrgType t, PsimagLite::String data)
	{
		const PsimagLite::String obsTc = (t == DmrgType::TYPE_0) ? "c'" : "c";
		const SizeType           tt    = (t == DmrgType::TYPE_0) ? 0 : 1;
		return data + "DynamicDmrgType=" + ttos(tt) + ";\n"
		    + "string TSPOp0:OperatorExpression=\"" + obsTc + "\";\n";
	}

	void doType(DmrgType t, PsimagLite::String data, SizeType mpiRank)
	{
		PsimagLite::String obs     = (t == DmrgType::TYPE_0) ? "c" : "c'";
		PsimagLite::String insitu2 = "<gs|" + obs + "|P2>,<gs|" + obs + "|P3>";

		PsimagLite::String data2 = addTypeAndObs(t, data);

		Matsubaras<RealType> matsubaras(params_.ficticiousBeta, params_.nMatsubaras);

		ManyOmegasType manyOmegas(
		    data2, params_.precision, matsubaras, runner_.application());

		const bool               dryrun   = false;
		const PsimagLite::String rootname = "dmftDynamics";
		manyOmegas.run(dryrun, rootname, insitu2);

		if (mpiRank != 0)
			return;

		const PsimagLite::String rootIname   = "input";
		const PsimagLite::String rootOname   = "OUTPUT";
		const bool               skipFourier = true;

		Dmrg::InputCheck                inputCheck;
		typename InputNgType::Writeable ioW(inputCheck, data2);
		typename InputNgType::Readable  io(ioW);
		ProcOmegasType                  procOmegas(
                    io, params_.precision, skipFourier, rootIname, rootOname, matsubaras);

		procOmegas.run();

		readGimp(rootOname, matsubaras, t);
	}

	void readGimp(PsimagLite::String filename, const MatsubarasType& matsubaras, DmrgType t)
	{
		std::ifstream fin(filename);
		if (!fin || !fin.good() || fin.bad())
			err("readGimp: Could not open " + filename + "\n");

		if (t == DmrgType::TYPE_0)
			gimp_.resize(matsubaras.total());

		SizeType ind = 0;
		while (!fin.eof()) {
			RealType val = 0;
			fin >> val;
			SizeType n = 0;
			fin >> n;

			bool     centerSeen = false;
			RealType val1       = 0;
			RealType val2       = 0;

			for (SizeType i = 0; i <= CENTER; ++i) {
				SizeType site = 0;
				fin >> site;

				fin >> val1;

				fin >> val2;

				if (site == CENTER) {
					centerSeen = true;
					break;
				}
			}

			if (!centerSeen)
				err("Internal error: center " + ttos(CENTER)
				    + " not seen, freq id = " + ttos(ind) + "\n");

			if (t == DmrgType::TYPE_0)
				gimp_[ind] = ComplexType(val2, val1);
			else
				gimp_[ind] += ComplexType(val2, val1);

			++ind;

			if (ind >= gimp_.size())
				break;

			if (n == 1)
				continue;

			const SizeType tmp = CENTER + 1;
			n -= tmp;

			for (SizeType i = 0; i < 3 * n; ++i)
				fin >> val1;
		}

		if (ind < gimp_.size())
			err("readGimp: Not all values computed\n");
	}

	void scaleGimp()
	{
		const SizeType n      = gimp_.size();
		const RealType factor = -2.0; /// M_PI;
		for (SizeType i = 0; i < n; ++i)
			gimp_[i] *= factor;
	}

	ComplexType density() const
	{
		const SizeType n   = gimp_.size();
		ComplexType    sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += gimp_[i];

		return sum;
	}

	void writeGimpForDebugOnly() const
	{
		const SizeType n = gimp_.size();
		std::ofstream  fout("gimp.debug");
		if (!fout || !fout.good())
			err("Could not write to gimp.debug\n");

		Matsubaras<RealType> matsubaras(params_.ficticiousBeta, params_.nMatsubaras);

		for (SizeType i = 0; i < n; ++i) {
			const ComplexType value = gimp_[i];
			const RealType    omega = matsubaras.omega(i);
			fout << omega << " " << PsimagLite::real(value) << " "
			     << PsimagLite::imag(value) << "\n";
		}
	}

	const ParamsDmftSolverType& params_;
	DmrgRunnerType              runner_;
	VectorComplexType           gimp_;
};
}
#endif // IMPURITYSOLVER_DMRG_H
