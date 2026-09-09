#ifndef MANYOMEGAS_H
#define MANYOMEGAS_H
#include "DmrgRunner.h"
#include "OmegaParams.h"
#include <PsimagLite/InputNg.h>
#include <PsimagLite/InterNode.h>
#include <PsimagLite/LanczosSolver.h>
#include <PsimagLite/PsiBase64.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename RealType, typename OmegaParamsType> class ManyOmegas {

public:

	using VectorRealType  = typename PsimagLite::Vector<RealType>::Type;
	using DmrgRunnerType  = DmrgRunner<RealType>;
	using InputNgType     = typename DmrgRunnerType::InputNgType;
	using ApplicationType = PsimagLite::PsiApp;

	ManyOmegas(PsimagLite::String     data,
	           const OmegaParamsType& omegaParams,
	           const ApplicationType& app)
	    : data_(data)
	    , omegaParams_(omegaParams)
	    , app_(app)
	{ }

	void run(bool dryRun, PsimagLite::String root, const CmdLineOptions& cmdline_options)
	{
		// lambda
		PsimagLite::InterNode<> internode(PsimagLite::MPI::COMM_WORLD);
		int                     completed = 0;

		internode.parallelFor(
		    omegaParams_.offset(),
		    omegaParams_.total(),
		    [this, root, dryRun, cmdline_options, &completed](SizeType i, SizeType)
		    {
			    const RealType     omega = omegaParams_.omega(i);
			    PsimagLite::String data2 = addOmega(omega);
			    PsimagLite::String outputfile
			        = "\nOutputFile=\"" + root + ttos(i) + "\";\n";
			    data2 += outputfile;

			    CmdLineOptions cmdline_options2 = cmdline_options;

			    cmdline_options2.logfile
			        = std::string("runForinput") + ttos(i) + ".cout";

			    std::cerr << "ManyOmegas.h:: omega = " << omega;
			    std::cerr << " output=" << outputfile;
			    std::cerr << " logfile=" << cmdline_options2.logfile << "\n";

			    if (dryRun) {
				    std::cerr << "NOT done because -d\n";
				    ++completed;
				    return;
			    }

			    DmrgRunnerType runner(app_, data2, cmdline_options2);
			    runner.doOneRun();
			    ++completed;
		    });

		const SizeType mpiSize = PsimagLite::MPI::commSize(PsimagLite::MPI::COMM_WORLD);
		PsimagLite::Vector<int>::Type completedByRank(mpiSize, 0);
		completedByRank[PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD)] = completed;
		PsimagLite::MPI::allReduce(completedByRank);

		int totalCompleted = 0;
		for (SizeType rank = 0; rank < mpiSize; ++rank)
			totalCompleted += completedByRank[rank];

		if (totalCompleted != static_cast<int>(omegaParams_.total()))
			err("ManyOmegas: not all frequency tasks completed\n");

		// When there are at least as many tasks as ranks, every rank must complete
		// one. This collective invariant is a deterministic proof of distributed work.
		if (omegaParams_.total() < mpiSize)
			return;

		for (SizeType rank = 0; rank < mpiSize; ++rank) {
			if (completedByRank[rank] == 0)
				err("ManyOmegas: MPI rank " + ttos(rank)
				    + " did not complete a frequency task\n");
		}
	}

	PsimagLite::String addOmega(RealType wn) const
	{
		const PsimagLite::String str = "CorrectionVectorOmega=" + ttos(wn) + ";\n";
		return data_ + str;
	}

	PsimagLite::String     data_;
	const OmegaParamsType& omegaParams_;
	const ApplicationType& app_;
};
}
#endif // MANYOMEGAS_H
