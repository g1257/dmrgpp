#ifndef MANYOMEGAS_H
#define MANYOMEGAS_H
#include "Engine/DmrgRunner.h"
#include "Engine/OmegaParams.h"
#include "InputNg.h"
#include "InterNode.h"
#include "LanczosSolver.h"
#include "PsiBase64.h"
#include "PsimagLite.h"
#include "Vector.h"

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

		internode.parallelFor(
		    omegaParams_.offset(),
		    omegaParams_.total(),
		    [this, root, dryRun, cmdline_options](SizeType i, SizeType)
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
			    std::cerr << " logfile=" << cmdline_options2.logfile << " MPI rank=";
			    std::cerr << PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD)
			              << "\n";

			    if (dryRun) {
				    std::cerr << "NOT done because -d\n";
				    return;
			    }

			    DmrgRunnerType runner(app_, data2, cmdline_options2);
			    runner.doOneRun();
		    });
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
