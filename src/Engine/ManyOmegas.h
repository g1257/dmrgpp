#ifndef MANYOMEGAS_H
#define MANYOMEGAS_H
#include "DmrgRunner.h"
#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "PsimagLite.h"
#include "InterNode.h"
#include "OmegaParams.h"

namespace Dmrg {

template<typename ComplexOrRealType, typename OmegaParamsType>
class ManyOmegas {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef DmrgRunner<ComplexOrRealType> DmrgRunnerType;
	typedef typename DmrgRunnerType::InputNgType InputNgType;
	typedef PsimagLite::PsiApp ApplicationType;

	ManyOmegas(PsimagLite::String data,
	           RealType precision,
	           const OmegaParamsType& omegaParams,
	           const ApplicationType& app)
	    : data_(data), runner_(precision, app), omegaParams_(omegaParams)
	{}

	void run(bool dryRun, PsimagLite::String root, PsimagLite::String insitu)
	{
		//lambda
		PsimagLite::InterNode<> internode(PsimagLite::MPI::COMM_WORLD);

		internode.parallelFor(omegaParams_.offset(),
		                      omegaParams_.total(),
		                      [this, root, dryRun, insitu](SizeType i, SizeType)
		{
			const RealType omega = omegaParams_.omega(i);
			PsimagLite::String data2 = addOmega(omega);
			PsimagLite::String outputfile = "\nOutputFile=\"" + root + ttos(i) + "\";\n";
			data2 += outputfile;

			PsimagLite::String logfile = "runForinput" + ttos(i) + ".cout";

			if (dryRun) {
				std::cerr<<"ManyOmegas.h:: omega = "<<omega;
				std::cerr<<" output="<<outputfile;
				std::cerr<<" logfile="<<logfile<<" NOT done because -d\n";
				return;
			}

			runner_.doOneRun(data2, insitu, logfile);
		});
	}

	PsimagLite::String addOmega(RealType wn) const
	{
		const PsimagLite::String str = "real CorrectionVectorOmega=" + ttos(wn) + ";\n";
		return data_ + str;
	}

	PsimagLite::String data_;
	DmrgRunnerType runner_;
	const OmegaParamsType& omegaParams_;
};
}
#endif // MANYOMEGAS_H
