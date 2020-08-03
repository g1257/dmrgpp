#ifndef MANYOMEGAS_H
#define MANYOMEGAS_H
#include "DmrgRunner.h"
#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "PsimagLite.h"

namespace Dmrg {

template<typename ComplexOrRealType>
class ManyOmegas {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef DmrgRunner<ComplexOrRealType> DmrgRunnerType;
	typedef typename DmrgRunnerType::InputNgType InputNgType;
	typedef PsimagLite::PsiApp ApplicationType;

	struct OmegaParams {

		OmegaParams(PsimagLite::String data)
		{
			Dmrg::InputCheck inputCheck;
			typename InputNgType::Writeable ioWriteable(inputCheck, data);
			typename InputNgType::Readable io(ioWriteable);
			io.readline(begin, "OmegaBegin=");
			io.readline(step, "OmegaStep=");
			io.readline(total, "OmegaTotal=");
			io.readline(offset, "OmegaOffset=");
			io.readline(obs, "Observable=");
		}

		RealType begin;
		RealType step;
		SizeType offset;
		SizeType total;
		PsimagLite::String obs;
	};

	ManyOmegas(PsimagLite::String inputFile, RealType precision, const ApplicationType& app)
	    : inputfile_(inputFile), runner_(precision, app), omegaParams_(0)
	{
		InputNgType::Writeable::readFile(data_, inputFile);
		omegaParams_ = new OmegaParams(data_);
	}

	~ManyOmegas()
	{
		delete omegaParams_;
		omegaParams_ = nullptr;
	}

	void run(bool dryRun, PsimagLite::String root)
	{
		const PsimagLite::String obs = omegaParams_->obs;
		const PsimagLite::String insitu = "<gs|" + obs + "|P1>,<gs|" +
		        obs + "|P2>,<gs|" + obs + "|P3>";

		for (SizeType i = omegaParams_->offset; i < omegaParams_->total; ++i) {
			const RealType omega = i*omegaParams_->step + omegaParams_->begin;
			PsimagLite::String data2 = modifyOmega(omega);
			PsimagLite::String outputfile = "\nOutputFile=\"" + root + ttos(i) + "\";\n";
			data2 += outputfile;

			PsimagLite::String logfile = "runForinput" + ttos(i) + ".cout";

			if (dryRun) {
				std::cerr<<"ManyOmegas.h:: omega = "<<omega;
				std::cerr<<" output="<<outputfile;
				std::cerr<<" logfile="<<logfile<<" NOT done because -d\n";
				continue;
			}

			runner_.doOneRun(data2, insitu, logfile);
		}
	}

	PsimagLite::String modifyOmega(RealType omega) const
	{
		PsimagLite::String data = data_;
		const PsimagLite::String label = "$omega";
		const size_t pos = data.find(label);
		if (pos == PsimagLite::String::npos)
			err("Could not find " + label + " in " + inputfile_ + "\n");

		PsimagLite::String str = data.substr(0, pos);
		str += ttos(omega);
		const SizeType len2 = data.length();
		assert(len2 > pos + label.length());
		str += data.substr(pos + label.length(), len2 - pos - label.length());
		return str;
	}

	PsimagLite::String inputfile_;
	DmrgRunnerType runner_;
	OmegaParams* omegaParams_;
	PsimagLite::String data_;
};
}
#endif // MANYOMEGAS_H
