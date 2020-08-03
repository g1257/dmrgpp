#ifndef MANYOMEGAS_H
#define MANYOMEGAS_H
#include "DmrgRunner.h"
#include "PsiBase64.h"
#include "InputNg.h"
#include "LanczosSolver.h"
#include "Vector.h"
#include "ParamsDmftSolver.h"
namespace Dmft {

template<typename ComplexOrRealType>
class ManyOmegas {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef DmrgRunner<ComplexOrRealType> DmrgRunnerType;
	typedef typename DmrgRunnerType::InputNgType InputNgType;

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
		}

		RealType begin;
		RealType step;
		SizeType offset;
		SizeType total;
	};

	ManyOmegas(PsimagLite::String inputFile,
	           RealType precision,
	           bool echoInput)
	    : inputfile_(inputFile), runner_(precision), echoInput_(echoInput), omegaParams_(0)
	{
		InputNgType::Writeable::readFile(data_, inputFile);
		omegaParams_ = new OmegaParams(data_);
	}

	~ManyOmegas()
	{
		delete omegaParams_;
		omegaParams_ = nullptr;
	}

	void run(bool dryRun)
	{
		for (SizeType i = omegaParams_->offset; i < omegaParams_->total; ++i) {
			RealType omega = i*omegaParams_->step + omegaParams_->begin;
			PsimagLite::String data2 = modifyOmega(omega);
			PsimagLite::String sOptions = "";

			if (dryRun) {
				std::cerr<<"ManyOmegas.h:: omega = "<<omega<<" NOT done because -d\n";
				continue;
			}

			runner_.doOneRun(data2, sOptions);
		}
	}

	PsimagLite::String modifyOmega(RealType omega) const
	{
		PsimagLite::String data = data_;
		const PsimagLite::String label = "$omega";
		size_t pos = data.find(label);
		if (pos == PsimagLite::String::npos)
			err("Could not find " + label + " in " + inputfile_ + "\n");

		PsimagLite::String str = data.substr(0, pos);
		SizeType len = str.length();
		str += ttos(omega);
		str += data.substr(pos + len, data.length() - len - label.length());
		return str;
	}

	PsimagLite::String inputfile_;
	DmrgRunnerType runner_;
	bool echoInput_;
	OmegaParams* omegaParams_;
	PsimagLite::String data_;
};
}
#endif // MANYOMEGAS_H
