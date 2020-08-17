#ifndef PROCOMEGAS_H
#define PROCOMEGAS_H
#include "OmegaParams.h"
#include "PsimagLite.h"
#include "OmegasFourier.h"

namespace Dmrg {

template<typename ComplexOrRealType>
class ProcOmegas {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::PsiApp ApplicationType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef OmegaParams<InputNgType, RealType> OmegaParamsType;
	typedef OmegasFourier<RealType> OmegasFourierType;

	static const SizeType MAX_LINE_SIZE = 2048;

	ProcOmegas(PsimagLite::String inputFile, SizeType precision, const ApplicationType& app)
	    : inputfile_(inputFile), omegaParams_(0), omegasFourier_()
	{
		// set precision here
		PsimagLite::String data;
		InputNgType::Writeable::readFile(data, inputFile);
		omegaParams_ = new OmegaParamsType(data);
		//const SizeType centralSite = getCentralSite();
	}

	~ProcOmegas()
	{
		delete omegaParams_;
		omegaParams_ = nullptr;
	}

	void run(bool skipFourier, PsimagLite::String rootname)
	{
		for (SizeType i = omegaParams_->offset; i < omegaParams_->total; ++i) {
			const RealType omega = i*omegaParams_->step + omegaParams_->begin;

			procCommon(i, omega, skipFourier, rootname);

			if (skipFourier) return; // <<===  EARLY EXIT HERE

//			if (array_.size() == 0)
//				err("procAllQs: array is empty\n");

//			printSpectrum(array_);
		}
	}

private:

	void procCommon(SizeType ind, RealType omega, bool skipFourier, PsimagLite::String rootname)
	{
		//my $n = $GlobalNumberOfSites;

		PsimagLite::String inFile("runFor");
		inFile += rootname + ttos(ind) + ".cout";
		VectorRealType values1;
		VectorRealType values2;
		//SizeType maxSite = correctionVectorRead(values1, values2, inFile);

		//print STDERR "$0: omega=$omega maxSite=$maxSite\n"; <== LOGFILEOUT

		writeSpaceValues(omega, values1, values2);

		if (skipFourier) return; // <<=== EARLY EXIT HERE

		VectorRealType qValues;
		omegasFourier_.fourier(qValues, values1, values2);
		//print LOGFILEOUT "$0: Number of k values ".scalar(@qValues)."\n";

		VectorRealType array;
		omegasFourier_.writeFourier(array, qValues);
	}

	void writeSpaceValues(RealType omega, VectorRealType& v1, VectorRealType& v2)
	{
		const SizeType n = v1.size();
		if (v2.size() != n)
			err("writeSpaceValues: v1.size != v2.size\n");

		printToSpaceOut(ttos(omega) + " " + ttos(n) + "\n");

		for (SizeType i = 0; i < n; ++i)
			printToSpaceOut(ttos(i) + ttos(v1[i]) + ttos(v2[i]) + "\n");
	}

	void printToSpaceOut(PsimagLite::String)
	{
		err("printToSpaceOut: unimplemented\n");
	}

	SizeType correctionVectorRead(VectorRealType& v1, VectorRealType& v2, PsimagLite::String inFile)
	{
		PsimagLite::String status("clear");
		SizeType maxSite = 0;
		std::ifstream fin(inFile);

		if (!fin || !fin.good() || fin.bad())
			err("correctionVectorRead: Cannot read " + inFile + "\n");

		VectorStringType labels{"P2", "P3"}; // ORDER IMPORTANT HERE!

		const SizeType ns = MAX_LINE_SIZE;
		char* ss = new char[ns];
		VectorStringType tokens;
		while (fin.getline(ss, ns)) {

			PsimagLite::String s(ss);

			if (s.find("PsiApp: CmdLine") != PsimagLite::String::npos)
				continue;

			bool skip = true;
			for (SizeType i = 0; i < labels.size(); ++i) {
				bool isGs = (s.find("gs") != PsimagLite::String::npos ||
				        s.find("X0") != PsimagLite::String::npos);

				if (s.find(labels[i]) == PsimagLite::String::npos || !isGs)
					continue;

				status = labels[i];
				skip = false;
			}

			if (skip) continue;

			PsimagLite::split(tokens, s, ",");
			if (tokens.size() != 5)
				err("correctionVectorRead: Not 5 tokens in line in " + inFile + "\n");

			SizeType site = tokens[0];
			SizeType c = 0;
			for (SizeType i = 0; i < labels.size(); ++i) {
				++c;
				if (status != labels[i])
					continue;
				if (c == 1)
					v1[site] = tokens[1];
				else if (c == 2)
					v2[site] = tokens[1];
				else
					err("correctionVectorRead: counter c wrong in " + inFile + "\n");
			}

			if (maxSite < site) maxSite = site;
			status = "clear";
		}

		delete[] ss;
		ss = 0;

		++maxSite;
		//print LOGFILEOUT "$0: correctionVectorRead maxsite= $maxSite\n";

		return maxSite;
	}

	PsimagLite::String inputfile_;
	OmegaParamsType* omegaParams_;
	OmegasFourierType omegasFourier_;
};
}
#endif // PROCOMEGAS_H
