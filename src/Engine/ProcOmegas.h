#ifndef PROCOMEGAS_H
#define PROCOMEGAS_H
#include "OmegaParams.h"
#include "PsimagLite.h"
#include "OmegasFourier.h"

namespace Dmrg {

template<typename ComplexOrRealType>
class ProcOmegas {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::PsiApp ApplicationType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;

	ProcOmegas(PsimagLite::String inputFile, RealType precision, const ApplicationType& app)
	    : inputfile_(inputFile), runner_(precision, app), omegaParams_(0)
	{
		InputNgType::Writeable::readFile(data_, inputFile);
		omegaParams_ = new OmegaParams(data_);
		//const SizeType centralSite = getCentralSite();
	}

	~ProcOmegas()
	{
		delete omegaParams_;
		omegaParams_ = nullptr;
	}

	void main()
	{
		for (SizeType i = omegaParams_->offset; i < omegaParams_->total; ++i) {
			const RealType omega = i*omegaParams_->step + omegaParams_->begin;

			procCommon(ind, omega);

			if (noFourier_) return; // <<===  EARLY EXIT HERE

			if (array_.size() == 0)
				err("procAllQs: array is empty\n");

			printSpectrum(array_);
		}
	}

	void procCommon(SizeType ind, RealType omega)
	{
		//my $n = $GlobalNumberOfSites;

		PsimagLite::String inFile("runFor");
		inFile += inputRoot_ + ttos(ind) + ".cout";
		SizeType maxSite = correctionVectorRead(values1, values2, inFile);

		//print STDERR "$0: omega=$omega maxSite=$maxSite\n"; <== LOGFILEOUT

		writeSpaceValues(omega, values1, value2);

		if (noFourier_) return; // <<=== EARLY EXIT HERE

		OmegaUtils::fourier(qValues, values1, value2);
		//print LOGFILEOUT "$0: Number of k values ".scalar(@qValues)."\n";
		OmegaUtils::writeFourier(array, qValues);
	}

	void writeSpaceValues(RealType omega, VectorRealType& v1, VectorRealType& v2)
	{
		const SizeType n = array.size();
		printToSpaceOut(ttos(omega) + " " + ttos(n) + "\n");
		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String vv1 = v1[i];
			PsimagLite::String vv2 = v2[i];
			printToSpaceOut(ttos(i) + vv1  + vv2 + "\n");
		}
	}

	SizeType correctionVectorRead(VectorRealType& v1, VectorRealType& v2, PsimagLite::String inFile)
	{
		PsimagLite::String status("clear");
		SizeType maxSite = 0;
		std::ifstream fin(inFile);

		if (!fin || !fin.good() || fin.bad())
			err("correctionVectorRead: Cannot read " + inFile + "\n");

		VectorStringType labels{"P2", "P3"}; // ORDER IMPORTANT HERE!

		while (fin.getline(s, ns)) {

			if (s.substr("PsiApp\: +CmdLine/") != PsimagLite::String::npos)
				continue;

			bool skip = true;
			for (SizeType i = 0; i < labels.size(); ++i) {
				bool isGs = (s.substr("gs") != PsimagLite::String::npos ||
				        s.substr("X0") != PsimagLite::String::npos);

				if (s.substr(labels[i]) == PsimagLite::String::npos || !isGs)
					continue;

				status = labels[i];
				skip = false;
			}

			if (skip) continue;

			PsimagLite::split(tokens, s, ",");
			if (tokens.size() != 5)
				err("correctionVectorRead: Not 5 tokens in line in " + inFile + "\n");

			SizeType site = tokens[0];
			SizeType time = 0;
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

		++maxSite;
		//print LOGFILEOUT "$0: correctionVectorRead maxsite= $maxSite\n";

		return maxSite;
	}
};
}
#endif // PROCOMEGAS_H
