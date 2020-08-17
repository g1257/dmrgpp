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
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef PsimagLite::PsiApp ApplicationType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef OmegaParams<InputNgType, RealType> OmegaParamsType;
	typedef OmegasFourier<ComplexOrRealType, InputNgType::Readable> OmegasFourierType;

	static const SizeType MAX_LINE_SIZE = 2048;

	ProcOmegas(typename InputNgType::Readable& io,
	           SizeType precision,
	           bool skipFourier,
	           PsimagLite::String rootIname,
	           PsimagLite::String rootOname)
	    : rootIname_(rootIname),
	      rootOname_(rootOname),
	      omegaParams_(io),
	      omegasFourier_(skipFourier, io),
	      numberOfSites_(0)
	{
		// set precision here FIXME TODO
		io.readline(numberOfSites_, "TotalNumberOfSites");
		//const SizeType centralSite = getCentralSite();
	}

	void run()
	{
		VectorRealType values1(numberOfSites_);
		VectorRealType values2(numberOfSites_);
		VectorBoolType defined(numberOfSites_);

		std::ofstream* fout = nullptr;

		if (rootOname_ != "") fout = new std::ofstream(rootOname_);
		if (!fout || !*fout || fout->bad() || !fout->good())
			err("writeSpaceValues: Cannot write to " + rootOname_ + "\n");

		for (SizeType i = omegaParams_.offset; i < omegaParams_.total; ++i) {
			const RealType omega = i*omegaParams_.step + omegaParams_.begin;

			procCommon(i, omega, values1, values2, defined, fout);
		}

		if (fout)
			fout->close();

		delete fout;
		fout = nullptr;
	}

private:

	void procCommon(SizeType ind,
	                RealType omega,
	                VectorRealType& values1,
	                VectorRealType& values2,
	                VectorBoolType& defined,
	                std::ofstream* fout)
	{
		PsimagLite::String inFile("runFor");
		inFile += rootIname_ + ttos(ind) + ".cout";

		correctionVectorRead(values1, values2, defined, inFile);

		//print STDERR "$0: omega=$omega maxSite=$maxSite\n"; <== LOGFILEOUT

		if (fout)
			writeSpaceValues(*fout, omega, values1, values2);

		omegasFourier_.fourier(values1, values2);
		//print LOGFILEOUT "$0: Number of k values ".scalar(@qValues)."\n";
	}

	void writeSpaceValues(std::ofstream& fout,
	                      RealType omega,
	                      const VectorRealType& v1,
	                      const VectorRealType& v2)
	{
		const SizeType n = v1.size();
		if (v2.size() != n)
			err("writeSpaceValues: v1.size != v2.size\n");

		printToSpaceOut(fout, ttos(omega) + " " + ttos(n) + "\n");

		for (SizeType i = 0; i < n; ++i)
			printToSpaceOut(fout,
			                ttos(i) + " " +
			                ttos(v1[i]) + " " +
			                ttos(v2[i]) + "\n");
	}

	static void printToSpaceOut(std::ofstream& fout, PsimagLite::String str)
	{
		fout<<str;
	}

	void correctionVectorRead(VectorRealType& v1,
	                          VectorRealType& v2,
	                          VectorBoolType& defined,
	                          PsimagLite::String inFile)
	{
		PsimagLite::String status("clear");
		std::ifstream fin(inFile);

		if (!fin || !fin.good() || fin.bad())
			err("correctionVectorRead: Cannot read " + inFile + "\n");

		VectorStringType labels{"P2", "P3"}; // ORDER IMPORTANT HERE!

		const SizeType ns = MAX_LINE_SIZE;
		char* ss = new char[ns];
		std::fill(defined.begin(), defined.end(), false);
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

			VectorStringType tokens;
			PsimagLite::split(tokens, s, " ");
			if (tokens.size() != 5)
				err("correctionVectorRead: Not 5 tokens in line " + s + "\nFile= " + inFile + "\n");

			SizeType site = PsimagLite::atoi(tokens[0]);
			SizeType c = 0;
			for (SizeType i = 0; i < labels.size(); ++i) {
				++c;
				if (status != labels[i])
					continue;

				if (site >= numberOfSites_)
					err("correctionVectorRead: Site " + ttos(site) + " is too big\n");

				if (c == 1)
					v1[site] = PsimagLite::atof(tokens[1]);
				else if (c == 2)
					v2[site] = PsimagLite::atof(tokens[1]);
				else
					err("correctionVectorRead: counter c wrong in " + inFile + "\n");

				defined[site] = true;
			}

			status = "clear";
		}

		delete[] ss;
		ss = 0;
		checkSites(defined);
		//print LOGFILEOUT "$0: correctionVectorRead maxsite= $maxSite\n";
	}

	void checkSites(const VectorBoolType& defined)
	{
		const SizeType n = defined.size();
		for (SizeType i = 0; i < n; ++i)
			if (!defined[i])
				err("Undefined value for site= " + ttos(i) + "\n");
	}

	PsimagLite::String inputfile_;
	PsimagLite::String rootIname_;
	PsimagLite::String rootOname_;
	OmegaParamsType omegaParams_;
	OmegasFourierType omegasFourier_;
	SizeType numberOfSites_;
};
}
#endif // PROCOMEGAS_H
