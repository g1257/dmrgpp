#ifndef PROCOMEGAS_H
#define PROCOMEGAS_H
#include "OmegaParams.h"
#include "OmegasFourier.h"
#include "PsimagLite.h"

/* Limitations and Missing Features
 *
 * These may or maynot be implemented in the future
 *
 * (1) Ainur only
 *
 * (2) Cheby not supported
 *
 * (3) Fourier: Many geometries unsupported yet (see OmegasFourier.h)
 *
 * (4) Gnuplot output not supported yet
 *
 */
namespace Dmrg
{

template <typename ComplexOrRealType, typename OmegaParamsType>
class ProcOmegas
{

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef PsimagLite::PsiApp ApplicationType;
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef OmegasFourier<ComplexOrRealType, InputNgType::Readable> OmegasFourierType;
	typedef typename OmegasFourierType::VectorComplexType VectorComplexType;

	static const SizeType MAX_LINE_SIZE = 409600;

	class Qdata
	{

	public:

		Qdata(SizeType n)
		    : data_(n, nullptr)
		{
		}

		Qdata() { }

		void resize(SizeType n)
		{
			if (data_.size() != 0)
				err("Qdata::resize only supported for empty objects\n");
			data_.resize(n, nullptr);
		}

		void set(SizeType ind, const VectorComplexType& v)
		{
			if (data_.size() <= ind)
				err("Qdata::set out of bounds\n");

			if (data_[ind] != nullptr)
				err("Qdata::set already set\n");

			data_[ind] = new VectorComplexType(v);
		}

		const VectorComplexType& get(SizeType ind) const
		{
			if (data_.size() <= ind)
				err("Qdata::get out of bounds\n");

			if (data_[ind] == nullptr)
				err("Qdata::get empty location\n");

			return *data_[ind];
		}

	private:

		typename PsimagLite::Vector<VectorComplexType*>::Type data_;
	};

	ProcOmegas(typename InputNgType::Readable& io,
	    SizeType precision,
	    bool skipFourier,
	    PsimagLite::String rootIname,
	    PsimagLite::String rootOname,
	    const OmegaParamsType& omegaParams)
	    : rootIname_(rootIname)
	    , rootOname_(rootOname)
	    , omegaParams_(omegaParams)
	    , omegasFourier_(skipFourier, io)
	    , numberOfSites_(0)
	{
		// set precision here FIXME TODO
		io.readline(numberOfSites_, "TotalNumberOfSites");
		// const SizeType centralSite = getCentralSite();
	}

	void run()
	{
		VectorRealType values1(numberOfSites_);
		VectorRealType values2(numberOfSites_);
		VectorBoolType defined(numberOfSites_);

		std::ofstream* fout = nullptr;

		if (rootOname_ != "")
			fout = new std::ofstream(rootOname_);
		if (!fout || !*fout || fout->bad() || !fout->good())
			err("writeSpaceValues: Cannot write to " + rootOname_ + "\n");

		qData_.resize(omegaParams_.total() - omegaParams_.offset());

		for (SizeType i = omegaParams_.offset(); i < omegaParams_.total(); ++i) {
			const RealType omega = omegaParams_.omega(i);

			procCommon(i, omega, values1, values2, defined, fout);
			qData_.set(i - omegaParams_.offset(), omegasFourier_.data());
		}

		if (fout)
			fout->close();

		delete fout;
		fout = nullptr;
	}

	void printPgfplots(PsimagLite::String foutname)
	{
		std::ofstream fout(foutname);
		if (!fout || fout.bad() || !fout.good())
			err("writeSpaceValues: Cannot write to " + foutname + "\n");

		SizeType numberOfQs = 0;
		for (SizeType i = omegaParams_.offset(); i < omegaParams_.total(); ++i) {
			const RealType omega = omegaParams_.omega(i);
			const VectorComplexType& v = qData_.get(i - omegaParams_.offset());

			if (i == omegaParams_.offset()) {
				assert(numberOfQs == 0);
				numberOfQs = v.size();
			} else if (numberOfQs != v.size()) {
				err("INTERNAL ERROR: Omega set with non equal number of q points\n");
			}

			for (SizeType m = 0; m < numberOfQs; ++m) {
				RealType q = omegasFourier_.q(m);
				fout << q << " " << omega << " " << PsimagLite::imag(v[m]) << "\n";
			}
		}
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

		// print STDERR "$0: omega=$omega maxSite=$maxSite\n"; <== LOGFILEOUT

		if (fout)
			writeSpaceValues(*fout, omega, values1, values2);

		omegasFourier_.fourier(values1, values2);
		// print LOGFILEOUT "$0: Number of k values ".scalar(@qValues)."\n";
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
			    ttos(i) + " " + ttos(v1[i]) + " " + ttos(v2[i]) + "\n");
	}

	static void printToSpaceOut(std::ofstream& fout, PsimagLite::String str)
	{
		fout << str;
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

		VectorStringType labels { "P2", "P3" }; // ORDER IMPORTANT HERE!

		const SizeType ns = MAX_LINE_SIZE;
		char* ss = new char[ns];
		std::fill(defined.begin(), defined.end(), false);
		while (fin.getline(ss, ns)) {

			PsimagLite::String s(ss);

			if (s.find("PsiApp: CmdLine") != PsimagLite::String::npos)
				continue;

			bool skip = true;
			for (SizeType i = 0; i < labels.size(); ++i) {
				bool isGs = (s.find("gs") != PsimagLite::String::npos || s.find("X0") != PsimagLite::String::npos);

				if (s.find(labels[i]) == PsimagLite::String::npos || !isGs)
					continue;

				status = labels[i];
				skip = false;
			}

			if (skip)
				continue;

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
		checkSites(defined, inFile);
		// print LOGFILEOUT "$0: correctionVectorRead maxsite= $maxSite\n";
	}

	void checkSites(const VectorBoolType& defined,
	    PsimagLite::String inFile)
	{
		const SizeType n = defined.size();
		for (SizeType i = 0; i < n; ++i)
			if (!defined[i])
				err("Undefined value for site= " + ttos(i) + " file= " + inFile + "\n");
	}

	PsimagLite::String inputfile_;
	PsimagLite::String rootIname_;
	PsimagLite::String rootOname_;
	const OmegaParamsType& omegaParams_;
	OmegasFourierType omegasFourier_;
	SizeType numberOfSites_;
	Qdata qData_;
};
}
#endif // PROCOMEGAS_H
