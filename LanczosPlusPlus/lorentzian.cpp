#include "Sort.h"
#include "TypeToString.h"
#include "Vector.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <unistd.h>

typedef double                             RealType;
typedef std::complex<RealType>             ComplexType;
typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
typedef PsimagLite::Vector<RealType>::Type VectorRealType;

void load(VectorRealType& e, VectorRealType& w, PsimagLite::String file)
{
	std::ifstream fin(file.c_str());
	RealType      tmp = 0;
	while (!fin.eof() && !fin.bad() && fin.good()) {
		fin >> tmp;
		if (fin.eof() || fin.bad() || !fin.good())
			break;
		e.push_back(tmp);
		fin >> tmp;
		w.push_back(tmp);
	}

	std::cerr << "load: " << e.size() << " values found in " << file << "\n";
}

void sort(VectorRealType& e, VectorRealType& w)
{
	assert(e.size() == w.size());
	if (e.size() == 0)
		return;

	PsimagLite::Sort<VectorRealType> sort;
	VectorSizeType                   iperm(e.size());
	sort.sort(e, iperm);
	VectorRealType ww(w.size());
	for (SizeType i = 0; i < w.size(); ++i) {
		ww[i] = w[iperm[i]];
	}

	w = ww;
}

void prune(VectorRealType& e, VectorRealType& w, RealType& emin, RealType& emax, RealType& wabsmax)
{
	assert(e.size() == w.size());
	if (e.size() == 0)
		return;

	SizeType i = 0;
	for (; i < w.size(); ++i) {
		if (fabs(w[i]) > 1e-6)
			break;
	}

	if (i > 0)
		i--;

	int j = w.size() - 1;
	for (; j >= 0; j--) {
		if (fabs(w[j]) > 1e-6)
			break;
	}

	SizeType final = static_cast<SizeType>(j + 1);

	SizeType       counter = 0;
	SizeType       total   = (final - i);
	VectorRealType ee(total);
	VectorRealType ww(total);
	for (SizeType k = i; k < final; ++k) {
		ee[counter]   = e[k];
		ww[counter++] = w[k];
		if (e[k] > emax)
			emax = e[k];
		if (e[k] < emin)
			emin = e[k];
		if (fabs(w[k]) > wabsmax)
			wabsmax = fabs(w[k]);
	}

	e = ee;
	w = ww;
	std::cerr << "prun: emax=" << emax << " emin=" << emin << " wabsmax=" << wabsmax << "\n";
	std::cerr << "prune: " << e.size() << " values remain after pruning\n";
}

ComplexType lorentzian(ComplexType z, const VectorRealType& e, const VectorRealType& w)
{
	ComplexType sum = 0;
	for (SizeType i = 0; i < e.size(); ++i) {
		ComplexType tmp = (z - e[i]);
		sum += w[i] / tmp;
	}

	return sum;
}

ComplexType findOmega(SizeType           ind,
                      SizeType           total,
                      RealType           omegaStep,
                      RealType           omegaInit,
                      RealType           eps,
                      RealType           beta,
                      PsimagLite::String mode)
{
	if (mode == "real")
		return ComplexType(ind * omegaStep + omegaInit, eps);
	if (mode == "matsubara") {
		SizeType totalOver2 = static_cast<SizeType>(total * 0.5);
		assert(beta > 0);
		RealType factor = 2.0 * M_PI / beta;
		if (ind < totalOver2) {
			RealType tmp = (totalOver2 - ind);
			return ComplexType(eps, -factor * tmp);
		}

		RealType tmp = (1 + ind) - totalOver2;
		return ComplexType(eps, factor * tmp);
	}

	PsimagLite::String str(__FILE__);
	str += " " + ttos(__LINE__) + "\n";
	str += "findOmega: Uknown mode " + mode + "\n";
	throw PsimagLite::RuntimeError(str);
}

void usage(char* name, PsimagLite::String msg = "")
{
	if (msg != "")
		std::cerr << name << ": " << msg << "\n";
	std::cerr << "USAGE: " << name
	          << " -f file -t total -m mode [-e eps] [-b beta] [-s step] [-S start]\n";
	std::cerr << "\tmode is either real or matsubara\n";
	std::cerr << "\tbeta is mandatory in matsubara mode\n";
}

int main(int argc, char** argv)
{
	int                opt = 0;
	PsimagLite::String file;
	PsimagLite::String mode;
	RealType           eps      = 0.1;
	SizeType           total    = 0;
	RealType           beta     = 0.0;
	RealType           start    = 0;
	RealType           step     = 0;
	bool               hasStart = false;
	bool               hasStep  = false;
	while ((opt = getopt(argc, argv, "f:t:m:e:b:s:S:")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 't':
			total = atoi(optarg);
			break;
		case 'm':
			mode = optarg;
			break;
		case 'e':
			eps = atof(optarg);
			break;
		case 'b':
			beta = atof(optarg);
			break;
		case 's':
			step    = atof(optarg);
			hasStep = true;
			break;
		case 'S':
			start    = atof(optarg);
			hasStart = true;
			break;
		default: /* '?' */
			usage(argv[0]);
			return 1;
		}
	}

	if (file == "" || total == 0 || mode == "") {
		usage(argv[0]);
		return 2;
	}

	if (mode == "matsubara" && beta == 0) {
		usage(argv[0], "beta cannot be zero in matsubara mode");
		return 2;
	}

	// load
	VectorRealType e;
	VectorRealType w;
	load(e, w, file);
	// sort
	sort(e, w);
	// min, max, prune
	RealType emin    = 1e10;
	RealType emax    = -emin;
	RealType wabsmax = 0;
	prune(e, w, emin, emax, wabsmax);

	RealType omegaInit = (hasStart) ? start : emin;
	RealType omegaStep = (emax - omegaInit) / (total - 1);
	if (hasStep)
		omegaStep = step;
	RealType factor = 1.0 / wabsmax;

	for (SizeType i = 0; i < total; ++i) {
		ComplexType z     = findOmega(i, total, omegaStep, omegaInit, eps, beta, mode);
		RealType    omega = (mode == "real") ? std::real(z) : std::imag(z);
		ComplexType val   = lorentzian(z, e, w);
		val *= factor;
		std::cout << omega << " " << std::real(val) << " " << std::imag(val) << "\n";
	}
}
