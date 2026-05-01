#include "Engine/OneSector.h"
#include "Io/IoSimple.h"
#include "PsimagLite.h"

typedef double                                          RealType;
typedef PsimagLite::IoSimple::In                        InputType;
typedef LanczosPlusPlus::OneSector<RealType, InputType> OneSectorType;
typedef PsimagLite::Vector<OneSectorType*>::Type        VectorOneSectorType;
typedef PsimagLite::Vector<RealType>::Type              VectorRealType;
typedef OneSectorType::VectorSizeType                   VectorSizeType;
typedef OneSectorType::MatrixType                       MatrixType;

struct ThermalOptions {
	ThermalOptions(PsimagLite::String    operatorName_,
	               RealType              beta_,
	               RealType              mu_,
	               RealType              constant_,
	               const VectorSizeType& sites_)
	    : operatorName(operatorName_)
	    , beta(beta_)
	    , mu(mu_)
	    , constant(constant_)
	    , sites(sites_)
	{ }

	PsimagLite::String operatorName;
	RealType           beta;
	RealType           mu;
	RealType           constant;
	VectorSizeType     sites;
};

// Compute X^(s,s')_{n,n'} = \sum_{t,t'}U^{s*}_{n,t}A_{t,t'}^(s,s')U^{s'}_{t',n'}
void computeX(MatrixType&          x,
              const MatrixType&    a,
              const OneSectorType& sectorSrc,
              const OneSectorType& sectorDest)
{
	MatrixType tmp(x.n_row(), x.n_col());
	sectorDest.multiplyRight(tmp, a);
	sectorSrc.multiplyLeft(x, tmp);
}

SizeType findJnd(const VectorOneSectorType& sectors, const VectorSizeType& jndVector)
{
	for (SizeType i = 0; i < sectors.size(); ++i) {
		if (sectors[i]->isSector(jndVector))
			return i;
	}

	PsimagLite::String str(__FILE__);
	str += " " + ttos(__LINE__) + "\n";
	str += "findJnd can't find sector index\n";
	throw PsimagLite::RuntimeError(str);
}

void findOperatorAndMatrix(MatrixType&                a,
                           SizeType&                  jnd,
                           SizeType                   siteIndex,
                           SizeType                   ind,
                           const ThermalOptions&      opt,
                           const VectorOneSectorType& sectors,
                           InputType&                 io)
{
	if (opt.operatorName == "i") {
		jnd        = ind;
		SizeType n = sectors[ind]->size();
		a.resize(n, n);
		a.setTo(0);
		for (SizeType i = 0; i < n; ++i)
			a(i, i) = 1.0;

		return;
	} else if (opt.operatorName != "c") {
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "findOperatorAndMatrix: unknown operator " + opt.operatorName + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	SizeType spin = 0;
	assert(opt.sites.size() > siteIndex);
	SizeType           site  = opt.sites[siteIndex];
	PsimagLite::String label = "#Operator_c_";
	label += ttos(spin) + "_" + ttos(site);
	io.advance(label, 0);
	VectorSizeType jndVector;
	io.read(jndVector, "#SectorDest");
	if (jndVector.size() == 0)
		return;
	jnd = findJnd(sectors, jndVector);

	io.read(a, "#Matrix");
}

RealType computePartialZ(SizeType                   ind,
                         const ThermalOptions&      opt,
                         const VectorOneSectorType& sectors,
                         RealType                   factor)
{
	SizeType n   = sectors[ind]->size();
	RealType sum = 0.0;
	for (SizeType i = 0; i < n; ++i) {
		RealType e1  = sectors[ind]->eig(i);
		RealType arg = opt.beta * (factor - e1);
		sum += exp(arg);
	}

	return sum;
}

RealType computePartialE(SizeType                   ind,
                         const ThermalOptions&      opt,
                         const VectorOneSectorType& sectors,
                         RealType                   factor)
{
	SizeType n   = sectors[ind]->size();
	RealType sum = 0.0;
	for (SizeType i = 0; i < n; ++i) {
		RealType e1  = sectors[ind]->eig(i);
		RealType arg = opt.beta * (factor - e1);
		sum += exp(arg) * e1;
	}

	return sum;
}

RealType computeThisSector(SizeType                   ind,
                           const ThermalOptions&      opt,
                           const VectorOneSectorType& sectors,
                           InputType&                 io,
                           RealType                   factor,
                           RealType                   zInverse)
{
	SizeType   jnd = 0;
	MatrixType a;
	findOperatorAndMatrix(a, jnd, 0, ind, opt, sectors, io);

	SizeType n = a.n_row();
	if (n == 0)
		return 0.0;
	assert(n == sectors[ind]->size());

	SizeType m = a.n_col();
	assert(m == sectors[jnd]->size());
	if (m == 0)
		return 0.0;

	// Read operator 2   --> B
	MatrixType b;
	assert(opt.sites.size() == 2);
	if (opt.sites[0] == opt.sites[1]) {
		b = a;
	} else {
		assert(opt.sites[0] < opt.sites[1]);
		SizeType knd = 0;
		findOperatorAndMatrix(b, knd, 1, ind, opt, sectors, io);

		if (jnd != knd) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "computeThisSector: too many destination sectors\n";
			throw PsimagLite::RuntimeError(str);
		}
	}

	// Compute X^(s,s')_{n,n'} = \sum_{t,t'}U^{s*}_{n,t}A_{t,t'}^(s,s')U^{s'}_{t',n'}
	MatrixType x(n, m);
	computeX(x, a, *(sectors[ind]), *(sectors[jnd]));
	// Same for Y from B
	MatrixType y(n, m);
	computeX(y, b, *(sectors[ind]), *(sectors[jnd]));
	// Result is
	// \sum_{n,n'} X_{n,n'} Y_{n',n} exp(-i(E_n'-E_n)t))exp(-beta E_n)
	RealType sum     = 0.0;
	SizeType counter = 0;
	for (SizeType i = 0; i < n; ++i) {
		for (SizeType j = 0; j < m; ++j) {
			RealType e1  = sectors[ind]->eig(i);
			RealType e2  = sectors[jnd]->eig(j);
			RealType arg = opt.beta * (factor - e1);
			RealType val = x(i, j) * PsimagLite::conj(y(i, j)) * exp(arg) * zInverse;
			if (opt.operatorName != "i" && fabs(val) > 1e-12) {
				std::cout << (e1 - e2 + opt.mu) << " " << val << "\n";
				counter++;
			}

			sum += val;
		}
	}

	std::cerr << "Sector " << ind << " found " << counter << " values, sum=" << sum << "\n";
	return sum;
}

void computeAverageFor(const ThermalOptions& opt, const VectorOneSectorType& sectors, InputType& io)
{
	ThermalOptions optZ = opt;
	optZ.operatorName   = "i";
	VectorSizeType nupAndDown;
	VectorRealType muFactors(sectors.size(), 0);

	RealType zPartition = 0.0;
	RealType numerator  = 0.0;
	RealType energy     = 0.0;
	for (SizeType i = 0; i < sectors.size(); ++i) {
		io.read(nupAndDown, "#SectorSource");
		if (nupAndDown.size() != 2) {
			throw PsimagLite::RuntimeError("#SectorSource\n");
		}

		muFactors[i] = opt.mu * (nupAndDown[0] + nupAndDown[1]) + opt.constant;
		RealType tmp = computePartialZ(i, optZ, sectors, muFactors[i]);
		numerator += tmp * (nupAndDown[0] + nupAndDown[1]);
		energy += computePartialE(i, optZ, sectors, muFactors[i]);
		zPartition += tmp;
	}

	RealType zInverse = 1.0 / zPartition;
	std::cerr << "density=" << (numerator * zInverse) << " zPartition=" << zPartition << "\n";
	std::cerr << "energy=" << (energy * zInverse) << " zPartition=" << zPartition << "\n";

	if (opt.sites.size() < 2)
		return;

	io.rewind();
	RealType sum = 0.0;
	for (SizeType i = 0; i < sectors.size(); ++i) {
		sum += computeThisSector(i, opt, sectors, io, muFactors[i], zInverse);
	}

	std::cerr << "operator=" << opt.operatorName;
	std::cerr << " beta=" << opt.beta << " mu=" << opt.mu;
	std::cerr << " partition=" << zPartition << " sum=" << sum << "\n";
}

void usage(char* name, PsimagLite::String msg = "")
{
	if (msg != "")
		std::cerr << name << ": " << msg << "\n";
	std::cerr << "USAGE: " << name << " -f file -c operator -b beta ";
	std::cerr << " -s site1[,site2] [-m mu] [-C constant]\n";
}

int main(int argc, char** argv)
{
	int                                          opt = 0;
	PsimagLite::String                           operatorName;
	PsimagLite::String                           file;
	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	VectorSizeType                               sites(2, 0);
	RealType                                     beta     = 0;
	RealType                                     mu       = 0;
	RealType                                     constant = 0;

	while ((opt = getopt(argc, argv, "f:c:b:s:m:C:")) != -1) {
		switch (opt) {
		case 'c':
			operatorName = optarg;
			break;
		case 'f':
			file = optarg;
			break;
		case 'b':
			beta = atof(optarg);
			break;
		case 's':
			PsimagLite::split(tokens, optarg, ",");
			break;
		case 'm':
			mu = atof(optarg);
			break;
		case 'C':
			constant = atof(optarg);
			break;
		default: /* '?' */
			usage(argv[0]);
			return 1;
		}
	}

	if (file == "" || operatorName == "") {
		usage(argv[0]);
		return 2;
	}

	if (tokens.size() > sites.size()) {
		usage(argv[0], "Too many sites");
		return 3;
	}

	if (tokens.size() == 0)
		sites.clear();

	for (SizeType i = 0; i < tokens.size(); ++i)
		sites[i] = atoi(tokens[i].c_str());

	if (sites.size() > 0 && sites[1] < sites[0]) {
		usage(argv[0], "site1 must be smaller than site2");
	}

	InputType io(file);
	SizeType  total = 0;
	io.readline(total, "#TotalSectors=");
	VectorOneSectorType sectors(total);

	for (SizeType i = 0; i < sectors.size(); ++i) {
		sectors[i] = new OneSectorType(io);
		// sectors[i]->info(std::cout);
	}

	ThermalOptions options(operatorName, beta, mu, constant, sites);
	io.rewind();
	computeAverageFor(options, sectors, io);

	for (SizeType i = 0; i < sectors.size(); ++i)
		delete sectors[i];
}
