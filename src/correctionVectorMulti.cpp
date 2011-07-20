typedef double RealType;

#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<RealType> ConcurrencyType;
#endif

#include "LineChangerLinear.h"
#include "Cloner.h"
#include "Runner.h"

typedef Dmrg::LineChangerLinear<RealType> LineChangerType;
typedef Dmrg::Cloner<LineChangerType> ClonerType;
typedef Dmrg::Runner RunnerType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -i infile  -d dataRoot";
	std::cerr<<" -o inputRoot\n";
}

int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	int opt;
	std::string infile = "";
	std::string dataRoot = "";
	std::string inputRoot = "";
	
	while ((opt = getopt(argc, argv,
		"i:d:o:")) != -1) {
		switch (opt) {
		case 'i':
			infile = optarg;
			break;
		case 'd':
			dataRoot = optarg;
			break;
		case 'o':
			inputRoot = optarg;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	if (infile=="" || dataRoot=="" || inputRoot == "") {
		usage(argv[0]);
		return 1;
	}
	size_t total = concurrency.nprocs();

	ClonerType cloner(infile,inputRoot,".inp");
	LineChangerType lcl1("CorrectionVectorOmega=",0.5,0.0,"","");
	cloner.push(lcl1);
	LineChangerType lcl2("OutputFile=",1.0,0.0,dataRoot,".txt");
	cloner.push(lcl2);

	RunnerType run("./dmrg",inputRoot,".inp");

	concurrency.loopCreate(total);
	size_t i = 0;
	while(concurrency.loop(i)) {
		cloner.createInputFile(i);
		run(i);
	}
}
