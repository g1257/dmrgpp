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

int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	ClonerType cloner("input507-45a.inp","inputRoot",".inp");
	LineChangerType lcl1("CorrectionVectorOmega=",0.5,0.0,"","");
	cloner.push(lcl1);
	LineChangerType lcl2("OutputFile=",1.0,0.0,"dataRoot",".txt");
	cloner.push(lcl2);
	
	RunnerType run("./dmrg","inputRoot",".inp");
	
	size_t total = 10;
	concurrency.loopCreate(total);
	size_t i = 0;
	while(concurrency.loop(i)) {
		cloner.createInputFile(i);
		run(i);
	}

}

