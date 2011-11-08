// BEGIN LICENSE BLOCK
/*
Copyright (c) 2011 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

------------------------------------------------------------------
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
------------------------------------------------------------------

*/
// END LICENSE BLOCK

#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "KernelPolynomial.h"
#include "IoSimple.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -f file  -b omega1";
	std::cerr<<" -e omega2 -s omegaStep \n";
	std::cerr<<"Conditions: omega1<omega2 omegaStep>0 \n";
}

int main(int argc,char *argv[])
{
	int opt = 0;
	std::string file="";
	RealType wbegin = 0;
	RealType wend = 0;
	RealType wstep = 0;
	while ((opt = getopt(argc, argv,
		"f:b:e:s")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'b':
			wbegin = atof(optarg);
			break;
		case 'e':
			wend = atof(optarg);
			break;
		case 's':
			wstep = atof(optarg);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	if (file=="" || wbegin>=wend || wstep<=0) {
		usage(argv[0]);
		return 1;
	}

	PsimagLite::IoSimple::In io(file);
	
	// need to get momentsEven, momentsOdd, a, and b
	size_t n = momentsEven.size();
	std::vector<RealType> moments(2*n);
	for (size_t i=0;i<n;i++) {
		moments[2*i] = momentsEven[i];
		moments[2*i+1] = momentsOdd[i];
	}
	
	typedef PsimagLite::KernelPolynomial<RealType> KernelPolynomialType;
	KernelPolynomialType kernelPolynomial(moments,a,b);
	
	std::vector<RealType> v;
	kernelPolynomial.jackson(v,wbegin,wend,wstep);
	
	RealType omega = wbegin;
	size_t i = 0;
	while(omega<wend) {
		std::cout<<omega<<" "<<v[i++];
		omega += wstep;
	}
}

