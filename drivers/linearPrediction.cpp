// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************

*/
// END LICENSE BLOCK

#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "IoSimple.h"
#include "LinearPrediction.h"

typedef double FieldType;
typedef PsimagLite::LinearPrediction<FieldType> LinearPredictionType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -f file -l label -p n\n";
}

int main(int argc,char *argv[])
{
	int opt = 0;
	std::string file="";
	std::string label="";
	size_t p = 0;
	while ((opt = getopt(argc, argv,
			"f:p:l:")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'l':
			label = optarg;
			break;
		case 'p':
			p = atoi(optarg);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	if (file=="" || label=="" || p==0) {
		usage(argv[0]);
		return 1;
	}

	PsimagLite::IoSimple::In io(file);
	std::vector<FieldType> y;
	io.read(y,label);
	size_t n = y.size();
	std::cout<<"#Found "<<n<<" points in file "<<file<<"\n";
	LinearPredictionType linearPrediction(y);
	linearPrediction.predict(p);
	for (size_t i=0;i<p+n;i++) {
		std::cout<<i<<" "<<linearPrediction(i)<<"\n";
	}
}



