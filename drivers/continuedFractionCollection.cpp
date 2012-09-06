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
#include "TridiagonalMatrix.h"
#include "ContinuedFraction.h"
#include "ContinuedFractionCollection.h"

using namespace PsimagLite;
typedef double RealType;
typedef TridiagonalMatrix<RealType> TridiagonalMatrixType;
typedef ContinuedFraction<RealType,TridiagonalMatrixType> ContinuedFractionType;
typedef ContinuedFractionCollection<ContinuedFractionType>
	ContinuedFractionCollectionType;
typedef ContinuedFractionType::PlotParamsType PlotParamsType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -f file  -b omega1";
	std::cerr<<" -e omega2 -s omegaStep -d delta\n";
	std::cerr<<"Conditions: omega1<omega2 omegaStep>0 delta>0\n";
}

void plotAll(const ContinuedFractionCollectionType& cfCollection,
		const PlotParamsType& params)
{
	ContinuedFractionCollectionType::PlotDataType v;
	cfCollection.plot(v,params);
	std::cout.precision(12);
	for (size_t x=0;x<v.size();x++) {
		std::cout<<v[x].first<<" "<<std::imag(v[x].second);
		std::cout<<" "<<std::real(v[x].second)<<"\n";
	}
}

void plotOneByOne(const ContinuedFractionCollectionType& cfCollection,
		const PlotParamsType& params)
{
	
	std::cout.precision(12);
	for (size_t i=0;i<cfCollection.size();i++) {
		ContinuedFractionCollectionType::PlotDataType v;
		cfCollection.plotOne(i,v,params);
		for (size_t x=0;x<v.size();x++) {
			std::cout<<v[x].first<<" "<<std::imag(v[x].second);
			std::cout<<" "<<std::real(v[x].second)<<"\n";
		}
	}
}

int main(int argc,char *argv[])
{
	int opt = 0;
	std::string file="";
	RealType wbegin = 0;
	RealType wend = 0;
	RealType wstep = 0;
	RealType delta = 0;
	bool oneByOne = false;
	while ((opt = getopt(argc, argv,
		"f:b:e:s:d:1")) != -1) {
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
		case 'd':
			delta = atof(optarg);
			break;
		case '1':
			oneByOne = true;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	if (file=="" || wbegin>=wend || wstep<=0 || delta<=0) {
		usage(argv[0]);
		return 1;
	}

	IoSimple::In io(file);
	ContinuedFractionCollectionType cfCollection(io);
	PlotParamsType params(wbegin,wend,wstep,delta);
	if (!oneByOne) plotAll(cfCollection,params);
	else plotOneByOne(cfCollection,params);
}

