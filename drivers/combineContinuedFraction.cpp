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
#include "String.h"

using namespace PsimagLite;
typedef double RealType;
typedef TridiagonalMatrix<RealType> TridiagonalMatrixType;
typedef ContinuedFraction<TridiagonalMatrixType> ContinuedFractionType;
typedef ContinuedFractionCollection<ContinuedFractionType> ContinuedFractionCollectionType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" file1 file2";
}

int main(int argc,char *argv[])
{
	if (argc<2) {
		usage(argv[0]);
		return 1;
	}

	ContinuedFractionCollectionType cfCollection(FREQ_REAL);

	String s = "#Avector";
	for (int x = 1;x<argc;x++) {
		IoSimple::In io(argv[x]);
		io.advance(s,IoSimple::In::LAST_INSTANCE);
		ContinuedFractionType cf(io);
		cfCollection.push(cf);
	}

	IoSimple::Out ioOut(std::cout);
	ioOut.setPrecision(12);
	cfCollection.save(ioOut);

}

