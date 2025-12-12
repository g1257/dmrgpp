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

#include "AkimaSpline.h"
#include "Vector.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using ComplexOrRealType = double;
using VectorType = std::vector<ComplexOrRealType>;
using AkimaSplineType = PsimagLite::AkimaSpline<VectorType>;

void readTwoColumnData(const std::string& file, VectorType& v0, VectorType& v1)
{
	std::ifstream fin(file.c_str());
	if (!fin || !fin.good() || fin.bad())
		throw std::runtime_error("Cannot open file\n");
	while (!fin.eof()) {
		std::string s;
		fin >> s;
		if (s[0] == '#')
			continue;
		ComplexOrRealType x = std::atof(s.c_str());
		SizeType size = v0.size();
		if (size > 1 && x < v0[size - 1])
			break;
		v0.push_back(x);
		fin >> s;
		if (s[0] == '#')
			continue;
		v1.push_back(atof(s.c_str()));
	}
	fin.close();
}

int main(int argc, char* argv[])
{
	VectorType x, s;
	readTwoColumnData(argv[1], x, s);

	AkimaSplineType akimaSpline(x, s);
	ComplexOrRealType xstart = std::atof(argv[2]);
	ComplexOrRealType xend = std::atof(argv[3]);
	SizeType total = std::atoi(argv[4]);
	ComplexOrRealType xstep = (xend - xstart) / total;

	for (ComplexOrRealType x = xstart; x < xend; x += xstep) {
		std::cout << x << " " << akimaSpline(x) << "\n";
	}
}
