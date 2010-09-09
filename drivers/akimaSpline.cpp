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
#include "Vector.h"
#include "AkimaSpline.h"

typedef double FieldType;
typedef Vector<FieldType> VectorType;
typedef AkimaSpline<VectorType> AkimaSplineType;
typedef AkimaSplineType::IntervalType IntervalType;

void readTwoColumnData(const std::string& file,VectorType& v0,VectorType& v1)
{
	std::ifstream fin(file.c_str());
	if (!fin || !fin.good() || fin.bad()) throw
		std::runtime_error("Cannot open file\n");
	while(!fin.eof()) {
		std::string s;
		fin>>s;
		if (s[0]=='#') continue;
		FieldType x = std::atof(s.c_str());
		size_t size=v0.size();
		if (size>1 && x<v0[size-1]) break;
		v0.push_back(x);
		fin>>s;
		if (s[0]=='#') continue;
		v1.push_back(atof(s.c_str()));
	}
	fin.close();
}

int main(int argc,char *argv[])
{
	VectorType x,s;
	readTwoColumnData(argv[1],x,s);

	AkimaSplineType akimaSpline(x,s);
	FieldType xstart = std::atof(argv[2]);
	FieldType xend = std::atof(argv[3]);
	size_t total = std::atoi(argv[4]);
	FieldType xstep = (xend-xstart)/total;

	for (FieldType x=xstart;x<xend;x+=xstep) {
		std::cout<<x<<" "<<akimaSpline(x)<<"\n";
	}
}



