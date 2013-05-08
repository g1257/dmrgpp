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

#include "LanczosSolver.h"
#include "DavidsonSolver.h"
#include "CrsMatrix.h"
#include "Random48.h"

using namespace PsimagLite;

typedef double RealType;
typedef double ComplexOrRealType;

struct SolverParameters {
	typedef  ::RealType RealType;

	SolverParameters()
		: steps(200),
		  tolerance(1e-10),
		  stepsForEnergyConvergence(100),
		  lotaMemory(true),
		  options("")
	{}

	size_t steps;
	RealType tolerance;
	size_t stepsForEnergyConvergence;
	bool lotaMemory;
	std::string options;

};

typedef Vector<ComplexOrRealType>::Type VectorType;
typedef CrsMatrix<ComplexOrRealType> SparseMatrixType;
typedef LanczosOrDavidsonBase<SolverParameters,SparseMatrixType,VectorType> SparseSolverType;
typedef LanczosSolver<SolverParameters,SparseMatrixType,VectorType> LanczosSolverType;
typedef DavidsonSolver<SolverParameters,SparseMatrixType,VectorType> DavidsonSolverType;

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -n rank [-x] [-d ] [-c max_columns] [-m max_value] [-r seed]\n";
	exit(1);
}

int main(int argc,char *argv[])
{
	int opt = 0;
	bool useDavidson = false;
	size_t n = 0;
	RealType maxValue = 0;
	size_t maxCol = 0;
	size_t seed = 0;
	bool lotaMemory = true;

	while ((opt = getopt(argc, argv,
		"n:c:m:r:dx")) != -1) {
		switch (opt) {
		case 'd':
			useDavidson=true;
			break;
		case 'n' :
			n = atoi(optarg);
			break;
		case 'c' :
			maxCol = atoi(optarg);
			break;
		case 'm' :
			maxValue = atof(optarg);
			break;
		case 'r':
			seed = atoi(optarg);
			break;
		case 'x':
			lotaMemory=false;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	// sanity checks
	if (n==0) usage(argv[0]);
	if (maxCol==0) maxCol = 1 + size_t(0.1*n);
	if (std::norm(maxValue)<1e-6) maxValue = 1.0;
	if (seed==0) seed = 3443331;

	// create a random matrix:
	Random48<RealType> random(seed);
	SparseMatrixType sparse(n,n);
	Vector<bool>::Type seenThisColumn(n);
	size_t counter = 0;
	for (size_t i=0;i<n;i++) {
		sparse.setRow(i,counter);
		// random vector:
		size_t x = 1+size_t(random()*maxCol);
		for (size_t j=0;j<seenThisColumn.size();j++) seenThisColumn[j]=false;
		for (size_t j=0;j<x;j++) {
			size_t col = size_t(random()*n);
			if (seenThisColumn[col]) continue;
			seenThisColumn[col]=true;
			ComplexOrRealType val = random()*maxValue;
		
			sparse.pushValue(val);
			sparse.pushCol(col);
			counter++;
		}
	}
	sparse.setRow(n,counter);
	sparse.checkValidity();
	// symmetrize:
	SparseMatrixType sparse2;
	transposeConjugate(sparse2,sparse);
	sparse += sparse2;
	sparse.checkValidity();
	assert(isHermitian(sparse));

	// sparse solver setup
	SolverParameters params;
	params.lotaMemory=lotaMemory;
	LanczosSolverType lanczosSolver(sparse,params);
	DavidsonSolverType davisonSolver(sparse,params);
	SparseSolverType* solver = 0;

	// select solver
	if (useDavidson) solver = &davisonSolver;
	else solver = &lanczosSolver;

	// diagonalize matrix
	RealType gsEnergy = 0;
	VectorType gsVector(n);
	solver->computeGroundState(gsEnergy,gsVector);

	std::cout<<"Energy="<<gsEnergy<<"\n";
}

