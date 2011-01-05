
#include "ConcurrencySerial.h"
#include "ContinuedFraction.h"
#include "HubbardLanczos.h"
#include "ParametersModelHubbard.h"
#include "Geometry.h"
#include "IoSimple.h"

using namespace Dmrg;
typedef double RealType;
typedef ConcurrencySerial<RealType> ConcurrencyType;
typedef Geometry<RealType> GeometryType;
typedef ParametersModelHubbard<RealType> ParametersModelType;
typedef IoSimple::In IoInputType;
typedef HubbardLanczos<RealType,ParametersModelType,GeometryType> ModelType;
typedef ContinuedFraction<ModelType,ConcurrencyType> ContinuedFractionType;
typedef typename ContinuedFractionType::TridiagonalMatrixType TridiagonalMatrixType;

RealType greenFunction(const TridiagonalMatrixType& ab,RealType z,RealType norma)
{
	//std::vector<RealType> eigs(T.n_row());
	//utils::diag(T,eigs,'V');
	RealType isign = 1.0;
	return norma * ab.computeContinuedFraction(z,isign);

}

int main(int argc,char *argv[])
{
	//! setup distributed parallelization
	ConcurrencyType concurrency(argc,argv);

	//Setup the Geometry
	IoInputType io(argv[1]);
	GeometryType geometry(io);

	// read model parameters
	ParametersModelType mp(io);

	// print license
	std::string license = "Copyright (c) 2009 , UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[DMRG++, Version 2.0.0]\n"
"\n"
"*********************************************************\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED. \n"
"\n"
"Please see full open source license included in file LICENSE.\n"
"*********************************************************\n"
"\n"
"\n"
"// END LICENSE BLOCK\n"
;
	if (concurrency.root()) std::cerr<<license;

	std::vector<RealType> qns;
	io.read(qns,"TargetQuantumNumbers");
	if (qns.size()<2) throw std::runtime_error("HubbardLanczos::ctor(...)\n");
	size_t nup=geometry.numberOfSites()*qns[0];
	size_t ndown=geometry.numberOfSites()*qns[1];
	//! Setup the Model
	ModelType model(nup,ndown,mp,geometry);
	size_t i = atoi(argv[2]);
	size_t j = atoi(argv[3]);

	ContinuedFractionType cf(model);

	//! get the g.s.:
	std::cout<<"Energy="<<cf.gsEnergy()<<"\n";

	std::cout<<"gf(i="<<i<<",j="<<j<<")\n";
	RealType norma;
	TridiagonalMatrixType ab;
	cf.getGreenFunction(ab,norma,i,j);
	for (int i=0;i<atoi(argv[4]);i++) {
		RealType omega = atof(argv[5])*i;
		std::cout<<omega<<" "<<greenFunction(ab,omega,norma)<<"\n";
	}
}

