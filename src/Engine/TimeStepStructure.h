// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."
 
*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file TimeStepStructure.h
 *
 *  FIXME DOCumentation
 */
#ifndef TIMESTEP_STRUCT_H
#define TIMESTEP_STRUCT_H

#include "Utils.h"
#include "TargetStructureParams.h"
#include "SimpleReader.h"

namespace Dmrg {
		
	template<typename OperatorType>
	struct TimeStepStructure {
		TimeStepStructure() :
			tau(0),
			timeSteps(0),
			advanceEach(0)
		{
		}
		
		std::string filename;
		typename OperatorType::RealType tau;
		size_t timeSteps;
		size_t advanceEach;
		std::vector<size_t> sites;
		std::vector<size_t> startingLoops;
		std::vector<OperatorType> aOperators;
		std::vector<size_t> electrons;
	};
	
// 	template<typename OperatorType,typename ModelType>
// 	inline TargetStructureParams<TimeStepStructure<OperatorType>,ModelType>&
// 	operator<=(TargetStructureParams<TimeStepStructure<OperatorType>,ModelType>& tsp,SimpleReader& reader)
// 	{
// 		typedef typename ModelType::RealType RealType;
// 		std::vector<size_t> sites,loops;
// 		std::string s;
// 		reader.read(s); // filename
// 		RealType tau=0;
// 		reader.read(tau);
// 		size_t timeSteps=0;
// 		reader.read(timeSteps);
// 		size_t advanceEach=0;
// 		reader.read(advanceEach);
// 		reader.read(sites);
// 		reader.read(loops);
// 		
// 		tsp.init(s,tau,timeSteps,advanceEach,sites,loops);
// 		
// 		for (size_t i=0;i<sites.size();i++) {
// 			//std::string s;
// 			reader.read(s);
// 			if (s == "cooked") {
// 				reader.read(s);
// 				std::vector<size_t> v;
// 				reader.read(v);
// 				tsp.setCookedData(i,s,v);
// 			} else {
// 				psimag::Matrix<RealType> m;
// 				reader.read(m);
// 				tsp.setRawData(i,m);
// 			}
// 			int fermiSign=0;
// 			reader.read(fermiSign);
// 			std::pair<size_t,size_t> jmValues;
// 			reader.read(jmValues);
// 			RealType angularFactor;
// 			reader.read(angularFactor);
// 			tsp.set(i,fermiSign,jmValues,angularFactor);
// 		}
// 		
// 		return tsp;
// 	}
// 	
	template<typename OperatorType>
	inline std::ostream&
	operator<<(std::ostream& os,const TimeStepStructure<OperatorType>& t)
	{
		os<<"#TimeStepStructure.operators"<<t.aOperators.size()<<"\n";
		for (size_t i=0;i<t.aOperators.size();i++) {
			os<<"#TimeStepStructure.operator "<<i<<"\n";
			os<<t.aOperators[i];
		}
		os<<"#TimeStepStructure.electrons\n";
		os<<t.electrons;
		os<<"#TimeStepStructure.site="<<t.sites;
		os<<"#TimeStepStructure.startingLoop="<<t.startingLoops<<"\n";
		os<<"#TimeStepStructure.filename="<<t.filename<<"\n";
		os<<"#TimeVectorsfilename.tau="<<t.tau<<"\n";
		os<<"#TimeVectorsfilename.timeSteps="<<t.timeSteps<<"\n";
		os<<"#TimeVectorsfilename.advanceEach="<<t.advanceEach<<"\n";
		return os;
	}
	
	
} // namespace Dmrg 

/*@}*/
#endif //TIMESTEP_STRUCT_H
