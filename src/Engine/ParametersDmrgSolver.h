// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

/*! \file ParametersDmrgSolver.h
 *
 *  Contains the parameters for the DmrgSolver class and implements functionality to
 *  read them from a JSON file
 *
 */
#ifndef PARAMETERSDMRGSOLVER_HEADER_H
#define PARAMETERSDMRGSOLVER_HEADER_H

#include "Utils.h"
#include "SimpleReader.h"

namespace Dmrg {
	struct FiniteLoop {
		int stepLength; // how much to go right (+) or left (-)
		size_t keptStates; // kept states
		int saveOption; // to save or not to save	
	};

	std::istream &operator>>(std::istream& is,FiniteLoop& fl)
	{
		is>>fl.stepLength;
		is>>fl.keptStates;
		is>>fl.saveOption;
		return is;
	}

	std::ostream &operator<<(std::ostream& os,const FiniteLoop& fl)
	{
		os<<fl.stepLength<<" ";
		os<<fl.keptStates<<" ";
		os<<fl.saveOption;
		return os;
	}

	struct DmrgCheckPoint {
		bool enabled;
		size_t index;
		std::string filename;
		
		template<typename IoInputType>
		void load(IoInputType& io)
		{
			io.readline(enabled,"CheckPointEnabled=");
			io.readline(index,"index=");
			io.readline(filename,"filename=");
		}
	};

	std::istream &operator>>(std::istream& is,DmrgCheckPoint& c)
	{
		is>>c.enabled;
		is>>c.index;
		is>>c.filename;
		return is;	
	}

	//! Structure that contains the Dmrg parameters
	template<typename FieldType>
	struct ParametersDmrgSolver {
		std::string filename; // filename to save observables and continued fractions
		size_t keptStatesInfinite; // number of states kept by the Dmrg Truncation procedure (infinite loop only)
		std::vector<FiniteLoop> finiteLoop; // number of states kept by the Dmrg Truncation procedure (finite loops only)
		std::string version;
		std::string options; // reserved for future use
		std::vector<FieldType> targetQuantumNumbers;
		DmrgCheckPoint checkpoint;
		size_t nthreads;
		
		//! Read Dmrg parameters from inp file
		template<typename IoInputType>
		ParametersDmrgSolver(IoInputType& io) 
		{
			io.readline(options,"SolverOptions="); 
			io.readline(version,"Version=");
			io.readline(filename,"OutputFile=");
			io.readline(keptStatesInfinite,"InfiniteLoopKeptStates=");
			io.read(finiteLoop,"FiniteLoops");
			if (options.find("hasQuantumNumbers")!=std::string::npos) 
				io.read(targetQuantumNumbers,"TargetQuantumNumbers");
			if (options.find("checkpoint")!=std::string::npos)
				checkpoint.load(io);
			nthreads=1; // provide a default value
			if (options.find("hasThreads")!=std::string::npos)
				io.readline(nthreads,"Threads=");
			
		} 

	};

	//! Read Dmrg parameters from JSON file
	/* ParametersDmrgSolver&
	operator <= (ParametersDmrgSolver& parameters, const dca::JsonReader& reader) 
	{
		const dca::JsonAccessor<std::string>& dmrg(reader["programSpecific"]["DMRG"]["solverDmrg"]);
		
		parameters.version <= dmrg["version"];
		parameters.filename <= dmrg["outputfile"];
		parameters.numberOfKeptStates <= dmrg["numberOfKeptStates"];
		
		return parameters;
	} */

	//! print dmrg parameters
	template<typename FieldType>
	std::ostream &operator<<(std::ostream &os,ParametersDmrgSolver<FieldType> const &parameters)
	{
		os<<"#This is DMRG++\n";
		os<<"parameters.version="<<parameters.version<<"\n";
		os<<"parameters.filename="<<parameters.filename<<"\n";
		os<<"parameters.options="<<parameters.options<<"\n";
		os<<"parameters.keptStatesInfinite="<<parameters.keptStatesInfinite<<"\n";
		utils::vectorPrint(parameters.finiteLoop,"finiteLoop",os);
		if (parameters.options.find("hasQuantumNumbers")!=std::string::npos) {
			os<<"parameters.targetQuantumNumbers=";
			for (size_t i=0;i<parameters.targetQuantumNumbers.size();i++) os<<parameters.targetQuantumNumbers[i]<<" ";
			os<<"\n";
		} else {
			os<<"parameters.targetQuantumNumbers=search\n";
		}

		os<<"parameters.nthreads="<<parameters.nthreads<<"\n";
		return os;
	}
} // namespace Dmrg
/*@}*/

#endif
