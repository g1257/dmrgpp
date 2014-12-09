/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

#include "TypeToString.h"
#include "Vector.h"
#include "Provenance.h"
#include "BoostSerializationHeaders.h"

namespace Dmrg {
/* PSIDOC FiniteLoop
\subsection{Enabling finite loops}

To enable finite loops make sure that the option `nofiniteloops` is <b>not</b>
present under `SolverOptions=`. Remember that the entry <b>FiniteLoops</b>
in the input file is a series of space-separated numbers. More than one space is allowed.
The first
number is the number of finite algorithm ``movements,'' followed by series of three
numbers for
each movement. Of the three numbers, the first is the number of sites to go forward
if positive
or backward if negative. The second number is the <i>m</i> for this movement and the
last number
is either 0 or 1, 0 will not save state data to disk and 1 will save all data to be
able to calculate
observables. The first movement starts from where the infinite loop left off, at the
 middle of the
lattice.

\subsection{Example of a Finite loops line in the input file}

\begin{verbatim}
	FiniteLoops 4 7 200 0 -7 200 0 7 200 1 7 200 1
\end{verbatim}

The number 4 implies 4 finite loops. The first fine loop is \verb!7 200 0!, meaning
go forward 7 steps, use \emph{m=200} for this finite sweep, and 0: do not store
transformation in disk.
The next is \verb!-7 200 0!, which goes backwards 7 sites, etc.
Remember that the finite loops start at the middle of the lattice, where the
infinite loop left off.
\todo{ADD FIGURE SHOWING WHAT THIS DOES.}

\subsection{The third number in the triplet}
The save option is a bitwise option where the
first bit means save or don't save, and the second bit
compute the g.s. or WFT it.
So there are 4 combinations (as of today):
\begin{table}
\begin{tabular}{ll}\toprule
Value & Description\\
0       & Don't save, compute the ground state\\
1       & Save, compute the ground state\\
2       & Don't save, WFT the ground state\\
3       & Save, WFT the ground state\\
\end{tabular}
\end{table}
\subsection{Caveats and Troubleshooting}

If `nofiniteloops` is an option in the options line of the input file then
the <b>FiniteLoops</b> line in the input file is ignored, and no finite loops are done.
In this case, DMRG++ stops when the infinite algorithm has finished.

Make sure the first number is the number of triplets that follow.

Make sure you don't fall off the lattice, by going forward or backwards too much.
Remember that at least one site must remain for the system part of the lattice.
So on a 16 site chain, when you start the finite loops you're at the middle, you
can go forward at most 7 sites, and backwards at most 7 sites.

There is some checking done to the finite loops input, see PTEXREF{139},
but you might find that it's not comprehensive.
*/
struct FiniteLoop {
	int stepLength; // how much to go right (+) or left (-)
	unsigned int keptStates; // kept states
	int saveOption; // to save or not to save

	FiniteLoop() : stepLength(0),keptStates(0),saveOption(0)
	{}

	FiniteLoop(int sl,unsigned int ks,int so)
	    : stepLength(sl),keptStates(ks),saveOption(so)
	{}

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & stepLength;
		ar & keptStates;
		ar & saveOption;
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "FiniteLoop";
		const char* start = (const char *)&stepLength;
		const char* end = (const char*)&keptStates;
		SizeType total = mres.memResolv(&stepLength,end-start,str + " stepLength");

		start = end;
		end = (const char*)&saveOption;
		total += mres.memResolv(&keptStates,end-start,str + " keptStates");
		total += mres.memResolv(&saveOption,sizeof(*this)-total,str + " saveOption");

		return total;
	}
};

//PTEX_LABEL{139}
inline void checkFiniteLoops(const PsimagLite::Vector<FiniteLoop>::Type& finiteLoop,
                             SizeType totalSites)
{
	PsimagLite::String s = "checkFiniteLoops: I'm falling out of the lattice ";
	PsimagLite::String loops = "";
	int x = totalSites/2-1; // must be signed
	if (totalSites & 1) x++;
	if (finiteLoop[0].stepLength<0) x++;
	int prevDeltaSign = 1;
	SizeType sopt = 0; // have we started saving yet?
	for (SizeType i=0;i<finiteLoop.size();i++)  {
		SizeType thisSaveOption = (finiteLoop[i].saveOption & 1);
		if (sopt == 1 && thisSaveOption ==0) {
			s = "Error for finite loop number " + ttos(i) + "\n";
			s += "Once you say 1 on a finite loop, then all";
			s += " finite loops that follow must have 1.";
			throw PsimagLite::RuntimeError(s.c_str());
		}
		if (sopt == 0 && thisSaveOption ==1) {
			sopt = 1;
			if (SizeType(x) != 1 && SizeType(x)!=totalSites-2) {
				s = __FILE__ + PsimagLite::String(": FATAL: for finite loop number ")
				        + ttos(i) + "\n";
				s += "Saving finite loops must start at the left or";
				s += " right end of the lattice\n";
				throw PsimagLite::RuntimeError(s.c_str());
			}
		}
		// naive location:
		int delta = finiteLoop[i].stepLength;
		x += delta;
		loops = loops + ttos(delta) + " ";

		// take care of bounces:
		if (i>0 && delta*prevDeltaSign < 0) x += prevDeltaSign;
		prevDeltaSign = 1;
		if (delta<0) prevDeltaSign = -1;

		// check that we don't fall out
		bool flag = false;
		if (x<=0) {
			s = s + "on the left end\n";
			flag = true;
		}
		if (SizeType(x)>=totalSites-1) {
			s = s + "on the right end\n";
			flag = true;
		}
		if (flag) {
			// complain and die if we fell out:
			s = s + "Loops so far: " + loops + "\n";
			s =s + "x=" + ttos(x) + " last delta=" +
			        ttos(delta);
			s =s + " sites=" + ttos(totalSites);
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

}

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

	DmrgCheckPoint() : enabled(false)
	{}

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & enabled;
		ar & filename;
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x = sizeof(DmrgCheckPoint),
	                   PsimagLite::String msg = "") const
	{
		assert(x == sizeof(DmrgCheckPoint));

		PsimagLite::String str = msg;
		msg += "DmrgCheckPoint";
		const char* start = (const char *)&enabled;
		const char* end = (const char*)&filename;
		SizeType total = mres.memResolv(&enabled,end-start,str + " enabled");

		total += mres.memResolv(&filename,x-total,str + " filename");

		return total;
	}

	bool enabled;
	PsimagLite::String filename;
};

std::istream &operator>>(std::istream& is,DmrgCheckPoint& c)
{
	is>>c.enabled;
	is>>c.filename;
	return is;
}

/* PSIDOC ParametersDmrgSolver
\begin{itemize}
\item[Model=string]
A string indicating the model, be it HubbardOneBand HeisenbergSpinOneHalf, etc.

PSIDOCCOPY dmrgSolverOptions

\item[version=string]  A mandatory string that is read and ignored. Usually contains the
result
of doing \verb!git rev-parse HEAD!.

\item[outputfile=string]  The output file. This file will be created if non-existent,
 and if it
exits it will be truncated.

\item[InfiniteLoopKeptStates=integer]  \emph{m} value for the infinite algorithm.

\item[FiniteLoops=vector]
A series of space-separated numbers. More than one space is allowed.
The first number is the number of finite algorithm movements, followed by series
of three numbers for each movement. Of the three numbers, the first
is the number of sites to go forward if positive or backward if negative.
The second number is the <i>m</i> for this movement and the last number is either 0 or 1,
0 will not save state data to disk and 1 will save all data to be able to calculate
observables.
The first movement starts from where the infinite loop left off, at the middle of the
 lattice.
See the below for more information and examples on Finite Loops.

\item[TargetElectronsUp=integer]

\item[TargetElectronsDown=integer]
\end{itemize}
*/
template<typename FieldType,typename InputValidatorType>
struct ParametersDmrgSolver {
	typedef ParametersDmrgSolver<FieldType, InputValidatorType> ThisType;

	SizeType electronsUp;
	SizeType electronsDown;
	SizeType nthreads;
	SizeType sitesPerBlock;
	SizeType maxMatrixRankStored;
	SizeType keptStatesInfinite;
	int useReflectionSymmetry;
	FieldType tolerance;
	std::pair<bool,FieldType> gsWeight;
	PsimagLite::String filename;
	PsimagLite::String version;
	PsimagLite::String options;
	PsimagLite::String model;
	PsimagLite::String insitu;
	PsimagLite::String fileForDensityMatrixEigs;
	DmrgCheckPoint checkpoint;
	typename PsimagLite::Vector<FieldType>::Type targetQuantumNumbers;
	typename PsimagLite::Vector<SizeType>::Type adjustQuantumNumbers;
	typename PsimagLite::Vector<FiniteLoop>::Type finiteLoop;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version1)
	{
		ar & electronsUp;
		ar & electronsDown;
		ar & nthreads;
		ar & sitesPerBlock;
		ar & maxMatrixRankStored;
		ar & keptStatesInfinite;
		ar & useReflectionSymmetry;
		ar & tolerance;
		ar & gsWeight;
		ar & filename;
		ar & version;
		ar & options;
		ar & model;
		ar & insitu;
		ar & fileForDensityMatrixEigs;
		ar & checkpoint;
		ar & targetQuantumNumbers;
		ar & adjustQuantumNumbers;
		ar & finiteLoop;
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg) const
	{
		PsimagLite::String str = msg;
		msg += "ParametersDmrgSolver";
		const char* start = (const char *)&electronsUp;
		const char* end = (const char*)&electronsDown;
		SizeType total = mres.memResolv(&electronsUp,end-start,str + "electronsUp");
		start = end;
		end = (const char*)&nthreads;
		total += mres.memResolv(&electronsDown,end-start,str + "electronsDown");
		start = end;
		end = (const char*)&sitesPerBlock;
		total += mres.memResolv(&nthreads,end-start,str + "nthreads");
		start = end;
		end = (const char*)&maxMatrixRankStored;
		total += mres.memResolv(&sitesPerBlock,end-start,str + "sitesPerBlock");
		start = end;
		end = (const char*)&keptStatesInfinite;
		total += mres.memResolv(&maxMatrixRankStored,end-start,
		                        str + "maxMatrixRankStored");
		start = end;
		end = (const char*)&useReflectionSymmetry;
		total += mres.memResolv(&keptStatesInfinite,end-start,
		                        str + "keptStatesInfinite");
		start = end;
		end = (const char*)&tolerance;
		total += mres.memResolv(&useReflectionSymmetry,end-start,
		                        str + "useReflectionSymmetry");
		start = end;
		end = (const char*)&gsWeight;
		total += mres.memResolv(&tolerance,end-start,str + "tolerance");
		start = end;
		end = (const char*)&filename;
		total += mres.memResolv(&gsWeight,end-start,str + "gsWeight");

		start = end;
		end = (const char*)&version;
		total += mres.memResolv(&filename,end-start,str + "filename");
		start = end;
		end = (const char*)&options;
		total += mres.memResolv(&version,end-start,str + "version");
		start = end;
		end = (const char*)&model;
		total += mres.memResolv(&options,end-start,str + "options");
		start = end;
		end = (const char*)&insitu;
		total += mres.memResolv(&model,end-start,str + "model");
		start = end;
		end = (const char*)&fileForDensityMatrixEigs;
		total += mres.memResolv(&insitu,end-start,str + "insitu");
		start = end;
		end = (const char*)&checkpoint;
		total += mres.memResolv(&fileForDensityMatrixEigs,end-start,
		                        str + "fileForDensityMatrixEigs");

		start = end;
		end = (const char*)&targetQuantumNumbers;
		total += mres.memResolv(&checkpoint,end-start, str + "checkpoint");

		start = end;
		end = (const char*)&adjustQuantumNumbers;
		total += mres.memResolv(&targetQuantumNumbers,end-start,
		                        str + "targetQuantumNumbers");
		start = end;
		end = (const char*)&finiteLoop;
		total += mres.memResolv(&adjustQuantumNumbers,end-start,
		                        str + "adjustQuantumNumbers");
		total += mres.memResolv(&finiteLoop,sizeof(*this)-total, str + "finiteLoop");

		return total;
	}

	ParametersDmrgSolver(PsimagLite::String filename)
	{
		std::ifstream ifs(filename.c_str());
		boost::archive::text_iarchive ia(ifs);
		ia >> (*this);
	}

	//! Read Dmrg parameters from inp file
	ParametersDmrgSolver(InputValidatorType& io)
	    : sitesPerBlock(1),
	      maxMatrixRankStored(0),
	      gsWeight(false,0.0)
	{
		io.readline(model,"Model=");
		io.readline(options,"SolverOptions=");
		io.readline(version,"Version=");
		io.readline(filename,"OutputFile=");
		io.readline(keptStatesInfinite,"InfiniteLoopKeptStates=");
		typename PsimagLite::Vector<FieldType>::Type tmpVec;
		io.read(tmpVec,"FiniteLoops");
		for (SizeType i=0;i<tmpVec.size();i+=3) {
			typename PsimagLite::Vector<int>::Type xTmp(3);
			for (SizeType j=0;j<xTmp.size();j++) xTmp[j]=int(tmpVec[i+j]);
			FiniteLoop fl(xTmp[0],xTmp[1],xTmp[2]);
			finiteLoop.push_back(fl);
		}

		SizeType repeat = 0;

		try {
			io.readline(repeat,"RepeatFiniteLoopsTimes=");
		}  catch (std::exception& e) {}

		SizeType fromFl = 0;
		try {
			io.readline(fromFl,"RepeatFiniteLoopsFrom=");
		}  catch (std::exception& e) {}

		SizeType upToFl = finiteLoop.size()-1;
		try {
			io.readline(upToFl,"RepeatFiniteLoopsTo=");
		}  catch (std::exception& e) {}

		if (upToFl>=finiteLoop.size()) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: RepeatFiniteLoopsTo=" + ttos(upToFl);
			s += " is larger than current finite loops\n";
			s += "\nMaximum is " + ttos(finiteLoop.size())+ "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
		if (fromFl>upToFl) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: RepeatFiniteLoopsFrom=" + ttos(fromFl);
			s += " is larger than RepeatFiniteLoopsTo\n";
			s += "\nMaximum is " + ttos(upToFl)+ "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}
		upToFl++;

		for (SizeType i=0;i<repeat;i++) {
			for (SizeType j=fromFl;j<upToFl;j++) {
				FiniteLoop fl = finiteLoop[j];
				finiteLoop.push_back(fl);
			}
		}

		if (options.find("hasQuantumNumbers")!=PsimagLite::String::npos) {
			PsimagLite::String s = "*** WARNING: hasQuantumNumbers ";
			s += "option is obsolete in input file\n";
			std::cerr<<s;
		}
		try {
			io.read(targetQuantumNumbers,"TargetQuantumNumbers");
			PsimagLite::String s = "*** WARNING: TargetQuantumNumbers ";
			s += "is deprecated in input file\n";
			std::cerr<<s;
		} catch (std::exception& e){}

		bool hasElectrons = false;
		try {
			io.readline(electronsUp,"TargetElectronsUp");
			io.readline(electronsDown,"TargetElectronsDown");
			hasElectrons = true;
		} catch (std::exception& e) {}

		if (hasElectrons && targetQuantumNumbers.size()>0) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: Specifying both TargetElectronsUp/Down ";
			s += "and TargetQuantumNumbers is an error.";
			s += "\nSpecify one or the other only.\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		if (!hasElectrons && targetQuantumNumbers.size()==0) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: Either TargetElectronsUp/Down or TargetQuantumNumbers ";
			s += "must be specified.\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		if (options.find("useSu2Symmetry")!=PsimagLite::String::npos && hasElectrons) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: TargetElectronsUp/Down cannot be specified while ";
			s += "using SU(2) symmetry\n";
			s += "\nTargetQuantumNumbers must be specified instead.\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		try {
			io.read(adjustQuantumNumbers,"AdjustQuantumNumbers");
		} catch (std::exception& e) {}

		tolerance = -1.0;
		try {
			io.readline(tolerance,"TruncationTolerance=");
		} catch (std::exception& e) {}

		if (options.find("checkpoint")!=PsimagLite::String::npos)
			io.readline(checkpoint.filename,"CheckpointFilename=");
		else if (options.find("restart")!=PsimagLite::String::npos)
			io.readline(checkpoint.filename,"RestartFilename=");

		nthreads=1; // provide a default value
		try {
			io.readline(nthreads,"Threads=");
		} catch (std::exception& e) {}

		if (nthreads==0) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: nthreads cannot be zero\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		useReflectionSymmetry=0;
		try {
			io.readline(useReflectionSymmetry,"UseReflectionSymmetry=");
		} catch (std::exception& e) {}
		fileForDensityMatrixEigs="";
		try {
			io.readline(fileForDensityMatrixEigs,"FileForDensityMatrixEigs=");
		} catch (std::exception& e) {}

		insitu = "";
		try {
			io.readline(insitu,"insitu=");
		} catch (std::exception& e) {}

		try {
			io.readline(sitesPerBlock,"SitesPerBlock=");
		} catch (std::exception& e) {}

		try {
			io.readline(gsWeight.second,"GsWeight=");
			gsWeight.first = true;
		} catch (std::exception& e) {}

		try {
			io.readline(maxMatrixRankStored,"MaxMatrixRankStored=");
		} catch (std::exception& e) {}
	}
};

// print dmrg parameters
template<typename FieldType,typename InputValidatorType>
std::ostream &operator<<(std::ostream &os,
                         const ParametersDmrgSolver<FieldType,InputValidatorType> & p)
{
	os<<"#This is DMRG++\n";
	Provenance provenance;
	os<<provenance;
	os<<"parameters.version="<<p.version<<"\n";
	os<<"parameters.model="<<p.model<<"\n";
	os<<"parameters.filename="<<p.filename<<"\n";
	os<<"parameters.options="<<p.options<<"\n";
	os<<"parameters.keptStatesInfinite="<<p.keptStatesInfinite<<"\n";
	os<<"finiteLoop\n";
	os<<p.finiteLoop;

	if (p.targetQuantumNumbers.size()>0) {
		os<<"parameters.targetQuantumNumbers=";
		for (SizeType i=0;i<p.targetQuantumNumbers.size();i++)
			os<<p.targetQuantumNumbers[i]<<" ";
		os<<"\n";
	} else {
		os<<"parameters.electronsUp="<<p.electronsUp<<"\n";
		os<<"parameters.electronsDown="<<p.electronsDown<<"\n";
	}
	if (p.tolerance>0)
		os<<"parameters.tolerance="<<p.tolerance<<"\n";
	os<<"parameters.nthreads="<<p.nthreads<<"\n";
	os<<"parameters.useReflectionSymmetry="<<p.useReflectionSymmetry<<"\n";
	if (p.checkpoint.filename!="")
		os<<"parameters.restartFilename="<<p.checkpoint.filename<<"\n";
	if (p.fileForDensityMatrixEigs!="")
		os<<"parameters.fileForDensityMatrixEigs="<<p.fileForDensityMatrixEigs<<"\n";

	if (p.gsWeight.first)
		os<<"GsWeight="<<p.gsWeight.second<<"\n";

	if (p.options.find("MatrixVectorStored")==PsimagLite::String::npos)
		os<<"MaxMatrixRankStored="<<p.maxMatrixRankStored<<"\n";

	return os;
}
} // namespace Dmrg
/*@}*/

#endif

