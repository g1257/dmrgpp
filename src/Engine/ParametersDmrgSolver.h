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
#include "IoSimple.h"

namespace Dmrg {
/* PSIDOC FiniteLoop
\subsection{Enabling finite loops}

To enable finite loops make sure that the option `nofiniteloops` is \emph{not}
present under `SolverOptions=`. Remember that the entry \verb!FiniteLoops!
in the input file is a series of space-separated numbers. More than one space is allowed.
The first
number is the number of finite algorithm ``movements,'' followed by series of three
numbers for
each movement. Of the three numbers, the first is the number of sites to go forward
if positive
or backward if negative. The second number is the $m$ for this movement and the
last number
is a bitwise option described in \S\ref{sec:thirdnumber}.
The first movement starts from where the infinite loop left off, at the
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

\subsection{The third number in the triplet}\label{sec:thirdnumber}
The third number in the triplet is a bitwise option where the
first bit means save or don't save, the second bit
compute the g.s. or WFT it without eigenvalue or eigenvector updates,
and the third bit compute the g.s. of WFT it updating eigenvalue and eigenvector.
It is a fatal error to have both bits 1 and 2 set.
\begin{table}
\begin{tabular}{ll}\toprule
Bit & Description\\
0       & save or don't save state for the observe code\\
1       & compute the g.s. or fast WFT it\\
2       & compute the g.s. or slowly WFT it\\
\end{tabular}
\caption{Meaning of each bit of the third number in the
finite loop triplet. It is a fatal error to have both bits 1 and 2 set.}
\end{table}
\subsection{Caveats and Troubleshooting}

If `nofiniteloops` is an option in the options line of the input file then
the \verb!FiniteLoops! line in the input file is ignored, and no finite loops are done.
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
A string indicating the model, be it HubbardOneBand Heisenberg, etc.

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
The second number is the \emph{m} for this movement and the last number is either 0 or 1,
0 will not save state data to disk and 1 will save all data to be able to calculate
observables.
The first movement starts from where the infinite loop left off, at the middle of the
 lattice.
See the below for more information and examples on Finite Loops.

\end{itemize}
*/
template<typename FieldType,typename InputValidatorType>
struct ParametersDmrgSolver {

	typedef ParametersDmrgSolver<FieldType, InputValidatorType> ThisType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<FieldType>::Type VectorFieldType;

	SizeType nthreads;
	SizeType sitesPerBlock;
	SizeType maxMatrixRankStored;
	SizeType keptStatesInfinite;
	SizeType excited;
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
	VectorSizeType adjustQuantumNumbers;
	typename PsimagLite::Vector<FiniteLoop>::Type finiteLoop;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version1)
	{
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
		ar & adjustQuantumNumbers;
		ar & finiteLoop;
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg) const
	{
		return 0;
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
	      excited(0),
	      gsWeight(false,0.0)
	{
		io.readline(model,"Model=");
		io.readline(options,"SolverOptions=");
		io.readline(version,"Version=");
		io.readline(filename,"OutputFile=");
		io.readline(keptStatesInfinite,"InfiniteLoopKeptStates=");
		readFiniteLoops(io,finiteLoop,0);

		if (options.find("hasQuantumNumbers")!=PsimagLite::String::npos) {
			PsimagLite::String s = "*** FATAL: hasQuantumNumbers ";
			s += "option is no longer allowed in input file\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		VectorSizeType tmpVector;
		try {
			io.read(tmpVector,"TargetQuantumNumbers");
		} catch (std::exception&){}


		if (tmpVector.size()>0) {
			PsimagLite::String s = "*** FATAL: TargetQuantumNumbers ";
			s += "is no longer allowed in input file\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		try {
			io.read(adjustQuantumNumbers,"AdjustQuantumNumbers");
		} catch (std::exception& e) {}

		tolerance = -1.0;
		try {
			io.readline(tolerance,"TruncationTolerance=");
		} catch (std::exception& e) {}

		if (options.find("checkpoint")!=PsimagLite::String::npos) {
			io.readline(checkpoint.filename,"CheckpointFilename=");
			checkRestart(filename,checkpoint.filename,options,"CheckpointFilename=");
		} else if (options.find("restart")!=PsimagLite::String::npos) {
			io.readline(checkpoint.filename,"RestartFilename=");
			checkRestart(filename,checkpoint.filename,options,"CheckpointFilename=");
		}

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

		try {
			io.readline(excited,"Excited=");
		} catch (std::exception& e) {}
	}

	template<typename SomeInputType>
	static void readFiniteLoops(SomeInputType& io,
	                            PsimagLite::Vector<FiniteLoop>::Type& vfl,
	                            SizeType offset)
	{
		VectorFieldType tmpVec;
		io.read(tmpVec,"FiniteLoops");
		readFiniteLoops_(io,vfl,tmpVec,offset);
	}

	template<typename SomeInputType>
	static void readFiniteLoops_(SomeInputType& io,
	                            PsimagLite::Vector<FiniteLoop>::Type& vfl,
	                            const VectorFieldType& tmpVec,
	                            SizeType offset)
	{
		for (SizeType i=0;i<tmpVec.size();i+=3) {
			typename PsimagLite::Vector<int>::Type xTmp(3);
			for (SizeType j=0;j<xTmp.size();j++) xTmp[j]=int(tmpVec[i+j]);
			FiniteLoop fl(xTmp[0],xTmp[1],xTmp[2]);
			vfl.push_back(fl);
		}

		SizeType repeat = 0;

		try {
			io.readline(repeat,"RepeatFiniteLoopsTimes=");
		}  catch (std::exception& e) {}

		SizeType fromFl = offset;
		try {
			io.readline(fromFl,"RepeatFiniteLoopsFrom=");
			fromFl += offset;
		}  catch (std::exception& e) {}

		SizeType upToFl = vfl.size()-1;
		try {
			io.readline(upToFl,"RepeatFiniteLoopsTo=");
			fromFl += offset;
		}  catch (std::exception&) {}

		if (upToFl >= vfl.size()) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: RepeatFiniteLoopsTo=" + ttos(upToFl);
			s += " is larger than current finite loops\n";
			s += "\nMaximum is " + ttos(vfl.size())+ "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		if (fromFl > upToFl) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: RepeatFiniteLoopsFrom=" + ttos(fromFl);
			s += " is larger than RepeatFiniteLoopsTo\n";
			s += "\nMaximum is " + ttos(upToFl)+ "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		upToFl++;

		for (SizeType i=0;i<repeat;i++) {
			for (SizeType j=fromFl;j<upToFl;j++) {
				FiniteLoop fl = vfl[j];
				vfl.push_back(fl);
			}
		}
	}

private:

	static void checkRestart(PsimagLite::String filename1,
	                         PsimagLite::String filename2,
	                         PsimagLite::String options,
	                         PsimagLite::String label)
	{
		checkFilesNotEqual(filename1,filename2,label);
		checkTwoSiteDmrg(filename2,options);
	}

	static void checkFilesNotEqual(PsimagLite::String filename1,
	                               PsimagLite::String filename2,
	                               PsimagLite::String label)
	{
		if (filename1 != filename2) return;
		PsimagLite::String s (__FILE__);
		s += "\nFATAL: " + filename1 + "is equal to " + filename2;
		s += " in label " + label + " in input file\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}

	static void checkTwoSiteDmrg(PsimagLite::String filename2,
	                             PsimagLite::String options)
	{
		PsimagLite::IoSimple::In io(filename2);
		PsimagLite::String optionsOld;
		io.readline(optionsOld,"parameters.options");
		bool bOld = (optionsOld.find("twositedmrg")!=PsimagLite::String::npos);
		bool b = (options.find("twositedmrg")!=PsimagLite::String::npos);
		if (bOld == b) return;
		PsimagLite::String s (__FILE__);
		PsimagLite::String s1 = ": Previous run used twositedmrg ";
		s1 += " but this one does not.";
		PsimagLite::String s2 = ": Previous run didn't use twositedmrg ";
		s2 += " but this one does.";
		if (bOld) s += s1;
		else s += s2;
		throw PsimagLite::RuntimeError(s + "\n");
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
	os<<"FiniteLoops ";
	os<<p.finiteLoop;

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

