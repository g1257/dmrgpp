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
#include "ArchiveFiles.h"
#include "PsimagLite.h"
#include "RestartStruct.h"
#include "FiniteLoop.h"

namespace Dmrg {

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
	typedef PsimagLite::Matrix<FieldType> MatrixFieldType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef std::pair<FieldType, SizeType> PairRealSizeType;
	typedef typename PsimagLite::Vector<FiniteLoop>::Type VectorFiniteLoopType;

	SizeType nthreads;
	SizeType sitesPerBlock;
	SizeType maxMatrixRankStored;
	SizeType keptStatesInfinite;
	SizeType excited;
	SizeType dumperBegin;
	SizeType dumperEnd;
	SizeType precision;
	int useReflectionSymmetry;
	PairRealSizeType truncationControl;
	PsimagLite::String filename;
	PsimagLite::String version;
	PsimagLite::String options;
	PsimagLite::String model;
	PsimagLite::String insitu;
	PsimagLite::String fileForDensityMatrixEigs;
	PsimagLite::String recoverySave;
	RestartStruct checkpoint;
	VectorSizeType adjustQuantumNumbers;
	VectorFiniteLoopType finiteLoop;
	FieldType degeneracyMax;
	FieldType denseSparseThreshold;

	template<class Archive>
	void serialize(Archive&, const unsigned int)
	{}

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
	ParametersDmrgSolver(InputValidatorType& io,
	                     PsimagLite::String sOptions,
	                     bool earlyExit = false,
	                     bool isObserveCode = false)
	    : sitesPerBlock(1),
	      maxMatrixRankStored(0),
	      excited(0),
	      dumperBegin(0),
	      dumperEnd(0),
	      precision(6),
	      recoverySave("0"),
	      degeneracyMax(1e-12),
	      denseSparseThreshold(0.1)
	{
		io.readline(model,"Model=");
		io.readline(options,"SolverOptions=");
		options += sOptions;
		io.readline(version,"Version=");
		io.readline(filename,"OutputFile=");
		if (earlyExit) return;

		io.readline(keptStatesInfinite,"InfiniteLoopKeptStates=");
		readFiniteLoops(io,finiteLoop);

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
		} catch (std::exception&) {}

		truncationControl = PairRealSizeType(-1.0,keptStatesInfinite);
		try {
			PsimagLite::String s("");
			VectorStringType tokens;
			io.readline(s,"TruncationTolerance=");
			PsimagLite::split(tokens, s, ",");
			truncationControl.first = atof(tokens[0].c_str());
			if (tokens.size() > 1)
				truncationControl.second = atoi(tokens[1].c_str());
			warnIfFiniteMlessThanMin(finiteLoop, truncationControl.second);
		} catch (std::exception&) {}

		nthreads=1; // provide a default value
		try {
			io.readline(nthreads,"Threads=");
		} catch (std::exception&) {}

		if (nthreads==0) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: nthreads cannot be zero\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		useReflectionSymmetry=0;
		try {
			io.readline(useReflectionSymmetry,"UseReflectionSymmetry=");
		} catch (std::exception&) {}
		fileForDensityMatrixEigs="";
		try {
			io.readline(fileForDensityMatrixEigs,"FileForDensityMatrixEigs=");
		} catch (std::exception&) {}

		insitu = "";
		try {
			io.readline(insitu,"insitu=");
		} catch (std::exception&) {}

		try {
			io.readline(sitesPerBlock,"SitesPerBlock=");
		} catch (std::exception&) {}

		try {
			io.readline(maxMatrixRankStored,"MaxMatrixRankStored=");
		} catch (std::exception&) {}

		try {
			io.readline(excited,"Excited=");
		} catch (std::exception&) {}

		try {
			io.readline(degeneracyMax,"DegeneracyMax=");
		} catch (std::exception&) {}

		try {
			io.readline(recoverySave,"RecoverySave=");
		} catch (std::exception&) {}

		try {
			io.readline(dumperBegin,"KroneckerDumperBegin=");
		} catch (std::exception&) {}

		try {
			io.readline(dumperEnd,"KroneckerDumperEnd=");
		} catch (std::exception&) {}

		if (options.find("KroneckerDumper") != PsimagLite::String::npos) {
			if (options.find("MatrixVectorStored") != PsimagLite::String::npos) {
				PsimagLite::String msg("FATAL: KroneckerDumper cannot run with ");
				throw PsimagLite::RuntimeError(msg + "MatrixVectorStored\n");
			}
		} else {
			if (dumperBegin > 0 || dumperEnd > 0) {
				PsimagLite::String msg("FATAL: KroneckerDumperBegin|End needs ");
				throw PsimagLite::RuntimeError(msg + "KroneckerDumper in SolverOptions\n");
			}
		}

		try {
			io.readline(precision,"Precision=");
		} catch (std::exception&) {}

		try {
			io.readline(denseSparseThreshold, "DenseSparseThreshold=");
		} catch (std::exception&) {}

		if (isObserveCode) return;
		bool hasRestart = false;
		if (options.find("restart")!=PsimagLite::String::npos) {
			io.readline(checkpoint.filename,"RestartFilename=");
			ArchiveFiles<ThisType>::unpackIfNeeded(checkpoint.filename);
			checkRestart(filename,checkpoint.filename,options,"RestartFilename=");
			hasRestart = true;
		} else {
			PsimagLite::String str;
			try {
				io.readline(str,"RestartFilename=");
				PsimagLite::String s = "WARNING: RestartFilename ignored in input ";
				s += "because restart option not present in SolverOptions.\n";
				std::cerr<<s;
			} catch (std::exception&) {}
		}

		if (hasRestart) {
			try {
				io.readline(checkpoint.into,"RestartInto=");
			} catch (std::exception&) {}

			if (checkpoint.into != "All" && checkpoint.into != "GroundState") {
				PsimagLite::String str = "FATAL: RestartInto=All | GroundState\n";
				throw PsimagLite::RuntimeError(str);
			}

			try {
				io.readline(checkpoint.labelForPsi,"RestartLabelForPsi=");
			} catch (std::exception&) {}

			try {
				io.readline(checkpoint.labelForEnergy,"RestartLabelForEnergy=");
			} catch (std::exception&) {}
		}
	}

	template<typename SomeInputType>
	static void readFiniteLoops(SomeInputType& io,
	                            VectorFiniteLoopType& vfl)
	{
		if (io.version() < io.versionAinur()) {
			VectorFieldType tmpVec;
			io.read(tmpVec,"FiniteLoops");
			readFiniteLoops_(io,vfl,tmpVec);
		} else {
			MatrixFieldType tmpMat;
			io.readMatrix(tmpMat,"FiniteLoops");
			readFiniteLoops_(io,vfl,tmpMat);
		}
	}

	template<typename SomeInputType>
	static void readFiniteLoops_(SomeInputType& io,
	                             VectorFiniteLoopType& vfl,
	                             const VectorFieldType& tmpVec)
	{
		for (SizeType i = 0; i < tmpVec.size(); i += 3) {
			typename PsimagLite::Vector<int>::Type xTmp(3);
			assert(2 + i < tmpVec.size());
			for (SizeType j = 0; j < xTmp.size(); ++j)
				xTmp[j] = static_cast<int>(tmpVec[i+j]);
			FiniteLoop fl(xTmp[0], xTmp[1], xTmp[2]);
			vfl.push_back(fl);
		}

		readFiniteLoops_(io, vfl);
	}

	template<typename SomeInputType>
	static void readFiniteLoops_(SomeInputType& io,
	                             VectorFiniteLoopType& vfl,
	                             const MatrixFieldType& tmpMat)
	{
		for (SizeType i = 0; i < tmpMat.rows(); ++i) {
			FiniteLoop fl(tmpMat(i,0), tmpMat(i,1), tmpMat(i,2));
			vfl.push_back(fl);
		}

		readFiniteLoops_(io, vfl);
	}

	template<typename SomeInputType>
	static void readFiniteLoops_(SomeInputType& io,
	                             VectorFiniteLoopType& vfl)
	{
		SizeType repeat = 0;

		try {
			io.readline(repeat,"RepeatFiniteLoopsTimes=");
		}  catch (std::exception& e) {}

		SizeType fromFl = 0;
		try {
			io.readline(fromFl,"RepeatFiniteLoopsFrom=");
		}  catch (std::exception& e) {}

		if (vfl.size() == 0) {
			std::cerr<<"WARNING: No finite loops found\n";
		}

		SizeType upToFl = vfl.size()-1;
		try {
			io.readline(upToFl,"RepeatFiniteLoopsTo=");
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

	static void warnIfFiniteMlessThanMin(const VectorFiniteLoopType& vfl, SizeType minM)
	{
		for (SizeType i = 0; i < vfl.size(); ++i) {
			if (vfl[i].keptStates >= minM) continue;
			std::cout<<"WARNING: Triplet number "<<i<<" has m= "<<vfl[i].keptStates;
			std::cout<<" which is less than minimum m = "<<minM;
			std::cout<<" as found in TruncationTolerance\n";
		}
	}

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
		bool doNotCheck = (options.find("doNotCheckTwoSiteDmrg") != PsimagLite::String::npos);

		if (bOld == b) return;

		PsimagLite::String s (__FILE__);
		PsimagLite::String s1 = ": Previous run used twositedmrg ";
		s1 += " but this one does not.";
		PsimagLite::String s2 = ": Previous run didn't use twositedmrg ";
		s2 += " but this one does.";
		if (bOld) s += s1;
		else s += s2;
		if (doNotCheck) {
			std::cerr<<"WARNING: "<<s<<"\n";
			return;
		}

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
	if (p.options.find("KroneckerDumper") != PsimagLite::String::npos) {
		os<<"parameters.dumperBegin="<<p.dumperBegin<<"\n";
		os<<"parameters.dumperEnd="<<p.dumperEnd<<"\n";
	}

	os<<"parameters.precision="<<p.precision<<"\n";
	os<<"parameters.keptStatesInfinite="<<p.keptStatesInfinite<<"\n";
	os<<"FiniteLoops ";
	os<<p.finiteLoop;

	os<<"RecoverySave="<<p.recoverySave<<"\n";
	if (p.truncationControl.first > 0) {
		os<<"parameters.tolerance="<<p.truncationControl.first<<",";
		os<<p.truncationControl.second<<"\n";
	}

	os<<"parameters.degeneracyMax="<<p.degeneracyMax<<"\n";
	os<<"parameters.denseSparseThreshold="<<p.denseSparseThreshold<<"\n";
	os<<"parameters.nthreads="<<p.nthreads<<"\n";
	os<<"parameters.useReflectionSymmetry="<<p.useReflectionSymmetry<<"\n";
	os<<p.checkpoint;

	if (p.fileForDensityMatrixEigs!="")
		os<<"parameters.fileForDensityMatrixEigs="<<p.fileForDensityMatrixEigs<<"\n";

	if (p.options.find("MatrixVectorStored")==PsimagLite::String::npos)
		os<<"MaxMatrixRankStored="<<p.maxMatrixRankStored<<"\n";

	return os;
}
} // namespace Dmrg
/*@}*/

#endif

