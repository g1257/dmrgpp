/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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
#include "Io/IoSelector.h"
#include "PsimagLite.h"
#include "RestartStruct.h"
#include "FiniteLoop.h"
#include "Io/IoSerializerStub.h"
#include "Recovery.h"
#include "ProgressIndicator.h"
#include <sstream>

namespace Dmrg {

/* PSIDOC ParametersDmrgSolver
\begin{itemize}
\item[Model=string]
A string indicating the model, be it HubbardOneBand Heisenberg, etc.

\item[version=string]  A mandatory string that is read and ignored. Usually contains the
result
of doing \verb!git rev-parse HEAD!.

\item[outputfile=string]  The root for the output file.
If IoNg is used the output file will be added the hdf5 extension.
If IoNg is not used then no extension will be added.
This file will be created if non-existent,
 and if it
exits it will be truncated.

\item[InfiniteLoopKeptStates] If an integer then this is
 the \emph{m} value for the infinite algorithm. Else it is the name
 of the filename to restart from. If the filename is numeric, use double
 quotes around it to have it interpreted as a string.

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
template<typename FieldType,typename InputValidatorType, typename QnType>
struct ParametersDmrgSolver {

	typedef ParametersDmrgSolver<FieldType, InputValidatorType, QnType> ThisType;
	typedef typename QnType::PairSizeType PairSizeType;
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
	bool autoRestart;
	PairRealSizeType truncationControl;
	PsimagLite::String filename;
	PsimagLite::String version;
	PsimagLite::String options;
	PsimagLite::String model;
	PsimagLite::String insitu;
	PsimagLite::String fileForDensityMatrixEigs;
	PsimagLite::String recoverySave;
	RestartStruct checkpoint;
	typename QnType::VectorQnType adjustQuantumNumbers;
	VectorFiniteLoopType finiteLoop;
	FieldType degeneracyMax;
	FieldType denseSparseThreshold;

	void write(PsimagLite::String label,
	           PsimagLite::IoSerializer& ioSerializer) const
	{
		PsimagLite::String root = label;

		ioSerializer.createGroup(root);

		ioSerializer.write(root + "/nthreads", nthreads);
		ioSerializer.write(root + "/sitesPerBlock", sitesPerBlock);
		ioSerializer.write(root + "/maxMatrixRankStored", maxMatrixRankStored);
		ioSerializer.write(root + "/keptStatesInfinite", keptStatesInfinite);
		ioSerializer.write(root + "/excited", excited);
		ioSerializer.write(root + "/dumperBegin", dumperBegin);
		ioSerializer.write(root + "/dumperEnd", dumperEnd);
		ioSerializer.write(root + "/precision", precision);
		ioSerializer.write(root + "/useReflectionSymmetry", useReflectionSymmetry);
		ioSerializer.write(root + "/truncationControl", truncationControl);
		ioSerializer.write(root + "/filename", filename);
		ioSerializer.write(root + "/version", version);
		ioSerializer.write(root + "/options", options);
		ioSerializer.write(root + "/model", model);
		ioSerializer.write(root + "/insitu", insitu);
		ioSerializer.write(root + "/fileForDensityMatrixEigs", fileForDensityMatrixEigs);
		ioSerializer.write(root + "/recoverySave", recoverySave);
		checkpoint.write(label + "/checkpoint", ioSerializer);
		ioSerializer.write(root + "/adjustQuantumNumbers", adjustQuantumNumbers);
		ioSerializer.write(root + "/finiteLoop", finiteLoop);
		ioSerializer.write(root + "/degeneracyMax", degeneracyMax);
		ioSerializer.write(root + "/denseSparseThreshold", denseSparseThreshold);
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg) const
	{
		return 0;
	}

	//! Read Dmrg parameters from inp file
	ParametersDmrgSolver(InputValidatorType& io,
	                     PsimagLite::String sOptions,
	                     bool earlyExit = false,
	                     bool isObserveCode = false)
	    : sitesPerBlock(1),
	      maxMatrixRankStored(0),
	      keptStatesInfinite(0),
	      excited(0),
	      dumperBegin(0),
	      dumperEnd(0),
	      precision(6),
	      autoRestart(false),
	      recoverySave("no"),
	      adjustQuantumNumbers(0, QnType(false, VectorSizeType(), PairSizeType(0, 0), 0)),
	      degeneracyMax(1e-12),
	      denseSparseThreshold(0.2)
	{
		io.readline(model,"Model=");
		io.readline(options,"SolverOptions=");
		options += sOptions;
		io.readline(version,"Version=");

		try {
			io.readline(filename,"OutputFile=");
		} catch (std::exception&) {
			filename = io.filename();
		}

		filename = filenameFromRootname(filename);

		if (earlyExit) return;

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

		tmpVector.clear();
		try {
			io.read(tmpVector,"AdjustQuantumNumbers");
		} catch (std::exception&) {}

		if (tmpVector.size() > 0)
			QnType::adjustQns(adjustQuantumNumbers,
			                  tmpVector,
			                  modeFromModel(model));

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
			if (options.find("twositedmrg") == PsimagLite::String::npos) {
				std::cerr<<"WARNING: TruncationTolerance used without twositedmrg\n";
				std::cout<<"WARNING: TruncationTolerance used without twositedmrg\n";
			}
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
		PsimagLite::String restartFrom;
		bool hasRestartFrom = getValueIfPresent(restartFrom, "RestartFilename=", io);
		PsimagLite::String infLoops = "0";
		bool infLoopsIsAnInt = true;
		try {
			io.readline(infLoops, "InfiniteLoopKeptStates=");
			std::istringstream iss(infLoops);
			iss >> keptStatesInfinite;
			infLoopsIsAnInt = (iss.eof());
		} catch (std::exception&) {
			keptStatesInfinite = 0;
		}

		if (options.find("restart")!=PsimagLite::String::npos) {
			if (!infLoopsIsAnInt and hasRestartFrom) {
				PsimagLite::String tmp = "FATAL: RestartFilename found in input ";
				err(tmp + "while InfiniteLoopKeptStates not an int\n");
			}

			checkpoint.setFilename("");
			if (!infLoopsIsAnInt) {
				checkpoint.setFilename(infLoops);
				// remove double quotes if present
				SizeType begin = (infLoops[0] == '"') ? 1 : 0;
				SizeType last = infLoops.length();
				assert(last > 0);
				--last;
				SizeType end = (infLoops[last] == '"') ? last : last + 1;
				checkpoint.setFilename(infLoops.substr(begin, end - begin));
			} else {
				if (keptStatesInfinite > 0) {
					PsimagLite::String tmp = "WARNING: The numeric value of ";
					tmp += "InfiniteLoopKeptStates will be ignored\n";
					std::cerr<<tmp;
					std::cout<<tmp;
				}
			}

			if (hasRestartFrom)
				checkpoint.setFilename(restartFrom);

			if (checkpoint.filename() == "") {
				PsimagLite::String tmp = "FATAL: RestartFilename NOT found in input ";
				err(tmp + "AND InfiniteLoopKeptStates is an int\n");
			}

			checkpoint.setFilename(filenameFromRootname(checkpoint.filename()));
			checkRestart(filename, checkpoint.filename(), options);
			hasRestart = true;
		} else {
			if (hasRestartFrom) {
				PsimagLite::String tmp = "FATAL: RestartFilename found in input ";
				err(tmp + "but no restart found in SolverOptions.\n");
			}

			if (!infLoopsIsAnInt) {
				PsimagLite::String tmp = "FATAL: InfiniteLoopKeptStates not an integer ";
				err(tmp + "but no restart found in SolverOptions.\n");
			}

			if (keptStatesInfinite == 0) {
				PsimagLite::String tmp = "FATAL: InfiniteLoopKeptStates must be ";
				err(tmp + "a positive integer for a NON restart run.\n");
			}
		}

		Recovery<ThisType, int>::checkOptions(recoverySave, options);
		Recovery<ThisType, int>::autoRestart(*this);

		if (hasRestart) {
			checkpoint.read(io);
		}
	}

	template<typename SomeInputType>
	static bool getValueIfPresent(PsimagLite::String& str,
	                              PsimagLite::String label,
	                              SomeInputType& io)
	{
		try {
			io.readline(str, label);
			return true;
		} catch (std::exception&) {
			return false;
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
			io.read(tmpMat, "FiniteLoops");
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

	static void checkRestart(PsimagLite::String filename1,
	                         PsimagLite::String filename2,
	                         PsimagLite::String options)
	{
		checkFilesNotEqual(filename1, filename2);
		checkTwoSiteDmrg(filename2, options);
	}

	// print dmrg parameters
	friend std::ostream& operator<<(std::ostream& os,
	                                const ParametersDmrgSolver& p)
	{
		os<<"This is DMRG++\n";
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

private:

	static SizeType modeFromModel(PsimagLite::String model)
	{
		if (model == "HubbardAncilla") return 3;
		if (model == "FeAsBasedSc") return 1;
		err("Not supported modeFromModel for model= " + model + "\n");
		return 0;
	}

	static void warnIfFiniteMlessThanMin(const VectorFiniteLoopType& vfl, SizeType minM)
	{
		for (SizeType i = 0; i < vfl.size(); ++i) {
			if (vfl[i].keptStates >= minM) continue;
			std::cout<<"WARNING: Triplet number "<<i<<" has m= "<<vfl[i].keptStates;
			std::cout<<" which is less than minimum m = "<<minM;
			std::cout<<" as found in TruncationTolerance\n";
		}
	}

	static void checkFilesNotEqual(PsimagLite::String filename1,
	                               PsimagLite::String filename2)
	{
		if (filename1 != filename2) return;
		PsimagLite::String s (__FILE__);
		s += "\nFATAL: " + filename1 + "is equal to " + filename2 + "\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}

	static void checkTwoSiteDmrg(PsimagLite::String filename2,
	                             PsimagLite::String options)
	{
		PsimagLite::IoSelector::In io(filename2);
		PsimagLite::String optionsOld;
		try {
			io.read(optionsOld,"PARAMETERS/options");
		} catch (...) {
			std::cerr<<"WARNING: could not read PARAMETERS/options from ";
			std::cerr<<io.filename()<<"\n";
			return;
		}

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

	PsimagLite::String filenameFromRootname(PsimagLite::String f) const
	{
		size_t findIndex = f.find(".txt");
		if (findIndex != PsimagLite::String::npos)
			f.replace(findIndex, PsimagLite::String(".txt").length(), ".hd5");

		findIndex = f.find(".inp");
		if (findIndex != PsimagLite::String::npos)
			f.replace(findIndex, PsimagLite::String(".inp").length(), ".hd5");

		findIndex = f.find(".ain");
		if (findIndex != PsimagLite::String::npos)
			f.replace(findIndex, PsimagLite::String(".ain").length(), ".hd5");

		findIndex = f.find(".hd5");
		if (findIndex == PsimagLite::String::npos)
			f += ".hd5";

		return f;
	}
};
} // namespace Dmrg
/*@}*/

#endif

