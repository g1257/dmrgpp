/*
Copyright (c) 2009-2014-2019, UT-Battelle, LLC
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
#include "Io/IoSelector.h"
#include "PsimagLite.h"
#include "RestartStruct.h"
#include "FiniteLoop.h"
#include "Io/IoSerializerStub.h"
#include "Recovery.h"
#include "ProgressIndicator.h"
#include <sstream>
#include "Options.h"
#include "TruncationControl.h"
#include "AlgebraicStringToNumber.h"
#include <sys/types.h>
#include <unistd.h>

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
template<typename FieldType, typename InputValidatorType, typename QnType>
struct ParametersDmrgSolver {

	typedef ParametersDmrgSolver<FieldType, InputValidatorType, QnType> ThisType;
	typedef typename QnType::PairSizeType PairSizeType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<FieldType>::Type VectorFieldType;
	typedef PsimagLite::Matrix<FieldType> MatrixFieldType;
	typedef PsimagLite::Matrix<PsimagLite::String> MatrixStringType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef std::pair<FieldType, SizeType> PairRealSizeType;
	using FiniteLoopType = FiniteLoop<FieldType>;
	typedef typename PsimagLite::Vector<FiniteLoopType>::Type VectorFiniteLoopType;
	typedef Options<InputValidatorType> OptionsType;
	using TruncationControlType = TruncationControl<FieldType>;

	SizeType nthreads;
	SizeType nthreads2;
	SizeType sitesPerBlock;
	SizeType maxMatrixRankStored;
	SizeType keptStatesInfinite;
	SizeType dumperBegin;
	SizeType dumperEnd;
	SizeType precision;
	SizeType numberOfExcited;
	SizeType gemmRnb;
	SizeType opOnSiteThreshold;
	bool autoRestart;
	TruncationControlType truncationControl;
	PsimagLite::String filename;
	PsimagLite::String version;
	OptionsType options;
	PsimagLite::String model;
	PsimagLite::String insitu;
	PsimagLite::String recoverySave;
	PsimagLite::String printHamiltonianAverage;
	PsimagLite::String saveDensityMatrixEigenvalues;
	PsimagLite::String findSymmetrySector;
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
		ioSerializer.write(root + "/nthreads2", nthreads2);
		ioSerializer.write(root + "/sitesPerBlock", sitesPerBlock);
		ioSerializer.write(root + "/maxMatrixRankStored", maxMatrixRankStored);
		ioSerializer.write(root + "/keptStatesInfinite", keptStatesInfinite);
		ioSerializer.write(root + "/numberOfExcited", numberOfExcited);
		ioSerializer.write(root + "/gemmRnb", gemmRnb);
		ioSerializer.write(root + "/dumperBegin", dumperBegin);
		ioSerializer.write(root + "/dumperEnd", dumperEnd);
		ioSerializer.write(root + "/precision", precision);
		truncationControl.write(root + "/truncationControl", ioSerializer);
		ioSerializer.write(root + "/filename", filename);
		ioSerializer.write(root + "/version", version);
		options.write(root + "/options", ioSerializer);
		ioSerializer.write(root + "/model", model);
		ioSerializer.write(root + "/insitu", insitu);
		ioSerializer.write(root + "/recoverySave", recoverySave);
		ioSerializer.write(root + "/printHamiltonianAverage", printHamiltonianAverage);
		ioSerializer.write(root + "/saveDensityMatrixEigenvalues", saveDensityMatrixEigenvalues);
		checkpoint.write(label + "/checkpoint", ioSerializer);
		ioSerializer.write(root + "/adjustQuantumNumbers", adjustQuantumNumbers);
		ioSerializer.write(root + "/finiteLoop", finiteLoop);
		ioSerializer.write(root + "/degeneracyMax", degeneracyMax);
		ioSerializer.write(root + "/denseSparseThreshold", denseSparseThreshold);
		ioSerializer.write(root + "/opOnSiteThreshold", opOnSiteThreshold);
		ioSerializer.write(root + "/findSymmetrySector", findSymmetrySector);
	}

	//! Read Dmrg parameters from inp file
	ParametersDmrgSolver(InputValidatorType& io,
	                     PsimagLite::String sOptions,
	                     bool earlyExit = false,
	                     bool isObserveCode = false)
	    : nthreads(1),
	      nthreads2(1),
	      sitesPerBlock(1),
	      maxMatrixRankStored(0),
	      keptStatesInfinite(0),
	      dumperBegin(0),
	      dumperEnd(0),
	      precision(6),
	      numberOfExcited(1),
	      gemmRnb(0),
	      opOnSiteThreshold(0),
	      autoRestart(false),
	      options("SolverOptions=", io),
	      recoverySave(""),
	      adjustQuantumNumbers(0, QnType(false, VectorSizeType(), PairSizeType(0, 0), 0)),
	      degeneracyMax(1e-12),
	      denseSparseThreshold(0.2)
	{
		if (options.isSet("minimizeDisk"))
			options += ",noSaveWft,noSaveStacks,noSaveData";

		io.readline(model,"Model=");
		options += sOptions;
		io.readline(version,"Version=");

		bool ciRun = options.isSet("ciRun");

		try {
			io.readline(filename,"OutputFile=");
			ciRun = false;
		} catch (std::exception&) {
			filename = io.filename();
		}

		filename = filenameFromRootname(filename,
		                                options.isSet("addPidToOutputName"),
		                                ciRun);

		if (earlyExit) return;

		PsimagLite::String infLoops = "0";
		bool infLoopsIsAnInt = true;
		readFiniteAndInfiniteLoops(infLoops, infLoopsIsAnInt, io);

		if (options.isSet("hasQuantumNumbers")) {
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
		try {
			io.readline(nthreads, "Threads=");
		} catch (std::exception&) {}


		try {
			io.readline(nthreads2, "ThreadsLevelTwo=");
		} catch (std::exception&) {}

		if (nthreads == 0 || nthreads2 == 0) {
			PsimagLite::String s (__FILE__);
			s += "\nFATAL: nthreads and nthreads2 cannot be zero\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		try {
			io.readline(gemmRnb, "GemmRnb=");
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
			io.readline(numberOfExcited, "NumberOfExcited=");
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

		if (options.isSet("KroneckerDumper")) {
			if (options.isSet("MatrixVectorStored")) {
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

		try {
			io.readline(printHamiltonianAverage, "PrintHamiltonianAverage=");
		} catch (std::exception&) {}

		try {
			io.readline(saveDensityMatrixEigenvalues, "SaveDensityMatrixEigenvalues=");
		} catch (std::exception&) {}

		try {
			io.readline(opOnSiteThreshold, "OpOnSiteThreshold=");
		} catch (std::exception&) {}

		if (options.isSet("findSymmetrySector"))
			findSymmetrySector="1==1";

		PsimagLite::String tmpString;
		try {
			io.readline(tmpString, "FindSymmetrySector=");
		} catch (std::exception&) {}

		if (findSymmetrySector != "" && tmpString != "")
			err(PsimagLite::String("Either findSymmetrySector in SolverOptions or ") +
			    PsimagLite::String("FindSymmetrySector= line or neither of them, but not both\n"));

		if (tmpString != "")
			findSymmetrySector = tmpString;

		if (options.isSet("restart")) {
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
			}

			if (hasRestartFrom)
				checkpoint.setFilename(restartFrom);

			if (checkpoint.filename() == "") {
				PsimagLite::String tmp = "FATAL: RestartFilename NOT found in input ";
				err(tmp + "AND InfiniteLoopKeptStates is an int\n");
			}

			checkpoint.setFilename(filenameFromRootname(checkpoint.filename(), false, false));
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

	static bool getValueIfPresent(PsimagLite::String& str,
	                              PsimagLite::String label,
	                              InputValidatorType& io)
	{
		try {
			io.readline(str, label);
			return true;
		} catch (std::exception&) {
			return false;
		}
	}

	void readFiniteAndInfiniteLoops(PsimagLite::String& infLoops,
	                                bool& infLoopsIsAnInt,
	                                InputValidatorType& io)
	{
		try {
			io.readline(infLoops, "InfiniteLoopKeptStates=");
			std::istringstream iss(infLoops);
			iss >> keptStatesInfinite;
			infLoopsIsAnInt = (iss.eof());
		} catch (std::exception&) {
			keptStatesInfinite = 0;
		}

		truncationControl.read(io, keptStatesInfinite, options.isSet("twositedmrg"));

		readFiniteLoops(io, finiteLoop, truncationControl, -1);
	}

	void readFiniteLoops(InputValidatorType& io,
	                     VectorFiniteLoopType& vfl,
	                     const TruncationControlType& truncationC,
	                     int lastSite) const
	{
		if (io.version() < io.versionAinur()) {
			VectorStringType tmpVec;
			io.read(tmpVec,"FiniteLoops");
			readFiniteLoops_(io, vfl, tmpVec, truncationC);
		} else {
			MatrixStringType tmpMat;
			io.read(tmpMat, "FiniteLoops");
			readFiniteLoops_(io, vfl, tmpMat, truncationC, lastSite);
		}
	}

	void readFiniteLoops_(InputValidatorType& io,
	                      VectorFiniteLoopType& vfl,
	                      const VectorStringType& tmpVec,
	                      const TruncationControlType& truncationC) const
	{
		for (SizeType i = 0; i < tmpVec.size(); i += 3) {
			typename PsimagLite::Vector<int>::Type xTmp(2);
			assert(2 + i < tmpVec.size());
			for (SizeType j = 0; j < xTmp.size(); ++j)
				xTmp[j] = PsimagLite::atoi(tmpVec[i+j]);

			FiniteLoopType fl(xTmp[0], xTmp[1], tmpVec[i + 2], truncationC);
			vfl.push_back(fl);
		}

		readFiniteLoopsRepeat(io, vfl);
	}

	void readFiniteLoops_(InputValidatorType& io,
	                      VectorFiniteLoopType& vfl,
	                      const MatrixStringType& tmpMat,
	                      const TruncationControlType& truncationC,
	                      int lastSite) const
	{
		SizeType numberOfSites = 0;
		io.readline(numberOfSites, "TotalNumberOfSites=");
		bool latticeIsOdd = (numberOfSites & 1);

		typedef AlgebraicStringToNumber<FieldType> AlgebraicStringToNumberType;
		AlgebraicStringToNumberType algebraicStringToNumber("FiniteLoops", numberOfSites);

		std::cout<<"FiniteLoopLengths=[";

		for (SizeType i = 0; i < tmpMat.rows(); ++i) {
			int length = (tmpMat(i, 0) == "@auto") ? autoNumber(i, numberOfSites, lastSite)
			                                       : algebraicStringToNumber.
			                                         procLength(tmpMat(i, 0));
			SizeType m = PsimagLite::atoi(tmpMat(i, 1));
			FiniteLoopType fl(length, m, tmpMat(i, 2), truncationC);
			vfl.push_back(fl);
			if (lastSite >= 0) {
				lastSite += length;
				if (latticeIsOdd && lastSite == 1)
					lastSite = 0;
			}

			std::cout<<length;
			if (i + 1 < tmpMat.rows()) {
				std::cout<<", ";
			}
		}

		std::cout<<"];\n";

		readFiniteLoopsRepeat(io, vfl);
	}

	int autoNumber(SizeType ind, SizeType numberOfSites, int lastSite) const
	{
		// if lastSite < 0 then this isn't called from checkpoint
		// so it doesn't matter as it will be overwritten later
		if (lastSite < 0) return lastSite;

		return (ind == 0) ? autoNumberFirst(numberOfSites, lastSite)
		                  : autoNumberAfterFirst(ind, numberOfSites, lastSite);
	}

	int autoNumberFirst(SizeType numberOfSites, SizeType lastSite) const
	{
		bool isRestart = this->options.isSet("restart");
		bool allinsystem = this->options.isSet("geometryallinsystem");
		if (isRestart || allinsystem) {

			SizeType oneOrTwo = (numberOfSites & 1) ? 1 : 2;
			if (lastSite != 1 && (lastSite + oneOrTwo) != numberOfSites)
				err("@auto: Internal error; please report this problem\n");

			int x = (numberOfSites - 2);
			if (lastSite + oneOrTwo == numberOfSites) {
				x *= (-1);
			}

			return x;
		} else {
			SizeType y = (numberOfSites - 2);
			if (y & 1) ++y;
			return y/2;
		}
	}

	int autoNumberAfterFirst(SizeType ind, SizeType numberOfSites, SizeType lastSite) const
	{
		assert(ind > 0);
		assert(numberOfSites >= 2);
		SizeType x = (numberOfSites - 2);
		if (x & 1) --x;
		if (lastSite != x && lastSite != 0)
			err("autoNumberAfterFirst: Internal error, please report to DMRG++ mailing list\n");

		return (lastSite == 0) ? x : -x;
	}

	static void readFiniteLoopsRepeat(InputValidatorType& io,
	                                  VectorFiniteLoopType& vfl)
	{
		SizeType repeat = 0;

		try {
			io.readline(repeat,"RepeatFiniteLoopsTimes=");
		}  catch (std::exception&) {}

		SizeType fromFl = 0;
		try {
			io.readline(fromFl,"RepeatFiniteLoopsFrom=");
		}  catch (std::exception&) {}

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
				FiniteLoopType fl = vfl[j];
				vfl.push_back(fl);
			}
		}
	}

	static void checkRestart(PsimagLite::String filename1,
	                         PsimagLite::String filename2,
	                         const OptionsType& options)
	{
		checkFilesNotEqual(filename1, filename2);
		checkTwoSiteDmrg(filename2, options);
	}

private:

	static SizeType modeFromModel(PsimagLite::String model)
	{
		if (model == "HubbardAncilla") return 3;
		if (model == "FeAsBasedSc" ||
		        model == "HolsteinThin") return 1;
		err("Not supported modeFromModel for model= " + model + "\n");
		return 0;
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
	                             const OptionsType& options)
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
		bool b = options.isSet("twositedmrg");
		bool doNotCheck = options.isSet("doNotCheckTwoSiteDmrg");

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

	static PsimagLite::String filenameFromRootname(PsimagLite::String f,
	                                               bool addPidToOutputName,
	                                               bool ciRun)
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

		if (ciRun) {
			PsimagLite::String strNumber = findNumber(f);
			if (strNumber != "")
				f = "data" + strNumber + ".hd5";
		}

		if (addPidToOutputName) f += "." + ttos(getpid());
		return f;
	}

	static PsimagLite::String findNumber(PsimagLite::String f)
	{
		size_t findIndex = f.find(".hd5");
		if (findIndex == PsimagLite::String::npos)
			return "";

		PsimagLite::String tmp;
		for (SizeType i = 0; i < findIndex; ++i) {
			SizeType j = findIndex - i;
			if (f[j] >= 48 && f[j] <= 57)
				tmp += f[j];
			else
				break;
		}

		const SizeType len = tmp.length();
		PsimagLite::String result = tmp;
		for (SizeType i = 0; i < len; ++i) {
			SizeType j = len - i;
			result[i] = tmp[j];
		}

		return result;
	}
};
} // namespace Dmrg
/*@}*/

#endif

