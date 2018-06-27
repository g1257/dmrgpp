/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

/*! \file Recovery.h
 *
 *
 */

#ifndef DMRG_RECOVER_H
#define DMRG_RECOVER_H

#include "Checkpoint.h"
#include "Vector.h"
#include "ProgramGlobals.h"
#include "ProgressIndicator.h"
#include <fstream>
#include <sys/types.h>
#include <dirent.h>

namespace Dmrg {

template<typename ParametersType,typename CheckpointType>
class Recovery  {

	enum OptionEnum {DISABLED, BY_DELTATIME, BY_LOOP};

	struct OptionSpec {

		OptionSpec() : optionEnum(BY_LOOP), value(1) {}

		OptionEnum optionEnum;
		SizeType value;
	};

	struct OpaqueRestart {

		OpaqueRestart() : loopIndex(0), stepCurrent(0)
		{}

		SizeType loopIndex;
		SizeType stepCurrent;
	};

public:

	enum {SYSTEM = ProgramGlobals::SYSTEM, ENVIRON = ProgramGlobals::ENVIRON};

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename CheckpointType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename CheckpointType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename CheckpointType::IoType IoType;
	typedef typename CheckpointType::TargetingType TargetingType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename CheckpointType::MemoryStackType MemoryStackType;
	typedef typename CheckpointType::DiskStackType DiskStackType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorBlockType;

	Recovery(const VectorBlockType& siteIndices,
	         typename IoType::Out& ioOut,
	         const CheckpointType& checkpoint,
	         const WaveFunctionTransfType& wft,
	         const BasisWithOperatorsType& pS,
	         const BasisWithOperatorsType& pE)
	    : progress_("Recovery"),
	      siteIndices_(siteIndices),
	      checkpoint_(checkpoint),
	      wft_(wft),
	      pS_(pS),
	      pE_(pE),
	      savedTime_(0),
	      counter_(0)
	{
		procOptions();

		if (checkpoint_.parameters().options.find("recoveryEnableRead") ==
		        PsimagLite::String::npos) return;

		if (!checkpoint_.parameters().autoRestart) return;

		readRecovery();

//		assert(ioOut.filename() == checkpoint_.parameters().filename);
//		ioOut.close();
//		copyFile(checkpoint_.parameters().filename,
//		         checkpoint_.parameters().checkpoint.filename);
//		ioOut.open(checkpoint_.parameters().filename, IoType::ACC_RDW);

		VectorStringType parts;
		makeThreeParts(parts, checkpoint_.parameters().checkpoint.filename);
		if (parts.size() == 3) counter_ = 1 + atoi(parts[1].c_str());
	}

	~Recovery()
	{
		// NEEDS TO CLEAN ALL Recovery + digit + filename files
		// FIXME TODO
		for (SizeType i = 0; i < files_.size(); ++i)
			unlink(files_[i].c_str());
	}

	SizeType indexOfFirstFiniteLoop() const
	{
		return opaqueRestart_.loopIndex;
	}

	SizeType stepCurrent(ProgramGlobals::DirectionEnum direction) const
	{
		return (checkpoint_.parameters().autoRestart) ? opaqueRestart_.stepCurrent :
		                                                nonRecoveryStepCurrent(direction);
	}

	// this function is called before the ctor
	static void autoRestart(ParametersType& params)
	{
		if (params.options.find("recoveryEnableRead") == PsimagLite::String::npos)
			return;

		// params.filename must have been corrected already if necessary
		PsimagLite::String recoveryFile = getRecoveryFile(params.filename);
		if (recoveryFile == "")
			return;

		// *  add the line RestartFilename= pointing to the data file of the
		// run to be restarted.
		params.checkpoint.filename = recoveryFile;
		//params.checkRestart(params.filename, recoveryFile, params.options, "INTERNAL=");

		params.autoRestart = true;
	}

	bool byLoop(SizeType loopIndex) const
	{
		return (optionSpec_.optionEnum == BY_LOOP &&
		        (loopIndex % optionSpec_.value) == 0);
	}

	bool byTime() const
	{
		if (optionSpec_.optionEnum != BY_DELTATIME) return false;

		bool firstCall = (savedTime_ == 0);
		SizeType time = PsimagLite::ProgressIndicator::time();
		SizeType deltaTime = time - savedTime_;
		savedTime_ = time;
		return (!firstCall && deltaTime > optionSpec_.value);
	}

	void write(const TargetingType& psi,
	           SizeType loopIndex,
	           SizeType stepCurrent,
	           int lastSign,
	           typename IoType::Out& ioOutCurrent) const
	{
		PsimagLite::String prefix(recoveryFilePrefix());
		prefix += ttos(counter_++);
		PsimagLite::String savedName(prefix + checkpoint_.parameters().filename);
		files_.push_back(savedName);
		ioOutCurrent.flush();

		//copyFile(savedName.c_str(), ioOutCurrent.filename());

		typename IoType::Out ioOut(savedName, IoType::ACC_TRUNC);

		writeEnergies(ioOut, ioOutCurrent.filename());

		writeRecovery(ioOut, loopIndex, stepCurrent);

		// taken from end of finiteDmrgLoops
		checkpoint_.write(pS_, pE_, ioOut);
		ioOut.createGroup("FinalPsi");
		psi.write(siteIndices_[stepCurrent], ioOut, "FinalPsi");
		ioOut.write(lastSign, "LastLoopSign");

		// wft dtor
		wft_.write(ioOut);

		ioOut.close();

		// checkpoint stacks
		checkpoint_.checkpointStacks(savedName);

		if (counter_ >= checkpoint_.parameters().recoveryMaxFiles)
			counter_ = 0;
	}

private:

	static PsimagLite::String recoveryFilePrefix() { return "Recovery"; }

	void procOptions()
	{
		PsimagLite::String str = checkpoint_.parameters().recoverySave;

		if (str == "") return;

		if (str == "no") {
			optionSpec_.optionEnum = DISABLED;
			return;
		}

		if (str.length() < 3) dieWithError(str);

		if (str[0] == 'l' && str[1] == '%') {
			PsimagLite::String each = str.substr(2, str.length() - 2);
			optionSpec_.optionEnum = BY_LOOP;
			optionSpec_.value = atoi(each.c_str());
			std::cerr<<"Recovery by loop every "<<optionSpec_.value<<" loops\n";
			return;
		}

		if (str.length() < 4) dieWithError(str);

		if (str[0] == 'd' && str[1] == 't' && str[2] == '>') {
			PsimagLite::String each = str.substr(3, str.length() - 3);
			optionSpec_.optionEnum = BY_DELTATIME;
			optionSpec_.value = atoi(each.c_str());
			std::cerr<<"Recovery by delta time greater than "<<optionSpec_.value<<" seconds\n";
			return;
		}

		dieWithError(str);
	}

	void dieWithError(PsimagLite::String str) const
	{
		err("Syntax error for RecoverySave expression " + str + "\n");
	}

	static PsimagLite::String getRecoveryFile(PsimagLite::String filename)
	{
		const PsimagLite::String prefix = recoveryFilePrefix();
		std::vector<PsimagLite::String> files;
		listFilesInDirectory(files, ".");

		if (files.size() == 0) return "";

		PsimagLite::String saved("");
		SizeType max = 0;

		for (SizeType i = 0; i < files.size(); ++i) {
			std::vector<PsimagLite::String> parts;
			makeThreeParts(parts, files[i]);
			if (parts.size() != 3 || parts[0] != prefix || parts[2] != filename)
				continue;
			SizeType counter = atoi(parts[1].c_str());
			if (counter >= max) {
				max = counter;
				saved = files[i];
			}
		}

		return saved;
	}

	static void makeThreeParts(std::vector<PsimagLite::String>& parts,
	                           PsimagLite::String filename)
	{
		const PsimagLite::String prefix = recoveryFilePrefix();
		const SizeType len = prefix.length();
		if (filename.substr(0, len) != prefix) return;
		parts.push_back(prefix);

		PsimagLite::String buffer("");
		for (SizeType i = len; i < filename.length(); ++i) {
			if (isAdigit(filename[i])) buffer += filename[i];
			break;
		}

		if (buffer == "") return;

		parts.push_back(buffer);

		SizeType lastPartLen = filename.length() - buffer.length() - len;
		parts.push_back(filename.substr(len + buffer.length(), lastPartLen));
	}

	static bool isAdigit(char c)
	{
		return (c > 47 && c < 58);
	}

	static void listFilesInDirectory(std::vector<PsimagLite::String>& files,
	                                 PsimagLite::String path)
	{
		DIR* dir = 0;
		dirent* ent = 0;
		if ((dir = opendir(path.c_str())) != 0) {
			while ((ent = readdir(dir)) != 0) {
				files.push_back(ent->d_name);
			}

			closedir(dir);
			return;
		}

		/* could not open directory */
		perror("");
	}

	static void copyFile(PsimagLite::String destName, PsimagLite::String sourceName)
	{
		std::ifstream source(sourceName.c_str(), std::ios::binary);
		std::ofstream dest(destName.c_str(), std::ios::binary);
		dest << source.rdbuf();
		source.close();
		dest.close();
	}

	void writeRecovery(typename IoType::Out& ioOut,
	                   SizeType loopIndex,
	                   SizeType stepCurrent) const
	{
		ioOut.createGroup("Recovery");

		ioOut.write(loopIndex, "Recovery/loopIndex")	;
		ioOut.write(stepCurrent, "Recovery/stepCurrent");
	}

	void readRecovery()
	{
		typename IoType::In ioIn2(checkpoint_.parameters().checkpoint.filename);

		ioIn2.read(opaqueRestart_.loopIndex, "Recovery/loopIndex");
		ioIn2.read(opaqueRestart_.stepCurrent, "Recovery/stepCurrent");
		ioIn2.close();
	}

	// set initial site to add to either system or environment:
	// this is a bit tricky and has been a source of endless bugs
	// basically we have pS on the left and pE on the right,
	// and we need to determine which site is to be added
	int nonRecoveryStepCurrent(ProgramGlobals::DirectionEnum direction) const
	{
		// all right, now we can get the actual site to add:
		SizeType sitesPerBlock = checkpoint_.parameters().sitesPerBlock;
		VectorSizeType siteToAdd(sitesPerBlock);
		// left-most site of pE
		for (SizeType j = 0; j < sitesPerBlock; ++j)
			siteToAdd[j] = pE_.block()[j];

		if (direction == ProgramGlobals::EXPAND_ENVIRON) {
			// right-most site of pS
			for (SizeType j = 0; j < sitesPerBlock; ++j)
				siteToAdd[j] = pS_.block()[pS_.block().size() - 1 - j];
		}

		// now stepCurrent_ is such that sitesIndices_[stepCurrent_] = siteToAdd
		// so:
		int sc = PsimagLite::isInVector(siteIndices_, siteToAdd);

		if (sc < 0)
			err("nonRecoveryStepCurrent(...): step current error\n");

		return sc; // phew!!, that's all folks, now bugs, go away!!
	}

	void writeEnergies(typename IoType::Out& ioOut,
	                   PsimagLite::String file) const
	{
		PsimagLite::String energyLabel = checkpoint_.parameters().checkpoint.labelForEnergy;
		ioOut.flush();
		typename IoType::In ioIn(file);
		SizeType total = 0;
		ioIn.read(total, energyLabel + "/Size");
		if (total == 0)
			err("writeEnergies: no energies?\n");

		ioOut.createGroup(energyLabel);

		ioOut.write(total, energyLabel + "/Size");

		for (SizeType i = 0; i < total; ++i) {
			typename CheckpointType::RealType x = 0.0;
			ioIn.read(x, energyLabel + "/" + ttos(i));
			ioOut.write(x, energyLabel + "/" + ttos(i));
		}

		ioIn.close();
	}

	PsimagLite::ProgressIndicator progress_;
	OptionSpec optionSpec_;
	OpaqueRestart opaqueRestart_;
	const VectorBlockType& siteIndices_;
	const CheckpointType& checkpoint_;
	const WaveFunctionTransfType& wft_;
	const BasisWithOperatorsType& pS_;
	const BasisWithOperatorsType& pE_;
	mutable SizeType savedTime_;
	mutable SizeType counter_;
	mutable VectorStringType files_;
};     //class Recovery

} // namespace Dmrg
/*@}*/
#endif

