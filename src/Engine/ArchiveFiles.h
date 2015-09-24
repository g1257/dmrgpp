/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

/*! \file ArchiveFiles.h
 *
 *
 */

#ifndef DMRG_ARCHIVE_FILES_H
#define DMRG_ARCHIVE_FILES_H

#include "Vector.h"
#include "ProgramGlobals.h"
#include "ProgressIndicator.h"
#include "TarPack.h"
#include "IoSimple.h"

namespace Dmrg {

template<typename ParametersType>
class ArchiveFiles  {

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	ArchiveFiles(const ParametersType& parameters,
	             PsimagLite::String filename,
	             bool addExtra,
	             PsimagLite::String extra)
	    : parameters_(parameters)
	{
		files_.push_back(filename);
		PsimagLite::String rootname = parameters.filename;
		appendToList(files_,rootname);
		if (!addExtra && extra != "")
			files_.push_back(extra);
	}

	void deletePackedFiles()
	{
		if (parameters_.options.find("tarEnable") == PsimagLite::String::npos)
			return;

		PsimagLite::String tarname = tarName(parameters_.filename);
		Dmrg::TarPack tarPack(tarname);
		for (SizeType i = 0; i < files_.size(); ++i)
			tarPack.add(files_[i]);

		if (parameters_.options.find("tarNoDelete") != PsimagLite::String::npos)
			return;

		for (SizeType i = 1; i < files_.size(); ++i)
			unlink(files_[i].c_str());
	}

	static void unpackIfNeeded(PsimagLite::String rootname)
	{
		PsimagLite::IoSimple::In io;
		try {
			io.open(rootname);
		} catch (std::exception&) {
			appendToList(filesToDelete_,rootname);
			UnTarPack untarPack(tarName(rootname));

			for (SizeType i = 0; i < filesToDelete_.size(); ++i)
				untarPack.extract<NormalExtract>(filesToDelete_[i],false,0);
		}
	}

	void listOrClear(PsimagLite::String filename, PsimagLite::String what) const
	{
		if (what == "DELETE")
			clearFiles(filename);
		else if (what == "list")
			listFiles(filename);
		else
			throw PsimagLite::RuntimeError("-F DELETE | list\n");
	}

	static void staticDelete()
	{
		for (SizeType i = 0; i < filesToDelete_.size(); ++i)
			unlink(filesToDelete_[i].c_str());
		filesToDelete_.clear();
	}

	static PsimagLite::String rootName(PsimagLite::String filename)
	{
		PsimagLite::String rootname = filename;
		size_t index =rootname.find(".", 0);
		if (index != PsimagLite::String::npos) {
			rootname.erase(index,filename.length());
		}

		return rootname;
	}

	static PsimagLite::String coutName(PsimagLite::String filename)
	{
		PsimagLite::String rootname = filename;
		size_t index =rootname.find(".", 0);
		if (index != PsimagLite::String::npos) {
			rootname.erase(index,filename.length());
		}

		return "runFor" + rootname + ".cout";
	}

private:

	void clearFiles(PsimagLite::String filename) const
	{
		for (SizeType i = 0; i < files_.size(); ++i) {
			if (files_[i] == filename) continue;
			unlink(files_[i].c_str());
		}

		PsimagLite::String coutname = coutName(filename);
		unlink(coutname.c_str());
	}

	void listFiles(PsimagLite::String filename) const
	{
		for (SizeType i = 0; i < files_.size(); ++i) {
			if (files_[i] == filename) continue;
			std::cout<<files_[i]<<"\n";
		}

		PsimagLite::String coutname = coutName(filename);
		std::cout<<coutname<<"\n";
	}

	static PsimagLite::String tarName(PsimagLite::String rootname)
	{
		PsimagLite::String tarname = rootname;
		size_t index = tarname.find(".");
		if (index != PsimagLite::String::npos)
			tarname.erase(index,rootname.length());

		return tarname + ".tar";
	}

	static void appendToList(VectorStringType& files,
	                         PsimagLite::String rootname)
	{
		files.push_back(rootname);
		files.push_back(ProgramGlobals::SYSTEM_STACK_STRING+rootname);
		files.push_back(ProgramGlobals::ENVIRON_STACK_STRING+rootname);
		files.push_back(ProgramGlobals::WFT_STRING+rootname);
	}

	static VectorStringType filesToDelete_;
	const ParametersType& parameters_;
	VectorStringType files_;
}; //class ArchiveFiles

template<typename T>
typename ArchiveFiles<T>::VectorStringType ArchiveFiles<T>::filesToDelete_;

} // namespace Dmrg
/*@}*/
#endif

