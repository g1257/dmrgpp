
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2011, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
/** \ingroup PsimagLite */
/*@{*/

/*! \file GitRevision.h
 *
 *
 */


#ifndef GIT_REVISION_H_
#define GIT_REVISION_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

namespace PsimagLite {

class GitRevision {

	static const size_t MAX_LINE  = 4096;

public:
	
	GitRevision()
	: tempFile_("version.txt"),gitCommand_("git")
	{
		backupTemporaryFile();
		id_ = getId();
		diff_ = hasDiff();
	}

	const std::string& id() const { return id_; }
	
	std::string diff() const { return diff_; }

private:

	void backupTemporaryFile()
	{
		std::ifstream fin(tempFile_.c_str());
		if (!fin) return;
		std::string cmd1 = "cp " + tempFile_ + "   " + tempFile_ + ".bak";
		system(cmd1.c_str());
		fin.close();
	}

	std::string getId()
	{
		std::string cmd = gitCommand_ +   "  rev-parse HEAD > " + tempFile_ + " 2>/dev/null ";
		system(cmd.c_str());
		std::string revision = "unknown";
		std::ifstream fin(tempFile_.c_str());
		if (!fin || fin.bad()) return revision;
		fin>>revision;
		fin.close();
		if (revision.find("command not found")!=std::string::npos)
		       return "unknown";
		return revision;
	}

	std::string hasDiff()
	{
		std::string cmd = gitCommand_ + "  diff --shortstat > " + tempFile_ + " 2>/dev/null ";
		system(cmd.c_str());
		std::ifstream fin(tempFile_.c_str());
		if (!fin || fin.bad()) return "unknown";
		char* diff = new char[MAX_LINE];
		fin.getline(diff,MAX_LINE);
		std::string sdiff(diff);
		delete [] diff;
		fin.close();
		if (sdiff.find("command not found")!=std::string::npos)
			return "unknown";

		return sdiff;
	}

	std::string tempFile_;
	std::string gitCommand_;
	std::string id_;
	std::string diff_;

}; // class GitRevision

std::ostream& operator<<(std::ostream& os,const GitRevision& gitRev)
{
	os<<"#revision="<<gitRev.id()<<"\n";
	os<<"#HasDiff="<<gitRev.diff()<<"\n";
}

//int main(int argc,char *argv[])
//{
//	GitRevision gitRev;
//	std::cout<<gitRev;
//}
} // namespace PsimagLite

/*@}*/
#endif  //GIT_REVISION_H_

