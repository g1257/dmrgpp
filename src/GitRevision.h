
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
#include <unistd.h>

namespace PsimagLite {

class GitRevision {

	static const size_t MAX_LINE  = 4096;

public:
	
	GitRevision(const std::string& dir,const std::string& prefix)
	: prefix_(prefix),tempFile_("version.txt"),gitCommand_("git"),ret_(0),retp_(0)
	{
		std::string currentPath = getCurrentPath();
		ret_ = chdir(dir.c_str());
		backupTemporaryFile();
		id_ = getId();
		diff_ = hasDiff();
		ret_ = chdir(currentPath.c_str());
	}

	const std::string& id() const { return id_; }
	
	std::string diff() const { return diff_; }

	std::string prefix() const { return prefix_; }

private:

	void backupTemporaryFile()
	{
		std::ifstream fin(tempFile_.c_str());
		if (!fin) return;
		std::string cmd1 = "cp " + tempFile_ + "   " + tempFile_ + ".bak";
		ret_ = system(cmd1.c_str());
		fin.close();
	}

	std::string getId()
	{
		std::string cmd = gitCommand_ +   "  rev-parse HEAD > " + tempFile_ + " 2>/dev/null ";
		//std::cerr<<"Command is "<<cmd<<"\n";
		ret_ = system(cmd.c_str());
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
		if (id_ == "unknown") return "unknown";
		std::string cmd = gitCommand_ + "  diff --shortstat > " + tempFile_ + " 2>/dev/null ";
		ret_ = system(cmd.c_str());
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

	std::string getCurrentPath()
	{
		size_t bufsize = 4096;
		char* buf = new char[bufsize];
		retp_ = getcwd(buf,bufsize);

		std::string result(buf);
		delete [] buf;
		return result;
	}

	std::string prefix_;
	std::string tempFile_;
	std::string gitCommand_;
	int ret_;
	char *retp_;
	std::string id_;
	std::string diff_;

}; // class GitRevision

std::ostream& operator<<(std::ostream& os,const GitRevision& gitRev)
{
	os<<"const static char *"<<gitRev.prefix()<<"Revision=\""<<gitRev.id()<<"\";\n";
	os<<"const static char *"<<gitRev.prefix()<<"Diff=\""<<gitRev.diff()<<"\";\n";
	return os;
}

//int main(int argc,char *argv[])
//{
//	GitRevision gitRev;
//	std::cout<<gitRev;
//}
} // namespace PsimagLite

/*@}*/
#endif  //GIT_REVISION_H_

