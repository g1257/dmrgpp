
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
#include <cstdlib>
#include <unistd.h>
#include "String.h"

namespace PsimagLite {

class GitRevision {

	static const SizeType MAX_LINE  = 4096;

public:

	GitRevision(const String& dir,const String& prefix)
	    : prefix_(prefix),tempFile_("version.txt"),gitCommand_("git"),ret_(0),retp_(0)
	{
		String currentPath = getCurrentPath();
		ret_ = chdir(dir.c_str());
		backupTemporaryFile();
		id_ = getId();
		diff_ = hasDiff();
		ret_ = chdir(currentPath.c_str());
	}

	const String& id() const { return id_; }

	String diff() const { return diff_; }

	String prefix() const { return prefix_; }

private:

	void backupTemporaryFile()
	{
		std::ifstream fin(tempFile_.c_str());
		if (!fin) return;
		String cmd1 = "cp " + tempFile_ + "   " + tempFile_ + ".bak";
		ret_ = system(cmd1.c_str());
		fin.close();
	}

	String getId()
	{
		String cmd = gitCommand_ + "  rev-parse HEAD > " + tempFile_ + " 2>/dev/null ";
		ret_ = system(cmd.c_str());
		String revision = "unknown";
		std::ifstream fin(tempFile_.c_str());
		if (!fin || fin.bad()) return revision;
		fin>>revision;
		fin.close();
		if (revision.find("command not found")!=String::npos)
			return "unknown";
		return revision;
	}

	String hasDiff()
	{
		if (id_ == "unknown") return "unknown";
		String cmd = gitCommand_ + "  diff --shortstat > " + tempFile_ + " 2>/dev/null ";
		ret_ = system(cmd.c_str());
		std::ifstream fin(tempFile_.c_str());
		if (!fin || fin.bad()) return "unknown";
		char* diff = new char[MAX_LINE];
		fin.getline(diff,MAX_LINE);
		String sdiff(diff);
		delete [] diff;
		fin.close();
		if (sdiff.find("command not found")!=String::npos)
			return "unknown";

		return sdiff;
	}

	String getCurrentPath()
	{
		SizeType bufsize = 4096;
		char* buf = new char[bufsize];
		retp_ = getcwd(buf,bufsize);

		String result(buf);
		delete [] buf;
		return result;
	}

	String prefix_;
	String tempFile_;
	String gitCommand_;
	int ret_;
	char *retp_;
	String id_;
	String diff_;

}; // class GitRevision

std::ostream& operator<<(std::ostream& os,const GitRevision& gitRev)
{
	os<<"const static char *"<<gitRev.prefix()<<"Revision=\""<<gitRev.id()<<"\";\n";
	os<<"const static char *"<<gitRev.prefix()<<"Diff=\""<<gitRev.diff()<<"\";\n";
	return os;
}

} // namespace PsimagLite

/*@}*/
#endif  //GIT_REVISION_H_

