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

/*! \file TarPack.h
 *
 * This class inspired from ideas at
 * https://github.com/lindenb/cclindenb/
 * blob/28f2164dbed87cbec51839d0cf1e2a2a8e563788/src/core/lindenb/io/tarball.h
 *
 * TarPack: creating a tar archive
#include "TarPack.h"

int main(int argc,char** argv)
{
	if (argc < 2) {
		std::cerr<<"USAGE: "<<argv[0]<<" file.tar\n";
		return 1;
	}

	Dmrg::UnTarPack untarPack(argv[1]);
	Dmrg::UnTarPack::VectorStringType files;
	untarPack.list(files);
	for (SizeType i = 0; i < files.size(); ++i)
		std::cout<<files[i]<<"\n";
}
 *
 *
 * TarPack: extracting a tar archive
 *
#include "TarPack.h"

int main(int argc,char** argv)
{
	if (argc < 2) {
		std::cerr<<"USAGE: "<<argv[0]<<" file.tar\n";
		return 1;
	}

	Dmrg::UnTarPack untarPack(argv[1]);
	Dmrg::UnTarPack::VectorStringType files;
	untarPack.list(files);
	for (SizeType i = 0; i < files.size(); ++i)
		std::cout<<files[i]<<"\n";
}
 *
 *
 */

#ifndef DMRG_TAR_PACK_H
#define DMRG_TAR_PACK_H
#include "Vector.h"
#include <iostream>
#include <fstream>
#include <cstring>

namespace Dmrg {

struct PosixTarHeader {

	typedef long long unsigned int LongType;

	PosixTarHeader(PsimagLite::String filename,
	               LongType fileSize)
	{
		std::memset(this,0,sizeof(PosixTarHeader));
		if (filename == "" && fileSize == 0) return;

		std::sprintf(magic,"ustar  ");
		std::sprintf(mtime,"%011lo",time(0));
		std::sprintf(mode,"%07o",0644);
		std::snprintf(uname,32,"nobody");
		std::sprintf(gname,"nobody");
		std::snprintf(name,100,"%s",filename.c_str());
		typeflag[0]=0;
		std::sprintf(size,"%011llo",static_cast<LongType>(fileSize));
	}

	PosixTarHeader (std::ifstream& fin)
	{
		fin.read((char*)(this),sizeof(*this));
	}

	char name[100];
	char mode[8];
	char uid[8];
	char gid[8];
	char size[12];
	char mtime[12];
	char checksum[8];
	char typeflag[1];
	char linkname[100];
	char magic[6];
	char version[2];
	char uname[32];
	char gname[32];
	char devmajor[8];
	char devminor[8];
	char prefix[155];
	char pad[12];
}; // struct PosixTarHeader

class TarHeader {

public :

	typedef PosixTarHeader::LongType LongType;

	enum RecordEnum {RECORD_CLOSED, RECORD_OPEN};

	TarHeader(PsimagLite::String filename, LongType len)
	    : header_(filename,len),len_(len),status_(RECORD_OPEN)
	{
		checksum();
	}

	TarHeader(std::ifstream& fin)
	    : header_(fin),len_(0),status_(RECORD_CLOSED)
	{}

	~TarHeader()
	{
		if (len_ == 0) return;
		statusCheck(RECORD_CLOSED,"dtor");
	}

	void writeTo(std::ofstream& fout) const
	{
		statusCheck(RECORD_OPEN,"writeTo");
		fout.write((const char*)(&header_),
		           sizeof(PosixTarHeader));
	}

	void endRecord(std::ofstream& fout)
	{
		statusCheck(RECORD_OPEN,"endRecord");
		if (len_ == 0) {
			status_ = RECORD_CLOSED;
			return;
		}

		char c='\0';
		LongType len = len_;
		while ((len%sizeof(PosixTarHeader)) != 0)
		{
			fout.write(&c,sizeof(char));
			++len;
		}

		status_ = RECORD_CLOSED;
	}

	const PosixTarHeader& header() const { return header_; }

private:

	void checksum()
	{
		unsigned int sum = 0;
		char *p = (char *)(&header_);
		char *q = p + sizeof(PosixTarHeader);
		while (p < header_.checksum) sum += *p++ & 0xff;
		for (int i = 0; i < 8; ++i)  {
			sum += ' ';
			++p;
		}

		while (p < q) sum += *p++ & 0xff;

		std::sprintf(header_.checksum,"%06o",sum);
	}

	void statusCheck(RecordEnum whatItShouldBe, PsimagLite::String func) const
	{
		if (status_ == whatItShouldBe) return;
		std::cerr<<"TarHeader: WARNING: status is "<<statusToString(status_);
		std::cerr<<" but it should have been "<<statusToString(whatItShouldBe);
		std::cerr<<" from function = "<<func<<"\n";
	}

	PsimagLite::String statusToString(RecordEnum st) const
	{
		if (st == RECORD_CLOSED) return "RECORD_CLOSED";
		return "RECORD_OPEN";
	}

	PosixTarHeader header_;
	LongType len_;
	RecordEnum status_;
}; // class TarHeader

class TarHelper {

	typedef TarHeader::LongType LongType;

public:

	static LongType fileSize(PsimagLite::String filename)
	{
	    std::ifstream fin(filename.c_str(),
		                 std::ifstream::ate | std::ifstream::binary);
	    return fin.tellg();
	}

	static void goodDescriptorOrThrow(std::ifstream& f,
	                                  PsimagLite::String filename)
	{
		if (!f || !f.good() || !f.is_open() || f.bad())
			throw PsimagLite::RuntimeError("Cannot read from " + filename + "\n");
	}

	static void goodDescriptorOrThrow(std::ofstream& f,
	                                  PsimagLite::String filename)
	{
		if (!f || !f.is_open())
			throw PsimagLite::RuntimeError("Cannot write to " + filename + "\n");
	}
};

class TarPack  {

	typedef TarHeader::LongType LongType;

public:

	TarPack(PsimagLite::String filename)
	    : fout_(filename.c_str(),std::ifstream::binary)
	{
		TarHelper::goodDescriptorOrThrow(fout_,filename);
	}

	~TarPack()
	{
	    fout_.close();
	}

	void add(PsimagLite::String filename, PsimagLite::String nameInArchive = "")
	{
		if (nameInArchive == "") nameInArchive = filename;

		std::ifstream fin(filename.c_str(),std::ifstream::binary);
		TarHelper::goodDescriptorOrThrow(fin, filename);

		LongType len = TarHelper::fileSize(filename);


		TarHeader tarHeader(nameInArchive,len);

		tarHeader.writeTo(fout_);

		std::copy(std::istreambuf_iterator<char>(fin),
		          std::istreambuf_iterator<char>(),
		          std::ostreambuf_iterator<char>(fout_));

		tarHeader.endRecord(fout_);
	}


private:

	std::ofstream fout_;

};     //class TarPack

class UnTarPack  {

	typedef TarHeader::LongType LongType;

public:

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	UnTarPack(PsimagLite::String filename)
	    : fin_(filename.c_str(),std::ifstream::binary)
	{
		TarHelper::goodDescriptorOrThrow(fin_,filename);
	}

	~UnTarPack()
	{
	    fin_.close();
	}

	void list(VectorStringType& files)
	{
		while (!fin_.eof()) {
			TarHeader tarHeader(fin_);
			std::stringstream ss;
			ss<<std::oct<<tarHeader.header().size;
			LongType len = 0;
			ss>>len;
			if (fin_.eof()) break;
			LongType currentPos = fin_.tellg();
			fin_.seekg(currentPos + len);
			readPadding(len);
			PsimagLite::String name(tarHeader.header().name);
			files.push_back(name);
			//std::cout<<name<<" "<<len<<"\n";
		}
	}

private:

	void readPadding(LongType len)
	{
		if (len == 0) return;
		SizeType padding = len % sizeof(PosixTarHeader);
		padding = sizeof(PosixTarHeader) - padding;
		if (padding == 0) return;
		SizeType allocation = (padding > 512) ? padding : 512;
		char *c = new char[allocation];
		fin_.read(c,padding);
		delete [] c;
	}

	std::ifstream fin_;

};     //class UnTarPack
} // namespace Dmrg
/*@}*/
#endif

