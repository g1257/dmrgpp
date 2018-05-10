/*
Copyright (c) 2009-2018, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.]
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

/** \ingroup PsimagLite */
/*@{*/

/*! \file IoSimple.h
 *
 *  This class handles Input/Output in a simple way
 */

#ifndef IOSIMPLE_HEADER_H
#define IOSIMPLE_HEADER_H

#include <iostream>
#include <vector>
#include <fstream>
#include "Matrix.h"
#include <cstdlib>
#include "Map.h"
#include "Concurrency.h"
#include "Stack.h"

namespace PsimagLite {
//! IoSimple class handles Input/Output (IO) for the Dmrg++ program
class IoSimple {

	template<typename T>
	struct PrintWithEqualSign {
		enum { True = Loki::TypeTraits<T>::isArith ||IsEnum<T>::True };
	};

public:

	class Out {

	public:

		Out()  : rank_(0),fout_(0) {}

		Out(std::ostream& os) : rank_(0), filename_("OSTREAM")
		{
			fout_=(std::ofstream *)&os;
		}

		Out(const String& fn)
		    : rank_(Concurrency::rank()),filename_(fn),fout_(0)
		{
			if (rank_!=0) return;
			if (!fout_) fout_=new std::ofstream;
#ifdef PSI_PUBSETBUF
			fout_->rdbuf()->pubsetbuf(0,0);
#endif
			fout_->open(fn.c_str());
			if (!(*fout_) || !fout_->good())
				throw RuntimeError("Out: error while opening file!\n");
		}

		~Out()
		{
			if (rank_!=0) return;
			if (filename_=="OSTREAM" || !fout_) return;
			fout_->close();
			delete fout_;
		}

		bool ng() const { return false; }

		const String& filename() const { return filename_; }

		void open(String const &fn,
		          std::ios_base::openmode mode)
		{
			if (rank_!=0) return;
			if (filename_=="OSTREAM")
				throw RuntimeError("open: not possible\n");
			filename_=fn;
			if (!fout_) fout_=new std::ofstream;
#ifdef PSI_PUBSETBUF
			fout_->rdbuf()->pubsetbuf(0,0);
#endif
			fout_->open(fn.c_str(),mode);
			if (!(*fout_) || !fout_->good())
				throw RuntimeError("Out: error while opening file!\n");
		}

		void close()
		{
			if (filename_=="OSTREAM")
				throw RuntimeError("close: not possible\n");
			fout_->close();
		}

		template<typename T>
		void writeline(T,
		               PsimagLite::String,
		               PsimagLite::OstringStream& msg,
		               SizeType)
		{
			printline(msg);
		}

		void printline(const String &s)
		{
			if (rank_!=0) return;
			(*fout_)<<s<<"\n";
		}

		void printline(OstringStream &s)
		{
			if (rank_!=0) return;
			(*fout_)<<s.str()<<"\n";
			s.flush();
			s.seekp(std::ios_base::beg);
		}

		template<typename T>
		void write(const T&x,
		           const String& label,
		           typename EnableIf<PrintWithEqualSign<T>::True, int>::Type = 0)
		{
			(*fout_)<<label<<"="<<x<<"\n";
		}

		template<class T>
		void write(const T& something,
		           const String& label,
		           typename EnableIf<!PrintWithEqualSign<T>::True, int>::Type = 0)
		{
			if (rank_!=0) return;
			if (!(*fout_) || !fout_->good())
				throw RuntimeError("Out: file not open!\n");
			(*fout_)<<label<<"\n";
			(*fout_)<<something<<"\n";
		}

		int rank() { return rank_; }

		void setPrecision(SizeType x)
		{
			if (!fout_) return;
			fout_->precision(x);
		}

		template<typename X>
		friend Out& operator<<(Out& io, const X& t)
		{
			if (io.rank_!=0) return io;
			(*(io.fout_))<<t;
			return io;
		}

	private:

		int rank_;
		String filename_;
		std::ofstream* fout_;
	};

	class In {

	public:

		typedef int long LongIntegerType;
		static const LongIntegerType LAST_INSTANCE = -1;
		static const LongIntegerType ONLY_INSTANCE = 0;
		typedef unsigned int long LongSizeType;

		In() { }

		In(String const &fn) : filename_(fn), fin_(fn.c_str())
		{
			if (!fin_ || !fin_.good() || fin_.bad()) {
				String s = "IoSimple::ctor(...): Can't open file "
				        +filename_+"\n";
				throw RuntimeError(s.c_str());
			}
		}

		~In()
		{
			fin_.close();
		}

		bool ng() const { return false; }

		void open(String const &fn)
		{
			filename_=fn;
			fin_.open(fn.c_str());
			if (!fin_ || !fin_.good() || fin_.bad()) {
				String s = "IoSimpleIn::open(...) failed for file "
				        + filename_ + "\n";
				throw RuntimeError(s.c_str());
			}
		}

		void close()
		{
			fin_.close();
		}

		template<typename X>
		SizeType readline(X &x,const String &s,LongIntegerType level=0)
		{
			String temp;
			bool found=false;
			bool foundOnce =false;
			LongSizeType counter=0;
			if (fin_.bad() || !fin_.good()) throw RuntimeError("Readline\n");
			while (!fin_.eof()) {
				fin_>>temp;
				if (fin_.eof()) break;
				if (temp.substr(0,s.size())==s) {
					foundOnce =true;
					IstringStream temp2(temp.substr(s.size(),
					                                temp.size()));
					temp2 >> x;
					if (level>=0 && counter==LongSizeType(level)) {
						found=true;
						break;
					}
					counter++;
				}
			}

			if (!foundOnce || (!found && level!=LAST_INSTANCE)) {
				String emessage =
				        "IoSimple::In::readline(): Not found "+s+
				        " in file "+filename_;
				throw RuntimeError(emessage.c_str());
			}

			if (level==LAST_INSTANCE) {
				fin_.close();
				fin_.open(filename_.c_str());
				readline(x,s,counter-1);
			}

			return counter;

		}

		template<typename X>
		typename EnableIf<IsVectorLike<X>::True,std::pair<String,SizeType> >::Type
		read(X &x,
		     String const &s,
		     LongIntegerType level=0,
		     bool beQuiet = false)
		{
			std::pair<String,SizeType> sc = advance(s,level,beQuiet);
			int xsize;
			fin_>>xsize;
			if (xsize == 0) return sc;
			x.resize(xsize);
			for (int i=0;i<xsize;i++) {
				typename X::value_type tmp;
				fin_>>tmp;
				x[i]=tmp;
			}
			return sc;
		}

		template<typename X>
		void read(X &mat,
		          String const &s,
		          LongIntegerType level= 0,
		          typename EnableIf<IsMatrixLike<X>::True, int>::Type  = 0)
		{
			advance(s,level);
			fin_>>mat;
		}

		template<typename X>
		typename EnableIf<IsStackLike<X>::True,std::pair<String,SizeType> >::Type
		read(X &x,
		     String const &s,
		     LongIntegerType level=0,
		     bool beQuiet = false)
		{
			std::pair<String,SizeType> sc = advance(s,level,beQuiet);
			fin_ >> x;
			return sc;
		}

		std::pair<String,SizeType> advance(String const &s,
		                                   LongIntegerType level=0,
		                                   bool beQuiet=false)
		{

			String temp="NOTFOUND";
			String tempSaved="NOTFOUND";
			LongSizeType counter=0;
			bool found=false;

			while (!fin_.eof()) {
				fin_>>temp;
				if (fin_.eof() || !fin_.good() || fin_.bad()) break;

				if (temp.substr(0,s.size())==s) {
					tempSaved = temp;
					if (level>=0 && counter==LongSizeType(level)) {
						found=true;
						break;
					}
					counter++;
				}
			}

			if (level==LAST_INSTANCE && tempSaved!="NOTFOUND") {
				fin_.close();
				fin_.open(filename_.c_str());
				if (counter>1) advance(s,counter-2);
				return std::pair<String,SizeType>(tempSaved,counter);
			}

			if (!found && tempSaved=="NOTFOUND") {
				if (!beQuiet) {
					std::cerr<<"Not found "<<s<<" in file "<<filename_;
					std::cerr<<" level="<<level<<" counter="<<counter<<"\n";
				}
				throw RuntimeError("IoSimple::In::read()\n");
			}

			return std::pair<String,SizeType>(tempSaved,counter);
		}

		SizeType count(const String& s)
		{
			SizeType i = 0;
			while (i<1000) {
				try {
					advance(s,0,true);
					i++;
				} catch (std::exception& e) {
					rewind();
					return i;
				}
			}

			String ss = "IoSimple::count(...): too many "
			        +s+" in file "+filename_+"\n";
			throw RuntimeError(s.c_str());

		}

		void rewind()
		{
			fin_.clear(); // forget we hit the end of file
			fin_.seekg(0, std::ios::beg); // move to the start of the file
		}

		bool eof() const { return fin_.eof(); }

		template<typename X>
		friend void operator>>(In& io, X& t)
		{
			(io.fin_)>>t;
		}

	private:

		String filename_;
		std::ifstream fin_;
	};
}; //class IoSimple

template<>
struct IsInputLike<IoSimple::In> {
	enum {True = true};
};

template<>
struct IsOutputLike<IoSimple::Out> {
	enum {True = true};
};

} // namespace PsimagLite

namespace Spf {

class IoSimpleIn : public PsimagLite::IoSimple::In {

public:

	IoSimpleIn(const char* fn) : PsimagLite::IoSimple::In(PsimagLite::String(fn))
	{}
};
}
/*@}*/
#endif

