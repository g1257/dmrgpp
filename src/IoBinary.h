/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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

/*! \file IoBinary.h
 *
 *  This class handles Input/Output in binary

Objectives:

* Support search (with labels)

* Robust, that is, consistency checks.

* Fast

* As small file size as possible.

------------------------ SPECIFICATION	 -------------------------------

The binary file is composed of records. Each record is composed of:

(1) label
(2) checksum
(3) total record size
(4) type (*)
(5) size (*)
(6) payload (*)

Parts marked with (*) are optional.

(1) label
Size of this part: variable
1.1 the 5 characters LABEL ( 5 bytes)
1.2 the size of the actual label (sizeof(size_t), usually 8 bytes)
1.3 the label (variable size)

(2) checksum including (3), (4), (5) and (6).
Size of this part: 1 byte if we use, let's say, CRC-8
Rationale: This is for consistency check.

(3) total record size including (4), (5) and (6).
Size of this part = sizeof(size_t) usually 8 bytes.
Rationale: This is for consistency check.

(4) type
Size of this part: 1 byte
4.1 The 4 lower bits describe the native type, for example,
int, size_t, double, float, char *, ..., we have space here for 16 possibilities
4.2 The 4 higher bits describe composite structures, for example,
vector, matrix, etc. We might reserve the last bit to indicate complex.

(5) The size of the payload
Size of this part: sizeof(size_t), usually 8 bytes

(6) The actual payload.
Size of this part: variable

----------------------------END OF SPECIFICATION -----------------------------------

Example 1: print the label "Hello, World!"
5 bytes: 4c 41 42 45 4c <-- The word LABEL
8 bytes: 00 00 00 00 00 00 00 0E  --> the number 14
12 bytes: 48 65 6c 6c 6f 20 57 6f 72 6c 64 21 <--- The phrase Hello World!
1 byte: ??  <-- checksum
8 bytes  00 00 00 00 00 00 00 00  <-- size of the record
There's no (4), (5), (6).

Total size = 34 bytes
If it were written in ascii = 12 bytes.


Example 2: print the number 12 with label "TotalSites"
5 bytes: 4c 41 42 45 4c <-- The word LABEL
8 bytes: 00 00 00 00 00 00 00 0A  --> the number 10
10 bytes: 54 6f 74 61 6c 53 69 74 65 73 <-- The word TotalSites
1 byte: ??  <-- checksum
8 bytes: 09 <-- size of the record
1 byte: 01 <--- type==integer
8 bytes: 00 00 00 00 00 00 00 0C <-- the number 12 (payload)

Total size = 41 bytes
If it were written in ascii = 12 bytes.


Example 3: print the vector of 6 doubles 0.3 0.55 0.786  -0.59 9.33 -1.0 with label "Energies"
5 bytes: 4c 41 42 45 4c <-- The word LABEL
8 bytes: 00 00 00 00 00 00 00 08  --> the number 8
8 bytes: 45 6e 65 72 67 69 65 73 <-- The word Energies
1 byte: ?? <-- checksum
8 bytes: 09 <-- size of the record
1 byte: 13 <--- type==19, which would mean double and vector
48 bytes:  ... <-- the numbers making up the vector (payload)

Total Size = 74 bytes.
If it were written in ascii = 43 bytes.
----------------------------------------------------------------------
 */
  
#ifndef IO_BINARY_H
#define IO_BINARY_H

#include <iostream>
#include <string>
#include <vector>
#include "Stack.h"
#include <utility>
#include "Matrix.h"
#include "TypeToString.h"
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

namespace PsimagLite {
	//! IoBinary class handles Input/Output (IO) in binary format
	class IoBinary {

		typedef unsigned char TypeType;
		typedef std::pair<size_t,size_t> PairType;

		enum {TYPE_INT=1,TYPE_SIZE_T=2,TYPE_FLOAT=3,TYPE_DOUBLE=4};

		enum {TYPE_VECTOR = 16,TYPE_MATRIX = 32, TYPE_PAIR = 48, TYPE_UNKNOWN = 64};

#ifdef HAVE_BINARY_IO
		static const size_t ENABLED=1;
#else
		static const size_t ENABLED=0;
#endif

	public:

		class Out {

		public:

			Out()  : rank_(0),fout_(-1) {}

			Out(const std::string& fn,int rank) : rank_(rank),filename_(fn),fout_(-1)
			{
				if (!ENABLED) return;
				if (rank_!=0) return;
				fout_ = ::open(fn.c_str(),O_WRONLY | O_CREAT | O_NOCTTY | O_TRUNC, S_IRUSR | S_IWUSR);
				if (fout_<0) {
					std::string s(__FILE__);
					s += " " + ttos(__LINE__) + "\n";
					s += "Cannot open file " + filename_ + " for writing\n";
					throw std::runtime_error(s.c_str());
				}
			}

			~Out()
			{
				if (!ENABLED) return;
				if (rank_!=0) return;
				if (fout_<0) return;
				::close(fout_);
			}

			void open(std::string const &fn,std::ios_base::openmode mode,int rank)
			{
				if (!ENABLED) return;
				rank_=rank;
				if (rank_!=0) return;
				if (filename_=="OSTREAM")
					throw std::runtime_error("open: not possible\n");
				filename_=fn;
				//if (!fout_) fout_=new std::ofstream;
				int flags = 0;
				if (mode==std::ios_base::trunc) flags = O_WRONLY | O_CREAT | O_NOCTTY | O_TRUNC;
				else if (mode==std::ios_base::app) flags = O_WRONLY | O_APPEND;
				fout_ = ::open(fn.c_str(),flags, S_IRUSR | S_IWUSR);
				if (fout_<0)
					throw std::runtime_error("Out: error while opening file!\n");
			}

			void close()
			{
				if (!ENABLED) return;
				if (fout_<0)
					throw std::runtime_error("close: not possible\n");
				filename_="FILE_IS_CLOSED";
				::close(fout_);
			}


//			void printline(std::ostringstream &s)
//			{
//				if (rank_!=0) return;
//				(*fout_)<<s.str()<<"\n";
//				s.flush();
//				s.seekp(std::ios_base::beg);
//			}



			template<typename X>
			void printVector(X const &x,std::string const &label)
			{
				if (!ENABLED) return;
				if (rank_!=0) return;
				printLabel(label);

				char check = 0;
				write(fout_, (const void *)&check,1);
				size_t total = 0;
				write(fout_, (const void *)&total,sizeof(total));

				typename X::value_type dummy;
				TypeType type = TYPE_VECTOR | charTypeOf(dummy);
				write(fout_,(const void *)&type,sizeof(type));
				size_t length = x.size();
				write(fout_,(const void *)&length,sizeof(length));

				const typename X::value_type* const ptr = &(x[0]);
				write(fout_, (const void *)ptr,length*sizeof(dummy));
			}

			//! Specialization for bool (FIXME: std::vector<bool> is evil, remove it!)
			void printVector(const std::vector<bool>& x,std::string const &label)
			{
				if (!ENABLED) return;
				if (rank_!=0) return;
				printLabel(label);

				char check = 0;
				write(fout_, (const void *)&check,1);
				size_t total = 0;
				write(fout_, (const void *)&total,sizeof(total));

				bool dummy;
				TypeType type = TYPE_VECTOR | charTypeOf(dummy);
				write(fout_,(const void *)&type,sizeof(type));
				std::vector<char> xx(x.size());
				convertFromBool(xx,x);
				size_t length = xx.size();
				write(fout_,(const void *)&length,sizeof(length));
				write(fout_,(const void *)&(xx[0]),length);
			}

			template<class T>
			void print(const std::string& s,const T&  something)
			{
				if (!ENABLED) return;
				if (rank_!=0) return;
				makeSureFileIsOpen();
				printLabel(s);
				char check = 0;
				write(fout_, (const void *)&check,1);
				size_t total = 1;
				write(fout_, (const void *)&total,sizeof(total));
				TypeType type = charTypeOf(something);
				write(fout_,(const void *)&type,sizeof(type));
				write(fout_,(const void*)&something,sizeof(something));
			}

			void print(const std::string& s)
			{
				if (!ENABLED) return;
				printLabel(s);
				char check = 0;
				write(fout_, (const void *)&check,1);
				size_t total = 0;
				write(fout_, (const void *)&total,sizeof(total));
			}

			template<typename SomeMatrixType>
			void printMatrix(const SomeMatrixType& mat,std::string const &s)
			{
				if (!ENABLED) return;
				if (rank_!=0) return;
				printLabel(s);

				char check = 0;
				write(fout_, (const void *)&check,1);
				size_t total = 0;
				write(fout_, (const void *)&total,sizeof(total));

				typename SomeMatrixType::value_type dummy;
				TypeType type = TYPE_MATRIX | charTypeOf(dummy);
				write(fout_,(const void *)&type,sizeof(type));

				mat.print(fout_);
			}

			template<typename SomeType>
			void printMatrix(const std::stack<SomeType>& mat,std::string const &s)
			{
				if (!ENABLED) return;
				if (rank_!=0) return;
				printLabel(s);

				char check = 0;
				write(fout_, (const void *)&check,1);
				size_t total = 0;
				write(fout_, (const void *)&total,sizeof(total));

				SomeType dummy;
				TypeType type = TYPE_MATRIX | charTypeOf(dummy);
				write(fout_,(const void *)&type,sizeof(type));
				std::print(fout_,mat);
			}

			int rank() { return rank_; }

//			template<typename X>
//			friend Out& operator<<(Out& io,const X& t);

		private:

			void printLabel(const std::string& s)
			{
				if (rank_!=0) return;
				makeSureFileIsOpen();

				char label1[] = {'L','A','B','E','L'};
				write(fout_,label1,5);
				size_t length = s.length();
				write(fout_,(const void *)&length,sizeof(length));
				write(fout_, (const void *)&(s[0]),length);
			}

			void makeSureFileIsOpen() const
			{
				if (fout_<0) {
					std::string s(__FILE__);
					s += ttos(__LINE__) + "\n";
					s += "File " + filename_ + " is not open\n";
					throw std::runtime_error(s.c_str());
				}
			}

			int rank_;
			std::string filename_;
			int fout_;
		};

		class In {

		public:

			static const int LAST_INSTANCE=-1;
			static const size_t MEMORY_POOL_SIZE = 10240;
			typedef int long long LongIntegerType;
			typedef unsigned int long long LongSizeType;

			In() : filename_(""),fin_(-1),memoryPool_(0)
			{}

			In(std::string const &fn) : filename_(fn), fin_(-1),memoryPool_(0)
			{
				if (!ENABLED) return;
				fin_= ::open(fn.c_str(),O_RDONLY);
				if (fin_<0) {
					std::string s = "IoBinary::ctor(...): Can't open file "
							+filename_+"\n";
					throw std::runtime_error(s.c_str());
				}
			}

			~In()
			{
				if (!ENABLED) return;
				if (memoryPool_) delete [] memoryPool_;
				if (fin_>0) ::close(fin_);
			}

			void open(std::string const &fn)
			{
				if (!ENABLED) return;
				filename_=fn;
				fin_= ::open(fn.c_str(),O_RDONLY);
				if (fin_<0) {
					std::string s = "IoBinaryIn::open(...) failed for file "
							+ filename_ + "\n";
					throw std::runtime_error(s.c_str());
				}
			}

			void close() 
			{
				if (!ENABLED) return;
				filename_="FILE_IS_CLOSED";
				assert(fin_>=0);
				::close(fin_);
			}

			template<typename X>
			size_t readline(X &x,const std::string &s,LongIntegerType level=0)
			{
				if (!ENABLED) return 0;
				std::pair<std::string,size_t> sc = advance(s,level,false);
				//std::cerr<<"FOUND--------->\n";
				char check = 0;
				size_t total = 0;
				PairType type;
				readCheckTotalAndType(check,total,type);

				if (total!=1) {
					std::string str(__FILE__);
					str += ttos(__FILE__) + "\n";
					str += "Expecting a single item\n";
					throw std::runtime_error(str.c_str());
				}

				::read(fin_,&x,sizeof(x));
				return sc.second;
				
			}

			template<typename X>
			std::pair<std::string,size_t> read(X &x,
							   std::string const &s,
							   LongIntegerType level=0,
							   bool beQuiet = false)
			{
				if (!ENABLED) return std::pair<std::string,size_t>("NOT_ENABLED",0);
				std::pair<std::string,size_t> sc = advance(s,level,beQuiet);
				//std::cerr<<"FOUND--------->\n";
				char check = 0;
				size_t total = 0;
				PairType type;
				readCheckTotalAndType(check,total,type);
				if (type.second == TYPE_VECTOR) {
					readVector(x);
					return sc;
				}

				throw std::runtime_error("Unimplemented read\n");
			}

			//! Assumes something of the form 
			//! label[key]=value
//			template<typename X>
//			void read(std::map<std::string,X> &x,
//			          std::string const &s,
//			          LongIntegerType level=0)
//			{
//				size_t counter=0;
//				bool beQuiet = true;
//				while(true) {
//					try {
//						std::pair<std::string,size_t> sc = advance(s,level,beQuiet);
//						// sc.first contains the full string and also value
//						std::string key;
//						X val=0;
//						getKey(key,val,sc.first);
				
//						x[key]=val;
//						counter++;
//					} catch (std::exception& e) {
//						break;
//					}
//				}
//				rewind();
//				if (counter==0) {
//					std::string s2 (__FILE__);
//					s2 += " No " + s + " found in the input file or ";
//					s2 += " could not parse it\n";
//					throw std::runtime_error(s2.c_str());
//				}
//			}
			
//			template<typename X>
//			std::pair<std::string,size_t> readKnownSize(X &x,
//			                                            std::string const &s,
//			                                            LongIntegerType level=0)
//			{
//				std::pair<std::string,size_t> sc = advance(s,level);
				
//				for (size_t i=0;i<x.size();i++) {
//					typename X::value_type tmp;
//					fin_>>tmp;
//					x[i]=tmp;
//				}
//				return sc;
//			}

			size_t count(const std::string& s)
			{
				if (!ENABLED) return 0;
				size_t counter = 0;
				while(true) {
					std::string temp = readNextLabel();
					if (temp=="NOTFOUND") break;
					counter++;
				}
				return counter;
			}

			std::pair<std::string,size_t> advance(std::string const &s,
							      LongIntegerType level=0,
							      bool beQuiet=false)
			{
				if (!ENABLED) return std::pair<std::string,size_t>("NOT_ENALBED",0);
				//std::string temp="NOTFOUND";
				std::string tempSaved="NOTFOUND";
				LongSizeType counter=0;
				bool found=false;
				//size_t c = 0;
				while(true) {

					std::string temp = readNextLabel();
					if (temp=="NOTFOUND") break;
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
					::close(fin_);
					fin_ = ::open(filename_.c_str(),O_RDONLY);
					if (counter>1) advance(s,counter-2);
					return std::pair<std::string,size_t>(tempSaved,counter);
				}
				
				//std::cerr<<"count="<<c<<"\n";
				if (!found && tempSaved=="NOTFOUND") {
					if (!beQuiet) {
						std::cerr<<"Not found "<<s<<" in file "<<filename_;
						std::cerr<<" level="<<level<<" counter="<<counter<<"\n";
					}
					throw std::runtime_error("IoBinary::In::read()\n");
				}
				//std::cerr<<"------------\n";
				return std::pair<std::string,size_t>(tempSaved,counter);
			}
			
//			size_t count(const std::string& s)
//			{
//				size_t i = 0;
//				while(i<1000) {
//					try {
//						advance(s,0,true);
//						i++;
//					} catch (std::exception& e) {
//						rewind();
//						return i;
//					}
//				}
//				std::string ss = "IoBinary::count(...): too many "
//					+s+" in file "+filename_+"\n";
//				throw std::runtime_error(s.c_str());
				
//			}

//			template<typename T>
//			void read(std::vector<std::pair<T,T> > &x,
//			          std::string const &s,
//			          LongIntegerType level=0)
//			{
//				advance(s,level);
//				int xsize;
//				fin_>>xsize;
//				x.resize(xsize);
//				T tmp1,tmp2;
//				for (int i=0;i<xsize;i++) {
//					fin_>>tmp1;
//					fin_>>tmp2;
//					x[i]=std::pair<T,T>(tmp1,tmp2);
//				}
//			}

//			template<typename X,template<typename> class SomeType>
//			void readSparseVector(SomeType<X> &x,
//			                      std::string const &s,
//			                      LongIntegerType level=0)
//			{
//				advance(s,level);
//				int xsize;
//				fin_>>xsize;
//				x.resize(xsize);
//				fin_>>xsize;
//				for (int i=0;i<xsize;i++) {
//					int index;
//					X value;
//					fin_>>index;
//					fin_>>value;
//					x[index]=value;
//				}
				
//			}

			template<typename X>
			void readMatrix(X &mat,
					std::string const &s,
					LongIntegerType level= 0)
			{
				if (!ENABLED) return;
				advance(s,level);
				std::cerr<<"FOUND--------->\n";
				char check = 0;
				size_t total = 0;
				PairType type;
				readCheckTotalAndType(check,total,type);
				if (type.second == TYPE_MATRIX) {
					readMatrix(mat,s);
					return;
				}

				throw std::runtime_error("Unimplemented read\n");
			}

//			template<
//				typename FieldType,
//				template <typename> class SparseMatrixTemplate,
//				template<typename,template<typename> class>
//			class X>
//			void readMatrix(X<FieldType,SparseMatrixTemplate>& op,
//			                const std::string& s,
//			                LongIntegerType level=0)
//			{
//				advance(s,level);
//				fin_>>op.data;
//				fin_>>op.fermionSign;
//				fin_>>op.j;
//			}

			void rewind()
			{
				if (!ENABLED) return;
				throw std::runtime_error("IoBinary::rewind(): Unimplemented\n");
//				fin_.clear(); // forget we hit the end of file
//				fin_.seekg(0, std::ios::beg); // move to the start of the file
			}

			const char* filename() const 
			{
				return filename_.c_str();
			}

//			template<typename X>
//			friend void operator>>(In& io,X& t);

			std::string readNextLabel()
			{
				if (!ENABLED) return "NOT_ENABLED";
				std::string magicLabel = "LABEL";
				char c = 0;
				size_t j= 0;
				while (true) {
					int x = ::read(fin_,&c,1);
					if (x!=1) {
						std::string s(__FILE__);
						s += " " + ttos(__LINE__);
						s += "\nCannot read next label\n";
						return "NOTFOUND";
						//throw std::runtime_error(s.c_str());
					}
					if (c==magicLabel[j]) {
						j++;
					} else {
						j=0;
					}
					//std::cerr<<"j="<<j<<" "<<magicLabel.length()<<"\n";
					if (j==magicLabel.length()) break;
				}

				size_t y = 0;
				int x = ::read(fin_,&y,sizeof(y));
				if (x!=sizeof(y)) {
					std::string s(__FILE__);
					s += " " + ttos(__LINE__);
					s += "\nCannot read next label (2)\n";
					throw std::runtime_error(s.c_str());
				}

				if (size_t(y)>=MEMORY_POOL_SIZE) {
					std::string s(__FILE__);
					s += " "  + ttos(__LINE__) + "\nMemory Pool " + ttos(MEMORY_POOL_SIZE);
					s +=  " is too small, needed " + ttos(y) + "\n";
				}
				if (!memoryPool_) memoryPool_ = new TypeType[MEMORY_POOL_SIZE];
				x = ::read(fin_,memoryPool_,y);
				if (size_t(x)!=y) {
					std::string s(__FILE__);
					s += " " + ttos(__LINE__);
					s += "\nCannot read next label (3)\n";
					throw std::runtime_error(s.c_str());
				}
				std::string ret="";
				for (size_t i=0;i<size_t(x);i++) ret += memoryPool_[i];
				std::cerr<<"FOUND!!! ret="<<ret<<"\n";
				return ret;

			}

			void readCheckTotalAndType(char& check,size_t& total,PairType& type)
			{
				if (!ENABLED) return;
				check = 0;
				::read(fin_,&check,1);
				total = 0;
				::read(fin_,&total,sizeof(total));
				type.first = 0;
				::read(fin_,&type.first,1);
				type.second = type.first;
				type.second >>= 4;
				type.first &= 15;
			}

			static std::string nameOfType(const PairType& type)
			{
				if (!ENABLED) return "NOT_ENABLED";
				std::string type1 = "UNKNOWN";
				if (type.first==TYPE_INT) type1 = "INT";
				else if (type.first==TYPE_SIZE_T) type1="SIZE_T";
				else if (type.first==TYPE_FLOAT) type1="FLOAT";
				else if (type.first==TYPE_DOUBLE) type1="DOUBLE";

				std::string type2 = "";
				if (type.second==TYPE_VECTOR) type2="VECTOR";
				if (type.second==TYPE_MATRIX) type2="MATRIX";
				if (type.second==TYPE_PAIR) type2="PAIR";
				if (type.second==TYPE_UNKNOWN) type2="UNKNOWN";

				return type1 + " " + type2;
			}

		private:

			template<typename X>
			void readVector(std::vector<X>& x)
			{
				size_t xsize = 0;
				::read(fin_,&xsize,sizeof(xsize));

				std::cerr<<"xsize="<<xsize<<"\n";
				x.resize(xsize);

				X tmp;
				for (size_t i=0;i<xsize;i++) {
					int l = ::read(fin_,&tmp,sizeof(tmp));
					if (l!=sizeof(tmp)) throw std::runtime_error("Mmm!\n");
					x[i]=tmp;
				}
			}

			template<typename X>
			void readMatrix(Matrix<X> &mat)
			{
				mat.read(fin_);
			}

			//! full contains label[key]=value
//			template<typename X>
//			void getKey(std::string& key,X& x,const std::string& full)
//			{
//				size_t i=0;
//				for (;i<full.length();i++) {
//					if (full[i]=='[') break;
//				}
//				key = "";
//				size_t j=i+1;
//				for (;j<full.length();j++) {
//					if (full[j]==']') break;
//					key += full[j];
//				}
//				j++;
//				if (full[j++]!='=') {
//					std::string s(__FILE__);
//					s += "Something failed while parsing line " + full;
//					s += " of input file\n";
//					throw std::runtime_error(s.c_str());
//				}
//				std::string val="";
//				for (size_t k=j;k<full.length();k++)
//					val += full[k];
//				x = atof(val.c_str());
//			}

			std::string filename_;
			int 	fin_;
			TypeType* memoryPool_;
		};

		template<typename T>
		static size_t charTypeOf(const T& dummy)
		{
			std::string s(__FILE__);
			s += " " + ttos(__LINE__) + "\nNo known type\n";
			std::cerr<<"ATTENTION!!!! "<<s;
			return TYPE_UNKNOWN;
		}

		static size_t charTypeOf(const double& dummy) { return TYPE_DOUBLE; }

		static size_t charTypeOf(const float& dummy) { return TYPE_FLOAT; }

		static size_t charTypeOf(const size_t& dummy) { return TYPE_SIZE_T; }

		template<typename T>
		static size_t charTypeOf(const std::pair<T,T>& dummy)
		{
			T dummy2;
			return TYPE_PAIR | charTypeOf(dummy2);
		}

		static void convertFromBool(std::vector<char>& xx,const std::vector<bool>& x)
		{
			for (size_t i=0;i<xx.size();i++) xx[i] = (x[i]) ? 48 : 49;
		}

	}; //class IoBinary

//	template<typename T>
//	IoBinary::Out& operator<<(IoBinary::Out& io,const T& t)
//	{
//		if (io.rank_!=0) return io;
//		io.print(t);
//		return io;
//	}

//	template<typename T>
//	void operator>>(IoBinary::In& io,T& t)
//	{
//		io.read(t);
//	}

//	template<typename T1,typename T2>
//	void printMap(std::ostream& os, const std::map<T1,T2>& x,const std::string& label)
//	{
//		typedef typename std::map<T1,T2>::const_iterator MapIteratorType;
//		for (MapIteratorType it = x.begin();it!=x.end();++it) {
//			os<<label<<"["<<it->first<<"]="<<it->second<<"\n";
//		}
//	}

} // namespace PsimagLite 


//namespace Spf {

//	class IoBinaryIn : public PsimagLite::IoBinary::In {
//	public:
//		IoBinaryIn(const char* fn) : PsimagLite::IoBinary::In(std::string(fn)) { }
//	};
//}
/*@}*/	
#endif
