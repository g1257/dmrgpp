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
#include "String.h"
#include "Concurrency.h"

namespace PsimagLite {
	//! IoSimple class handles Input/Output (IO) for the Dmrg++ program 
	class IoSimple {

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

			void open(String const &fn,
					std::ios_base::openmode mode)
			{
				if (rank_!=0) return;
				if (filename_=="OSTREAM") 
					throw RuntimeError("open: not possible\n");
				filename_=fn;
				if (!fout_) fout_=new std::ofstream;
				fout_->open(fn.c_str(),mode);
				if (!(*fout_) || !fout_->good())
					throw RuntimeError("Out: error while opening file!\n");
			}

			void close()
			{
				if (filename_=="OSTREAM")
					throw RuntimeError("close: not possible\n");
				filename_="FILE_IS_CLOSED";
				fout_->close();
			}

			void printline(const String &s)
			{
				if (rank_!=0) return;
				(*fout_)<<s<<"\n";
			}

			void printline(PsimagLite::OstringStream &s)
			{
				if (rank_!=0) return;
				(*fout_)<<s.str()<<"\n";
				s.flush();
				s.seekp(std::ios_base::beg);
			}

			template<typename X>
			void printVector(X const &x,String const &label)
			{
				if (rank_!=0) return;
				(*fout_)<<label<<"\n";
				(*fout_)<<x.size()<<"\n";
				for (SizeType i=0;i<x.size();i++) (*fout_)<<x[i]<<"\n";
			}

			template<class T>
			void print(const String& label,const T&  something)
			{
				if (rank_!=0) return;
				if (!(*fout_) || !fout_->good())
					throw RuntimeError("Out: file not open!\n");
				(*fout_)<<label;
				(*fout_)<<something<<"\n";
			}

			void print(const String& something)
			{
				if (rank_!=0) return;
				(*fout_)<<something;
				SizeType last = something.length();
				if (last>0) last--;
				if (something[last]!='\n') (*fout_)<<"\n";
			}

			template<typename SomePrintableType>
			void print(const SomePrintableType& something)
			{
				if (rank_!=0) return;
				(*fout_)<<something;
			}

			template<typename X>
			void printMatrix(Matrix<X> const &mat,String const &s)
			{
				if (rank_!=0) return;
				(*fout_)<<s<<"\n";
				(*fout_)<<mat;
			}

			template<typename X>
			void printMatrix(X const &mat,String const &s)
			{
				if (rank_!=0) return;
				(*fout_)<<s<<"\n";
				(*fout_)<<mat;
			}

			int rank() { return rank_; }

			void flush()
			{
				if (rank_!=0 || !fout_) return;
				fout_->flush();
			}

			void setPrecision(SizeType x)
			{
				if (!fout_) return;
				fout_->precision(x);
			}

			template<typename X>
			friend Out& operator<<(Out& io,const X& t);

		private:

			int rank_;
			String filename_;
			std::ofstream* fout_;
		};

		class In {

		public:

			static const int LAST_INSTANCE=-1;
			typedef int long long LongIntegerType;
			typedef unsigned int long long LongSizeType;

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
				filename_="FILE_IS_CLOSED"; 
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
				while(!fin_.eof()) {
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
				x.resize(xsize);
				for (int i=0;i<xsize;i++) {
					typename X::value_type tmp;
					fin_>>tmp;
					x[i]=tmp;
				}
				return sc;
			}

//			template<typename VectorLikeType>
//			void read(VectorLikeType &x,
//			          String const &s,
//			          LongIntegerType level=0)
//			{
//				advance(s,level);
//				int xsize;
//				fin_>>xsize;
//				x.resize(xsize);
//				typename VectorLikeType::value_type::first_type tmp1;
//				typename VectorLikeType::value_type::second_type tmp2;
//				for (int i=0;i<xsize;i++) {
//					fin_>>tmp1;
//					fin_>>tmp2;
//					x[i]=typename VectorLikeType::value_type(tmp1,tmp2);
//				}
//			}

			//! Assumes something of the form 
			//! label[key]=value
			template<typename MapType>
			typename EnableIf<IsMapLike<MapType>::True,void>::Type
			read(MapType& x,
			     String const &s,
			     LongIntegerType level=0)
			{
				SizeType counter=0;
				bool beQuiet = true;
				while(true) {
					try {
						std::pair<String,SizeType> sc = advance(s,level,beQuiet);
						// sc.first contains the full string and also value
						String key;
						typename MapType::mapped_type val=0;
						getKey(key,val,sc.first);
				
						x[key]=val;
						counter++;
					} catch (std::exception& e) {
						break;
					}
				}
				rewind();
				if (counter==0) {
					String s2 (__FILE__);
					s2 += " No " + s + " found in the input file or ";
					s2 += " could not parse it\n";
					throw RuntimeError(s2.c_str());
				}
			}
			
			template<typename X>
			std::pair<String,SizeType> readKnownSize(X &x,
			                                            String const &s,
			                                            LongIntegerType level=0)
			{
				std::pair<String,SizeType> sc = advance(s,level);
				
				for (SizeType i=0;i<x.size();i++) {
					typename X::value_type tmp;
					fin_>>tmp;
					x[i]=tmp;
				}
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
				//SizeType c = 0;
				while(!fin_.eof()) {
					fin_>>temp;
					//c++;
					//std::cerr<<"Line="<<temp<<" target="<<s<<"\n";
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
				
				//std::cerr<<"count="<<c<<"\n";
				if (!found && tempSaved=="NOTFOUND") {
					if (!beQuiet) {
						std::cerr<<"Not found "<<s<<" in file "<<filename_;
						std::cerr<<" level="<<level<<" counter="<<counter<<"\n";
					}
					throw RuntimeError("IoSimple::In::read()\n");
				}
				//std::cerr<<"------------\n";
				return std::pair<String,SizeType>(tempSaved,counter);
			}
			
			SizeType count(const String& s)
			{
				SizeType i = 0;
				while(i<1000) {
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

			template<typename X,template<typename> class SomeType>
			void readSparseVector(SomeType<X> &x,
			                      String const &s,
			                      LongIntegerType level=0)
			{
				advance(s,level);
				int xsize;
				fin_>>xsize;
				x.resize(xsize);
				fin_>>xsize;
				for (int i=0;i<xsize;i++) {
					int index;
					X value;
					fin_>>index;
					fin_>>value;
					x[index]=value;
				}
				
			}

			template<typename X>
			void readMatrix(X &mat,
			                String const &s,
			                LongIntegerType level= 0)
			{
				advance(s,level);
				fin_>>mat;
			}

			template<
				typename FieldType,
				template <typename> class SparseMatrixTemplate,
				template<typename,template<typename> class>
			class X>
			void readMatrix(X<FieldType,SparseMatrixTemplate>& op,
			                const String& s,
			                LongIntegerType level=0)
			{
				advance(s,level);
				fin_>>op.data;
				fin_>>op.fermionSign;
				fin_>>op.j;
			}

			void rewind()
			{
				fin_.clear(); // forget we hit the end of file
				fin_.seekg(0, std::ios::beg); // move to the start of the file
			}

			const char* filename() const 
			{
				return filename_.c_str();
			}

			template<typename X>
			friend void operator>>(In& io,X& t);

		private:

			//! full contains label[key]=value
			template<typename X>
			void getKey(String& key,X& x,const String& full)
			{
				SizeType i=0;
				for (;i<full.length();i++) {
					if (full[i]=='[') break;
				}
				key = "";
				SizeType j=i+1;
				for (;j<full.length();j++) {
					if (full[j]==']') break;
					key += full[j];
				}
				j++;
				if (full[j++]!='=') {
					String s(__FILE__);
					s += "Something failed while parsing line " + full;
					s += " of input file\n";
					throw RuntimeError(s.c_str());
				}
				String val="";
				for (SizeType k=j;k<full.length();k++)
					val += full[k];
				x = atof(val.c_str());
			}

			String filename_;
			std::ifstream fin_;
		};
	}; //class IoSimple

	template<typename T>
	IoSimple::Out& operator<<(IoSimple::Out& io,const T& t)
	{
		if (io.rank_!=0) return io;
		(*(io.fout_))<<t;
		return io;
	}

	template<typename T>
	void operator>>(IoSimple::In& io,T& t)
	{
		io.fin_>>t;
	}

} // namespace PsimagLite 


namespace Spf {

	class IoSimpleIn : public PsimagLite::IoSimple::In {
	public:
		IoSimpleIn(const char* fn) : PsimagLite::IoSimple::In(PsimagLite::String(fn)) { }
	};
}
/*@}*/	
#endif
