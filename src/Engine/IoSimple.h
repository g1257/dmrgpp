// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file IoSimple.h
 *
 *  This class handles Input/Output for the Dmrg++ program in a simple way
 */
  
#ifndef IOSIMPLE_HEADER_H
#define IOSIMPLE_HEADER_H

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include "Matrix.h"

namespace Dmrg {
	//! IoSimple class handles Input/Output (IO) for the Dmrg++ program 
	class IoSimple {
	public:
		class Out {
			public:
				Out()  : rank_(0) {}
				Out(const std::string& fn,int rank) : rank_(rank),filename_(fn)
				{
					if (rank_!=0) return;
					fout_.open(fn.c_str());
					if (!fout_ || !fout_.good()) throw std::runtime_error("Out: error while opening file!\n");
				}

				void open(std::string const &fn,std::ios_base::openmode mode,int rank)
				{
					rank_=rank;
					if (rank_!=0) return;
					filename_=fn;
					fout_.open(fn.c_str(),mode);
					if (!fout_ || !fout_.good()) throw std::runtime_error("Out: error while opening file!\n");
				}

				void close() 
				{
					filename_="FILE_IS_CLOSED"; 
					fout_.close(); 
				}

				void printline(const std::string &s) 
				{
					if (rank_!=0) return;
					fout_<<s<<"\n";
				}

				void printline(std::ostringstream &s) 
				{
					if (rank_!=0) return;
					fout_<<s.str()<<"\n";
					s.flush();
					s.seekp(std::ios_base::beg);
				}

				

				template<typename X>
				void printVector(X const &x,std::string const &label)
				{
					if (rank_!=0) return;
					fout_<<label<<"\n";
					fout_<<x.size()<<"\n";
					for (size_t i=0;i<x.size();i++) fout_<<x[i]<<"\n";
				}

// 				template<typename X>
// 				void printSparseVector(X const &x,std::string const &label)
// 				{
// 					if (rank_!=0) return;
// 					x.print(fout_,label);
// 					
// 				}
				
				template<typename X>
				void printStack(const std::stack<X>& st,std::string const &label)
				{
					if (rank_!=0) return;
					std::stack<X> st2=st;
					
					fout_<<label<<"\n";
					fout_<<st.size()<<"\n";
					while(st2.size()>0) {
						X x = st2.top();
						st2.pop();
						fout_<<x<<"\n";
					}
				}

				template<class T>
				void print(const T&  something)
				{
					if (rank_!=0) return;
					if (!fout_ || !fout_.good()) throw std::runtime_error("Out: file not open!\n");
					fout_<<something;
				}

				void print(const std::string& something)
				{
					if (rank_!=0) return;
					fout_<<something;
				}

				template<typename X>
				void printMatrix(psimag::Matrix<X> const &mat,std::string const &s)
				{
					if (rank_!=0) return;
					fout_<<s<<"\n";
					fout_<<mat.n_row()<<" "<<mat.n_col()<<"\n";
					
					for (size_t i=0;i<mat.n_row();i++) {
						for (size_t j=0;j<mat.n_col();j++) fout_<<mat(i,j)<<" ";
						fout_<<"\n";
					}
				}

				template<typename X>
				void printMatrix(X const &mat,std::string const &s)
				{
					if (rank_!=0) return;
					fout_<<s<<"\n";
					fout_<<mat;
					
				}

				int rank() { return rank_; }

				~Out() 
				{
					if (rank_!=0) return;
					fout_.close();
				}

			private:
				int rank_;
				std::string filename_;
				std::ofstream 	fout_;
		};

		class In {
			public:
				static const size_t LAST_INSTANCE=4096;

				In() { }

				In(std::string const &fn) : filename_(fn), fin_(fn.c_str())
				{
					if (!fin_ || !fin_.good() || fin_.bad()) 
						throw std::runtime_error("IoSimple::ctor(...): Can't open file\n");
				}

				void open(std::string const &fn)
				{
					filename_=fn;
					fin_.open(fn.c_str());
					if (!fin_ || !fin_.good() || fin_.bad()) {
						throw std::runtime_error("IoSimpleIn::open(...) failed\n");
					}
				}

				void close() 
				{
					filename_="FILE_IS_CLOSED"; 
					fin_.close(); 
				}

				template<typename X>
				size_t readline(X &x,const std::string &s,size_t level=0) 
				{
					std::string temp;
					bool found=false;
					bool foundOnce =false;
					size_t counter=0;
					while(!fin_.eof()) {
						fin_>>temp;
						if (fin_.eof()) break;
						if (temp.substr(0,s.size())==s) {
							foundOnce =true;
							std::istringstream temp2(temp.substr(s.size(),temp.size()));
							temp2 >> x;
							if (counter==level) {
								found=true;
								break;
							}
							counter++;
						}
					}
					if (!foundOnce || (!found && level!=LAST_INSTANCE)) {
						std::string emessage = "IoSimple::In::readline(): Not found "+s+" in file "+filename_;
						throw std::runtime_error(s);
					}

					if (level==LAST_INSTANCE) {
						fin_.close();
						fin_.open(filename_.c_str());
						readline(x,s,counter-1);
					}	
					return counter;
					
				}

				template<typename X>
				std::pair<std::string,size_t> read(X &x,std::string const &s,int level=0)
				{
					std::pair<std::string,size_t> sc = advance(s,level);
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

				template<typename X>
				void read(std::stack<X>& x,const std::string& s,size_t level=0)
				{
					std::vector<X> tmpVec;
					std::pair<std::string,size_t> sc = read(tmpVec,s,level);
					if (sc.second!=level) {
						fin_.close();
						fin_.open(filename_.c_str());
						read(tmpVec,s,sc.second-1);
					}
					x.empty();
					for (int i=tmpVec.size()-1;i>=0;i--)
					{
						x.push(tmpVec[i]);
					}
				}

				std::pair<std::string,size_t> advance(std::string const &s,int level=0)
				{
					std::string temp="NOTFOUND";
					std::string tempSaved="NOTFOUND";
					int counter=0;
					bool found=false;
					//size_t c = 0;
					while(!fin_.eof()) {
						fin_>>temp;
						//c++;
						//if (temp[0]=='#') std::cerr<<"Line="<<temp<<" target="<<s<<" count="<<c<<"\n";
						if (fin_.eof() || !fin_.good() || fin_.bad()) break;
						
						if (temp.substr(0,s.size())==s) {
							tempSaved = temp;
							if (counter==level) {
								found=true;
								break;
							}
							counter++;
						}
					}
					//std::cerr<<"count="<<c<<"\n";
					if (!found && tempSaved=="NOTFOUND") {
						std::cerr<<"Not found "<<s<<" in file "<<filename_;
						std::cerr<<" level="<<level<<" counter="<<counter<<"\n";
						throw std::runtime_error("IoSimple::In::read()\n");
					}
					//std::cerr<<"------------\n";
					return std::pair<std::string,size_t>(tempSaved,counter);
				}

				template<typename T>
				void read(std::vector<std::pair<T,T> > &x,std::string const &s,int level=0)
				{
					advance(s,level);
					int xsize;
					fin_>>xsize;
					x.resize(xsize);
					T tmp1,tmp2;
					for (int i=0;i<xsize;i++) {
						fin_>>tmp1;
						fin_>>tmp2;
						x[i]=std::pair<T,T>(tmp1,tmp2);
					}
				}

				template<typename X,template<typename> class SomeType>
				void readSparseVector(SomeType<X> &x,std::string const &s,int level=0)
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
				void readMatrix(psimag::Matrix<X> &mat,std::string const &s,int level= 0)
				{
					advance(s,level);
					int nrow,ncol;

					fin_>>nrow;
					fin_>>ncol;
					mat.resize(nrow,ncol);

					for (size_t i=0;i<mat.n_row();i++) for (size_t j=0;j<mat.n_col();j++) fin_>>mat(i,j);
				}

				template<typename FieldType,template <typename> class SparseMatrixTemplate,
    				template<typename,template<typename> class> class X>
				void readMatrix(X<FieldType,SparseMatrixTemplate>& op,const std::string& s,int level=0)
				{
					advance(s,level);
					fin_>>op.data;
					fin_>>op.fermionSign;
					fin_>>op.j;
				}

				template<typename X>
				void readMatrix(X& mat,const std::string& s,int level=0)
				{
					advance(s,level);
					fin_>>mat;
				}

				void rewind()
				{
					fin_.clear();              // forget we hit the end of file
					fin_.seekg(0, std::ios::beg);   // move to the start of the file
				}

				~In() { fin_.close(); }

			private:
				std::string filename_;
				std::ifstream 	fin_;
		};
	}; //class IoSimple
} // namespace Dmrg

/*@}*/	
#endif
