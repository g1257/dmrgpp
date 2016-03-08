/*
Copyright (c) 2012, UT-Battelle, LLC
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

/*! \file InputNg.h
 *
 *
 */

#ifndef INPUT_NG_H
#define INPUT_NG_H

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cassert>
#include "Vector.h"
#include <cstdlib>
#include <algorithm>
#include "Map.h"
#include "Matrix.h"
#include "loki/TypeTraits.h"

namespace PsimagLite {

template<typename InputCheckType>
class InputNg {

	class MyCompare {

		enum {FIRST,SECOND};

		typedef std::pair<String,SizeType> PairType;

	public:

		bool operator()(const String& x1,const String& x2)
		{
			PairType p1 = mysplit(x1);
			PairType p2 = mysplit(x2);
			if (p1.second==0 && p2.second==0)
				return (x1<x2);
			if (p1.second==0) return true;
			if (p2.second==0) return false;
			if (p1.first!=p2.first) return (x1<x2);
			return (p1.second<p2.second);
		}

	private:

		PairType mysplit(const String& x) const
		{
			SizeType mode=FIRST;
			String xfirst = "";
			String xsecond = "";
			for (SizeType i=0;i<x.length();i++) {
				if (x[i]=='@') {
					mode=SECOND;
					continue;
				}
				if (mode==FIRST) xfirst += x[i];
				else xsecond += x[i];
			}
			if (xsecond=="") return PairType(xfirst,0);
			return PairType(xfirst,atoi(xsecond.c_str()));
		}
	};

	typedef MyCompare MyCompareType;

	typedef typename Map<String,String,MyCompareType>::Type MapStrStrType;
	typedef typename Map<String,Vector<String>::Type,MyCompareType>::Type MapStrVecType;

public:

	class Writeable {

		enum {WHITESPACE,ENDOFLINE,EQUALSIGN,ALPHA_CHAR,NUMERIC_CHAR,COMMENT};

		enum {IN_LABEL,IN_VALUE_OR_LABEL,IN_VALUE_TEXT,IN_VALUE_NUMERIC,IN_COMMENT};

	public:

		Writeable(const String& file,const InputCheckType& inputCheck)
		    : data_(""),
		      line_(0),
		      state_(IN_LABEL),
		      numericVector_(0),
		      lastLabel_(""),
		      file_(file),
		      inputCheck_(inputCheck),
		      verbose_(false)
		{
			std::ifstream fin(file.c_str());
			if (!fin || !fin.good() || fin.bad()) {
				String s(__FILE__);
				s += " Cannot open file " + file + "\n";
				throw RuntimeError(s.c_str());
			}

			char c=0;
			while (!fin.eof()) {
				fin.get(c);
				data_ += c;
			}

			fin.close();
			check();

			if (verbose_) {
				std::cout<<"START\n";
				printMap(mapStrStr_,"StrStr");
				std::cout<<"END\nSTART\n";
				printMap(mapStrVec_,"StrVec");
				std::cout<<"END\n";
			}
		}

		Writeable(const String& data, const InputCheckType& inputCheck, int)
		    : data_(data),
		      line_(0),
		      state_(IN_LABEL),
		      numericVector_(0),
		      lastLabel_(""),
		      file_("-"),
		      inputCheck_(inputCheck),
		      verbose_(false)
		{
			check();

			if (verbose_) {
				std::cout<<"START\n";
				printMap(mapStrStr_,"StrStr");
				std::cout<<"END\nSTART\n";
				printMap(mapStrVec_,"StrVec");
				std::cout<<"END\n";
			}
		}

		void set(MapStrStrType& mapStrStr,
		         MapStrVecType& mapStrVec,
		         Vector<String>::Type& labelsForRemoval) const
		{
			mapStrStr=mapStrStr_;
			mapStrVec=mapStrVec_,
			        labelsForRemoval=labelsForRemoval_;
		}

		const String& filename() const { return file_; }

	private:

		void check()
		{
			String buffer="";
			for (SizeType i=0;i<data_.length();i++) {
				SizeType type = findTypeOf(data_.at(i));
				if (state_ == IN_COMMENT && type != ENDOFLINE)
					continue;
				switch(type) {
				case ENDOFLINE:
					line_++;
					if (state_==IN_COMMENT) {
						state_ = IN_LABEL;
						break;
					}
					if (buffer=="") break;
					saveBuffer(buffer,ENDOFLINE);
					buffer="";
					break;
				case WHITESPACE:
					if (buffer=="" || state_==IN_COMMENT) break;
					saveBuffer(buffer,WHITESPACE);
					buffer="";
					break;
				case EQUALSIGN:
					if (buffer=="" || state_==IN_COMMENT) break;
					saveBuffer(buffer,EQUALSIGN);
					buffer="";
					break;
				case COMMENT:
					state_=IN_COMMENT;
					break;
				default:
					if (state_==IN_VALUE_OR_LABEL) {
						if (type==ALPHA_CHAR) {
							checkNumbers();
							numericVector_.clear();
							state_=IN_LABEL;
						} else {
							state_=IN_VALUE_NUMERIC;
						}
					}
					buffer += data_.at(i);
					break;
				}
			}
			if (numericVector_.size()>0) checkNumbers();
		}

		void saveBuffer(const String& buffer,SizeType whatchar)
		{
			String s(__FILE__);
			String adjLabel="";
			switch(state_) {
			case IN_LABEL:
				if (verbose_) std::cout<<"Read label="<<buffer<<"\n";
				lastLabel_=buffer;
				if (whatchar==EQUALSIGN) state_=IN_VALUE_TEXT;
				else state_=IN_VALUE_NUMERIC;
				break;
			case IN_VALUE_OR_LABEL:
				std::cerr<<"Line="<<line_<<"\n";
				s += "Error while buffer=" + buffer;
				s += String(" and current line=") + String("\n");
				break;
			case IN_VALUE_TEXT:
				if (verbose_) std::cout<<"Read text value="<<buffer<<"\n";
				adjLabel = adjLabelForDuplicates(lastLabel_,mapStrStr_);
				mapStrStr_[adjLabel] = buffer;
				state_=IN_LABEL;
				inputCheck_.check(adjLabel,buffer,line_);
				break;
			case IN_VALUE_NUMERIC:
				if (verbose_) std::cout<<"Read numeric value="<<buffer<<"\n";
				numericVector_.push_back(buffer);
				state_=IN_VALUE_OR_LABEL;
				break;
			}
		}

		SizeType findTypeOf(char c) const
		{
			if (c=='\n') return ENDOFLINE;
			if (c==' ' || c=='\t') return WHITESPACE;
			if (c=='=') return EQUALSIGN;
			if (c>=48 && c<=58) return NUMERIC_CHAR;
			if (c=='.' || c=='+' || c=='-') return NUMERIC_CHAR;
			if (c=='(' || c==')' || c==',') return NUMERIC_CHAR;
			if (c=='#') return COMMENT;
			return ALPHA_CHAR;
		}

		void checkNumbers()
		{
			if (numericVector_.size()==1) {
				String s(__FILE__);
				s += " use equal sign instead of space in line "+ttos(line_) + "\n";
				throw RuntimeError(s.c_str());
				return;
			}
			String s(__FILE__);
			if (numericVector_.size()==0) {
				std::cerr<<"Line="<<line_<<"\n";
				throw RuntimeError(s.c_str());
			}
			SizeType adjExpected = atoi(numericVector_[0].c_str());

			if (!inputCheck_.check(lastLabel_,numericVector_,line_) &&
			    numericVector_.size()!=adjExpected+1) {
				std::cout<<" Number of numbers to follow is wrong, expected ";
				std::cout<<adjExpected<<" got ";
				std::cout<<(numericVector_.size()-1)<<"\n";
				std::cerr<<"Line="<<line_<<"\n";
				throw RuntimeError(s.c_str());
			}

			String adjLabel=adjLabelForDuplicates(lastLabel_,mapStrVec_);
			mapStrVec_[adjLabel]=numericVector_;

		}

		template<typename SomeMapType>
		String adjLabelForDuplicates(const String& label,SomeMapType& mymap)
		{
			String rootLabel = findRootLabel(label);
			int x = findLastOccurrence(rootLabel,mymap);
			if (x<0) return label;
			labelsForRemoval_.push_back(rootLabel);
			x++;
			String newlabel = rootLabel + "@" + ttos(x);
			if (verbose_) std::cerr<<"NEWLABEL=*"<<newlabel<<"*\n";
			return newlabel;
		}

		template<typename SomeMapType>
		int findLastOccurrence(const String& root1,SomeMapType& mymap)
		{
			int x = -1;
			for (typename SomeMapType::iterator it=mymap.begin();it!=mymap.end();++it) {
				String root2 = findRootLabel(it->first);
				if (root1==root2) x++;
			}
			return x;
		}

		template<typename MapType>
		typename EnableIf<IsMapLike<MapType>::True,void>::Type
		printMap(MapType& mp,const String& label)
		{
//			typedef typename MapType::key_type KeyType;
//			typedef typename MapType::mapped_type MappedType;
			std::cout<<label<<"\n";
			typename MapType::iterator it;
			for (it=mp.begin();it!=mp.end();++it) {
				std::cout<<it->first<<" "<<it->second<<"\n";
			}
		}

		String data_;
		SizeType line_;
		SizeType state_;
		Vector<String>::Type numericVector_;
		String lastLabel_;
		String file_;
		InputCheckType inputCheck_;
		bool verbose_;
		typename Map<String,String,MyCompareType>::Type mapStrStr_;
		typename Map<String,Vector<String>::Type,MyCompareType>::Type mapStrVec_;
		Vector<String>::Type labelsForRemoval_;
	}; // class Writeable

	class Readable {

		typedef typename Map<String,String,MyCompareType>::Type::iterator MapStringIteratorType;
		typedef typename Map<String,Vector<String>::Type,MyCompareType>::Type::iterator
		                 MapStringVectorIteratorType;

	public:

		Readable(const Writeable& inputWriteable) : file_(inputWriteable.filename())
		{
			inputWriteable.set(mapStrStr_,mapStrVec_,labelsForRemoval_);
		}

		template<typename SomeMemResolvType>
		SizeType memResolv(SomeMemResolvType& mres,
		                   SizeType,
		                   String msg = "") const
		{
			String str = msg;
			str += "InputNgReadable";

			const char* start = reinterpret_cast<const char *>(this);
			const char* end = reinterpret_cast<const char *>(&mapStrStr_);
			SizeType total = mres.memResolv(&file_, end-start, str + " file");

			start = end;
			end = reinterpret_cast<const char *>(&mapStrVec_);
			total += mres.memResolv(&mapStrStr_, end-start, str + " mapStrStr");

			start = end;
			end = reinterpret_cast<const char *>(&labelsForRemoval_);
			total += sizeof(mapStrVec_); //mres.memResolv(&mapStrVec_, end-start, str + " mapStrVec");

			total += mres.memResolv(&labelsForRemoval_,
			                        sizeof(*this) - total,
			                        str + " labelsForRemoval");

			return total;
		}

		void readline(String& val,
		              const String& label,
		              bool clean = true,
		              bool forceRemoval = false)
		{
			String label2 = label2label(label);
			MapStringIteratorType it = findFirstValueForLabel(label2,mapStrStr_);
			if (it==mapStrStr_.end()) throwWithMessage(label,label2);

			val= it->second.c_str();

			if (clean) cleanLabelsIfNeeded(label2,mapStrStr_,it,forceRemoval);
		}

		template<typename FloatingType>
		typename EnableIf<Loki::TypeTraits<FloatingType>::isFloat,void>::Type
		readline(FloatingType& val,const String& label)
		{
			String label2 = label2label(label);
			MapStringIteratorType it = findFirstValueForLabel(label2,mapStrStr_);
			if (it==mapStrStr_.end()) throwWithMessage(label,label2);

			val= atof(it->second.c_str());

			cleanLabelsIfNeeded(label2,mapStrStr_,it);
		}

		void readline(long long int& val,const String& label)
		{
			String label2 = label2label(label);
			MapStringIteratorType it = findFirstValueForLabel(label2,mapStrStr_);
			if (it==mapStrStr_.end()) throwWithMessage(label,label2);

			val= atoi(it->second.c_str());

			cleanLabelsIfNeeded(label2,mapStrStr_,it);
		}

		void readline(SizeType& val,const String& label)
		{
			String label2 = label2label(label);
			MapStringIteratorType it = findFirstValueForLabel(label2,mapStrStr_);
			if (it==mapStrStr_.end()) throwWithMessage(label,label2);

			val= atoi(it->second.c_str());

			cleanLabelsIfNeeded(label2,mapStrStr_,it);
		}

		void readline(int& val,const String& label)
		{
			String label2 = label2label(label);
			MapStringIteratorType it = findFirstValueForLabel(label2,mapStrStr_);
			if (it==mapStrStr_.end()) throwWithMessage(label,label2);

			val= atoi(it->second.c_str());

			cleanLabelsIfNeeded(label2,mapStrStr_,it);
		}

		void read(SizeType& val,const String& label)
		{
			String label2 = label2label(label);

			MapStringIteratorType it =  findFirstValueForLabel(label2,mapStrStr_);
			if (it==mapStrStr_.end()) throwWithMessage(label,label2);

			val= atoi(it->second.c_str());

			cleanLabelsIfNeeded(label2,mapStrStr_,it);
		}

		template<typename MapLikeType>
		typename EnableIf<IsMapLike<MapLikeType>::True,void>::Type
		read(MapLikeType& val,const String& label)
		{
			String label2 = label2label(label);

			typedef typename Map<String,String,MyCompareType>::Type::iterator MyIteratorType;
			for (MyIteratorType it=mapStrStr_.begin();it!=mapStrStr_.end();++it) {
				String mystr = it->first;
				long unsigned int it0 = mystr.find(label2);
				if (it0 == String::npos) continue;
				String::iterator it1 = find(mystr.begin(),mystr.end(),'[');
				if (it1 == mystr.end()) continue;
				String::iterator it2 = find(mystr.begin(),mystr.end(),']');
				if (it2 == mystr.end()) {
					String str("Malformed ");
					str += mystr + " entry\n";
					throw RuntimeError(str);
				}
				String mystr2 = mystr;
				mystr2.erase(0,label2.length()+1);
				mystr2.erase(mystr2.length()-1,1);
				val[mystr2] = atof(it->second.c_str());
			}
		}

		template<typename VectorLikeType>
		typename EnableIf<IsVectorLike<VectorLikeType>::True,void>::Type
		read(VectorLikeType& val,const String& label)
		{
			String label2 = label2label(label);
			typedef typename VectorLikeType::value_type NumericType;
			MapStringVectorIteratorType it = findFirstValueForLabel(label2,mapStrVec_);
			if (it==mapStrVec_.end()) throwWithMessage(label,label2);

			SizeType len =  it->second.size();
			assert(len>1);
			val.resize(len-1);
			for (SizeType i=0;i<len-1;i++) {
				val[i]=stringToComplexOrReal<NumericType>(it->second[i+1].c_str());
			}
			cleanLabelsIfNeeded(label2,mapStrVec_,it);
		}

		template<typename VectorLikeType>
		void readKnownSize(VectorLikeType& val,const String& label)
		{
			String label2 = label2label(label);
			typedef typename VectorLikeType::value_type NumericType;
			MapStringVectorIteratorType it =  findFirstValueForLabel(label2,mapStrVec_);
			if (it==mapStrVec_.end()) throwWithMessage(label,label2);

			SizeType len =  it->second.size();
			val.resize(len);
			for (SizeType i=0;i<len;i++) {
				val[i]=static_cast<NumericType>(atof(it->second[i].c_str()));
			}
			cleanLabelsIfNeeded(label2,mapStrVec_,it);
		}

		template<typename FloatingType>
		typename EnableIf<Loki::TypeTraits<FloatingType>::isFloat,void>::Type
		readMatrix(Matrix<FloatingType>& m,const String& label)
		{
			String label2 = label2label(label);

			MapStringVectorIteratorType it =  findFirstValueForLabel(label2,mapStrVec_);
			if (it==mapStrVec_.end()) throwWithMessage(label,label2);

			if (it->second.size()<2 || atoi(it->second[0].c_str())<=0 ||
			    atoi(it->second[1].c_str())<=0) {
				String s(__FILE__);
				s += " readMatrix: \n";
				throw RuntimeError(s.c_str());
			}

			SizeType nrow = SizeType(atoi(it->second[0].c_str()));
			SizeType ncol = SizeType(atoi(it->second[1].c_str()));
			m.resize(nrow,ncol);
			if (it->second.size()<2+nrow*ncol) {
				String s(__FILE__);
				s += " readMatrix: \n";
				throw RuntimeError(s.c_str());
			}
			SizeType k = 2;
			for (SizeType i=0;i<m.n_row();i++)
				for (SizeType j=0;j<m.n_col();j++)
					m(i,j) = atof(it->second[k++].c_str());

			cleanLabelsIfNeeded(label2,mapStrVec_,it);
		}

		template<typename FloatingType>
		typename EnableIf<Loki::TypeTraits<FloatingType>::isFloat,void>::Type
		readMatrix(Matrix<std::complex<FloatingType> >& m,
		           const String& label)
		{
			String label2 = label2label(label);

			MapStringVectorIteratorType it =  findFirstValueForLabel(label2,mapStrVec_);
			if (it==mapStrVec_.end()) throwWithMessage(label,label2);

			if (it->second.size()<2 || atoi(it->second[0].c_str())<=0 ||
			    atoi(it->second[1].c_str())<=0) {
				String s(__FILE__);
				s += " readMatrix: \n";
				throw RuntimeError(s.c_str());
			}

			SizeType nrow = SizeType(atoi(it->second[0].c_str()));
			SizeType ncol = SizeType(atoi(it->second[1].c_str()));
			m.resize(nrow,ncol);
			if (it->second.size()<2+nrow*ncol) {
				String s(__FILE__);
				s += " readMatrix: \n";
				throw RuntimeError(s.c_str());
			}
			SizeType k = 2;
			for (SizeType i=0;i<m.n_row();i++) {
				for (SizeType j=0;j<m.n_row();j++) {
					IstringStream is(it->second[k++]);
					is >> m(i,j);
				}
			}
			cleanLabelsIfNeeded(label2,mapStrVec_,it);
		}

		const String& filename() const
		{
			return file_;
		}

	private:

		template<typename SomeMapType>
		void cleanLabelsIfNeeded(const String& label,
		                         SomeMapType& mymap,
		                         typename SomeMapType::iterator& it,
		                         bool forceRemoval = false)
		{
			Vector<String>::Type::iterator it2 = find(labelsForRemoval_.begin(),
			                                          labelsForRemoval_.end(),
			                                          label);
			if (it2!=labelsForRemoval_.end() || forceRemoval) mymap.erase(it);
		}

		String label2label(const String& label)
		{
			SizeType len = label.length();
			if (len==0) {
				String s(__FILE__);
				s += " readline: label cannot be null\n";
				throw RuntimeError(s.c_str());
			}
			if (label.at(len-1)=='=') len--;
			return label.substr(0,len);
		}

		template<typename SomeMapType>
		typename SomeMapType::iterator findFirstValueForLabel(const String& label,
		                                                      SomeMapType& mymap)
		{
			for (typename SomeMapType::iterator it=mymap.begin();it!=mymap.end();++it) {
				String root2 = findRootLabel(it->first);
				if (label==root2) {
					return it;
				}
			}
			return mymap.end();
		}

		template<typename ComplexOrRealType>
		typename EnableIf<IsComplexNumber<ComplexOrRealType>::True,ComplexOrRealType>::Type
		stringToComplexOrReal(const String& s) const
		{
			typedef typename Real<ComplexOrRealType>::Type RealType;

			if (s[0]!='(') {
				return stringToReal(s.c_str());
			}

			String buffer("");
			SizeType start = 0;
			for (SizeType i = 1; i < s.length(); ++i) {
				start = i;
				if (s[i] == ',') break;
				buffer += s[i];
			}

			RealType r = stringToReal(buffer);

			start++;
			buffer = "";
			for (SizeType i = start; i < s.length(); ++i) {
				if (s[i] == ')') break;
				buffer += s[i];
			}

			RealType img = stringToReal(buffer);

			return ComplexOrRealType(r,img);
		}

		template<typename ComplexOrRealType>
		typename EnableIf<!IsComplexNumber<ComplexOrRealType>::True,
		typename Real<ComplexOrRealType>::Type>::Type
		stringToComplexOrReal(const String& s) const
		{
			return static_cast<typename Real<ComplexOrRealType>::Type>(stringToReal(s));
		}

		double stringToReal(const String& s) const
		{
			for (SizeType i = 0; i < s.length(); ++i) {
				char c = s[i];
				bool b1 = (c<48 || c>57);
				bool b2 = (c != '.' && c != '-' && c != '+' && c != 'e');
				if (b1 && b2) {
					String str = s +" is not a real number\n";
					str += "Suggestion: Add -DUSE_COMPLEX to your Makefile.\n";
					throw RuntimeError(str);
				}
			}

			return atof(s.c_str());
		}

		void throwWithMessage(const String& label,const String& label2="")
		{
			String s("Message issued by: ");
			s += String(__FILE__) + "\n";
			s += "ATTENTION: ERROR MESSAGE, PLEASE READ: ";
			s += " The (probably) mandatory label: " + label;
			if (label2.length()>0 && label2!=label) s += " (a.k.a. " + label2 +")";
			s += " was not found in the input file.\n";
			throw RuntimeError(s.c_str());
		}

		//serializr start class InputNgReadable
		//serializr normal file_
		String file_;
		//serializr normal mapStrStr_
		typename Map<String,String,MyCompareType>::Type mapStrStr_;
		//serializr normal mapStrVec_
		typename Map<String,Vector<String>::Type,MyCompareType>::Type mapStrVec_;
		//serializr normal labelsForRemoval_
		Vector<String>::Type labelsForRemoval_;
	}; // class Readable

	static String findRootLabel(const String& label)
	{
		String buffer="";
		SizeType len = label.length();
		for (SizeType i=0;i<len;i++) {
			if (label.at(i)=='@') break;
			buffer += label.at(i);
		}
		return buffer;
	}

}; //InputNg

class InputEmptyCheck {

};

} // namespace PsimagLite

/*@}*/

#endif // INPUT_NG_H

