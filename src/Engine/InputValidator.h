/*
Copyright (c) 2012, UT-Battelle, LLC
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

/*! \file InputValidator.h
 *
 *
 */

#ifndef INPUT_VALIDATOR_H
#define INPUT_VALIDATOR_H

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include <cstdlib>
#include "Matrix.h"
#include "Utils.h"

namespace PsimagLite {

class InputValidator {
	
	enum {WHITESPACE,ENDOFLINE,EQUALSIGN,ALPHA_CHAR,NUMERIC_CHAR};
	
	enum {IN_LABEL,IN_VALUE_OR_LABEL,IN_VALUE_TEXT,IN_VALUE_NUMERIC};
	
public:
	
	InputValidator(const std::string& file)
	: file_(file),
	  data_(""),
	  line_(0),
	  stage_(IN_LABEL),
	  numericVector_(0),
	  lastLabel_(""),
	  MagicLabel_("FiniteLoops"),
	  labelsWithKnownSize_(3),
	  verbose_(false)
	{
		labelsWithKnownSize_[0] = "JMVALUES";
		labelsWithKnownSize_[1] = "RAW_MATRIX";
		labelsWithKnownSize_[2] = "Connectors";

		std::ifstream fin(file.c_str());
		if (!fin || !fin.good() || fin.bad()) {
			std::string s(__FILE__);
			s += " Cannot open file " + file + "\n";
			throw std::runtime_error(s.c_str());
		}
		
		char c=0;
		while(!fin.eof()) {
			fin.get(c);
			data_ += c;
		}
		fin.close();
		check();
		if (verbose_) {
			printMap(mapStrStr_,"StrStr");
			printMap(mapStrVec_,"StrVec");
		}
	}

	void readline(std::string& val,const std::string& label)
	{
		std::string label2 = label2label(label);
		std::map<std::string,std::string>::iterator it =  findFirstValueForLabel(label2,mapStrStr_);
		if (it==mapStrStr_.end()) throw std::runtime_error("InputValidator");

		val= it->second.c_str();

		cleanLabelsIfNeeded(label2,mapStrStr_,it);
	}

	void readline(double& val,const std::string& label)
	{
		std::string label2 = label2label(label);
		std::map<std::string,std::string>::iterator it =  findFirstValueForLabel(label2,mapStrStr_);
		if (it==mapStrStr_.end()) throw std::runtime_error("InputValidator");

		val= atof(it->second.c_str());

		cleanLabelsIfNeeded(label2,mapStrStr_,it);
	}

	void readline(long long int& val,const std::string& label)
	{
		std::string label2 = label2label(label);
		std::map<std::string,std::string>::iterator it =  findFirstValueForLabel(label2,mapStrStr_);
		if (it==mapStrStr_.end()) throw std::runtime_error("InputValidator");

		val= atoi(it->second.c_str());

		cleanLabelsIfNeeded(label2,mapStrStr_,it);
	}

	void readline(size_t& val,const std::string& label)
	{
		std::string label2 = label2label(label);
		std::map<std::string,std::string>::iterator it =  findFirstValueForLabel(label2,mapStrStr_);
		if (it==mapStrStr_.end()) throw std::runtime_error("InputValidator");

		val= atoi(it->second.c_str());

		cleanLabelsIfNeeded(label2,mapStrStr_,it);
	}

	void readline(int& val,const std::string& label)
	{
		std::string label2 = label2label(label);
		std::map<std::string,std::string>::iterator it =  findFirstValueForLabel(label2,mapStrStr_);
		if (it==mapStrStr_.end()) throw std::runtime_error("InputValidator");

		val= atoi(it->second.c_str());

		cleanLabelsIfNeeded(label2,mapStrStr_,it);
	}

	void read(long unsigned int& val,const std::string& label)
	{
		std::string label2 = label2label(label);

		std::map<std::string,std::string>::iterator it =  findFirstValueForLabel(label2,mapStrStr_);
		if (it==mapStrStr_.end()) throw std::runtime_error("InputValidator");


		val= atoi(it->second.c_str());

		cleanLabelsIfNeeded(label2,mapStrStr_,it);
	}

	template<typename T>
	void read(std::vector<T>& val,const std::string& label)
	{
		std::string label2 = label2label(label);

		std::map<std::string,std::vector<double> >::iterator it =  findFirstValueForLabel(label2,mapStrVec_);
		if (it==mapStrVec_.end()) throw std::runtime_error("InputValidator");

		size_t len =  it->second.size();
		val.resize(len-1);
		for (size_t i=0;i<len-1;i++) {
			val[i]=static_cast<T>(it->second[i+1]);
		}
		cleanLabelsIfNeeded(label2,mapStrVec_,it);
	}

	template<typename T>
	void readKnownSize(std::vector<T>& val,const std::string& label)
	{
		std::string label2 = label2label(label);

		std::map<std::string,std::vector<double> >::iterator it =  findFirstValueForLabel(label2,mapStrVec_);
		if (it==mapStrVec_.end()) throw std::runtime_error("InputValidator");

		size_t len =  it->second.size();
		val.resize(len);
		for (size_t i=0;i<len;i++) {
			val[i]=static_cast<T>(it->second[i+1]);
		}
		cleanLabelsIfNeeded(label2,mapStrVec_,it);
	}

	template<typename T>
	void readMatrix(PsimagLite::Matrix<T>& m,const std::string& label)
	{
		std::string label2 = label2label(label);

		std::map<std::string,std::vector<double> >::iterator it =  findFirstValueForLabel(label2,mapStrVec_);
		if (it==mapStrVec_.end()) throw std::runtime_error("InputValidator");


		std::vector<double>& v = it->second;
		if (v.size()<2 || v[0]<=0 || v[1]<=0) {
			std::string s(__FILE__);
			s += " readMatrix: \n";
			throw std::runtime_error(s.c_str());
		}
		size_t nrow = size_t(v[0]);
		size_t ncol = size_t(v[1]);
		m.resize(nrow,ncol);
		if (v.size()<2+nrow*ncol) {
			std::string s(__FILE__);
			s += " readMatrix: \n";
			throw std::runtime_error(s.c_str());
		}
		size_t k = 2;
		for (size_t i=0;i<m.n_row();i++)
			for (size_t j=0;j<m.n_row();j++)
				m(i,j) = v[k++];

		cleanLabelsIfNeeded(label2,mapStrVec_,it);
	}

private:

	template<typename SomeMapType>
	void cleanLabelsIfNeeded(const std::string& label,SomeMapType& mymap,typename SomeMapType::iterator& it)
	{
		std::vector<std::string>::iterator it2 = find(labelsForRemoval_.begin(),labelsForRemoval_.end(),label);
		if (it2!=labelsForRemoval_.end()) mymap.erase(it);
	}

	template<typename T1,typename T2>
	void printMap(std::map<T1,T2>& mp,const std::string& label)
	{
		std::cout<<label<<"\n";
		typename  std::map<T1,T2>::iterator it;
		for (it=mp.begin();it!=mp.end();++it) {
			std::cout<<it->first<<" "<<it->second<<"\n";
		}
	}

	std::string label2label(const std::string& label)
	{
		size_t len = label.length();
		if (len==0) {
			std::string s(__FILE__);
			s += " readline: label cannot be null\n";
			throw std::runtime_error(s.c_str());
		}
		if (label.at(len-1)=='=') len--;
		return label.substr(0,len);
	}

	void check()
	{
		std::string buffer="";
		for (size_t i=0;i<data_.length();i++) {
			size_t type = findTypeOf(data_.at(i));
			switch(type) {
			case ENDOFLINE:
				line_++;
				if (buffer=="") break;
				saveBuffer(buffer,ENDOFLINE);
				buffer="";
				break;
			case WHITESPACE:
				if (buffer=="") break;
				saveBuffer(buffer,WHITESPACE);
				buffer="";
				break;
			case EQUALSIGN:
				if (buffer=="") break;
				saveBuffer(buffer,EQUALSIGN);
				buffer="";
				break;
			default:
				if (stage_==IN_VALUE_OR_LABEL) {
					if (type==ALPHA_CHAR) {
						checkNumbers();
						numericVector_.clear();
						stage_=IN_LABEL;
					} else {
						stage_=IN_VALUE_NUMERIC;
					}
				}
				buffer += data_.at(i);
				break;
			}
		}
		if (numericVector_.size()>0) checkNumbers();
	}
	
	void saveBuffer(const std::string& buffer,size_t whatchar)
	{
		std::string s(__FILE__);
		std::string adjLabel="";
		switch(stage_) {
		case IN_LABEL:
			if (verbose_) std::cout<<"Read label="<<buffer<<"\n";
			lastLabel_=buffer;
			if (whatchar==EQUALSIGN) stage_=IN_VALUE_TEXT;
			else stage_=IN_VALUE_NUMERIC;
			break;
		case IN_VALUE_OR_LABEL:
			std::cerr<<"Line="<<line_<<"\n";
			s += "Error while buffer=" + buffer;
			s += std::string(" and current line=") + std::string("\n");
			break;
		case IN_VALUE_TEXT:
			if (verbose_) std::cout<<"Read text value="<<buffer<<"\n";
			adjLabel = adjLabelForDuplicates(lastLabel_,mapStrStr_);
			mapStrStr_[adjLabel] = buffer;
			stage_=IN_LABEL;
			break;
		case IN_VALUE_NUMERIC:
			if (verbose_) std::cout<<"Read numeric value="<<buffer<<"\n";
			numericVector_.push_back(atof(buffer.c_str()));
			stage_=IN_VALUE_OR_LABEL;
			break;
		}
	}
	
	size_t findTypeOf(char c) const
	{
		if (c=='\n') return ENDOFLINE;
		if (c==' ' || c=='\t') return WHITESPACE;
		if (c=='=') return EQUALSIGN;
		if (c>=48 && c<=58) return NUMERIC_CHAR;
		if (c=='.' || c=='+' || c=='-') return NUMERIC_CHAR;
		return ALPHA_CHAR;
	}
	
	void checkNumbers()
	{
		if (numericVector_.size()==1) {
			std::string s(__FILE__);
			s += " use equal sign instead of space in line "+ttos(line_) + "\n";
			throw std::runtime_error(s.c_str());
			return;
		}
		std::string s(__FILE__);
		if (numericVector_.size()==0) {
			std::cerr<<"Line="<<line_<<"\n";
			throw std::runtime_error(s.c_str());
		}
		size_t adjExpected = size_t(numericVector_[0]);
		if (lastLabel_==MagicLabel_) adjExpected *= 3;

		int x = PsimagLite::isInVector(labelsWithKnownSize_,lastLabel_);
		if (x>=0) adjExpected = numericVector_.size()-1;
		if (numericVector_.size()==adjExpected+1) {
			std::string adjLabel=adjLabelForDuplicates(lastLabel_,mapStrVec_);
			mapStrVec_[adjLabel]=numericVector_;
			return;
		}
			
		std::cout<<" Number of numbers to follow is wrong, expected ";
		std::cout<<(numericVector_.size()-1)<<" got "<<adjExpected<<"\n";
		std::cerr<<"Line="<<line_<<"\n";
		throw std::runtime_error(s.c_str());
	}

	template<typename SomeMapType>
	std::string adjLabelForDuplicates(const std::string& label,SomeMapType& mymap)
	{
		std::string rootLabel = findRootLabel(label);
		int x = findLastOccurrence(rootLabel,mymap);
		if (x<0) return label;
		labelsForRemoval_.push_back(rootLabel);
		x++;
		std::string newlabel = rootLabel + "@" + ttos(x);
		if (verbose_) std::cerr<<"NEWLABEL=*"<<newlabel<<"*\n";
		return newlabel;
	}

	template<typename SomeMapType>
	int findLastOccurrence(const std::string& root1,SomeMapType& mymap)
	{
		int x = -1;
		for (typename SomeMapType::iterator it=mymap.begin();it!=mymap.end();++it) {
			std::string root2 = findRootLabel(it->first);
			if (root1==root2) x++;
		}
		return x;
	}

	std::string findRootLabel(const std::string& label) const
	{
		std::string buffer="";
		size_t len = label.length();
		for (size_t i=0;i<len;i++) {
			if (label.at(i)=='@') break;
			buffer += label.at(i);
		}
		return buffer;
	}

	template<typename SomeMapType>
	typename SomeMapType::iterator findFirstValueForLabel(const std::string& label,SomeMapType& mymap)
	{
		for (typename SomeMapType::iterator it=mymap.begin();it!=mymap.end();++it) {
			std::string root2 = findRootLabel(it->first);
			if (label==root2) {
				return it;
			}
		}
		return mymap.end();
	}
	
	std::string file_,data_;
	size_t line_;
	size_t stage_;
	std::vector<double> numericVector_;
	std::string lastLabel_;
	const std::string MagicLabel_;
	std::vector<std::string> labelsWithKnownSize_;
	bool verbose_;
	std::map<std::string,std::string> mapStrStr_;
	std::map<std::string,std::vector<double> > mapStrVec_;
	std::vector<std::string> labelsForRemoval_;
}; //InputValidator
} // namespace PsimagLite

/*@}*/

#endif // INPUT_VALIDATOR_H
