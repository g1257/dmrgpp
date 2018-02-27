/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/** \ingroup DMRG */
/*@{*/

/*! \file VerySparseMatrix.h
 *
 *  A class to represent a sparse matrix in trivial format
 *
 */

#ifndef VERY_SPARSE_MATRIX_HEADER_H
#define VERY_SPARSE_MATRIX_HEADER_H

#include "Sort.h" // in PsimagLite

namespace Dmrg {
// Yet another sparse matrix class
template<class ComplexOrRealType>
class VerySparseMatrix {

	typedef std::pair<SizeType,SizeType> PairType;
	typedef PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;

public:

	typedef ComplexOrRealType value_type;

	VerySparseMatrix(SizeType rank) : rank_(rank),sorted_(true)
	{}

	template<typename CrsMatrixType>
	VerySparseMatrix(const CrsMatrixType& crs)
	    : rank_(crs.rows()),
	      values_(crs.nonZeros()),
	      coordinates_(crs.nonZeros()),
	      sorted_(true)
	{
		assert(crs.rows()==crs.cols());
		SizeType counter = 0;
		for (SizeType i=0;i<rank_;++i) {
			for (int k=crs.getRowPtr(i);k<crs.getRowPtr(i+1);++k) {
				// (i,crs.getCol(k)) --> coordinate
				coordinates_[counter] = PairType(i,crs.getCol(k));
				values_[counter++] = crs.getValue(k);
			}
		}
	}

	ComplexOrRealType& operator()(SizeType row,SizeType col)
	{
		std::pair<SizeType,SizeType> coordinate(row,col);
		int x=PsimagLite::isInVector(coordinates_,coordinate);
		if (x<0) {
			coordinates_.push_back(coordinate);
			values_.push_back(0);
			x=values_.size()-1;
			sorted_=false;
		}
		return values_[x];
	}

	bool operator!=(const VerySparseMatrix<ComplexOrRealType>& other) const
	{
		return !(*this == other);
	}

	void clear()
	{
		rank_=0;
		values_.clear();
		coordinates_.clear();
		sorted_=false;
	}

	template<typename CrsMatrixType>
	void operator=(const CrsMatrixType& crs)
	{
		clear();
		rank_=crs.rank();
		sorted_=true;
		for (SizeType i=0;i<rank_;i++)
			for (int k=crs.getRowPtr(i);k<crs.getRowPtr(i+1);k++)
				// (i,crs.getCol(k)) --> coordinate
				this->operator()(i,crs.getCol(k))=crs.getValue(k);
		// crs.getValue(k) --> value
	}

	//! same as T& operator() but doesn't check for dupes
	//! faster than T& operator() but use with care.
	void set(SizeType row,SizeType col,const ComplexOrRealType& value)
	{
		std::pair<SizeType,SizeType> coordinate(row,col);
		coordinates_.push_back(coordinate);
		values_.push_back(value);
	}

	ComplexOrRealType operator()(SizeType row,SizeType col) const
	{
		std::pair<SizeType,SizeType> coordinate(row,col);
		int x=PsimagLite::isInVector(coordinates_,coordinate);
		if (x<0) {
			//T value=0;
			return 0;
		}
		return values_[x];
	}

	const ComplexOrRealType& operator()(SizeType i) const
	{
		return values_[i];
	}

	ComplexOrRealType& operator()(SizeType i)
	{
		return values_[i];
	}

	void operator+=(VerySparseMatrix<ComplexOrRealType>& other)
	{
		if (sorted_ && other.sorted())
			plusEqualOrd(other);
		else
			throw PsimagLite::RuntimeError("VerySparseMatrix::operator+=(): unsorted\n");
	}

	void sort()
	{
		VectorSizeType perm;
		sort(perm);
	}

	void sort(VectorSizeType& perm)
	{
		VectorSizeType rows(coordinates_.size());
		VectorSizeType cols(coordinates_.size());
		VectorComplexOrRealType values = values_;
		for (SizeType i=0;i<coordinates_.size();i++) {
			rows[i]=coordinates_[i].first;
			cols[i]=coordinates_[i].second;

		}

		perm.resize(rows.size());
		PsimagLite::Sort<VectorSizeType> sort;
		sort.sort(rows,perm);

		for (SizeType i=0;i<coordinates_.size();i++) {
			coordinates_[i].first = rows[i];
			coordinates_[i].second = cols[perm[i]];
			values_[i] = values[perm[i]];
		}
		sorted_=true;
	}

	SizeType rows() const { return rank_; }

	SizeType cols() const { return rank_; }

	SizeType getRow(SizeType i) const
	{
		assert(i < coordinates_.size());
		return coordinates_[i].first;
	}

	void getRow(VectorSizeType& cols,
	            SizeType row,
	            SizeType startIndex=0) const
	{
		cols.clear();
		for (SizeType i=startIndex;i<coordinates_.size();i++) {
			if (coordinates_[i].first == row) cols.push_back(i);
			if (sorted_ && coordinates_[i].first > row) break;
		}
	}

	SizeType getColumn(SizeType i) const
	{
		assert(i < coordinates_.size());
		return coordinates_[i].second;
	}

	void getColumn(VectorSizeType& rows,
	               SizeType col) const
	{
		rows.clear();
		for (SizeType i=0;i<coordinates_.size();i++)
			if (coordinates_[i].second == col) rows.push_back(i);
	}

	SizeType nonZeros() const { return values_.size(); }

	ComplexOrRealType getValue(SizeType i) const
	{
		assert(i < values_.size());
		return values_[i];
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const VerySparseMatrix<ComplexOrRealType>& m)
	{
		os<<m.rank_<<" "<<m.rank_;
		if (m.rank_ == 0) return os;
		for (SizeType i=0;i<m.values_.size();i++) {
			os<<m.coordinates_[i].first<<" ";
			os<<m.coordinates_[i].second<<" "<<m.values_[i]<<"\n";
		}

		return os;
	}

	friend std::istream& operator>>(std::istream& is,
	                                VerySparseMatrix<ComplexOrRealType>& m)
	{
		is>>m.rank_;
		if (m.rank_ == 0) return is;
		PsimagLite::String temp;
		is>>temp;
		assert(static_cast<SizeType>(atoi(temp.c_str())) == m.rank_);
		PairType coordinate;

		while (true) {
			is>>temp;
			if (temp[0] == '#') break;

			coordinate.first = atoi(temp.c_str());
			is>>temp;
			coordinate.second = atoi(temp.c_str());
			m.coordinates_.push_back(coordinate);

			is>>temp;
			m.values_.push_back(atof(temp.c_str()));
		}

		return is;
	}

	template<typename IoType>
	void saveToDisk(IoType& outHandle)
	{
		PsimagLite::String s = "rank="+ttos(rank_);
		outHandle.printline(s);
		outHandle.write(coordinates_,"coordinates");
		outHandle.write(values_,"values");
		s="#######\n";
		outHandle.printline(s);
	}

	template<typename IoType>
	void loadFromDisk(IoType& inHandle,bool check=false)
	{
		clear();
		PsimagLite::String s = "rank=";
		inHandle.readline(rank_,s);
		inHandle.read(coordinates_,"coordinates");
		if (check) checkCoordinates();
		inHandle.read(values_,"values");
		sorted_=true;
	}

	//! for debuggin only:
	void checkCoordinates() const
	{
		SizeType flag=0;
		for (SizeType i=0;i<coordinates_.size();i++) {
			if (coordinates_[i].first<0 || 	coordinates_[i].first>=rank_) {
				std::cerr<<"coordinates["<<i<<"].first="<< coordinates_[i].first<<"\n";
				flag=1;
				break;
			}
			if (coordinates_[i].second<0 || 	coordinates_[i].second>=rank_) {
				std::cerr<<"coordinates["<<i<<"].second="<< coordinates_[i].second<<"\n";
				flag=2;
				break;
			}
		}
		if (flag==0) return;
		std::cerr<<"rank="<<rank_<<"\n";
		throw PsimagLite::RuntimeError("VerySparseMatrix::checkCoordinates()\n");
	}

	bool sorted() const { return sorted_; }


private:

	void plusEqualOrd(VerySparseMatrix<ComplexOrRealType>& other)
	{
		SizeType j = 0;
		SizeType i = 0;
		// pre-alloc memory:
		VectorPairType newcoord(coordinates_.size()+other.coordinates_.size());
		VectorComplexOrRealType newvals(coordinates_.size()+other.coordinates_.size());
		SizeType counter = 0;
		while(i<coordinates_.size() && j<other.coordinates_.size()) {
			PairType row = coordinates_[i];
			PairType row2 = other.coordinates_[j];
			//std::cerr<<row<<" "<<row2<<"\n";
			if (row2 < row) {
				//newcoord.push_back(j);
				newcoord[counter]=row2;
				newvals[counter]=other.values_[j];
				counter++;
				j++;
			} else if (row2 == row) {
				newcoord[counter]=row2;
				newvals[counter]=values_[i] + other.values_[j];
				j++;
				i++;
				counter++;
			} else {
				newcoord[counter]=row;
				newvals[counter]=values_[i];
				i++;
				counter++;
			}
		}
		for (SizeType j1=j;j1<other.coordinates_.size();j1++) {
			newcoord[counter]=other.coordinates_[j1];
			newvals[counter]=other.values_[j1];
			counter++;
		}
		for (SizeType j1=i;j1<coordinates_.size();j1++) {
			newcoord[counter]=coordinates_[j1];
			newvals[counter]=values_[j1];
			counter++;
		}
		coordinates_.resize(counter);
		values_.resize(counter);
		for (SizeType j1=0;j1<counter;j1++) {
			coordinates_[j1] = newcoord[j1];
			values_[j1] = newvals[j1];
		}
		// this is sorted
	}

	bool notEqual(const char *s) const
	{
		std::cerr<<"notEqual="<<s<<"\n";
		return false;
	}

	SizeType rank_;
	VectorComplexOrRealType values_;
	VectorPairType coordinates_;
	bool sorted_;
}; // VerySparseMatrix

template<typename T>
bool isHermitian(const VerySparseMatrix<T>& m)
{
	T eps =1e-6;
	for (SizeType i=0;i<m.nonZeros();i++) {
		SizeType row = m.getRow(i);
		SizeType col = m.getColumn(i);
		if (fabs(m.getValue(i)-m(col,row))>eps) {
			return false;
		}
	}
	return true;
}
} // namespace Dmrg

namespace PsimagLite {
template<typename T>
struct IsMatrixLike<Dmrg::VerySparseMatrix<T> >
{
	enum {True = true};
};
}
/*@}*/
#endif
