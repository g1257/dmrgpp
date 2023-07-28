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

#include "CrsMatrix.h"
#include "Sort.h" // in PsimagLite

namespace Dmrg
{
// Yet another sparse matrix class
template <class ComplexOrRealType>
class VerySparseMatrix
{

	typedef std::pair<SizeType, SizeType> PairType;
	typedef PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;

public:

	typedef ComplexOrRealType value_type;

	explicit VerySparseMatrix(SizeType rows, SizeType cols)
	    : rows_(rows)
	    , cols_(cols)
	    , sorted_(true)
	{
	}

	explicit VerySparseMatrix(const PsimagLite::CrsMatrix<ComplexOrRealType>& crs)
	    : rows_(crs.rows())
	    , cols_(crs.cols())
	    , values_(crs.nonZeros())
	    , coordinates_(crs.nonZeros())
	    , sorted_(true)
	{
		SizeType counter = 0;
		for (SizeType i = 0; i < rows_; ++i) {
			for (int k = crs.getRowPtr(i); k < crs.getRowPtr(i + 1); ++k) {
				// (i,crs.getCol(k)) --> coordinate
				coordinates_[counter] = PairType(i, crs.getCol(k));
				values_[counter++] = crs.getValue(k);
			}
		}

		sort();
	}

	void makeDiagonal(SizeType n, ComplexOrRealType value)
	{
		rows_ = cols_ = n;
		sorted_ = true;
		coordinates_.resize(n);
		values_.resize(n);
		for (SizeType i = 0; i < n; ++i) {
			coordinates_[i] = PairType(i, i);
			values_[i] = value;
		}
	}

	void resize(SizeType rows, SizeType cols)
	{
		clear();
		rows_ = rows;
		cols_ = cols;
	}

	void resize(SizeType rows, SizeType cols, SizeType nonzeros)
	{
		coordinates_.resize(nonzeros);
		values_.resize(nonzeros);
		sorted_ = false;
		rows_ = rows;
		cols_ = cols;
	}

	ComplexOrRealType& operator()(SizeType row, SizeType col)
	{
		std::pair<SizeType, SizeType> coordinate(row, col);
		int x = PsimagLite::indexOrMinusOne(coordinates_, coordinate);
		if (x < 0) {
			coordinates_.push_back(coordinate);
			values_.push_back(0);
			x = values_.size() - 1;
			sorted_ = false;
		}

		return values_[x];
	}

	bool operator!=(const VerySparseMatrix<ComplexOrRealType>& other) const
	{
		return !(*this == other);
	}

	void clear()
	{
		rows_ = cols_ = 0;
		values_.clear();
		coordinates_.clear();
		sorted_ = false;
	}

	void operator=(const PsimagLite::CrsMatrix<ComplexOrRealType>& crs)
	{
		clear();
		rows_ = crs.rows();
		cols_ = crs.cols();
		sorted_ = true;
		for (SizeType i = 0; i < rows_; i++)
			for (int k = crs.getRowPtr(i); k < crs.getRowPtr(i + 1); k++)
				// (i,crs.getCol(k)) --> coordinate
				this->operator()(i, crs.getCol(k)) = crs.getValue(k);
		// crs.getValue(k) --> value
	}

	//! same as T& operator() but doesn't check for dupes
	//! faster than T& operator() but use with care.
	void push(SizeType row, SizeType col, const ComplexOrRealType& value)
	{
		std::pair<SizeType, SizeType> coordinate(row, col);
		coordinates_.push_back(coordinate);
		values_.push_back(value);
	}

	void set(SizeType counter,
	    SizeType row,
	    SizeType col,
	    const ComplexOrRealType& value)
	{
		assert(values_.size() == coordinates_.size());
		assert(counter < coordinates_.size());
		coordinates_[counter] = PairType(row, col);
		values_[counter] = value;
	}

	const ComplexOrRealType& operator()(SizeType i) const
	{
		return values_[i];
	}

	ComplexOrRealType& operator()(SizeType i)
	{
		return values_[i];
	}

	void transposeConjugate()
	{
		SizeType n = values_.size();
		assert(n == coordinates_.size());
		std::swap(rows_, cols_);
		for (SizeType i = 0; i < n; ++i) {
			values_[i] = PsimagLite::conj(values_[i]);
			std::swap(coordinates_[i].first, coordinates_[i].second);
		}
	}

	void operator+=(const VerySparseMatrix<ComplexOrRealType>& other)
	{
		if (sorted_ && other.sorted())
			plusEqualOrd(other);
		else
			throw PsimagLite::RuntimeError("VerySparseMatrix::operator+=(): unsorted\n");
	}

	void sort()
	{
		sorted_ = true;

		SizeType n = coordinates_.size();

		if (n == 0)
			return;

		VectorSizeType perm(n, 0);
		PsimagLite::Sort<VectorPairType> sort;
		sort.sort(coordinates_, perm);

		VectorComplexOrRealType values(n, 0.0);
		for (SizeType i = 0; i < n; ++i)
			values[i] = values_[perm[i]];

		// uniqueness
		PairType prevCoord = coordinates_[0];
		PsimagLite::Vector<bool>::Type mark(n, false);
		SizeType ind = 0;
		ComplexOrRealType val = values[ind];
		SizeType newsize = 0;
		for (SizeType i = 1; i < n; ++i) {
			if (prevCoord != coordinates_[i]) {
				prevCoord = coordinates_[i];
				values[ind] = val;
				val = values[i];
				ind = i;
				++newsize;
				continue;
			}

			mark[i] = true;
			val += values[i];
		}

		values[ind] = val;
		++newsize;

		VectorPairType coords(newsize);
		values_.resize(newsize);
		SizeType j = 0;
		for (SizeType i = 0; i < n; ++i) {
			if (mark[i])
				continue;
			coords[j] = coordinates_[i];
			values_[j++] = values[i];
		}

		assert(j == newsize);
		coordinates_.swap(coords);
	}

	SizeType rows() const { return rows_; }

	SizeType cols() const { return cols_; }

	SizeType getRow(SizeType i) const
	{
		assert(i < coordinates_.size());
		return coordinates_[i].first;
	}

	void getRow(VectorSizeType& cols,
	    SizeType row,
	    SizeType startIndex = 0) const
	{
		cols.clear();
		for (SizeType i = startIndex; i < coordinates_.size(); i++) {
			if (coordinates_[i].first == row)
				cols.push_back(i);
			if (sorted_ && coordinates_[i].first > row)
				break;
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
		for (SizeType i = 0; i < coordinates_.size(); i++)
			if (coordinates_[i].second == col)
				rows.push_back(i);
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
		os << m.rows_ << " " << m.cols_;
		if (m.rows_ == 0 || m.cols_ == 0)
			return os;
		for (SizeType i = 0; i < m.values_.size(); i++) {
			os << m.coordinates_[i].first << " ";
			os << m.coordinates_[i].second << " " << m.values_[i] << "\n";
		}

		return os;
	}

	friend std::istream& operator>>(std::istream& is,
	    VerySparseMatrix<ComplexOrRealType>& m)
	{
		is >> m.rows_;
		is >> m.cols_;
		if (m.rows_ == 0 || m.cols_ == 0)
			return is;
		PsimagLite::String temp;
		PairType coordinate;

		while (true) {
			is >> temp;
			if (temp[0] == '#')
				break;

			coordinate.first = atoi(temp.c_str());
			is >> temp;
			coordinate.second = atoi(temp.c_str());
			m.coordinates_.push_back(coordinate);

			is >> temp;
			m.values_.push_back(atof(temp.c_str()));
		}

		return is;
	}

	template <typename IoType>
	void saveToDisk(IoType& outHandle)
	{
		PsimagLite::String s = "rows=" + ttos(rows_);
		outHandle.printline(s);
		outHandle.printline("cols=" + ttos(cols_));
		outHandle.write(coordinates_, "coordinates");
		outHandle.write(values_, "values");
		s = "######\n";
		outHandle.printline(s);
	}

	template <typename IoType>
	void loadFromDisk(IoType& inHandle, bool check = false)
	{
		clear();
		inHandle.readline(rows_, "rows=");
		inHandle.readline(cols_, "cols=");
		inHandle.read(coordinates_, "coordinates");
		if (check)
			checkCoordinates();
		inHandle.read(values_, "values");
		sorted_ = true;
	}

	//! for debuggin only:
	void checkCoordinates() const
	{
		SizeType flag = 0;
		for (SizeType i = 0; i < coordinates_.size(); i++) {
			if (coordinates_[i].first < 0 || coordinates_[i].first >= rows_) {
				std::cerr << "coordinates[" << i << "].first=" << coordinates_[i].first << "\n";
				flag = 1;
				break;
			}
			if (coordinates_[i].second < 0 || coordinates_[i].second >= cols_) {
				std::cerr << "coordinates[" << i << "].second=" << coordinates_[i].second << "\n";
				flag = 2;
				break;
			}
		}

		if (flag == 0)
			return;
		std::cerr << "rank=" << rows_ << " " << cols_ << "\n";
		throw PsimagLite::RuntimeError("VerySparseMatrix::checkCoordinates()\n");
	}

	bool sorted() const { return sorted_; }

private:

	void plusEqualOrd(const VerySparseMatrix<ComplexOrRealType>& other)
	{
		if (coordinates_.size() == 0) {
			*this = other;
			return;
		}

		if (other.coordinates_.size() == 0)
			return;

		// pre-alloc memory:
		VectorPairType newcoord(coordinates_.size());
		VectorComplexOrRealType newvals(coordinates_.size());

		newcoord.reserve(coordinates_.size() + other.coordinates_.size());
		newvals.reserve(coordinates_.size() + other.coordinates_.size());

		// ----------------------------
		// check coordinates are sorted
		// ----------------------------

		SizeType n = coordinates_.size();
		if (n > 0)
			--n;
		for (SizeType i = 0; i < n; ++i) {
			assert(coordinates_[i] <= coordinates_[i + 1]);
		}

		SizeType m = other.coordinates_.size();
		if (m > 0)
			--m;
		for (SizeType i = 0; i < m; ++i) {
			assert(other.coordinates_[i] <= other.coordinates_[i + 1]);
		};

		// ---------------------------------------------------------
		// initialization to place something reasonable in newcoord
		// ---------------------------------------------------------

		// ------------------------------------------------------
		// note use "maxrc" as sentinel end marker to simplify coding
		// ------------------------------------------------------
		const unsigned int bigval = 1024 * 1024 * 1024;
		const PairType maxrc(bigval, bigval);

		SizeType jp = 0;
		SizeType ip = 0;
		SizeType counter = 0;

		PairType irc = (coordinates_.size() > 0) ? coordinates_[ip] : maxrc;
		PairType jrc = (other.coordinates_.size() > 0) ? other.coordinates_[jp] : maxrc;
		if (irc < jrc) {
			newcoord[0] = irc;
			newvals[0] = values_[0];
			ip = 1;
		} else {
			newcoord[0] = jrc;
			newvals[0] = other.values_[0];
			jp = 1;
		};

		bool has_work_ip = (ip < coordinates_.size());
		bool has_work_jp = (jp < other.coordinates_.size());
		bool has_work = has_work_ip || has_work_jp;
		while (has_work) {

			has_work_ip = (ip < coordinates_.size());
			has_work_jp = (jp < other.coordinates_.size());
			has_work = has_work_ip || has_work_jp;
			if (!has_work) {
				break;
			};

			// ------------------------------------------------------
			// unified treatment of  coordinate_ or other.coordinate_
			// ------------------------------------------------------
			PairType irc = (has_work_ip) ? coordinates_[ip] : maxrc;
			ComplexOrRealType val_irc = (has_work_ip) ? values_[ip] : 0;

			PairType jrc = (has_work_jp) ? other.coordinates_[jp] : maxrc;
			ComplexOrRealType val_jrc = (has_work_jp) ? other.values_[jp] : 0;

			PairType ijrc = (irc < jrc) ? irc : jrc;
			ComplexOrRealType val_ijrc = (irc < jrc) ? val_irc : val_jrc;

			if (irc < jrc) {
				// ----------------------------
				// pick coordinate_, advance ip
				// ----------------------------
				ip++;
			} else {
				// ----------------------------------
				// pick other.coordinate_, advance jp
				// ----------------------------------
				jp++;
			};

			PairType newrc = newcoord[counter];
			if (newrc == ijrc) {
				// ------------------------
				// common entry: just add accumulate value
				// ------------------------
				newvals[counter] += val_ijrc;
			} else {
				// -------------------------
				// advance to next new entry
				// -------------------------
				counter++;
				newvals[counter] = val_ijrc;
				newcoord[counter] = ijrc;
			};
		};

		// update vector
		// TODO: add swap method to just swap pointers

		SizeType nonzeros = counter + 1;
		coordinates_.resize(nonzeros);
		values_.resize(nonzeros);

		for (SizeType i = 0; i < nonzeros; i++) {
			coordinates_[i] = newcoord[i];
			values_[i] = newvals[i];
		}
	}

	bool notEqual(const char* s) const
	{
		std::cerr << "notEqual=" << s << "\n";
		return false;
	}

	SizeType rows_;
	SizeType cols_;
	VectorComplexOrRealType values_;
	VectorPairType coordinates_;
	bool sorted_;
}; // VerySparseMatrix

template <typename T>
bool isHermitian(const VerySparseMatrix<T>& m)
{
	T eps = 1e-6;
	for (SizeType i = 0; i < m.nonZeros(); i++) {
		SizeType row = m.getRow(i);
		SizeType col = m.getColumn(i);
		if (fabs(m.getValue(i) - m(col, row)) > eps) {
			return false;
		}
	}
	return true;
}

template <typename T>
void verySparseMatrixToDenseMatrix(PsimagLite::Matrix<T>& m,
    const VerySparseMatrix<T>& vsm)
{
	m.resize(vsm.rows(), vsm.cols());
	m.setTo(0.0);
	SizeType nonzeroes = vsm.nonZeros();
	for (SizeType x = 0; x < nonzeroes; ++x)
		m(vsm.getRow(x), vsm.getColumn(x)) = vsm.getValue(x);
}

template <typename T>
void fullMatrixToVerySparseMatrix(VerySparseMatrix<T>& vsm,
    const PsimagLite::Matrix<T>& m)
{
	vsm.resize(m.rows(), m.cols());
	for (SizeType i = 0; i < m.rows(); ++i) {
		for (SizeType j = 0; j < m.cols(); ++j) {
			const T& val = m(i, j);
			if (val == 0.0)
				continue;
			vsm.push(i, j, val);
		}
	}
}
} // namespace Dmrg

namespace PsimagLite
{
template <typename T>
struct IsMatrixLike<Dmrg::VerySparseMatrix<T>> {
	enum { True = true };
};
}
/*@}*/
#endif
