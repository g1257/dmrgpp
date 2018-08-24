/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file BlockDiagonalMatrix.h
 *
 *  A class to represent a block diagonal matrix
 *
 */
#ifndef BLOCK_DIAGONAL_MATRIX_H
#define BLOCK_DIAGONAL_MATRIX_H
#include <vector>
#include <iostream>
#include "Matrix.h" // in PsimagLite
#include "ProgramGlobals.h"
#include "Concurrency.h"
#include "NoPthreadsNg.h"
#include "CrsMatrix.h"
#include "PsimagLite.h"
#include "EnforcePhase.h"
#include "Io/IoSelector.h"

namespace Dmrg {

template<typename T>
struct IsBasisType {
	enum {True = false};
};

// A block matrix class
// Blocks can be of any type and are templated with the type MatrixInBlockTemplate
template<typename MatrixInBlockTemplate>
class BlockDiagonalMatrix {

public:

	enum SaveEnum {SAVE_ALL, SAVE_PARTIAL};

	typedef MatrixInBlockTemplate BuildingBlockType;
	typedef typename BuildingBlockType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;
	typedef PsimagLite::IoSelector::In IoInType;

	BlockDiagonalMatrix()
	    : isSquare_(true)
	{}

	BlockDiagonalMatrix(IoInType& io, PsimagLite::String label)
	{
		io.read(isSquare_, label + "/isSquare_");
		io.read(offsetsRows_, label + "/offsetRows_");
		io.read(offsetsCols_, label + "/offsetCols_");
		io.read(data_, label + "/data_");
	}

	template<typename SomeBasisType>
	BlockDiagonalMatrix(const SomeBasisType& basis,
	                    typename PsimagLite::EnableIf<
	                    IsBasisType<SomeBasisType>::True, int>::Type = 0)
	    : isSquare_(true),
	      offsetsRows_(basis.partition()),
	      offsetsCols_(basis.partition()),
	      data_(basis.partition() - 1)
	{
		SizeType n = offsetsRows_.size();
		assert(n == offsetsCols_.size());
		for (SizeType i = 0; i < n; ++i)
			offsetsRows_[i] = offsetsCols_[i] = basis.partition(i);
	}

	template<typename IoOutputType>
	void write(PsimagLite::String label1, IoOutputType& io) const
	{
		io.createGroup(label1);
		io.write(isSquare_, label1 + "/isSquare_");
		io.write(offsetsRows_, label1 + "/offsetRows_");
		io.write(offsetsCols_, label1 + "/offsetCols_");
		io.write(data_, label1 + "/data_");
	}

	void setTo(ComplexOrRealType value)
	{
		SizeType n = data_.size();
		if (n == 0)
			err("BlockDiagonalMatrix::setTo(...): cannot be called without structure\n");

		for (SizeType i = 0; i < n; ++i)
			data_[i].setTo(value);
	}

	void operator+=(const BlockDiagonalMatrix& m)
	{
		mustBeSquare("operator+=");
		BlockDiagonalMatrix c;
		if (offsetsRows_.size() < m.blocks())
			operatorPlus(c,*this,m);
		else operatorPlus(c,m,*this);
		*this = c;
	}

	void setBlock(SizeType i,int offset,MatrixInBlockTemplate const &m)
	{
		mustBeSquare("setBlock");
		assert(i < data_.size());
		data_[i]=m;
		assert(i < offsetsRows_.size() && i < offsetsCols_.size());
		offsetsRows_[i] = offsetsCols_[i] = offset;
	}

	void sumBlock(SizeType i,MatrixInBlockTemplate const &m)
	{
		assert(i < data_.size());
		data_[i] += m;
	}

	void enforcePhase()
	{
		SizeType n = data_.size();
		for (SizeType i = 0; i < n; ++i)
			EnforcePhase<ComplexOrRealType>::enforcePhase(data_[i]);
	}

	// rows aren't affected, columns may be truncated
	void truncate(const VectorSizeType& removedIndices2)
	{
		if (removedIndices2.size() == 0) return;

		mustBeSquare("truncate");
		SizeType n = data_.size();
		VectorIntType remap(cols(), -1);
		computeRemap(remap, removedIndices2);
		VectorSizeType offsetsOld = offsetsCols_;
		for (SizeType i = 0; i < n; ++i)
			truncate(i, remap, offsetsOld);
		isSquare_ = (rows() == cols());
		assert(cols() == rows() - removedIndices2.size());
	}

	SizeType rows() const
	{
		SizeType n = offsetsRows_.size();
		return (n == 0) ? 0 : offsetsRows_[n - 1];
	}

	SizeType cols() const
	{
		SizeType n = offsetsCols_.size();
		return (n == 0) ? 0 : offsetsCols_[n - 1];
	}

	const VectorSizeType& offsetsRows() const { return offsetsRows_; }

	const VectorSizeType& offsetsCols() const { return offsetsCols_; }

	SizeType offsetsRows(SizeType i) const
	{
		assert(i < offsetsRows_.size());
		return offsetsRows_[i];
	}

	SizeType offsetsCols(SizeType i) const
	{
		assert(i < offsetsCols_.size());
		return offsetsCols_[i];
	}

	SizeType blocks() const { return data_.size(); }

	void toSparse(PsimagLite::CrsMatrix<ComplexOrRealType>& fm) const
	{
		SizeType r = rows();
		SizeType c = cols();
		fm.resize(r, c);
		SizeType counter = 0;
		SizeType k = 0;
		for (SizeType i = 0; i < r; ++i) {
			fm.setRow(i, counter);
			if (k + 1 < offsetsRows_.size() && offsetsRows_[k + 1] <= i)
				++k;
			SizeType end = (k + 1 < offsetsCols_.size()) ? offsetsCols_[k + 1] : c;
			if (data_[k].rows() == 0 || data_[k].cols() == 0) continue;
			for (SizeType j = offsetsCols_[k]; j < end; ++j) {
				ComplexOrRealType val = data_[k](i - offsetsRows_[k], j - offsetsCols_[k]);
				if (PsimagLite::norm(val) == 0)
					continue;
				fm.pushValue(val);
				fm.pushCol(j);
				++counter;
			}
		}

		fm.setRow(r, counter);
		fm.checkValidity();
	}

	void toDense(PsimagLite::Matrix<ComplexOrRealType>& fm) const
	{
		SizeType r = rows();
		SizeType c = cols();
		fm.resize(r, c);
		fm.setTo(static_cast<ComplexOrRealType>(0));
		SizeType k = 0;
		for (SizeType i = 0; i < r; ++i) {
			if (k + 1 < offsetsRows_.size() && offsetsRows_[k + 1] <= i)
				++k;
			SizeType end = (k + 1 < offsetsCols_.size()) ? offsetsCols_[k + 1] : c;
			if (data_[k].rows() == 0 || data_[k].cols() == 0) continue;
			for (SizeType j = offsetsCols_[k]; j < end; ++j) {
				ComplexOrRealType val = data_[k](i - offsetsRows_[k], j - offsetsCols_[k]);
				fm(i, j) = val;
			}
		}
	}

	void diagAndEnforcePhase(SizeType m, VectorRealType& eigsTmp, char option)
	{
		assert(m < data_.size());
		PsimagLite::diag(data_[m], eigsTmp, option);
		EnforcePhase<ComplexOrRealType>::enforcePhase(data_[m]);
	}

	const MatrixInBlockTemplate& operator()(SizeType i) const
	{
		assert(i < data_.size());
		return data_[i];
	}

	void clear()
	{
		offsetsRows_.clear();
		offsetsCols_.clear();
		data_.clear();
	}

	void read(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer)
	{
		ioSerializer.read(isSquare_, label + "/isSquare");
		ioSerializer.read(offsetsRows_, label + "/offsetsRows");
		ioSerializer.read(offsetsCols_, label + "/offsetsCols");
		ioSerializer.read(data_, label + "/data");
	}

	void write(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);
		ioSerializer.write(label + "/isSquare", isSquare_);
		ioSerializer.write(label + "/offsetsRows", offsetsRows_);
		ioSerializer.write(label + "/offsetsCols", offsetsCols_);
		ioSerializer.write(label + "/data", data_);
	}

	friend std::ostream& operator<<(std::ostream& os, const BlockDiagonalMatrix& m)
	{
		PsimagLite::String str = (m.isSquare_) ? "1" : "0";
		os<<str<<"\n";
		os<<m.offsetsRows_;
		os<<m.offsetsCols_;
		os<<m.data_.size()<<"\n";
		for (SizeType i = 0; i < m.data_.size(); ++i)
			os<<m.data_[i];
		return os;
	}

	friend std::istream& operator>>(std::istream& is, BlockDiagonalMatrix& m)
	{
		int x = -1;
		PsimagLite::String temp;
		is>>temp;
		if (temp == "NAME=")
			is>>x;
		else
			x = atoi(temp.c_str());

		if (x != 0 && x != 1)
			err("std::istream& operator>> BlockDiagonalMatrix(1)\n");
		m.isSquare_ = (x == 1);
		is>>m.offsetsRows_;
		is>>m.offsetsCols_;
		int total = 0;
		is>>total;
		if (total < 0)
			err("std::istream& operator>> BlockDiagonalMatrix(2)\n");
		if (total == 0)
			return is;
		m.data_.resize(total);
		for (SizeType i = 0; i < m.data_.size(); ++i)
			is>>m.data_[i];
		return is;
	}

private:

	void computeRemap(VectorIntType& remap,
	                  const VectorSizeType& removedIndices2) const
	{
		VectorSizeType removedIndices = removedIndices2;
		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(removedIndices.size(), 0);
		sort.sort(removedIndices, iperm);
		SizeType c = cols();
		SizeType k = 0;
		VectorSizeType::iterator b = removedIndices.begin();
		for (SizeType j = 0; j < c; ++j) {
			if (std::find(b, removedIndices.end(), j) != removedIndices.end()) {
				++b;
				continue;
			}

			remap[j] = k++;
		}
	}

	// rows aren't affected, columns may be truncated,
	void truncate(SizeType ind, // <------ block to truncate
	              const VectorIntType& remap,
	              const VectorSizeType& offsetsOld)
	{
		assert(ind < data_.size());
		MatrixInBlockTemplate& m = data_[ind];
		assert(ind < offsetsOld.size());
		SizeType offsetOld = offsetsOld[ind];
		assert(ind < offsetsCols_.size());
		SizeType offsetNew = offsetsCols_[ind];
		SizeType c = m.cols();
		SizeType counter = 0;
		for (SizeType j = 0; j < c; ++j) {
			assert(j + offsetOld < remap.size());
			if (remap[j + offsetOld] >= 0)
				continue;
			++counter;
		}

		if (counter == 0) {
			offsetsCols_[ind + 1] = offsetNew + c;
			return;
		}

		if (counter == c) {
			deleteThisColBlock(ind);
			return;
		}

		assert(counter < c);

		SizeType r = m.rows();
		SizeType newCols = c - counter;
		assert(newCols < c);
		MatrixInBlockTemplate m2(r, newCols);
		for (SizeType i = 0; i < r; ++i) {
			for (SizeType j = 0; j < c; ++j) {
				if (remap[j + offsetOld] < 0) continue;
				SizeType cPrime = remap[j + offsetOld];
				assert(offsetNew <= cPrime);
				m2(i, cPrime - offsetNew) = m(i, j);
			}
		}

		m = m2;
		assert(ind + 1 < offsetsCols_.size());
		offsetsCols_[ind + 1] = offsetNew + newCols;
	}

	void deleteThisColBlock(SizeType ind)
	{
		assert(ind < data_.size());
		data_[ind].clear();
		assert(ind + 1 < offsetsCols_.size());
		offsetsCols_[ind + 1] = offsetsCols_[ind];
	}

	void mustBeSquare(PsimagLite::String msg) const
	{
		if (isSquare_) return;
		err("BlockDiagonalMatrix::" + msg + " must be square\n");
	}

	bool isSquare_;
	VectorSizeType offsetsRows_;
	VectorSizeType offsetsCols_;
	typename PsimagLite::Vector<MatrixInBlockTemplate>::Type data_;
}; // class BlockDiagonalMatrix

// Companion Functions

template<class MatrixInBlockTemplate>
bool isUnitary(const BlockDiagonalMatrix<MatrixInBlockTemplate>& B)
{
	bool flag=true;
	MatrixInBlockTemplate matrixTmp;

	for (SizeType m=0;m<B.blocks();m++) {
		matrixTmp = B(m);
		if (!isUnitary(matrixTmp)) flag=false;
	}
	return flag;
}

} // namespace Dmrg

namespace PsimagLite {
template<typename T>
struct IsMatrixLike<Dmrg::BlockDiagonalMatrix<T> > {
	enum {True=true};
};
} // namespace PsimagLite
/*@}*/

#endif

