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

namespace Dmrg {

// A block matrix class
// Blocks can be of any type and are templated with the type MatrixInBlockTemplate
//
// Note: In reality, Parallelization is disabled here because a LAPACK call
//        is needed and LAPACK is not necessarily thread safe.
template<typename MatrixInBlockTemplate>
class BlockDiagonalMatrix {

public:

	typedef MatrixInBlockTemplate BuildingBlockType;
	typedef typename BuildingBlockType::value_type FieldType;
	typedef BlockDiagonalMatrix<MatrixInBlockTemplate> BlockDiagonalMatrixType;
	typedef typename PsimagLite::Real<FieldType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	class LoopForDiag {

		typedef PsimagLite::Concurrency ConcurrencyType;

	public:

		LoopForDiag(BlockDiagonalMatrixType& C1,
		            VectorRealType& eigs1,
		            char option1)
		    : C(C1),
		      eigs(eigs1),
		      option(option1),
		      eigsForGather(C.blocks()),
		      weights(C.blocks())
		{

			for (SizeType m=0;m<C.blocks();m++) {
				eigsForGather[m].resize(C.offsets(m+1)-C.offsets(m));
				weights[m] =  C.offsets(m+1)-C.offsets(m);
			}

			eigs.resize(C.rank());
		}

		SizeType tasks() const { return C.blocks(); }

		void doTask(SizeType taskNumber, SizeType)
		{
			SizeType m = taskNumber;
			VectorRealType eigsTmp;
			PsimagLite::diag(C.data_[m],eigsTmp,option);
			enforcePhase(C.data_[m]);
			for (int j=C.offsets(m);j< C.offsets(m+1);j++)
				eigsForGather[m][j-C.offsets(m)] = eigsTmp[j-C.offsets(m)];

		}

		void gather()
		{
			for (SizeType m=0;m<C.blocks();m++) {
				for (int j=C.offsets(m);j< C.offsets(m+1);j++)
					eigs[j]=eigsForGather[m][j-C.offsets(m)];
			}
		}


		static void enforcePhase(PsimagLite::Matrix<FieldType>& a)
		{
			SizeType cols = a.cols();
			for (SizeType i = 0; i < cols; ++i) {
				enforcePhase(a, i);
			}
		}

	private:

		static void enforcePhase(PsimagLite::Matrix<FieldType>& a, SizeType col)
		{
			RealType sign1 = 0;
			SizeType rows = a.rows();
			for (SizeType j = 0; j < rows; ++j) {
				RealType val = PsimagLite::norm(a(j, col));
				if (val < 1e-6) continue;
				sign1 = (val > 0) ? 1 : -1;
				break;
			}

			assert(sign1 != 0);
			// get a consistent phase
			for (SizeType j = 0; j < rows; ++j)
				a(j, col) *= sign1;
		}

		BlockDiagonalMatrixType& C;
		VectorRealType& eigs;
		char option;
		typename PsimagLite::Vector<VectorRealType>::Type eigsForGather;
		typename PsimagLite::Vector<SizeType>::Type weights;
	};

	BlockDiagonalMatrix(int rank,int blocks) : rank_(rank),offsets_(blocks+1),data_(blocks)
	{
		offsets_[blocks]=rank;
	}

	template<typename SomeBasisType>
	BlockDiagonalMatrix(const SomeBasisType& basis)
	    : rank_(basis.size()), offsets_(basis.partition()), data_(offsets_.size() - 1)
	{
		SizeType n = offsets_.size();
		for (SizeType i = 0; i < n; ++i)
			offsets_[i] = basis.partition(i);
	}

	void setTo(FieldType value)
	{
		SizeType n = data_.size();
		if (n == 0)
			err("BlockDiagonalMatrix::setTo(...): cannot be called without structure\n");

		for (SizeType i = 0; i < n; ++i)
			data_[i].setTo(value);
	}

	void operator+=(const BlockDiagonalMatrixType& m)
	{
		BlockDiagonalMatrixType c;
		if (offsets_.size()<m.blocks()) operatorPlus(c,*this,m);
		else operatorPlus(c,m,*this);
		*this = c;
	}

	void setBlock(SizeType i,int offset,MatrixInBlockTemplate const &m)
	{
		assert(i<data_.size());
		data_[i]=m;
		assert(i<offsets_.size());
		offsets_[i]=offset;
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
			LoopForDiag::enforcePhase(data_[i]);
	}

	int rank() const { return rank_; }

	int offsets(int i) const { return offsets_[i]; }

	SizeType blocks() const { return data_.size(); }

	void toSparse(PsimagLite::CrsMatrix<FieldType>& fm) const
	{
		fm.resize(rank_, rank_);
		SizeType counter=0;
		SizeType k = 0;
		for (int i = 0; i < rank_; ++i) {
			fm.setRow(i,counter);
			if (k+1 < offsets_.size() && offsets_[k+1] <= i)
				++k;
			SizeType end = (k + 1 < offsets_.size()) ? offsets_[k + 1] : rank_;
			if (data_[k].rows() == 0 || data_[k].cols() == 0) continue;
			for (SizeType j = offsets_[k]; j < end; ++j) {
				FieldType val = data_[k](i-offsets_[k],j-offsets_[k]);
				if (PsimagLite::norm(val) == 0)
					continue;
				fm.pushValue(val);
				fm.pushCol(j);
				counter++;
			}
		}

		fm.setRow(rank_, counter);
		fm.checkValidity();
	}

	const MatrixInBlockTemplate& operator()(SizeType i) const
	{
		assert(i < data_.size());
		return data_[i];
	}

	template<class MatrixInBlockTemplate2>
	friend void operatorPlus(BlockDiagonalMatrix<MatrixInBlockTemplate2>& C,
	                         const BlockDiagonalMatrix<MatrixInBlockTemplate2>& A,
	                         const BlockDiagonalMatrix<MatrixInBlockTemplate2>& B)
	{
		SizeType i;
		int counter=0;
		C.rank_ = A.rank_;
		C.offsets_=A.offsets_;
		C.data_.resize(A.data_.size());
		for (i=0;i<A.data_.size();i++) {
			C.data_[i] = A.data_[i];

			while (B.offsets_[counter] <A.offsets_[i+1]) {
				accumulate(C.data_[i],B.data_[counter++]);
				if (counter>=B.offsets_.size()) break;
			}

			if (counter>=B.offsets_.size() && i<A.offsets_.size()-1)
				throw PsimagLite::RuntimeError("operatorPlus: restriction not met.\n");
		}
	}

	template<class MatrixInBlockTemplate2>
	friend void operatorPlus(BlockDiagonalMatrix<MatrixInBlockTemplate2>& A,
	                         const BlockDiagonalMatrix<MatrixInBlockTemplate2>& B)
	{
		SizeType i;
		int counter=0;

		for (i=0;i<A.data_.size();i++) {
			while (B.offsets_[counter] < A.offsets_[i+1]) {
				accumulate(A.data_[i],B.data_[counter++]);
				if (counter>=B.offsets_.size()) break;
			}

			if (counter>=B.offsets_.size() && i<A.offsets_.size()-1)
				throw PsimagLite::RuntimeError("operatorPlus: restriction not met.\n");
		}
	}

	template<class MatrixInBlockTemplate2>
	friend std::ostream &operator<<(std::ostream &s,
	                                const BlockDiagonalMatrix<MatrixInBlockTemplate2>& A)
	{
		for (SizeType m=0;m<A.blocks();m++) {
			int nrank = A.offsets(m+1)-A.offsets(m);
			s<<"block number "<<m<<" has rank "<<nrank<<"\n";
			s<<A.data_[m];
		}

		return s;
	}

private:

	int rank_; //the rank of this matrix
	typename PsimagLite::Vector<int>::Type offsets_; //starting of diagonal offsets for each block
	typename PsimagLite::Vector<MatrixInBlockTemplate>::Type data_; // data on each block

}; // class BlockDiagonalMatrix

// Companion Functions
// Parallel version of the diagonalization of a block diagonal matrix
template<typename SomeVectorType,typename SomeFieldType>
typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,
void>::Type
diagonalise(BlockDiagonalMatrix<PsimagLite::Matrix<SomeFieldType> >& C,
            SomeVectorType& eigs,
            char option)
{
	typedef typename BlockDiagonalMatrix<PsimagLite::Matrix<SomeFieldType> >::LoopForDiag LoopForDiagType;
	typedef PsimagLite::NoPthreadsNg<LoopForDiagType> ParallelizerType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	SizeType savedNpthreads = ConcurrencyType::npthreads;
	ConcurrencyType::npthreads = 1;
	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD,
	                              false);

	LoopForDiagType helper(C,eigs,option);

	threadObject.loopCreate(helper); // FIXME: needs weights

	helper.gather();

	ConcurrencyType::npthreads = savedNpthreads;
}

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
/*@}*/

#endif

