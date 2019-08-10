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

/*! \file ArrayOfMatStruct.h
 *
 *
 */

#ifndef ARRAY_OF_MAT_STRUCT_H
#define ARRAY_OF_MAT_STRUCT_H
#include "GenIjPatch.h"
#include "CrsMatrix.h"
#include "../KronUtil/MatrixDenseOrSparse.h"
#include "Profiling.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class ArrayOfMatStruct {

public:

	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorStorageType OperatorStorageType;
	typedef MatrixDenseOrSparse<SparseMatrixType> MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::RealType RealType;
	typedef GenIjPatch<LeftRightSuperType> GenIjPatchType;
	typedef typename GenIjPatchType::VectorSizeType VectorSizeType;
	typedef typename GenIjPatchType::BasisType BasisType;
	typedef typename MatrixDenseOrSparseType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	ArrayOfMatStruct(const OperatorStorageType& sparse1,
	                 const GenIjPatchType& patchOld,
	                 const GenIjPatchType& patchNew,
	                 typename GenIjPatchType::LeftOrRightEnumType leftOrRight,
	                 RealType threshold,
	                 bool useLowerPart)
	    : data_(patchNew(leftOrRight).size(), patchOld(leftOrRight).size())
	{
		const SparseMatrixType& sparse = sparse1.getCRS();
		const BasisType& basisOld = (leftOrRight == GenIjPatchType::LEFT) ?
		            patchOld.lrs().left() : patchOld.lrs().right();
		const BasisType& basisNew = (leftOrRight == GenIjPatchType::LEFT) ?
		            patchNew.lrs().left() : patchNew.lrs().right();
		const SizeType npatchOld = patchOld(leftOrRight).size();
		const SizeType npatchNew = patchNew(leftOrRight).size();

		const SizeType ipatchSize = npatchNew;
		const SizeType jpatchSize = npatchOld;

#ifdef NDEBUG
		const int idebug = 0;
#else
		// ---------------------------
		// turn on more debug checking
		// ---------------------------
		const int idebug = 1;
#endif



		// ---------------------------------
		// precompute the size of each patch
		// ---------------------------------

		std::vector<SizeType> ipatch_Size( ipatchSize );
		std::vector<SizeType> jpatch_Size( jpatchSize );

		// ----------------------------------
		// setup mapping from index to ipatch
		// ----------------------------------
		const SizeType ncols = sparse.cols();

		const SizeType invalid_jpatch_number = jpatchSize + 1;
		// ----------------------------------------
		// initiall fill array with invalid values
		// ----------------------------------------
		std::vector<SizeType> index_to_jpatch(ncols,  invalid_jpatch_number);

		std::vector<SizeType> offset_ipatch(ipatchSize);
		std::vector<SizeType> offset_jpatch(jpatchSize);

		for(SizeType ipatch=0; ipatch < ipatchSize; ++ipatch) {
			const SizeType igroup = patchNew(leftOrRight)[ipatch];
			const SizeType i1 = basisNew.partition(igroup);
			const SizeType i2 = basisNew.partition(igroup+1);
			assert( (0 <= i1) && (i1 <= i2) &&
			        (i2 <= sparse.rows()));


			offset_ipatch[ ipatch ] = i1;

			const SizeType isize = i2 - i1;

			ipatch_Size[ ipatch ] = isize;
			assert( ipatch_Size[ ipatch ] >= 1 );

		};

		for(SizeType jpatch=0; jpatch < jpatchSize; ++jpatch) {
			const SizeType jgroup = patchOld(leftOrRight)[jpatch];
			const SizeType j1 = basisOld.partition(jgroup);
			const SizeType j2 = basisOld.partition(jgroup+1);

			assert( (0 <= j1) && (j1 <= j2) &&
			        (j2 <= ncols) );

			offset_jpatch[ jpatch ] = j1;

			const SizeType jsize = j2 - j1;

			jpatch_Size[ jpatch ] = jsize;
			assert( jpatch_Size[ jpatch ] >= 1 );

			for(SizeType j=j1; j < j2; j++) {
				index_to_jpatch[ j ] = jpatch;
			};
		};


		// ------------
		// double check
		// ------------
		if (idebug >= 1) {

			SizeType sum_ipatch_Size = 0;
			for(SizeType ipatch=0; ipatch < ipatchSize; ipatch++) {
				sum_ipatch_Size += ipatch_Size[ ipatch ];
			};
			assert( (0 <= sum_ipatch_Size) && (sum_ipatch_Size <= sparse.rows()) );

			SizeType sum_jpatch_Size = 0;
			for(SizeType jpatch=0; jpatch < jpatchSize; jpatch++) {
				sum_jpatch_Size += jpatch_Size[ jpatch ];
			};
			assert( (0 <= sum_jpatch_Size) && (sum_jpatch_Size <= ncols) );

		};

		// -----------------------------------------
		// allocate data structure outside main loop
		// to avoid repeated allocation and deallocation
		// -----------------------------------------
		SizeType max_ipatchSize = 0;
		for(SizeType ipatch=0; ipatch < ipatchSize; ipatch++) {
			max_ipatchSize = (ipatch_Size[ ipatch ] > max_ipatchSize) ?
			            ipatch_Size[ ipatch ] :
			            max_ipatchSize;
		};

		std::vector<   std::vector<SizeType> > rowPtr1D(jpatchSize);
		for(SizeType jpatch=0; jpatch < jpatchSize; jpatch++) {
			rowPtr1D[jpatch].reserve( max_ipatchSize );
		};

		std::vector< SizeType > total_nz(jpatchSize,0);
		std::vector<bool> is_dense1D( jpatchSize, false );

		for(SizeType ipatch=0; ipatch < ipatchSize; ipatch++) {

			// ------------------------------------------------------
			// initialize  data structure to count number of nonzeros
			// per row in sparse matrix of  data_(ipatch,jpatch)
			// ------------------------------------------------------
			for(SizeType jpatch=0; jpatch < jpatchSize; jpatch++) {
				is_dense1D[ jpatch ] = false;
			};

			for(SizeType jpatch=0; jpatch < jpatchSize; jpatch++) {
				total_nz[ jpatch ] = 0;
			};


			const SizeType local_nrows = ipatch_Size[ipatch];
			for(SizeType jpatch=0; jpatch < jpatchSize; jpatch++) {
				rowPtr1D[ jpatch ].clear();
				rowPtr1D[ jpatch ].resize( local_nrows,0 );
			};

			// --------------------------------------
			// first pass to count number of nonzeros
			// --------------------------------------
			SizeType i1 = offset_ipatch[ ipatch ];
			SizeType i2 = i1 +  ipatch_Size[ipatch];

			for(SizeType irow=i1; irow < i2; irow++) {
				const SizeType istart = sparse.getRowPtr(irow);
				const SizeType iend = sparse.getRowPtr(irow+1);


				for(SizeType k=istart; k < iend; k++) {
					const SizeType jcol = sparse.getCol(k);
					const SizeType jpatch = index_to_jpatch[ jcol ];

					bool is_valid_jpatch = (0 <= jpatch) && (jpatch < jpatchSize);
					if (!is_valid_jpatch) continue;


					if (useLowerPart && (ipatch < jpatch))  continue;

					const SizeType indx = jpatch;

					const SizeType i1 = offset_ipatch[ ipatch ];

					const SizeType  local_irow = (irow - i1 );
					assert( (0 <= local_irow) && (local_irow < ipatch_Size[ ipatch ]) );

					(rowPtr1D[indx])[local_irow]++;
					total_nz[ indx ]++;
				};
			};

			// ---------------------
			// setup sparse matrices
			// ---------------------
			// ------------------------------------------------------
			// This option tries to directly use data_(ipatch,jpatch)
			// but requires some changes to the interface to explicitly
			// construct a fully dense or sparse matrix and extract back out
			// the dense matrix or sparse matrix
			// ------------------------------------------------------



			for(SizeType jpatch=0; jpatch < jpatchSize; ++jpatch) {

				data_(ipatch,jpatch) = 0;
				if (useLowerPart && (ipatch < jpatch))  continue;

				const SizeType lnrows = ipatch_Size[ipatch];
				const SizeType lncols = jpatch_Size[jpatch];

				const SizeType indx = jpatch;
				const SizeType nnz = total_nz[ indx ];

				const bool isDense = (nnz   >= threshold * lnrows * lncols);
				assert( nnz <= lnrows * lncols );

				is_dense1D[ indx ] = isDense;

				if (nnz == 0) continue; // <--- ATTENTION: EARLY EXIT

				data_(ipatch, jpatch) = new MatrixDenseOrSparseType(lnrows,
				                                                    lncols,
				                                                    isDense);

				MatrixDenseOrSparseType* pmat = data_(ipatch, jpatch);

				if (isDense) {
					// ---------------------------
					// store as fully dense matrix
					// ---------------------------
					MatrixType& dense_mat = pmat->getDense();
					assert( dense_mat.rows() == lnrows );
					assert( dense_mat.cols() == lncols );

					bool const need_zero_out = true;
					if (need_zero_out) {

						for( SizeType j=0; j < lncols; j++) {
							for( SizeType i=0; i < lnrows; i++) {
								dense_mat(i,j) = 0;
							};
						};
					};
				}
				else {
					// ----------------------
					// store as sparse matrix
					// ----------------------
					SparseMatrixType& sparse_mat = pmat->getSparse();
					assert( sparse_mat.rows() == lnrows );
					assert( sparse_mat.cols() == lncols );

					// ------------------------------
					// preallocate sufficient storage
					// ------------------------------
					sparse_mat.resize( lnrows, lncols, nnz );

					SizeType ip = 0;
					for(SizeType irow = 0; irow < lnrows; irow++) {
						SizeType icount = (rowPtr1D[ indx ])[irow];
						(rowPtr1D[ indx ])[irow] = ip;

						sparse_mat.setRow( irow, ip );

						ip += icount;
					};
					sparse_mat.setRow( lnrows, ip );
				};

			}; // for jpatch

			// ---------------------------------------
			// second pass to fill in numerical values
			// ---------------------------------------

			for(SizeType irow=i1; irow < i2; irow++) {
				const SizeType istart = sparse.getRowPtr(irow);
				const SizeType iend = sparse.getRowPtr(irow+1);


				for(SizeType k=istart; k < iend; k++) {
					const SizeType jcol = sparse.getCol(k);
					const SizeType jpatch = index_to_jpatch[ jcol ];

					const bool is_valid_jpatch = (0 <= jpatch) && (jpatch < jpatchSize);
					if (!is_valid_jpatch) continue;

					if (useLowerPart && (ipatch < jpatch))  continue;

					if (!data_(ipatch,jpatch)) continue; // <--- ATTENTION: EARLY EXIT

					const SizeType indx = jpatch;

					const SizeType i1 = offset_ipatch[ ipatch ];
					const SizeType local_irow = (irow - i1 );

					const SizeType j1 = offset_jpatch[ jpatch ];
					const SizeType local_jcol = (jcol - j1);

					const ComplexOrRealType aij = sparse.getValue(k);

					if (is_dense1D[ indx ]) {

						MatrixType& dense_mat = data_(ipatch,jpatch)->getDense();
						dense_mat(local_irow,local_jcol) = aij;
					}
					else {
						SparseMatrixType& sparse_mat = data_(ipatch,jpatch)->getSparse();

						const SizeType ip = (rowPtr1D[ indx ])[ local_irow ];

						sparse_mat.setCol( ip, local_jcol );
						sparse_mat.setValues(ip, aij );

						(rowPtr1D[ indx ])[ local_irow ]++;
					};
				};
			}; // for irow

			// --------------------------------
			// check data_(ipatch,jpatch)
			// --------------------------------
			if (idebug >= 1) {
				for(SizeType jpatch=0; jpatch < jpatchSize; ++jpatch) {

					if (useLowerPart && (ipatch < jpatch))  continue;

					SizeType indx = jpatch;

					if (is_dense1D[ indx ]) continue;

					if (data_(ipatch,jpatch))
						(data_(ipatch,jpatch)->getSparse()).checkValidity();
				};
			};
		}; // for ipatch
	}

	const MatrixDenseOrSparseType* operator()(SizeType i, SizeType j) const
	{
		assert(i<data_.n_row() && j<data_.n_col());
		return data_(i, j);
	}

	~ArrayOfMatStruct()
	{
		const SizeType cols = data_.cols();
		const SizeType rows = data_.rows();
		for (SizeType j = 0; j < cols; ++j)
			for (SizeType i = 0; i < rows; ++i)
				delete data_(i, j);
	}

private:

	ArrayOfMatStruct(const ArrayOfMatStruct&) = delete;

	ArrayOfMatStruct& operator=(const ArrayOfMatStruct&) = delete;

	PsimagLite::Matrix<MatrixDenseOrSparseType*> data_;
}; //class ArrayOfMatStruct
} // namespace Dmrg

/*@}*/

#endif // ARRAY_OF_MAT_STRUCT_H

