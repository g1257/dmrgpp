#include "dmrg_types.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "dmrg_vbatch.h"
#include "setup_nC.h"

template <typename T>
void setup_nC(SizeType noperator,
              SizeType npatches,
              /*
               * ---------
               * Intent(in)
               * ---------
               */
              const std::vector<T*>& Amatrix_,
              const VectorSizeType&  ld_Amatrix_,
              const std::vector<T*>& Bmatrix_,
              const VectorSizeType&  ld_Bmatrix_,
              const VectorSizeType&  left_patch_size_, /* sizes of left patches (INPUT) */
              const VectorSizeType&  right_patch_size_, /* sizes of right patches (INPUT) */
              /*
               * --------------------------------
               * Intent(out) but assume storage already allocated
               * --------------------------------
               */
              VectorSizeType& nC_,
              VectorSizeType& gnnz_A_,
              VectorSizeType& gnnz_B_)
{
	const IntegerType idebug = 1;

	SizeType ipatch    = 0;
	SizeType jpatch    = 0;
	SizeType ioperator = 0;

	for (jpatch = 1; jpatch <= npatches; jpatch++) {
		for (ipatch = 1; ipatch <= npatches; ipatch++) {
			nC_[indx2f(ipatch, jpatch, npatches)] = 0;
		};
	};

	for (ioperator = 1; ioperator <= noperator; ioperator++) {
		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			for (ipatch = 1; ipatch <= npatches; ipatch++) {

				T* Amat_ = Amatrix_[index3(ipatch, jpatch, ioperator, npatches)];
				T* Bmat_ = Bmatrix_[index3(ipatch, jpatch, ioperator, npatches)];
				SizeType ld_Amat
				    = ld_Amatrix_[index3(ipatch, jpatch, ioperator, npatches)];
				SizeType ld_Bmat
				    = ld_Bmatrix_[index3(ipatch, jpatch, ioperator, npatches)];

				SizeType nrowA = left_patch_size_[ipatch - 1];
				SizeType ncolA = left_patch_size_[jpatch - 1];

				SizeType nrowB = right_patch_size_[ipatch - 1];
				SizeType ncolB = right_patch_size_[jpatch - 1];

				SizeType nnz_Amat = 0;
				SizeType nnz_Bmat = 0;

				/*
	   --------------------------------------------------------------
	   Note: treat a NULL poIntegerTypeer as poIntegerTypeer to a matrix of all zeros
	   --------------------------------------------------------------
	   */
				if (Amat_ == NULL) {
					nnz_Amat = 0;
				} else {
					SizeType ia = 0;
					SizeType ja = 0;
					for (ja = 1; ja <= ncolA; ja++) {
						for (ia = 1; ia <= nrowA; ia++) {

							nnz_Amat += ((Amat_[indx2f(ia, ja, ld_Amat)]
							              == T(0))
							                 ? 0
							                 : 1);
						};
					};
				};

				if (Bmat_ == NULL) {
					nnz_Bmat = 0;
				} else {
					SizeType ib = 0;
					SizeType jb = 0;
					for (jb = 1; jb <= ncolB; jb++) {
						for (ib = 1; ib <= nrowB; ib++) {
							nnz_Bmat += ((Bmat_[indx2f(ib, jb, ld_Bmat)]
							              == T(0))
							                 ? 0
							                 : 1);
						};
					};
				};
				gnnz_A_[index3(ipatch, jpatch, ioperator, npatches)] = nnz_Amat;
				gnnz_B_[index3(ipatch, jpatch, ioperator, npatches)] = nnz_Bmat;

				SizeType is_zero_Amat = (nnz_Amat == 0);
				SizeType is_zero_Bmat = (nnz_Bmat == 0);
				if (is_zero_Amat || is_zero_Bmat) {
					/*
		  ---------------------------------------------------------
		  ignore this operator since kron(Amat,Bmat) is zero matrix
		  ---------------------------------------------------------
		  */
				} else {
					/*
		  ------------------------------------------
		  pair of non-zero Amat(k), Bmat(k) matrices
		  ------------------------------------------
		  */
					nC_[indx2f(ipatch, jpatch, npatches)]
					    = nC_[indx2f(ipatch, jpatch, npatches)] + 1;
				};
			}; /* end for jpatch */
		}; /* end for ipatch */
	}; /* end for ioperator */

	if (idebug >= 1) {
		/*
		 * ----------------------------------
		 * count total number of IntegerTypeeractions
		 * ----------------------------------
		 */
		SizeType total_number = 0;
		SizeType ipatch       = 0;
		SizeType jpatch       = 0;
		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			for (ipatch = 1; ipatch <= npatches; ipatch++) {
				total_number += nC_[indx2f(ipatch, jpatch, npatches)];
			};
		};
		printf("setup_nC:npatches=%lu,noperator=%lu,total_number=%lu\n",
		       npatches,
		       noperator,
		       total_number);
	};
}

template void setup_nC<MYTYPE>(SizeType,
                               SizeType,
                               /*
                                * ---------
                                * Intent(in)
                                * ---------
                                */
                               const std::vector<MYTYPE*>&,
                               const VectorSizeType&,
                               const std::vector<MYTYPE*>&,
                               const VectorSizeType&,
                               const VectorSizeType&,
                               const VectorSizeType&,
                               /*
                                * --------------------------------
                                * Intent(out) but assume storage already allocated
                                * --------------------------------
                                */
                               VectorSizeType&,
                               VectorSizeType&,
                               VectorSizeType&);
