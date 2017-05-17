#include "util.h"

template<typename ComplexOrRealType>
void csr_eye(const int nrow_B,
             const int ncol_B,
             PsimagLite::CrsMatrix<ComplexOrRealType>& b)
{
	/*
 * ---------------------------------------------------------------------------
 * Out:     sparse identity matrix in compressed sparse row format
 * ---------------------------------------------------------------------------
 */


	/*
  * ----------------------------------------------------
  * setup boolean array for fast mapping of column index
  * ----------------------------------------------------
  */

	int* nnz = new int[nrow_B];

	int ib = 0;



	/*
   * -------------------------------------------------------
   * first pass to setup number of nonzeros per row in B
   * -------------------------------------------------------
   */

	for(ib=0; ib < nrow_B; ib++) {
		nnz[ib] = 1;
	};


	/*
  * ---------------------------------
  * prefix sum to setup row pointers
  * ---------------------------------
  */

	int max_nnz = (nrow_B >= ncol_B) ? nrow_B : ncol_B;

	b.resize(nrow_B, ncol_B, max_nnz);
	b.setRow(0,0);
	for(ib=0; ib < nrow_B; ib++)
		b.setRow(ib+1,b.getRowPtr(ib) + nnz[ib]);

	/*
   * ------------------------
   * reset array for 2nd pass
   * ------------------------
   */
	for(ib=0; ib < nrow_B; ib++) {
		nnz[ib] = 0;
	};


	/*
   * --------------------------------------
   * second pass to fill in compress sparse row
   * data structure
   * --------------------------------------
   */
	for(ib=0; ib < nrow_B; ib++) {
		int jb = ib;
		int isvalid =  (0 <= jb) && (jb < ncol_B);
		if (isvalid) {
			int ipos = b.getRowPtr(ib) + nnz[ib];

			b.setValues(ipos,1);
			b.setCol(ipos,jb);
			nnz[ib] += 1;
		};
	};

	b.checkValidity();
	delete[] nnz;
}
