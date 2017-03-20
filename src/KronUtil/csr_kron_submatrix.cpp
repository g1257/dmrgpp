#include "util.h"
void csr_kron_submatrix( 
         const PsimagLite::CrsMatrix<double>& a,

      const PsimagLite::CrsMatrix<double>& b,
         
         const int nrindex, 
         const int ncindex, 
         const int max_nnz,
         const PsimagLite::Vector<int>::Type& rindex,
         const PsimagLite::Vector<int>::Type& cindex,

         PsimagLite::CrsMatrix<double>& h)
{
/*
 * -------------------------------------------------
 * extract a submatrix out of kronecker product
 * equivalent to 
 * C = kron(A,B), then H = C( rindex(:), cindex(:) )
 *
 * assume A, B are in sparse compressed ROW format
 * -------------------------------------------------
 */

	const int ncol_A = a.col();
	const int nrow_B = b.row();
	const int ncol_B = b.col();

    const int ncol_C = ncol_A * ncol_B;
#ifndef NDEBUG
	const int nrow_A = a.row();
    const int nrow_C = nrow_A * nrow_B;
#endif
    const int nrow_H = nrindex;
    const int ncol_H = ncindex;



 /*
  * -----------------------------
  * split cindex(:) into [jb,ja]
  * -----------------------------
  */
  int* ialist = new int[nrindex];
  int* iblist = new int[nrindex];

  int k = 0;
  for(k=0; k < nrindex; k++) {
     int ic = rindex[k];
     int ib = MOD( ic, nrow_B );
     int ia = (ic - ib)/nrow_B;

     assert((0 <= ia) && (ia < nrow_A));
     assert((0 <= ib) && (ib < nrow_B));
     assert((0 <= ic) && (ic < nrow_C));

     ialist[k] = ia;
     iblist[k] = ib;
     };

 /*
  * ------------------------------
  * setup mapping for column index
  * ------------------------------
  */
  int* cmap = new int[ncol_C];
  int jc = 0;
  for(jc=0; jc < ncol_C; jc++) {
      cmap[jc] = -1;
      };

  for(k=0; k < ncindex; k++) {
      int jc = cindex[k];

      assert((0 <= jc) && (jc < ncol_C));

      cmap[jc] = k;
      };


  int ih = 0;
  int ifree = 0;
  h.resize(nrow_H, ncol_H, max_nnz);
  for(ih=0; ih < nrow_H; ih++) {
     h.setRow(ih,ifree);

     int ia = ialist[ih];
     int ib = iblist[ih];

     int istarta = a.getRowPtr(ia);
     int ienda = a.getRowPtr(ia + 1);

     int istartb = b.getRowPtr(ib);
     int iendb = b.getRowPtr(ib + 1);

     int ka = 0;
     int kb = 0;
     for(ka=istarta; ka < ienda; ka++) {
     for(kb=istartb; kb < iendb; kb++) {
         int ja = a.getCol(ka);
         int jb = b.getCol(kb);
         double aij = a.getValue(ka);
         double bij = b.getValue(kb);

         int jc = jb + ja*ncol_B;
     

         
         int jh = cmap[jc];
         int isvalid = (0 <= jh) && (jh < ncol_H);
         if (isvalid) {
            double cij = aij * bij;

            assert((ifree < max_nnz));

            h.setCol(ifree, jh);
            h.setValues(ifree,cij);
            ifree++;
            };
         };
         };
     };
  h.setRow(nrow_H, ifree);
h.checkValidity();

  delete[] ialist;
  delete[] iblist;
  delete[] cmap;
}
