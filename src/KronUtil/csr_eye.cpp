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

        const int nnz_B = std::min( nrow_B, ncol_B );
        b.resize( nrow_B, ncol_B );
        b.reserve( nnz_B );


        int ip = 0;
        for( int irow=0; irow < nrow_B; irow++) {
                b.setRow(irow, ip);
                const int jcol = irow;
                const bool is_valid = (jcol < ncol_B);
                if (is_valid) {
                        b.pushValue( 1 );
                        b.pushCol( jcol );
                        ip++;
                };
        };

        b.setRow(nrow_B, ip);



	b.checkValidity();
}
