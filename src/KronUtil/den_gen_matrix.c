#include "util.h"
void den_gen_matrix( const int nrow_A, 
                     const int ncol_A, 
                     const double threshold, 
                           double a_[] )
{
#define A(ia,ja) a_[ (ia) + (ja)*nrow_A ]
/* 
 * -------------------------------
 * generate a random matix in (0,1)
 * accept only if   aij < threshold
 * full matrix if threshold > 1
 * sparse matrix if threshold << 1
 * -------------------------------
 */

 const double dzero = 0;

 int ia = 0;
 int ja = 0;

 for(ja=0; ja < ncol_A; ja++) {
 for(ia=0; ia < nrow_A; ia++) {
    double drand = ((double) rand())/( (double) RAND_MAX );
    double aij   = ((double) rand())/( (double) RAND_MAX );

    int is_accept = (drand <= threshold);

    A(ia,ja) = (is_accept) ? aij : dzero;

    };
    };
}
#undef A
