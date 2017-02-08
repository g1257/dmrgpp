#ifndef BLAS_H
#define BLAS_H

#ifdef __cplusplus
extern "C" {
#endif

extern
void dgemm_( const char *trans1, const char *trans2,
             const int *mm, const int *nn, const int *kk,
             const double *alpha,   const double *a, const int *ld1,
                                    const double *b, const int *ld2,
             const double *beta,          double *c, const int *ld3 );


#ifdef __cplusplus
}
#endif

#endif
