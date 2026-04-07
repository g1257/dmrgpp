#include "test_vbatch.h"

template<typename T>
void get_total_memory( SizeType noperator,
                       SizeType npatches,
                       const std::vector<T*>& Amatrix_,
                       const VectorSizeType& ld_Amatrix_,
                       const std::vector<T*>& Bmatrix_,
                       const VectorSizeType& ld_Bmatrix_,
                       const VectorSizeType& left_patch_size_,
                       const VectorSizeType& right_patch_size_,
                       SizeType& ptotal_memory_in_nbytes)
{
        double total_gflops = 0;
        double gmemA = 0;
        double gmemB = 0;
        double gmemBX = 0;
        double gmemXY = 0;

	VectorSizeType nC_(npatches*npatches,0);
        VectorSizeType gnnz_A_(npatches*npatches*noperator,0);
        VectorSizeType gnnz_B_(npatches*npatches*noperator,0);
	
        setup_nC<T>(noperator,
                  npatches,
                  Amatrix_, 
		  ld_Amatrix_,
                  Bmatrix_, 
		  ld_Bmatrix_,
                  left_patch_size_,
                  right_patch_size_,
                  nC_,
                  gnnz_A_,
                  gnnz_B_
                );

        estimate_work( npatches,
                       left_patch_size_,
                       right_patch_size_,
                       nC_,
                       &total_gflops,
                       &gmemA,
                       &gmemB,
                       &gmemBX,
                       &gmemXY);

        ptotal_memory_in_nbytes = sizeof(T) * (gmemA + gmemB + gmemBX + gmemXY);
}


template void get_total_memory<MYTYPE>( SizeType noperator,
                       SizeType npatches,
                       const std::vector<MYTYPE*>& Amatrix_,
                       const VectorSizeType& ld_Amatrix_,
                       const std::vector<MYTYPE*>& Bmatrix_,
                       const VectorSizeType& ld_Bmatrix_,
                       const VectorSizeType& left_patch_size_,
                       const VectorSizeType& right_patch_size_,
                       SizeType& ptotal_memory_in_nbytes
                       );

void estimate_work( SizeType npatches,
                    const VectorSizeType& left_patch_size_,
                    const VectorSizeType& right_patch_size_,
                    VectorSizeType& nC_,
                    double *ptotal_gflops,
                    double *pgmemA,
                    double *pgmemB,
                    double *pgmemBX,
                    double *pgmemXY
                    )
{
/*
 -------------------
 estimate total work
 -------------------
 */
 const SizeType ialign = 32;
 const double giga = 1000.0 * 1000.0 * 1000.0;

 assert( ptotal_gflops != NULL );
 assert( pgmemA != NULL );
 assert( pgmemB != NULL );
 assert( pgmemBX != NULL );
  
 double gmemA = 0.0;
 double gmemB = 0.0;
 double gmemBX = 0.0;
 double gmemXY = 0.0; 

 double total_flops = 0.0;
 {
 SizeType ipatch = 0;
 SizeType jpatch = 0;

 for(ipatch=1; ipatch <= npatches; ipatch++) {
         gmemXY += left_patch_size_[ipatch-1] * 
                    right_patch_size_[ipatch-1];
         };

 /* 
  * -----------------------
  * count both X, Y vectors 
  * -----------------------
  */
 gmemXY *= 2;  

 for(jpatch=1; jpatch <= npatches; jpatch++) {
 for(ipatch=1; ipatch <= npatches; ipatch++) {
    SizeType nop = nC_[((ipatch)-1)+((jpatch)-1)*npatches];
    if (nop <= 0) continue;

    double flops_total = 0.0;
    double flops_method1 = 0.0;
    double flops_method2 = 0.0;

    /*
     --------------------------------------
     Note: evaluate (B * X ) * transpose(A)
     --------------------------------------
     */
     
    
    SizeType nrowA = left_patch_size_[ipatch-1];
    SizeType ncolA = left_patch_size_[jpatch-1];
    SizeType nrowB = right_patch_size_[ipatch-1];
    SizeType ncolB = right_patch_size_[jpatch-1];
    SizeType ncolX  = ncolA;

    SizeType ldA = ialign * (( nrowA + (ialign-1))/ialign );
    SizeType ldB = ialign * (( nrowB + (ialign-1))/ialign );
    SizeType ldBX = ldB;
    gmemA += nop * ldA * ncolA;
    gmemB += nop * ldB * ncolB;
    gmemBX += nop * ldBX * ncolX;

    cal_kron_flops( nrowA, nrowB, ncolA, ncolB, 
            &flops_total, &flops_method1,   &flops_method2);

    total_flops += flops_method1*nop;
    };
    };
  };
  
 double total_gflops = total_flops/(giga);
 *ptotal_gflops = total_gflops;

 *pgmemA = gmemA;
 *pgmemB = gmemB;
 *pgmemBX = gmemBX;
 *pgmemXY = gmemXY;

}
