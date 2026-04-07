#include "dmrg_vbatch.h"
#include "setup_vbatch.h"
#include "setup_sparse_batch.h"
#include "BatchedGemm.h"

template<typename T>
BatchedGemm<T>::BatchedGemm(SizeType noperator,
                            SizeType npatches,
                            VectorSizeType& left_patch_size,
                            VectorSizeType& right_patch_size,
                            std::vector<T*>& Amatrix,
                            VectorSizeType& ld_Amatrix,
                            std::vector<T*>& Bmatrix,
                            VectorSizeType& ld_Bmatrix)
    : noperator_(noperator),
      npatches_(npatches),
      ld_Abatch_(0),
      ld_Bbatch_(0),
      ld_gAbatch_(npatches+1),
      ld_gBbatch_(npatches+1)
{
	dmrg_init();

	std::cout << "noperator=" << noperator << "\n";
	std::cout << "npatches=" << npatches << "\n";
	std::cout << "ld_Abatch=" << ld_Abatch_ << "\n";
	std::cout << "ld_Bbatch=" << ld_Bbatch_ << "\n";
	std::cout << "ld_Amatrix=" << ld_Amatrix << "\n";
	std::cout << "ld_Bmatrix=" << ld_Amatrix << "\n";
	
	const bool use_sparse =  true;
	if (use_sparse) {
		setup_sparse_batch<T>(noperator,
		                      npatches,
		                      left_patch_size,
		                      right_patch_size,
		                      Amatrix,
		                      ld_Amatrix,
		                      Bmatrix,
		                      ld_Bmatrix,
		                      left_patch_start_,
		                      right_patch_start_,
		                      xy_patch_start_,
		                      nC_,
		                      &gAbatch_,
		                      ld_gAbatch_,
		                      &gBbatch_,
		                      ld_gBbatch_);
		
	} else {
		
		setup_vbatch<T>(
		            noperator,
		            npatches,
		            left_patch_size,
		            right_patch_size,
		            left_patch_start_,
		            right_patch_start_,
		            xy_patch_start_,
		            Abatch_,
		            ld_Abatch_,
		            Bbatch_,
		            ld_Bbatch_,
		            Amatrix,
		            ld_Amatrix,
		            Bmatrix,
		            ld_Bmatrix);
	}
}

template<typename T>
BatchedGemm<T>::~BatchedGemm() {

   const bool use_sparse =  true;	
   if (use_sparse) {
       unsetup_sparse_batch<T>(&gAbatch_,&gBbatch_);
   }

}

template<typename T>
void BatchedGemm<T>::apply_Htarget(T* vin, T* vout)
{
	const bool use_sparse =  true;
	if (use_sparse) {
		
		apply_Htarget_sparse<T>(
		            noperator_,
		            npatches_,
		            left_patch_start_,
		            right_patch_start_,
		            xy_patch_start_,
		            nC_,
		            gAbatch_,
		            ld_gAbatch_,
		            gBbatch_,
		            ld_gBbatch_,
		            vin,
		            vout);
	}
	else {
		/*
		apply_Htarget_vbatch(
					noperator_,
					npatches_,
					left_patch_start_,
					right_patch_start_,
					xy_patch_start_,

					Abatch_, ld_Abatch_,
					Bbatch_, ld_Bbatch_,
					vin,
					vout);
		*/
		
	}

}

template class BatchedGemm<MYTYPE>;


// FOR DMRG COMPILATION BELOW

#if defined(USE_COMPLEX_Z)

template<>
BatchedGemm<double>::BatchedGemm(SizeType noperator,
                            SizeType npatches,
                            VectorSizeType& left_patch_size,
                            VectorSizeType& right_patch_size,
                            std::vector<double*>& Amatrix,
                            VectorSizeType& ld_Amatrix,
                            std::vector<double*>& Bmatrix,
                            VectorSizeType& ld_Bmatrix)
    : noperator_(noperator),
      npatches_(npatches),
      ld_Abatch_(0),
      ld_Bbatch_(0),
      ld_gAbatch_(npatches+1),
      ld_gBbatch_(npatches+1)	
{}

template<>
void BatchedGemm<double>::apply_Htarget(double* vin, double* vout)
{}

template<> BatchedGemm<double>::~BatchedGemm() {}

#else

#include<complex>

template<>
BatchedGemm<std::complex<double> >::BatchedGemm(SizeType noperator,
                            SizeType npatches,
                            VectorSizeType& left_patch_size,
                            VectorSizeType& right_patch_size,
                            std::vector<std::complex<double>* >& Amatrix,
                            VectorSizeType& ld_Amatrix,
                            std::vector<std::complex<double>* >& Bmatrix,
                            VectorSizeType& ld_Bmatrix)
    : noperator_(noperator),
      npatches_(npatches),
      ld_Abatch_(0),
      ld_Bbatch_(0),
      ld_gAbatch_(npatches+1),
      ld_gBbatch_(npatches+1)	
{}

template<>
void BatchedGemm<std::complex<double> >::apply_Htarget(std::complex<double>* vin, std::complex<double>* vout)
{}

template<> BatchedGemm<std::complex<double> >::~BatchedGemm() {}

#endif

