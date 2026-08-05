#include "dmrg_vbatch.h"
#include "setup_sparse_batch.h"
#include <cassert>
#include <complex>

template <typename T>
void unsetup_sparse_batch(std::vector<T*>& gAbatch_, std::vector<T*>& gBbatch_)
{

	T* pAmem = gAbatch_[0];
	T* pBmem = gBbatch_[0];
	assert(pAmem != NULL);
	assert(pBmem != NULL);

	delete pAmem;
	delete pBmem;
}

template void unsetup_sparse_batch<double>(std::vector<double*>&, std::vector<double*>&);

template void unsetup_sparse_batch<std::complex<double>>(std::vector<std::complex<double>*>&,
                                                         std::vector<std::complex<double>*>&);
