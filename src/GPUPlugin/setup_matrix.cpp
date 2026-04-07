#include "setup_matrix.h"
#include "dmrg_types.h"
#include <assert.h>

#if defined(USE_COMPLEX_Z)
std::complex<double> makeFloat(double zr, double zi) { return std::complex<double>(zr, zi); }
#else
double makeFloat(double zr, double) { return zr; }
#endif

SizeType f1(const SizeType& ipatch,
            const SizeType& jpatch,
            const SizeType& ioperator,
            const SizeType& nrow_Abatch,
            const SizeType& ncol_Abatch);

SizeType f2(const SizeType& ipatch,
            const SizeType& jpatch,
            const SizeType& ioperator,
            const SizeType& npatches);

SizeType f3(const SizeType& ipatch,
            const SizeType& jpatch,
            const SizeType& ioperator,
            const SizeType& npatches);

template <typename T>
void setup_matrix(SizeType              noperator,
                  SizeType              npatches,
                  std::vector<T>&       Abatch,
                  const VectorSizeType& left_patch_size_,
                  std::vector<T*>&      Amatrix_,
                  VectorSizeType&       ld_Amatrix_) // INTRODUCE Abatch_
{
	SizeType ioperator       = 0;
	SizeType ipatch          = 0;
	SizeType jpatch          = 0;
	SizeType total_left_size = 0;

	VectorSizeType left_patch_start_(npatches + 1, 0);

	SizeType istart = 1;
	for (ipatch = 1; ipatch <= npatches; ipatch++) {
		total_left_size += left_patch_size_[ipatch - 1];

		left_patch_start_[ipatch - 1] = istart;
		istart += left_patch_size_[ipatch - 1];
	};

	SizeType nrow_Abatch = total_left_size;
	SizeType ncol_Abatch = total_left_size;

	//   T* Abatch = new T[nrow_Abatch*ncol_Abatch*noperator];
	Abatch.resize(nrow_Abatch * ncol_Abatch * noperator);

	for (ioperator = 1; ioperator <= noperator; ioperator++) {
		for (jpatch = 1; jpatch <= npatches; jpatch++) {
			for (ipatch = 1; ipatch <= npatches; ipatch++) {

				T* Amat_ = &(Abatch[f1(
				    ipatch, jpatch, ioperator, nrow_Abatch, ncol_Abatch)]);

				SizeType ia_start = left_patch_start_[ipatch - 1];
				SizeType ja_start = left_patch_start_[jpatch - 1];

				SizeType nrow    = left_patch_size_[ipatch - 1];
				SizeType ncol    = left_patch_size_[jpatch - 1];
				SizeType ld_Amat = nrow;

				ld_Amatrix_[f3(ipatch, jpatch, ioperator, npatches)]
				    = left_patch_size_[ipatch - 1];
				Amatrix_[f2(ipatch, jpatch, ioperator, npatches)] = Amat_;

				SizeType i = 0;
				SizeType j = 0;
				for (j = 1; j <= ncol; j++) {
					for (i = 1; i <= nrow; i++) {
						SizeType ia   = (ia_start - 1) + i;
						SizeType ja   = (ja_start - 1) + j;
						double   dval = (ia + (ja - 1) * total_left_size);
						dval += ((double)(ioperator - 1))
						    * ((double)total_left_size)
						    * ((double)total_left_size);
						Amat_[((i)-1) + ((j)-1) * ld_Amat]
						    = makeFloat(dval, -dval);
					};
				};
			};
		};
	};
}

SizeType f1(const SizeType& ipatch,
            const SizeType& jpatch,
            const SizeType& ioperator,
            const SizeType& nrow_Abatch,
            const SizeType& ncol_Abatch)
{
	SizeType mySizeType = (((ipatch)-1) + ((jpatch)-1) * nrow_Abatch)
	    + ((long)((ioperator)-1) * (nrow_Abatch) * (ncol_Abatch));
	return mySizeType;
}

SizeType f2(const SizeType& ipatch,
            const SizeType& jpatch,
            const SizeType& ioperator,
            const SizeType& npatches)
{
	SizeType mySizeType
	    = (((ipatch)-1) + ((jpatch)-1) * npatches) + ((ioperator)-1) * (npatches) * (npatches);
	return mySizeType;
}

SizeType f3(const SizeType& ipatch,
            const SizeType& jpatch,
            const SizeType& ioperator,
            const SizeType& npatches)
{
	SizeType mySizeType
	    = (((ipatch)-1) + ((jpatch)-1) * npatches) + ((ioperator)-1) * (npatches) * (npatches);
	return mySizeType;
}

template void setup_matrix<MYTYPE>(SizeType,
                                   SizeType,
                                   std::vector<MYTYPE>& Abatch,
                                   const VectorSizeType&,
                                   std::vector<MYTYPE*>&,
                                   VectorSizeType&);
