#include "DMRGConfig.h"
#include "analysis.h"
#include "dmrg_types.h"
#include "dmrg_vbatch.h"
#include "estimate_work.h"
#include "setup_matrix.h"
#include "setup_nC.h"
#include "setup_sparse_batch.h"
#include "setup_vbatch.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MAGMA
#include "cuda.h"
#include "magma_v2.h"
#endif

/*
 ---------------------------------------
 simple program to test vbatch
 ---------------------------------------
*/

IntegerType main(IntegerType argc, char* argv[])
{
	const IntegerType ialign          = 32;
	const IntegerType idebug          = 1;
	SizeType          noperator       = 2;
	IntegerType       left_size       = 8;
	IntegerType       max_keep_states = 2000;

	const double giga = 1000.0 * 1000.0 * 1000.0;

	size_t nbytes_X = 0;
	size_t nbytes_Y = 0;

	/*
	 -----------------------------------------------------------
	 ./test_vbatch -o noperator -n left_sites -m max_keep_states
	 -----------------------------------------------------------
	 */
#ifdef USE_MAGMA
	const char* getopt_string = "g:o:n:m:";
#else
	const char* getopt_string = "o:n:m:";
#endif

	IntegerType opt = 0;
	while ((opt = getopt(argc, argv, getopt_string)) != -1) {
		switch (opt) {
		case 'n':
			left_size = atoi(optarg);
			break;
		case 'o':
			noperator = atoi(optarg);
			break;
		case 'm':
			max_keep_states = atoi(optarg);
			break;
#ifdef USE_MAGMA
		case 'g':
		{
			IntegerType max_gpus = atoi(optarg);
			if (max_gpus <= 0) {
				max_gpus = 1;
			};
			printf("calling dmrg_set_max_gpus(%d)\n", max_gpus);
			dmrg_set_max_gpus(max_gpus);
			break;
		}
#endif
		default: /* '?' */
#ifdef USE_MAGMA
			fprintf(stderr,
			        "Usage: %s [-g max_gpus] [-o noperator] [-n left_size] [-m "
			        "max_states]\n",
			        argv[0]);
#else
			fprintf(stderr,
			        "Usage: %s [-o noperator] [-n left_size] [-m max_states]\n",
			        argv[0]);
#endif
			exit(EXIT_FAILURE);
		};
	};

	dmrg_init();

	IntegerType right_size  = left_size;
	IntegerType target_up   = (left_size + right_size) / 2;
	IntegerType target_down = target_up;
	IntegerType max_patches = (1 + target_up) * (1 + target_down);
	/*
	 ---------------------------
	 assume left part is growing
	 ---------------------------
	*/
	const IntegerType keep_left_states  = 4 * max_keep_states;
	const IntegerType keep_right_states = max_keep_states;

	const IntegerType max_patches_dim = ialign * ICEIL(max_patches + 1, ialign);

	fprintf(stderr, "max_patches_dim=%i\n", max_patches_dim);
#ifdef USE_STACK
	IntegerType left_patch_size_[max_patches_dim];
	IntegerType left_patch_up_[max_patches_dim];
	IntegerType left_patch_down_[max_patches_dim];

	IntegerType right_patch_size_[max_patches_dim];
	IntegerType right_patch_up_[max_patches_dim];
	IntegerType right_patch_down_[max_patches_dim];
#else

	std::vector<SizeType> left_patch_size_(max_patches_dim + 1);
	std::vector<SizeType> left_patch_up_(max_patches_dim);
	std::vector<SizeType> left_patch_down_(max_patches_dim);

	std::vector<SizeType> right_patch_size_(max_patches_dim + 1);
	std::vector<SizeType> right_patch_up_(max_patches_dim);
	std::vector<SizeType> right_patch_down_(max_patches_dim);

#endif

	// #define left_patch_size(i) left_patch_size_[(i)-1]
	// #define right_patch_size(i) right_patch_size_[(i)-1]

	/*
	  {
	   IntegerType nthreads = 1;
	#ifdef _OPENMP
	   #pragma omp parallel
	   #pragma omp master
	   { nthreads =  omp_get_num_threads(); }
	#endif
	   printf("using %d threads\n", nthreads );
	   }
	*/

	SizeType npatches = gen_patches_comb(left_size,
	                                     right_size,
	                                     target_up,
	                                     target_down,
	                                     keep_left_states,
	                                     keep_right_states,

	                                     left_patch_size_,
	                                     right_patch_size_,

	                                     left_patch_up_,
	                                     left_patch_down_,

	                                     right_patch_up_,
	                                     right_patch_down_);

	/*
	 ---------------------------
	 estimate length of X vector
	 ---------------------------
	 */
	IntegerType xy_size = 0;
	{
		for (SizeType ipatch = 1; ipatch <= npatches; ipatch++) {
			IntegerType nrowX = right_patch_size_[ipatch - 1];
			IntegerType ncolX = left_patch_size_[ipatch - 1];
			xy_size += (nrowX * ncolX);
		};
	}

	printf("left_size=%d, right_size=%d, target_up=%d, target_down=%d\n",
	       left_size,
	       right_size,
	       target_up,
	       target_down);
	printf("keep_left_states=%d, keep_right_states=%d, noperator=%lu, xy_size=%d\n",
	       keep_left_states,
	       keep_right_states,
	       noperator,
	       xy_size);

	printf("npatches=%lu\n", npatches);

	if (idebug >= 1) {
		for (SizeType ipatch = 1; ipatch <= npatches; ipatch++) {
			printf("ipatch=%lu, left_patch_size=%lu, right_patch_size=%lu\n",
			       ipatch,
			       left_patch_size_[ipatch - 1],
			       right_patch_size_[ipatch - 1]);
		};
		fflush(stdout);
	}

	std::vector<FpType>      Abatch_;
	std::vector<FpType>      Bbatch_;
	FpType**                 pAbatch_ = NULL;
	FpType**                 pBbatch_ = NULL;
	std::vector<IntegerType> ld_pAbatch_(npatches, 0);
	std::vector<IntegerType> ld_pBbatch_(npatches, 0);

	std::vector<SizeType> left_patch_start_(npatches, 0);
	std::vector<SizeType> right_patch_start_(npatches, 0);
	std::vector<SizeType> xy_patch_start_(npatches, 0);
	std::vector<SizeType> nC_(npatches * npatches, 0);

	SizeType ld_Abatch = 0;
	SizeType ld_Bbatch = 0;
#define xy_patch_start(i) xy_patch_start_[(i) - 1]

	std::vector<FpType*>  Amatrix(npatches * npatches * noperator, nullptr);
	std::vector<FpType*>  Bmatrix(npatches * npatches * noperator, nullptr);
	std::vector<SizeType> ld_Amatrix(npatches * npatches * noperator, 0);
	std::vector<SizeType> ld_Bmatrix(npatches * npatches * noperator, 0);

	setup_matrix(noperator, npatches, Abatch_, left_patch_size_, Amatrix, ld_Amatrix);
	// assert( Amatrix != NULL );
	// assert( ld_Amatrix != NULL );

	setup_matrix(noperator, npatches, Bbatch_, right_patch_size_, Bmatrix, ld_Bmatrix);
	// assert( Bmatrix != NULL );
	// assert( ld_Bmatrix != NULL );

	const IntegerType use_sparse = (1 == 1);

	if (use_sparse) {
		setup_sparse_batch(noperator,
		                   npatches,
		                   left_patch_size_,
		                   right_patch_size_,

		                   Amatrix,
		                   ld_Amatrix,
		                   Bmatrix,
		                   ld_Bmatrix,

		                   left_patch_start_,
		                   right_patch_start_,
		                   xy_patch_start_,

		                   nC_,
		                   &pAbatch_,
		                   ld_pAbatch_,
		                   &pBbatch_,
		                   ld_pBbatch_);

	} else {
		nC_.resize(npatches * npatches, noperator);

		setup_vbatch(noperator,
		             npatches,
		             left_patch_size_,
		             right_patch_size_,
		             left_patch_start_,
		             right_patch_start_,
		             xy_patch_start_,
		             Abatch_,
		             ld_Abatch,
		             Bbatch_,
		             ld_Bbatch,
		             Amatrix,
		             ld_Amatrix,
		             Bmatrix,
		             ld_Bmatrix);
	};

	double total_gflops = 0;
	{
		double gmemA  = 0;
		double gmemB  = 0;
		double gmemBX = 0;
		double gmemXY = 0;

		estimate_work(npatches,
		              left_patch_size_,
		              right_patch_size_,
		              nC_,
		              &total_gflops,
		              &gmemA,
		              &gmemB,
		              &gmemBX,
		              &gmemXY);

		double gmemA_gbytes  = gmemA * sizeof(FpType) / giga;
		double gmemB_gbytes  = gmemB * sizeof(FpType) / giga;
		double gmemBX_gbytes = gmemBX * sizeof(FpType) / giga;
		double gmemXY_gbytes = gmemXY * sizeof(FpType) / giga;

		printf("total_gflops=%lf \n", total_gflops);
		printf("test_vbatch:estimated memory for Amat=%lf GBytes\n", gmemA_gbytes);
		printf("test_vbatch:estimated memory for Bmat=%lf GBytes\n", gmemB_gbytes);
		printf("test_vbatch:estimated memory for BXmat=%lf GBytes\n", gmemBX_gbytes);
		printf("test_vbatch:estimated memory for X and Y =%lf GBytes\n", gmemXY_gbytes);

		size_t total_memory_in_nbytes = 0;
		get_total_memory(noperator,
		                 npatches,
		                 Amatrix,
		                 ld_Amatrix,
		                 Bmatrix,
		                 ld_Bmatrix,
		                 left_patch_size_,
		                 right_patch_size_,
		                 total_memory_in_nbytes);
		printf("total_memory_in_nbytes = %lf GBytes\n",
		       (double)total_memory_in_nbytes / giga);
	}

#define Abatch(i, j) Abatch_[indx2f(i, j, ld_Abatch)]
#define Bbatch(i, j) Bbatch_[indx2f(i, j, ld_Bbatch)]

	size_t xy_size_dim = ialign * ICEIL(xy_size, ialign);

	nbytes_X   = (sizeof(FpType) * xy_size_dim);
	FpType* X_ = (FpType*)dmrg_malloc<FpType>(nbytes_X, nbytes_X);
	assert(X_ != nullptr);

	nbytes_Y   = (sizeof(FpType) * xy_size_dim);
	FpType* Y_ = (FpType*)dmrg_malloc<FpType>(nbytes_Y, nbytes_Y);
	assert(Y_ != nullptr);

#define X(i) X_[(i) - 1]
#define Y(i) Y_[(i) - 1]

#define hX(i) hX_[(i) - 1]
#define hY(i) hY_[(i) - 1]

	{
		FpType* hX_ = X_;

		IntegerType i = 0;
		for (i = 1; i <= xy_size; i++) {
			hX(i) = ((FpType)i) / ((FpType)xy_size);
		}
	}

	{
		IntegerType       itimes = 0;
		const IntegerType ntimes = 3;

		double total_time = -dmrg_get_wtime();
		for (itimes = 1; itimes <= ntimes; itimes++) {
			double ttime = -dmrg_get_wtime();
			if (use_sparse) {

				apply_Htarget_sparse(noperator,
				                     npatches,
				                     left_patch_start_,
				                     right_patch_start_,
				                     xy_patch_start_,
				                     nC_,
				                     pAbatch_,
				                     ld_pAbatch_,
				                     pBbatch_,
				                     ld_pBbatch_,
				                     X_,
				                     Y_);

			} else {
				apply_Htarget_vbatch(noperator,
				                     npatches,
				                     left_patch_start_,
				                     right_patch_start_,
				                     xy_patch_start_,
				                     Abatch_,
				                     ld_Abatch,
				                     Bbatch_,
				                     ld_Bbatch,
				                     X_,
				                     Y_);
			};
			ttime += dmrg_get_wtime();
			printf("itimes=%d, time=%lf sec, gflops/sec=%lf\n",
			       itimes,
			       ttime,
			       total_gflops / ttime);
		};
		total_time += dmrg_get_wtime();
		printf("total_time = %lf sec, ntimes = %d, averaged gflops = %lf\n",
		       total_time,
		       ntimes,
		       ntimes * total_gflops / total_time);
		printf("test_vbatch: memory X (%f GBytes) Y (%f GBytes)\n",
		       (double)nbytes_X / (giga),
		       (double)nbytes_Y / (giga));
	}

	FpType* hY_ = Y_;
#
	/*
	 * ---------------------------------
	 * generate summary statistics for Y
	 * ---------------------------------
	 */
	{
		IntegerType i     = 0;
		FpType      Y_avg = 0;
		double      Y_max = std::abs(hY(1));
		double      Y_min = std::abs(hY(1));

		for (i = 1; i <= xy_size; i++) {
			Y_avg += hY(i);
			Y_max = std::max(Y_max, std::abs(hY(i)));
			Y_min = std::min(Y_min, std::abs(hY(i)));
		};
		Y_avg = Y_avg / ((double)xy_size);
		printf("std::abs(Y_avg) = %le, Y_max = %le Y_min = %le \n",
		       std::abs(Y_avg),
		       Y_max,
		       Y_min);
	};

	if (use_sparse) {

		unsetup_sparse_batch(&pAbatch_, &pBbatch_);

	} else {
		/*unsetup_vbatch(
		                left_patch_start_,
		                right_patch_start_,
		                xy_patch_start_,
		                Abatch_,
		                Bbatch_
		              );*/
	};
	assert(X_ != nullptr);
	dmrg_free(X_);
	assert(Y_ != nullptr);
	dmrg_free(Y_);
	dmrg_finalize();

#ifdef USE_MAGMA
	magma_finalize();
#endif
	// free(ld_pAbatch_);
	// free(pBbatch_);
	exit(0);
	return (0);
}
