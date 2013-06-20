/*
Copyright (c) 2009,-2012 UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
/** \ingroup DMRG */
/*@{*/

/*! \file Parallel4PointDs.h
 *
 * DOC TBW FIXME
 */
#ifndef PARALLEL_4POINT_DS_H
#define PARALLEL_4POINT_DS_H

#include "Matrix.h"
#include "Mpi.h"

namespace Dmrg {

template<typename ModelType,typename FourPointCorrelationsType>
class Parallel4PointDs {

	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename FourPointCorrelationsType::MatrixType MatrixType;
	typedef typename MatrixType::value_type FieldType;

public:

	typedef typename ModelType::RealType RealType;

	Parallel4PointDs(MatrixType& fpd,
					 const FourPointCorrelationsType& fourpoint,
					 const ModelType& model,
					 const typename PsimagLite::Vector<SizeType>::Type& gammas,
					 const typename PsimagLite::Vector<PairType>::Type& pairs)
		: fpd_(fpd),fourpoint_(fourpoint),model_(model),gammas_(gammas),pairs_(pairs)
	{}

	void thread_function_(SizeType threadNum,SizeType blockSize,SizeType total,pthread_mutex_t* myMutex)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = PsimagLite::Concurrency::npthreads;

		for (SizeType p=0;p<blockSize;p++) {
			SizeType px = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (px>=total) continue;

			SizeType i = pairs_[px].first;
			SizeType j = pairs_[px].second;

			fpd_(i,j) = fourPointDelta(2*i,2*j,gammas_,model_,threadNum);
		}
	}

	//			template<typename SomeConcurrencyType,typename SomeOtherConcurrencyType>
	//			void sync(SomeConcurrencyType& conc,SomeOtherConcurrencyType& conc2)
	//			{
	//				conc.reduce(x_,conc2);
	//			}

private:

	template<typename SomeModelType>
	FieldType fourPointDelta(SizeType i,SizeType j,const typename PsimagLite::Vector<SizeType>::Type& gammas,const SomeModelType& model,SizeType threadId) const
	{
		SizeType hs = model.hilbertSize(0);
		SizeType nx = 0;
		while(hs) {
			hs>>=1;
			nx++;
		}
		nx /= 2;
		SizeType site = 0;
		const MatrixType& opC0 = model.naturalOperator("c",site,gammas[0] + 0*nx); // C_{gamma0,up}
		const MatrixType& opC1 = model.naturalOperator("c",site,gammas[1] + 1*nx); // C_{gamma1,down}
		const MatrixType& opC2 = model.naturalOperator("c",site,gammas[2] + 1*nx); // C_{gamma2,down}
		const MatrixType& opC3 = model.naturalOperator("c",site,gammas[3] + 0*nx); // C_{gamma3,up}

		return fourpoint_(
				'C',i,opC0,
			   'C',i+1,opC1,
			   'N',j,opC2,
			   'N',j+1,opC3,-1,threadId);
	}

	MatrixType& fpd_;
	const FourPointCorrelationsType& fourpoint_;
	const ModelType& model_;
	const typename PsimagLite::Vector<SizeType>::Type& gammas_;
	const typename PsimagLite::Vector<PairType>::Type& pairs_;
}; // class Parallel4PointDs
} // namespace Dmrg 

/*@}*/
#endif // PARALLEL_4POINT_DS_H
