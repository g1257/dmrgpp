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

/*! \file ParallelWft.h
 *
 * DOC TBW FIXME
 */
#ifndef PARALLEL_WFT_H
#define PARALLEL_WFT_H



namespace Dmrg {

template<typename RealType_,typename VectorWithOffsetType,typename WaveFunctionTransfType,typename LeftRightSuperType>
class ParallelWft {

public:

	typedef RealType_ RealType;

	ParallelWft(std::vector<VectorWithOffsetType>& targetVectors,
				size_t nk,
				const WaveFunctionTransfType& wft,
				const LeftRightSuperType& lrs)
		: targetVectors_(targetVectors),
		  nk_(nk),
		  wft_(wft),
		  lrs_(lrs)
	{}

	void thread_function_(size_t threadNum,size_t blockSize,size_t total,pthread_mutex_t* myMutex)
	{
		size_t nk = nk_;
		for (size_t p=0;p<blockSize;p++) {
			size_t ix = threadNum * blockSize + p + 1;
			if (ix>=targetVectors_.size()) break;
			VectorWithOffsetType phiNew = targetVectors_[0];
			wft_.setInitialVector(phiNew,targetVectors_[ix],lrs_,nk);
			targetVectors_[ix] = phiNew;
		}
	}

	//			template<typename SomeConcurrencyType,typename SomeOtherConcurrencyType>
	//			void sync(SomeConcurrencyType& conc,SomeOtherConcurrencyType& conc2)
	//			{
	//				conc.reduce(x_,conc2);
	//			}

private:

	std::vector<VectorWithOffsetType>& targetVectors_;
	size_t nk_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
}; // class ParallelWft
} // namespace Dmrg 

/*@}*/
#endif // PARALLEL_WFT_H
