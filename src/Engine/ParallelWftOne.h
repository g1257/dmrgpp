/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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
/** \file ParallelWftOne.h
*/

#ifndef DMRG_PARALLEL_WFT_ONE_H
#define DMRG_PARALLEL_WFT_ONE_H

#include "Vector.h"
#include "Concurrency.h"
#include "PackIndices.h"

namespace Dmrg {

template<typename VectorWithOffsetType,
         typename DmrgWaveStructType,
         typename LeftRightSuperType>
class ParallelWftOne {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename DmrgWaveStructType::VectorSizeType VectorSizeType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;

public:

	typedef typename VectorWithOffsetType::value_type VectorElementType;
	typedef typename PsimagLite::Real<VectorElementType>::Type RealType;

	ParallelWftOne(VectorWithOffsetType& psiDest,
	               const VectorWithOffsetType& psiSrc,
	               const LeftRightSuperType& lrs,
	               SizeType i0,
	               const VectorSizeType& nk,
	               const DmrgWaveStructType& dmrgWaveStruct,
	               typename ProgramGlobals::DirectionEnum dir)
	    : psiDest_(psiDest),
	      psiSrc_(psiSrc),
	      lrs_(lrs),
	      i0_(i0),
	      nk_(nk),
	      dmrgWaveStruct_(dmrgWaveStruct),
	      dir_(dir),
	      pack1_(0),
	      pack2_(0)
	{
		dmrgWaveStruct_.we.toSparse(we_);
        dmrgWaveStruct_.ws.toSparse(ws_);
		transposeConjugate(wsT_,ws_);
		transposeConjugate(weT_,we_);
		SizeType vOfNk = DmrgWaveStructType::volumeOf(nk);
		if (dir_ == ProgramGlobals::EXPAND_SYSTEM) {
			assert(dmrgWaveStruct_.lrs.right().permutationInverse().size()==
			       dmrgWaveStruct_.we.rows());
			assert(lrs_.left().permutationInverse().size()/vOfNk==
			       dmrgWaveStruct_.ws.cols());
			pack1_ = new PackIndicesType(lrs.left().permutationInverse().size());
			pack2_ = new PackIndicesType(lrs.left().permutationInverse().size()/vOfNk);
		} else {
			assert(dmrgWaveStruct_.lrs.left().permutationInverse().size()==
			       dmrgWaveStruct_.ws.rows());
			assert(lrs_.right().permutationInverse().size()/vOfNk==
			       dmrgWaveStruct_.we.cols());
			pack1_ = new PackIndicesType(lrs.super().permutationInverse().size()/
			                             lrs.right().permutationInverse().size());
			pack2_ = new PackIndicesType(vOfNk);
		}
	}

	~ParallelWftOne()
	{
		delete pack1_;
		delete pack2_;
	}

	SizeType tasks() const { return psiDest_.effectiveSize(i0_); }

	void doTask(SizeType taskNumber, SizeType)
	{
		SizeType start = psiDest_.offset(i0_);

		if (dir_ == ProgramGlobals::EXPAND_SYSTEM) {
			SizeType ip = 0;
			SizeType alpha = 0;
			SizeType kp = 0;
			SizeType jp = 0;
			pack1_->unpack(alpha,jp,(SizeType)lrs_.super().permutation(taskNumber+start));
			pack2_->unpack(ip,kp,(SizeType)lrs_.left().permutation(alpha));
			psiDest_.fastAccess(i0_,taskNumber) = createAux2b(psiSrc_,ip,kp,jp,wsT_,we_,nk_);
		} else {
			SizeType ip = 0;
			SizeType beta = 0;
			SizeType kp = 0;
			SizeType jp = 0;
			pack1_->unpack(ip,beta,(SizeType)lrs_.super().permutation(taskNumber+start));
			pack2_->unpack(kp,jp,(SizeType)lrs_.right().permutation(beta));
			psiDest_.fastAccess(i0_,taskNumber)=createAux1b(psiSrc_,ip,kp,jp,ws_,weT_,nk_);
		}
	}

private:

	// This class has pointers, disallow copy ctor and assignment
	template<typename T1, typename T2, typename T3>
	ParallelWftOne(const ParallelWftOne<T1,T2,T3>&);

	template<typename T1, typename T2, typename T3>
	ParallelWftOne& operator=(const ParallelWftOne<T1,T2,T3>&);

	template<typename SomeVectorType>
	SparseElementType createAux2b(const SomeVectorType& psiSrc,
	                              SizeType ip,
	                              SizeType kp,
	                              SizeType jp,
	                              const SparseMatrixType& wsT,
	                              const SparseMatrixType& we,
	                              const VectorSizeType& nk) const
	{
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		assert(nalpha==wsT.cols());

		SparseElementType sum=0;
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType beta = dmrgWaveStruct_.lrs.right().permutationInverse(kp+jp*volumeOfNk);

		for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
			SizeType alpha = wsT.getCol(k);
			SizeType begink = we.getRowPtr(beta);
			SizeType endk = we.getRowPtr(beta+1);
			for (SizeType k2=begink;k2<endk;++k2) {
				SizeType j = we.getCol(k2);
				SizeType x = dmrgWaveStruct_.lrs.super().
				        permutationInverse(alpha+j*nalpha);
				sum += wsT.getValue(k)*we.getValue(k2)*psiSrc.slowAccess(x);
			}
		}

		return sum;
	}

	template<typename SomeVectorType>
	SparseElementType createAux1b(const SomeVectorType& psiSrc,
	                              SizeType ip,
	                              SizeType kp,
	                              SizeType jp,
	                              const SparseMatrixType& ws,
	                              const SparseMatrixType& weT,
	                              const typename PsimagLite::Vector<SizeType>::Type& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType ni=dmrgWaveStruct_.ws.cols();
		SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip+kp*nip);

		SparseElementType sum=0;

		for (int k = ws.getRowPtr(alpha);k<ws.getRowPtr(alpha+1);k++) {
			SizeType i = ws.getCol(k);
			for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
				SizeType j = weT.getCol(k2);
				SizeType x = dmrgWaveStruct_.lrs.super().permutationInverse(i+j*ni);
				sum += ws.getValue(k)*weT.getValue(k2)*psiSrc.slowAccess(x);
			}
		}

		return sum;
	}

	VectorWithOffsetType& psiDest_;
	const VectorWithOffsetType& psiSrc_;
	const LeftRightSuperType& lrs_;
	SizeType i0_;
	const VectorSizeType& nk_;
	const DmrgWaveStructType& dmrgWaveStruct_;
	typename ProgramGlobals::DirectionEnum dir_;
	SparseMatrixType we_;
	SparseMatrixType ws_;
	PackIndicesType* pack1_;
	PackIndicesType* pack2_;
	SparseMatrixType wsT_;
	SparseMatrixType weT_;
}; // class ParallelWftOne
} // namespace Dmrg

/*@}*/
#endif // DMRG_PARALLEL_WFT_ONE_H

