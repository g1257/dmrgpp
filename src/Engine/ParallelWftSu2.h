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
/** \file ParallelWftSu2.h
*/

#ifndef DMRG_PARALLEL_WFT_SU2_H
#define DMRG_PARALLEL_WFT_SU2_H

#include "Vector.h"
#include "Concurrency.h"
#include "PackIndices.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename VectorWithOffsetType,
         typename DmrgWaveStructType,
         typename LeftRightSuperType>
class ParallelWftSu2 {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type
	VectorVectorWithOffsetType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename BasisWithOperatorsType::FactorsType FactorsType;

public:

	typedef typename VectorWithOffsetType::value_type VectorElementType;
	typedef typename PsimagLite::Real<VectorElementType>::Type RealType;

	ParallelWftSu2(VectorWithOffsetType& psiDest,
	               const VectorWithOffsetType& psiSrc,
	               const LeftRightSuperType& lrs,
	               SizeType i0,
	               const VectorSizeType& nk,
	               const DmrgWaveStructType& dmrgWaveStruct,
	               ProgramGlobals::DirectionEnum dir)
	    : psiDest_(psiDest),
	      psiSrc_(psiSrc),
	      lrs_(lrs),
	      i0_(i0),
	      nk_(nk),
	      dmrgWaveStruct_(dmrgWaveStruct),
	      dir_(dir),
	      we_(dmrgWaveStruct_.we),
	      ws_(dmrgWaveStruct_.ws),
	      pack1_(0),
	      pack2_(0)
	{
		transposeConjugate(wsT_,ws_);
		transposeConjugate(weT_,we_);

		if (dir_ == ProgramGlobals::EXPAND_SYSTEM) {
			assert(dmrgWaveStruct_.lrs.right().permutationInverse().size()==
			       dmrgWaveStruct_.we.row());
			assert(lrs_.left().permutationInverse().size()/volumeOf(nk)==
			       dmrgWaveStruct_.ws.col());
			pack1_ = new PackIndicesType(lrs.left().permutationInverse().size());
			pack2_ = new PackIndicesType(lrs.left().permutationInverse().size()/
			                             volumeOf(nk));
		} else {
			assert(dmrgWaveStruct_.lrs.left().permutationInverse().size()==
			       dmrgWaveStruct_.ws.row());
			assert(lrs_.right().permutationInverse().size()/volumeOf(nk)==
			       dmrgWaveStruct_.we.col());
			pack1_ = new PackIndicesType(lrs.super().permutationInverse().size()/
			                             lrs.right().permutationInverse().size());
			pack2_ = new PackIndicesType(volumeOf(nk));
		}

		const FactorsType& factorsSE = lrs_.super().getFactors();
		const FactorsType& factorsS = lrs_.left().getFactors();
		const FactorsType& factorsE = lrs_.right().getFactors();

		transposeConjugate(factorsInvSE_,factorsSE);
		transposeConjugate(factorsInvS_,factorsS);
		transposeConjugate(factorsInvE_,factorsE);
	}

	~ParallelWftSu2()
	{
		delete pack1_;
		delete pack2_;
	}

	static SizeType volumeOf(const VectorSizeType& v)
	{
		assert(v.size()>0);
		SizeType ret = v[0];
		for (SizeType i=1;i<v.size();i++) ret *= v[i];
		return ret;
	}

	SizeType tasks() const { return psiDest_.effectiveSize(i0_); }

	void doTask(SizeType taskNumber, SizeType)
	{
		SizeType start = psiDest_.offset(i0_);
		psiDest_.fastAccess(i0_,taskNumber) = 0.0;
		SizeType xx = taskNumber + start;

		for (int kI = factorsInvSE_.getRowPtr(xx);
		     kI < factorsInvSE_.getRowPtr(xx+1);
		     kI++) {
			if (dir_ == ProgramGlobals::EXPAND_SYSTEM) {
				SizeType ip = 0;
				SizeType alpha = 0;
				SizeType kp = 0;
				SizeType jp = 0;
				for (int kI=factorsInvSE_.getRowPtr(xx);
				     kI<factorsInvSE_.getRowPtr(xx+1);
				     kI++) {
					pack1_->unpack(alpha,jp,static_cast<SizeType>(factorsInvSE_.getCol(kI)));
					for (int k2I=factorsInvS_.getRowPtr(alpha);
					     k2I<factorsInvS_.getRowPtr(alpha+1);
					     k2I++) {
						pack2_->unpack(ip,kp,static_cast<SizeType>(factorsInvS_.getCol(k2I)));
						psiDest_.fastAccess(i0_,taskNumber) += factorsInvSE_.getValue(kI)*
						        factorsInvS_.getValue(k2I)*
						        createAux2b(psiSrc_,ip,kp,jp,wsT_,we_,nk_);
					}
				}

			} else {
				SizeType ip = 0;
				SizeType beta = 0;
				SizeType kp = 0;
				SizeType jp = 0;
				pack1_->unpack(ip,beta,static_cast<SizeType>(factorsInvSE_.getCol(kI)));
				for (int k2I = factorsInvE_.getRowPtr(beta);
				     k2I < factorsInvE_.getRowPtr(beta+1);
				     k2I++) {
					pack2_->unpack(kp,jp,static_cast<SizeType>(factorsInvE_.getCol(k2I)));
					psiDest_.fastAccess(i0_,taskNumber) += factorsInvSE_.getValue(kI)*
					        factorsInvE_.getValue(k2I)*
					        createAux1b(psiSrc_,ip,kp,jp,ws_,weT_,nk_);
				}
			}
		}
	}

private:

	// This class has pointers, disallow copy ctor and assignment
	template<typename T1, typename T2, typename T3>
	ParallelWftSu2(const ParallelWftSu2<T1,T2,T3>&);

	template<typename T1, typename T2, typename T3>
	ParallelWftSu2& operator=(const ParallelWftSu2<T1,T2,T3>&);

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
		assert(nalpha==wsT.col());

		const FactorsType& factorsE = dmrgWaveStruct_.lrs.right().getFactors();
		const FactorsType& factorsSE = dmrgWaveStruct_.lrs.super().getFactors();
		SizeType volumeOfNk = this->volumeOf(nk);
		SparseElementType sum=0;

		SizeType kpjp = kp+jp*volumeOfNk;
		assert(kpjp<dmrgWaveStruct_.lrs.right().permutationInverse().size());
		SizeType kpjpx = dmrgWaveStruct_.lrs.right().permutationInverse(kpjp);

		for (int k2I=factorsE.getRowPtr(kpjpx);k2I<factorsE.getRowPtr(kpjpx+1);k2I++) {
			SizeType beta = factorsE.getCol(k2I);
			for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
				SizeType alpha = wsT.getCol(k);
				for (int k2=we.getRowPtr(beta);k2<we.getRowPtr(beta+1);k2++) {
					SizeType j = we.getCol(k2);
					SizeType r = alpha + j*nalpha;
					for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
						SizeType x = factorsSE.getCol(kI);
						sum += wsT.getValue(k)*we.getValue(k2)*psiSrc.slowAccess(x)*
						        factorsSE.getValue(kI)*factorsE.getValue(k2I);
					}
				}
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
	                              const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = volumeOf(nk);
		SizeType ni=dmrgWaveStruct_.ws.col();
		SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
		const FactorsType& factorsS = dmrgWaveStruct_.lrs.left().getFactors();
		const FactorsType& factorsSE = dmrgWaveStruct_.lrs.super().getFactors();
		SparseElementType sum=0;

		SizeType ipkp=ip+kp*nip;
		for (int k2I=factorsS.getRowPtr(ipkp);k2I<factorsS.getRowPtr(ipkp+1);k2I++) {
			SizeType alpha = factorsS.getCol(k2I);
			for (int k=ws.getRowPtr(alpha);k<ws.getRowPtr(alpha+1);k++) {
				SizeType i = ws.getCol(k);
				for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
					SizeType j = weT.getCol(k2);
					SizeType r = i+j*ni;
					for (int kI=factorsSE.getRowPtr(r);kI<factorsSE.getRowPtr(r+1);kI++) {
						SizeType x = factorsSE.getCol(kI);
						sum += ws.getValue(k)*weT.getValue(k2)*psiSrc.slowAccess(x)*
						        factorsSE.getValue(kI)*factorsS.getValue(k2I);
					}
				}
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
	ProgramGlobals::DirectionEnum dir_;
	const SparseMatrixType& we_;
	const SparseMatrixType& ws_;
	PackIndicesType* pack1_;
	PackIndicesType* pack2_;
	SparseMatrixType wsT_;
	SparseMatrixType weT_;
	FactorsType factorsInvSE_;
	FactorsType factorsInvS_;
	FactorsType factorsInvE_;
}; // class ParallelWftSu2
} // namespace Dmrg

/*@}*/
#endif // DMRG_PARALLEL_WFT_SU2_H

