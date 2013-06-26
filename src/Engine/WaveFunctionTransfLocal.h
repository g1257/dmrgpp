/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file WaveFunctionTransformation.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_LOCAL_HEADER_H
#define WFT_LOCAL_HEADER_H

#include "PackIndices.h"
#include "ProgressIndicator.h"
#include "VectorWithOffsets.h" // so that std::norm() becomes visible here
#include "VectorWithOffset.h" // so that std::norm() becomes visible here
#include "WaveFunctionTransfBase.h"
#include "MatrixOrIdentity.h"

namespace Dmrg {

	template<typename DmrgWaveStructType,typename VectorWithOffsetType>
	class WaveFunctionTransfLocal : public
		WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType> {

		typedef PsimagLite::PackIndices PackIndicesType;

	public:
		typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
		typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef typename BasisType::FactorsType FactorsType;
		typedef typename DmrgWaveStructType::LeftRightSuperType
					LeftRightSuperType;
		typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;

		static const SizeType INFINITE = ProgramGlobals::INFINITE;
		static const SizeType EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
		static const SizeType EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;
		
		WaveFunctionTransfLocal(const SizeType& stage,
		                        const bool& firstCall,
		                        const SizeType& counter,
								const DmrgWaveStructType& dmrgWaveStruct,
								bool twoSiteDmrg)
		: stage_(stage),
		  firstCall_(firstCall),
		  counter_(counter),
		  dmrgWaveStruct_(dmrgWaveStruct),
		  twoSiteDmrg_(twoSiteDmrg),
		  progress_("WaveFunctionTransfLocal")
		{
			PsimagLite::OstringStream msg;
			msg<<"Constructing...";
			progress_.printline(msg,std::cout);
		}

		virtual void transformVector(VectorWithOffsetType& psiDest,
		                             const VectorWithOffsetType& psiSrc,
		                             const LeftRightSuperType& lrs,
		                             const typename PsimagLite::Vector<SizeType>::Type& nk) const

		{
			if (stage_==EXPAND_ENVIRON) {
				if (firstCall_) {
					transformVector1FromInfinite(psiDest,psiSrc,lrs,nk);
				} else if (counter_==0) {
					transformVector1bounce(psiDest,psiSrc,lrs,nk);
				} else {
					transformVector1(psiDest,psiSrc,lrs,nk);
				}
			}
			if (stage_==EXPAND_SYSTEM) {
				if (firstCall_)
					transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);
 				else if (counter_==0) 
					transformVector2bounce(psiDest,psiSrc,lrs,nk);
				else transformVector2(psiDest,psiSrc,lrs,nk);
			}
		}

	private:

		template<typename SomeVectorType>
		void transformVector1(SomeVectorType& psiDest,
		                      const SomeVectorType& psiSrc,
		                      const LeftRightSuperType& lrs,
							  const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			if (twoSiteDmrg_)
				return transformVector1FromInfinite(psiDest,psiSrc,lrs,nk);

			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				transformVector1(psiDest,psiSrc,lrs,i0,nk);
			}
		}

		template<typename SomeVectorType>
		void transformVector1(SomeVectorType& psiDest,
		                      const SomeVectorType& psiSrc,
		                      const LeftRightSuperType& lrs,
		                      SizeType i0,
		                      const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType nip = lrs.super().permutationInverse().size()/
					lrs.right().permutationInverse().size();

			assert(dmrgWaveStruct_.lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.row());
			assert(lrs.right().permutationInverse().size()/volumeOfNk==dmrgWaveStruct_.we.col());

			SizeType start = psiDest.offset(i0);
			SizeType final = psiDest.effectiveSize(i0)+start;
			
			const SparseMatrixType& ws = dmrgWaveStruct_.ws;
			const SparseMatrixType& we = dmrgWaveStruct_.we;
			SparseMatrixType weT;
			transposeConjugate(weT,we);
			
			PackIndicesType pack1(nip);
			PackIndicesType pack2(volumeOfNk);
			for (SizeType x=start;x<final;x++) {
				SizeType ip,beta,kp,jp;
				pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x));
				pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
				psiDest[x]=createAux1b(psiSrc,ip,kp,jp,ws,weT,nk);
			}
		}

		template<typename SomeVectorType>
		SparseElementType createAux1b(
				const SomeVectorType& psiSrc,
				SizeType ip,
				SizeType kp,
				SizeType jp,
				const SparseMatrixType& ws,
				const SparseMatrixType& weT,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType ni=dmrgWaveStruct_.ws.col();
			SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
			SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip+kp*nip);
			
			SparseElementType sum=0;

			for (int k = ws.getRowPtr(alpha);k<ws.getRowPtr(alpha+1);k++) {
				SizeType i = ws.getCol(k);
				for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
					SizeType j = weT.getCol(k2);
					SizeType x = dmrgWaveStruct_.lrs.super().permutationInverse(i+j*ni);
					sum += ws.getValue(k)*weT.getValue(k2)*psiSrc[x];
				}
			}

			return sum;
		}

		template<typename SomeVectorType>
		void transformVector1FromInfinite(SomeVectorType& psiDest,
										  const SomeVectorType& psiSrc,
										  const LeftRightSuperType& lrs,
										  const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				transformVector1FromInfinite(psiDest,psiSrc,lrs,i0,nk);
			}
		}

		template<typename SomeVectorType>
		void transformVector1FromInfinite(SomeVectorType& psiDest,
										  const SomeVectorType& psiSrc,
										  const LeftRightSuperType& lrs,
										  SizeType i0,
										  const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType nip = lrs.super().permutationInverse().size()/
					lrs.right().permutationInverse().size();

			assert(lrs.left().permutationInverse().size()==volumeOfNk ||
				   lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.row());
			assert(lrs.right().permutationInverse().size()/volumeOfNk==dmrgWaveStruct_.we.col());

			SizeType start = psiDest.offset(i0);
			SizeType final = psiDest.effectiveSize(i0)+start;

			SparseMatrixType ws(dmrgWaveStruct_.ws);
			SparseMatrixType we(dmrgWaveStruct_.we);
			SparseMatrixType weT;
			transposeConjugate(weT,we);

			PackIndicesType pack1(nip);
			PackIndicesType pack2(volumeOfNk);
			for (SizeType x=start;x<final;x++) {
				SizeType ip,beta,kp,jp;
				pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x));
				pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
				psiDest[x]=createAux1bFromInfinite(psiSrc,ip,kp,jp,ws,weT,nk);
			}
		}

		template<typename SomeVectorType>
		SparseElementType createAux1bFromInfinite(const SomeVectorType& psiSrc,
												  SizeType ip,
												  SizeType kp,
												  SizeType jp,
												  const SparseMatrixType& ws,
												  const SparseMatrixType& weT,
												  const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType ni=dmrgWaveStruct_.lrs.left().size(); //dmrgWaveStruct_.ws.col();
			SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
			MatrixOrIdentityType wsRef2(twoSiteDmrg_ && nip>volumeOfNk,ws);


			SparseElementType sum=0;
			for (SizeType k3=wsRef2.getRowPtr(ip);k3<wsRef2.getRowPtr(ip+1);k3++) {
				SizeType ip2 = wsRef2.getCol(k3);
				SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip2+kp*nip);

				for (int k = weT.getRowPtr(jp);k<weT.getRowPtr(jp+1);k++) {
					SizeType jp2 = weT.getCol(k);
					SizeType x = dmrgWaveStruct_.lrs.super().permutationInverse(alpha+jp2*ni);
					sum += weT.getValue(k)*psiSrc[x]*wsRef2.getValue(k3);

				}
			}
			return sum;
		}

		template<typename SomeVectorType>
		void transformVector2(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			if (twoSiteDmrg_)
				return transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);

			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				transformVector2(psiDest,psiSrc,lrs,i0,nk);
			}
		}

		template<typename SomeVectorType>
		void transformVector2(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,SizeType i0,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
			SizeType nalpha = lrs.left().permutationInverse().size();

			assert(dmrgWaveStruct_.lrs.right().permutationInverse().size()==dmrgWaveStruct_.we.row());
			assert(nip==dmrgWaveStruct_.ws.col());

			SizeType start = psiDest.offset(i0);
			SizeType final = psiDest.effectiveSize(i0)+start;
			
			const SparseMatrixType& we = dmrgWaveStruct_.we;
			const SparseMatrixType& ws = dmrgWaveStruct_.ws;
			SparseMatrixType wsT;
			transposeConjugate(wsT,ws);
			
			PackIndicesType pack1(nalpha);
			PackIndicesType pack2(nip);
			for (SizeType x=start;x<final;x++) {
				SizeType ip,alpha,kp,jp;
				pack1.unpack(alpha,jp,(SizeType)lrs.super().permutation(x));
				pack2.unpack(ip,kp,(SizeType)lrs.left().permutation(alpha));
				psiDest[x]=createAux2b(psiSrc,ip,kp,jp,wsT,we,nk);
			}
		}

		template<typename SomeVectorType>
		SparseElementType createAux2b(
				const SomeVectorType& psiSrc,
				SizeType ip,
				SizeType kp,
				SizeType jp,
				const SparseMatrixType& wsT,
				const SparseMatrixType& we,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
			assert(nalpha==wsT.col());

			SparseElementType sum=0;
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType beta = dmrgWaveStruct_.lrs.right().permutationInverse(kp+jp*volumeOfNk);

			for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
				SizeType alpha = wsT.getCol(k);
				SizeType begink = we.getRowPtr(beta);
				SizeType endk = we.getRowPtr(beta+1);
				for (SizeType k2=begink;k2<endk;++k2) {
					SizeType j = we.getCol(k2);
					SizeType x = dmrgWaveStruct_.lrs.super().permutationInverse(alpha+j*nalpha);
					sum += wsT.getValue(k)*we.getValue(k2)*psiSrc[x]; //*weRef.getValue(k3);
				}
			}
			return sum;
		}
		
		template<typename SomeVectorType>
		void transformVector2FromInfinite(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				transformVector2FromInfinite(psiDest,psiSrc,lrs,i0,nk);
			}
		}
		

		template<typename SomeVectorType>
		void transformVector2FromInfinite(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				SizeType i0,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
			SizeType nalpha = lrs.left().permutationInverse().size();

			assert(nip==dmrgWaveStruct_.ws.col());
			assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

			SizeType start = psiDest.offset(i0);
			SizeType final = psiDest.effectiveSize(i0)+start;
			
			const SparseMatrixType& we = dmrgWaveStruct_.we;
			const SparseMatrixType& ws = dmrgWaveStruct_.ws;
			SparseMatrixType wsT;
			transposeConjugate(wsT,ws);
			
			PackIndicesType pack1(nalpha);
			PackIndicesType pack2(nip);
			for (SizeType x=start;x<final;x++) {
				SizeType isn,jen;
				pack1.unpack(isn,jen,(SizeType)lrs.super().permutation(x));
				SizeType is,jpl;
				pack2.unpack(is,jpl,(SizeType)lrs.left().permutation(isn));
				psiDest[x]=createAux2bFromInfinite(psiSrc,is,jpl,jen,wsT,we,nk);
			}
		}

		// FIXME: INCOMING jen needs to be 4 times as big!!
		template<typename SomeVectorType>
		SparseElementType createAux2bFromInfinite(
				const SomeVectorType& psiSrc,
				SizeType is,
				SizeType jpl,
				SizeType jen,
				const SparseMatrixType& wsT,
				const SparseMatrixType& we,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
			SparseElementType sum=0;
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk;

			MatrixOrIdentityType weRef(twoSiteDmrg_ && ni>volumeOfNk,we);

			for (int k=wsT.getRowPtr(is);k<wsT.getRowPtr(is+1);k++) {
				SizeType ip = wsT.getCol(k);
				SparseElementType sum2 = 0;
				for (SizeType k2=weRef.getRowPtr(jen);k2<weRef.getRowPtr(jen+1);k2++) {
					SizeType jpr = weRef.getCol(k2);
					SizeType jp = dmrgWaveStruct_.lrs.right().permutationInverse(jpl + jpr*volumeOfNk);
					SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ip + jp*nalpha);
					sum2 += wsT.getValue(k)*psiSrc[y]*weRef.getValue(k2);

				}
				sum += sum2;
			}
			return sum;
		}
		
		template<typename SomeVectorType>
		void transformVector1bounce(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				transformVector1bounce(psiDest,psiSrc,lrs,i0,nk);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector1bounce(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				SizeType i0,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType nip = lrs.super().permutationInverse().size()/lrs.right().permutationInverse().size();
			PsimagLite::OstringStream msg;
			msg<<" We're bouncing on the right, so buckle up!";
			progress_.printline(msg,std::cout);
			
			assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());
			
			SizeType start = psiDest.offset(i0);
			SizeType final = psiDest.effectiveSize(i0)+start;
			
			SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
			PackIndicesType pack1(nip);
			PackIndicesType pack2(volumeOfNk);
			MatrixOrIdentityType wsRef(twoSiteDmrg_,dmrgWaveStruct_.ws);
			SizeType nip2 = (twoSiteDmrg_) ? dmrgWaveStruct_.ws.col() : nip;

			for (SizeType x=start;x<final;x++) {
				psiDest[x] = 0.0;
				SizeType ip,beta,kp,jp;
				pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x));
				for (SizeType k=wsRef.getRowPtr(ip);k<wsRef.getRowPtr(ip+1);k++) {
					SizeType ip2 = wsRef.getCol(k);
					pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
					SizeType ipkp = dmrgWaveStruct_.lrs.left().permutationInverse(ip2 + kp*nip2);
					SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ipkp + jp*nalpha);
					psiDest[x] += psiSrc[y]*wsRef.getValue(k);
				}
			}
		}
		
		template<typename SomeVectorType>
		void transformVector2bounce(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				transformVector2bounce(psiDest,psiSrc,lrs,i0,nk);
			}
		}
		

		template<typename SomeVectorType>
		void transformVector2bounce(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				SizeType i0,
				const typename PsimagLite::Vector<SizeType>::Type& nk) const
		{
			SizeType volumeOfNk = this->volumeOf(nk);
			SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
			SizeType nalpha = lrs.left().permutationInverse().size();
			
			PsimagLite::OstringStream msg;
			msg<<" We're bouncing on the left, so buckle up!";
			progress_.printline(msg,std::cout);
			
			assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

			SizeType start = psiDest.offset(i0);
			SizeType final = psiDest.effectiveSize(i0)+start;
			PackIndicesType pack1(nalpha);
			PackIndicesType pack2(nip);
			MatrixOrIdentityType weRef(twoSiteDmrg_,dmrgWaveStruct_.we);

			for (SizeType x=start;x<final;x++) {
				psiDest[x] = 0.0;

				SizeType ip,alpha,kp,jp;
				pack1.unpack(alpha,jp,(SizeType)lrs.super().permutation(x));
				pack2.unpack(ip,kp,(SizeType)lrs.left().permutation(alpha));

				for (SizeType k=weRef.getRowPtr(jp);k<weRef.getRowPtr(jp+1);k++) {
					SizeType jp2 = weRef.getCol(k);
					SizeType kpjp = dmrgWaveStruct_.lrs.right().permutationInverse(kp + jp2*volumeOfNk);

					SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ip + kpjp*nip);
					psiDest[x] += psiSrc[y] * weRef.getValue(k);
				}
			}	
		}

		const SizeType& stage_;
		const bool& firstCall_;
		const SizeType& counter_;
		const DmrgWaveStructType& dmrgWaveStruct_;
		bool twoSiteDmrg_;
		PsimagLite::ProgressIndicator progress_;
	}; // class WaveFunctionTransfLocal
} // namespace Dmrg

/*@}*/
#endif
