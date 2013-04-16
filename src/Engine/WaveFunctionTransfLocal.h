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
		typedef std::vector<SparseElementType> VectorType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef typename BasisType::FactorsType FactorsType;
		typedef typename DmrgWaveStructType::LeftRightSuperType
					LeftRightSuperType;
		typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;

		static const size_t INFINITE = ProgramGlobals::INFINITE;
		static const size_t EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
		static const size_t EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;
		
		WaveFunctionTransfLocal(const size_t& stage,
		                        const bool& firstCall,
		                        const size_t& counter,
								const DmrgWaveStructType& dmrgWaveStruct,
								bool twoSiteDmrg)
		: stage_(stage),
		  firstCall_(firstCall),
		  counter_(counter),
		  dmrgWaveStruct_(dmrgWaveStruct),
		  twoSiteDmrg_(twoSiteDmrg),
		  progress_("WaveFunctionTransfLocal",0)
		{
			std::ostringstream msg;
			msg<<"Constructing...";
			progress_.printline(msg,std::cout);
		}

		virtual void transformVector(VectorWithOffsetType& psiDest,
		                             const VectorWithOffsetType& psiSrc,
		                             const LeftRightSuperType& lrs,
		                             size_t nk) const

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
							  size_t nk) const
		{
			if (twoSiteDmrg_)
				return transformVector1FromInfinite(psiDest,psiSrc,lrs,nk);

			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector1(psiDest,psiSrc,lrs,i0,nk);
			}
		}

		template<typename SomeVectorType>
		void transformVector1(SomeVectorType& psiDest,
		                      const SomeVectorType& psiSrc,
		                      const LeftRightSuperType& lrs,
		                      size_t i0,
		                      size_t nk) const
		{
			size_t nip = lrs.super().permutationInverse().size()/
					lrs.right().permutationInverse().size();

			assert(dmrgWaveStruct_.lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.row());
			assert(lrs.right().permutationInverse().size()/nk==dmrgWaveStruct_.we.col());

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			const SparseMatrixType& ws = dmrgWaveStruct_.ws;
			const SparseMatrixType& we = dmrgWaveStruct_.we;
			SparseMatrixType weT;
			transposeConjugate(weT,we);
			
			PackIndicesType pack1(nip);
			PackIndicesType pack2(nk);
			for (size_t x=start;x<final;x++) {
				size_t ip,beta,kp,jp;
				pack1.unpack(ip,beta,(size_t)lrs.super().permutation(x));
				pack2.unpack(kp,jp,(size_t)lrs.right().permutation(beta));
				psiDest[x]=createAux1b(psiSrc,ip,kp,jp,ws,weT,nk);
			}
		}

		template<typename SomeVectorType>
		SparseElementType createAux1b(
				const SomeVectorType& psiSrc,
				size_t ip,
				size_t kp,
				size_t jp,
				const SparseMatrixType& ws,
				const SparseMatrixType& weT,
				size_t nk) const
		{
			size_t ni=dmrgWaveStruct_.ws.col();
			size_t nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/nk;
			size_t alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip+kp*nip);
			
			SparseElementType sum=0;

			for (int k = ws.getRowPtr(alpha);k<ws.getRowPtr(alpha+1);k++) {
				size_t i = ws.getCol(k);
				for (int k2=weT.getRowPtr(jp);k2<weT.getRowPtr(jp+1);k2++) {
					size_t j = weT.getCol(k2);
					size_t x = dmrgWaveStruct_.lrs.super().permutationInverse(i+j*ni);
					sum += ws.getValue(k)*weT.getValue(k2)*psiSrc[x];
				}
			}

			return sum;
		}

		template<typename SomeVectorType>
		void transformVector1FromInfinite(SomeVectorType& psiDest,
										  const SomeVectorType& psiSrc,
										  const LeftRightSuperType& lrs,
										  size_t nk) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector1FromInfinite(psiDest,psiSrc,lrs,i0,nk);
			}
		}

		template<typename SomeVectorType>
		void transformVector1FromInfinite(SomeVectorType& psiDest,
										  const SomeVectorType& psiSrc,
										  const LeftRightSuperType& lrs,
										  size_t i0,
										  size_t nk) const
		{
			size_t nip = lrs.super().permutationInverse().size()/
					lrs.right().permutationInverse().size();

			assert(lrs.left().permutationInverse().size()==nk ||
				   lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.row());
			assert(lrs.right().permutationInverse().size()/nk==dmrgWaveStruct_.we.col());

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;

			SparseMatrixType ws(dmrgWaveStruct_.ws);
			SparseMatrixType we(dmrgWaveStruct_.we);
			SparseMatrixType weT;
			transposeConjugate(weT,we);

			PackIndicesType pack1(nip);
			PackIndicesType pack2(nk);
			for (size_t x=start;x<final;x++) {
				size_t ip,beta,kp,jp;
				pack1.unpack(ip,beta,(size_t)lrs.super().permutation(x));
				pack2.unpack(kp,jp,(size_t)lrs.right().permutation(beta));
				psiDest[x]=createAux1bFromInfinite(psiSrc,ip,kp,jp,ws,weT,nk);
			}
		}

		template<typename SomeVectorType>
		SparseElementType createAux1bFromInfinite(const SomeVectorType& psiSrc,
												  size_t ip,
												  size_t kp,
												  size_t jp,
												  const SparseMatrixType& ws,
												  const SparseMatrixType& weT,
												  size_t nk) const
		{
			size_t ni=dmrgWaveStruct_.lrs.left().size(); //dmrgWaveStruct_.ws.col();
			size_t nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/nk;
			MatrixOrIdentityType wsRef2(twoSiteDmrg_ && nip>nk,ws);


			SparseElementType sum=0;
			for (size_t k3=wsRef2.getRowPtr(ip);k3<wsRef2.getRowPtr(ip+1);k3++) {
				size_t ip2 = wsRef2.getCol(k3);
				size_t alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip2+kp*nip);

				for (int k = weT.getRowPtr(jp);k<weT.getRowPtr(jp+1);k++) {
					size_t jp2 = weT.getCol(k);
					size_t x = dmrgWaveStruct_.lrs.super().permutationInverse(alpha+jp2*ni);
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
				size_t nk) const
		{
			if (twoSiteDmrg_)
				return transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);

			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector2(psiDest,psiSrc,lrs,i0,nk);
			}
		}

		template<typename SomeVectorType>
		void transformVector2(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,size_t i0,
				size_t nk) const
		{
			size_t nip = lrs.left().permutationInverse().size()/nk;
			size_t nalpha = lrs.left().permutationInverse().size();

			assert(dmrgWaveStruct_.lrs.right().permutationInverse().size()==dmrgWaveStruct_.we.row());
			assert(nip==dmrgWaveStruct_.ws.col());

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			const SparseMatrixType& we = dmrgWaveStruct_.we;
			const SparseMatrixType& ws = dmrgWaveStruct_.ws;
			SparseMatrixType wsT;
			transposeConjugate(wsT,ws);
			
			PackIndicesType pack1(nalpha);
			PackIndicesType pack2(nip);
			for (size_t x=start;x<final;x++) {
				size_t ip,alpha,kp,jp;
				pack1.unpack(alpha,jp,(size_t)lrs.super().permutation(x));
				pack2.unpack(ip,kp,(size_t)lrs.left().permutation(alpha));
				psiDest[x]=createAux2b(psiSrc,ip,kp,jp,wsT,we,nk);
			}
		}

		template<typename SomeVectorType>
		SparseElementType createAux2b(
				const SomeVectorType& psiSrc,
				size_t ip,
				size_t kp,
				size_t jp,
				const SparseMatrixType& wsT,
				const SparseMatrixType& we,
				size_t nk) const
		{
			size_t nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
			assert(nalpha==wsT.col());

			SparseElementType sum=0;
			size_t beta = dmrgWaveStruct_.lrs.right().permutationInverse(kp+jp*nk);

			for (int k=wsT.getRowPtr(ip);k<wsT.getRowPtr(ip+1);k++) {
				size_t alpha = wsT.getCol(k);
				size_t begink = we.getRowPtr(beta);
				size_t endk = we.getRowPtr(beta+1);
				for (size_t k2=begink;k2<endk;++k2) {
					size_t j = we.getCol(k2);
					size_t x = dmrgWaveStruct_.lrs.super().permutationInverse(alpha+j*nalpha);
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
				size_t nk) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector2FromInfinite(psiDest,psiSrc,lrs,i0,nk);
			}
		}
		

		template<typename SomeVectorType>
		void transformVector2FromInfinite(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				size_t i0,
				size_t nk) const
		{
			size_t nip = lrs.left().permutationInverse().size()/nk;
			size_t nalpha = lrs.left().permutationInverse().size();
			
			std::ostringstream msg;
			msg<<" We're moving to the finite loop, bumpy ride ahead!";
			progress_.printline(msg,std::cout);
			

			assert(nip==dmrgWaveStruct_.ws.col());
			assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			const SparseMatrixType& we = dmrgWaveStruct_.we;
			const SparseMatrixType& ws = dmrgWaveStruct_.ws;
			SparseMatrixType wsT;
			transposeConjugate(wsT,ws);
			
			PackIndicesType pack1(nalpha);
			PackIndicesType pack2(nip);
			for (size_t x=start;x<final;x++) {
				size_t isn,jen;
				pack1.unpack(isn,jen,(size_t)lrs.super().permutation(x));
				size_t is,jpl;
				pack2.unpack(is,jpl,(size_t)lrs.left().permutation(isn));
				psiDest[x]=createAux2bFromInfinite(psiSrc,is,jpl,jen,wsT,we,nk);
			}
		}

		// FIXME: INCOMING jen needs to be 4 times as big!!
		template<typename SomeVectorType>
		SparseElementType createAux2bFromInfinite(
				const SomeVectorType& psiSrc,
				size_t is,
				size_t jpl,
				size_t jen,
				const SparseMatrixType& wsT,
				const SparseMatrixType& we,
				size_t nk) const
		{
			size_t nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
			SparseElementType sum=0;
			size_t ni = dmrgWaveStruct_.lrs.right().size()/nk;

			MatrixOrIdentityType weRef(twoSiteDmrg_ && ni>nk,we);

			for (int k=wsT.getRowPtr(is);k<wsT.getRowPtr(is+1);k++) {
				size_t ip = wsT.getCol(k);
				SparseElementType sum2 = 0;
				for (size_t k2=weRef.getRowPtr(jen);k2<weRef.getRowPtr(jen+1);k2++) {
					size_t jpr = weRef.getCol(k2);
					size_t jp = dmrgWaveStruct_.lrs.right().permutationInverse(jpl + jpr*nk);
					size_t y = dmrgWaveStruct_.lrs.super().permutationInverse(ip + jp*nalpha);
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
				size_t nk) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector1bounce(psiDest,psiSrc,lrs,i0,nk);
			}
		}
		
		template<typename SomeVectorType>
		void transformVector1bounce(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				size_t i0,
				size_t nk) const
		{
			size_t nip = lrs.super().permutationInverse().size()/lrs.right().permutationInverse().size();
			std::ostringstream msg;
			msg<<" We're bouncing on the right, so buckle up!";
			progress_.printline(msg,std::cout);
			
			assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());
			
			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			
			size_t nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
			PackIndicesType pack1(nip);
			PackIndicesType pack2(nk);
			MatrixOrIdentityType wsRef(twoSiteDmrg_,dmrgWaveStruct_.ws);
			size_t nip2 = (twoSiteDmrg_) ? dmrgWaveStruct_.ws.col() : nip;

			for (size_t x=start;x<final;x++) {
				psiDest[x] = 0.0;
				size_t ip,beta,kp,jp;
				pack1.unpack(ip,beta,(size_t)lrs.super().permutation(x));
				for (size_t k=wsRef.getRowPtr(ip);k<wsRef.getRowPtr(ip+1);k++) {
					size_t ip2 = wsRef.getCol(k);
					pack2.unpack(kp,jp,(size_t)lrs.right().permutation(beta));
					size_t ipkp = dmrgWaveStruct_.lrs.left().permutationInverse(ip2 + kp*nip2);
					size_t y = dmrgWaveStruct_.lrs.super().permutationInverse(ipkp + jp*nalpha);
					psiDest[x] += psiSrc[y]*wsRef.getValue(k);
				}
			}
		}
		
		template<typename SomeVectorType>
		void transformVector2bounce(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				size_t nk) const
		{
			for (size_t ii=0;ii<psiDest.sectors();ii++) {
				size_t i0 = psiDest.sector(ii);
				transformVector2bounce(psiDest,psiSrc,lrs,i0,nk);
			}
		}
		

		template<typename SomeVectorType>
		void transformVector2bounce(
				SomeVectorType& psiDest,
				const SomeVectorType& psiSrc,
				const LeftRightSuperType& lrs,
				size_t i0,
				size_t nk) const
		{
			size_t nip = lrs.left().permutationInverse().size()/nk;
			size_t nalpha = lrs.left().permutationInverse().size();
			
			std::ostringstream msg;
			msg<<" We're bouncing on the left, so buckle up!";
			progress_.printline(msg,std::cout);
			
			assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

			size_t start = psiDest.offset(i0);
			size_t final = psiDest.effectiveSize(i0)+start;
			PackIndicesType pack1(nalpha);
			PackIndicesType pack2(nip);
			MatrixOrIdentityType weRef(twoSiteDmrg_,dmrgWaveStruct_.we);

			for (size_t x=start;x<final;x++) {
				psiDest[x] = 0.0;

				size_t ip,alpha,kp,jp;
				pack1.unpack(alpha,jp,(size_t)lrs.super().permutation(x));
				pack2.unpack(ip,kp,(size_t)lrs.left().permutation(alpha));

				for (size_t k=weRef.getRowPtr(jp);k<weRef.getRowPtr(jp+1);k++) {
					size_t jp2 = weRef.getCol(k);
					size_t kpjp = dmrgWaveStruct_.lrs.right().permutationInverse(kp + jp2*nk);

					size_t y = dmrgWaveStruct_.lrs.super().permutationInverse(ip + kpjp*nip);
					psiDest[x] += psiSrc[y] * weRef.getValue(k);
				}
			}
			
		}

		const size_t& stage_;
		const bool& firstCall_;
		const size_t& counter_;
		const DmrgWaveStructType& dmrgWaveStruct_;
		bool twoSiteDmrg_;
		PsimagLite::ProgressIndicator progress_;
	}; // class WaveFunctionTransfLocal
} // namespace Dmrg

/*@}*/
#endif
