/*
Copyright (c) 2009-2013, 2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file WaveFunctionTransfLocal.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *  this is for when NOT using SU(2)
 */

#ifndef WFT_LOCAL_HEADER_H
#define WFT_LOCAL_HEADER_H

#include "ProgressIndicator.h"
#include "VectorWithOffsets.h" // so that PsimagLite::norm() becomes visible here
#include "VectorWithOffset.h" // so that PsimagLite::norm() becomes visible here
#include "WaveFunctionTransfBase.h"
#include "MatrixOrIdentity.h"
#include "ParallelWftOne.h"
#include "Parallelizer.h"
#include "MatrixVectorKron/InitKronWft.h"
#include "MatrixVectorKron/KronMatrix.h"
#include "WftAccelBlocks.h"

namespace Dmrg {

template<typename DmrgWaveStructType,typename VectorWithOffsetType>
class WaveFunctionTransfLocal : public
        WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType> {


	typedef WaveFunctionTransfBase<DmrgWaveStructType,VectorWithOffsetType> BaseType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename BaseType::PackIndicesType PackIndicesType;

public:

	typedef typename BaseType::WftOptions WftOptions;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;
	typedef ParallelWftOne<VectorWithOffsetType,
	DmrgWaveStructType,
	LeftRightSuperType> ParallelWftType;
	typedef InitKronWft<LeftRightSuperType, WftOptions, DmrgWaveStructType> InitKronType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef WftAccelBlocks<BaseType> WftAccelBlocksType;

	WaveFunctionTransfLocal(const DmrgWaveStructType& dmrgWaveStruct,
	                        const WftOptions& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct),
	      wftOptions_(wftOptions),
	      wftAccelBlocks_(dmrgWaveStruct, wftOptions),
	      progress_("WaveFunctionTransfLocal")
	{
		PsimagLite::OstringStream msg;
		msg<<"Constructing Local";
		progress_.printline(msg,std::cout);
	}

	virtual void transformVector(VectorWithOffsetType& psiDest,
	                             const VectorWithOffsetType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             const VectorSizeType& nk) const

	{
		if (wftOptions_.dir == ProgramGlobals::EXPAND_ENVIRON) {
			if (wftOptions_.firstCall) {
				transformVector1FromInfinite(psiDest,psiSrc,lrs,nk);
			} else if (wftOptions_.counter == 0) {
				transformVector1bounce(psiDest,psiSrc,lrs,nk);
			} else {
				transformVector1(psiDest,psiSrc,lrs,nk);
			}

			return;
		}

		if (wftOptions_.dir == ProgramGlobals::EXPAND_SYSTEM) {
			if (wftOptions_.firstCall)
				transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);
			else if (wftOptions_.counter == 0)
				transformVector2bounce(psiDest,psiSrc,lrs,nk);
			else transformVector2(psiDest,psiSrc,lrs,nk);

			return;
		}

		err("WFT Local: Stage is not EXPAND_ENVIRON or EXPAND_SYSTEM\n");
	}

private:

	void transformVector1(VectorWithOffsetType& psiDest,
	                      const VectorWithOffsetType& psiSrc,
	                      const LeftRightSuperType& lrs,
	                      const VectorSizeType& nk) const
	{
		if (wftOptions_.twoSiteDmrg && wftOptions_.accel != WftOptions::ACCEL_PATCHES)
			return transformVector1FromInfinite(psiDest,psiSrc,lrs,nk);

		typename ProgramGlobals::DirectionEnum dir1 = ProgramGlobals::EXPAND_ENVIRON;

		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVectorParallel(psiDest,psiSrc,lrs,i0,nk,dir1);
		}
	}

	void transformVectorParallel(VectorWithOffsetType& psiDest,
	                             const VectorWithOffsetType& psiSrc,
	                             const LeftRightSuperType& lrs,
	                             SizeType i0,
	                             const VectorSizeType& nk,
	                             typename ProgramGlobals::DirectionEnum dir) const
	{
		if (wftOptions_.accel == WftOptions::ACCEL_PATCHES)
			return transformVectorParallelPatched(psiDest, psiSrc, lrs, i0, nk, dir);

		typedef PsimagLite::Parallelizer<ParallelWftType> ParallelizerType;
		ParallelizerType threadedWft(PsimagLite::Concurrency::npthreads,
		                             PsimagLite::MPI::COMM_WORLD);

		ParallelWftType helperWft(psiDest,
		                          psiSrc,
		                          lrs,
		                          i0,
		                          nk,
		                          dmrgWaveStruct_,
		                          dir);

		threadedWft.loopCreate(helperWft);
	}

	void transformVectorParallelPatched(VectorWithOffsetType& psiDest,
	                                    const VectorWithOffsetType& psiSrc,
	                                    const LeftRightSuperType& lrs,
	                                    SizeType i0,
	                                    const VectorSizeType& nk,
	                                    typename ProgramGlobals::DirectionEnum dir) const
	{
		SizeType qn = psiDest.qn(i0);
		SizeType iOld = findIold(psiSrc, psiDest, i0);
		InitKronType initKron(lrs, i0, qn, wftOptions_, dmrgWaveStruct_, iOld);
		KronMatrix<InitKronType> kronMatrix(initKron, "WFT");
		VectorType psiDestOneSector;
		psiDest.extract(psiDestOneSector, i0);

		VectorType psiSrcOneSector;
		psiSrc.extract(psiSrcOneSector, iOld);

		kronMatrix.matrixVectorProduct(psiDestOneSector, psiSrcOneSector);
		psiDest.setDataInSector(psiDestOneSector, i0);
	}

	void transformVector1FromInfinite(VectorWithOffsetType& psiDest,
	                                  const VectorWithOffsetType& psiSrc,
	                                  const LeftRightSuperType& lrs,
	                                  const VectorSizeType& nk) const
	{
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			SizeType iOld = findIold(psiSrc, psiDest, i0);
			tVector1FromInfinite(psiDest, i0, psiSrc, iOld, lrs, nk);
		}
	}

	void tVector1FromInfinite(VectorWithOffsetType& psiDest,
	                          SizeType i0,
	                          const VectorWithOffsetType& psiSrc,
	                          SizeType iOld,
	                          const LeftRightSuperType& lrs,
	                          const VectorSizeType& nk) const
	{
		if (wftOptions_.accel == WftOptions::ACCEL_TEMP)
			return transformTemp1FromInfinite(psiDest, psiSrc, lrs, i0, nk);

		if (wftOptions_.accel == WftOptions::ACCEL_BLOCKS)
			return wftAccelBlocks_.environFromInfinite(psiDest, i0, psiSrc, iOld, lrs, nk);

		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();

		assert(lrs.left().permutationInverse().size()==volumeOfNk ||
		       lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.rows());
		assert(lrs.right().permutationInverse().size()/volumeOfNk==dmrgWaveStruct_.we.cols());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);

		SparseMatrixType ws;
		dmrgWaveStruct_.ws.toSparse(ws);
		SparseMatrixType we;
		dmrgWaveStruct_.we.toSparse(we);
		SparseMatrixType weT;
		transposeConjugate(weT,we);

		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		for (SizeType x=0;x<total;x++) {
			SizeType ip,beta,kp,jp;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
			psiDest.fastAccess(i0,x)=createAux1bFromInfinite(psiSrc,ip,kp,jp,ws,weT,nk);
		}
	}

	SparseElementType createAux1bFromInfinite(const VectorWithOffsetType& psiSrc,
	                                          SizeType ip,
	                                          SizeType kp,
	                                          SizeType jp,
	                                          const SparseMatrixType& ws,
	                                          const SparseMatrixType& weT,
	                                          const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType ni=dmrgWaveStruct_.lrs.left().size();
		SizeType nip = dmrgWaveStruct_.lrs.left().permutationInverse().size()/volumeOfNk;
		MatrixOrIdentityType wsRef2(wftOptions_.twoSiteDmrg && nip>volumeOfNk,ws);
		SizeType start = weT.getRowPtr(jp);
		SizeType end = weT.getRowPtr(jp+1);
		SparseElementType sum=0;
		for (SizeType k3=wsRef2.getRowPtr(ip);k3<wsRef2.getRowPtr(ip+1);k3++) {
			int ip2 = wsRef2.getColOrExit(k3);
			if (ip2<0) continue;
			SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip2+kp*nip);

			for (SizeType k = start; k < end; k++) {
				SizeType jp2 = weT.getCol(k);
				SizeType x = dmrgWaveStruct_.lrs.super().permutationInverse(alpha+jp2*ni);
				sum += weT.getValue(k)*psiSrc.slowAccess(x)*wsRef2.getValue(k3);

			}
		}
		return sum;
	}

	void transformTemp1FromInfinite(VectorWithOffsetType& psiDest,
	                                const VectorWithOffsetType& psiSrc,
	                                const LeftRightSuperType& lrs,
	                                SizeType i0,
	                                const VectorSizeType& nk) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();

		SizeType nip2 = dmrgWaveStruct_.lrs.left().size()/volumeOfNk;

		assert(lrs.left().permutationInverse().size()==volumeOfNk ||
		       lrs.left().permutationInverse().size()==dmrgWaveStruct_.ws.rows());
		assert(lrs.right().permutationInverse().size()/volumeOfNk==dmrgWaveStruct_.we.cols());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);

		SparseMatrixType ws;
		dmrgWaveStruct_.ws.toSparse(ws);
		SparseMatrixType we;
		dmrgWaveStruct_.we.toSparse(we);
		SparseMatrixType weT;
		transposeConjugate(weT,we);

		MatrixOrIdentityType wsRef2(wftOptions_.twoSiteDmrg && nip>volumeOfNk,ws);

		MatrixType temporary;
		buildTemporaryEnviron(temporary, psiSrc, i0, we);

		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		for (SizeType x=0;x<total;x++) {
			SizeType ip = 0;
			SizeType beta = 0;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			SizeType kp = 0;
			SizeType jp = 0;
			pack2.unpack(kp, jp, lrs.right().permutation(beta));
			for (SizeType k3=wsRef2.getRowPtr(ip);k3<wsRef2.getRowPtr(ip+1);k3++) {
				int ip2 = wsRef2.getColOrExit(k3);
				if (ip2<0) continue;
				SizeType alpha = dmrgWaveStruct_.lrs.left().permutationInverse(ip2+kp*nip2);
				psiDest.fastAccess(i0,x) += temporary(jp, alpha)*wsRef2.getValue(k3);
			}
		}
	}

	void buildTemporaryEnviron(MatrixType& temporary,
	                           const VectorWithOffsetType& psiV,
	                           SizeType i0,
	                           const SparseMatrixType& we) const
	{
		SizeType ni=dmrgWaveStruct_.lrs.left().size();
		temporary.resize(we.cols(), ni);
		temporary.setTo(0.0);

		PackIndicesType packSuper(ni);
		SizeType total = psiV.effectiveSize(i0);
		SizeType offset = psiV.offset(i0);

		for (SizeType x = 0; x < total; ++x) {
			SizeType alpha = 0;
			SizeType jp2 = 0;
			packSuper.unpack(alpha, jp2, dmrgWaveStruct_.lrs.super().permutation(x + offset));
			SizeType start = we.getRowPtr(jp2);
			SizeType end = we.getRowPtr(jp2+ 1);
			for (SizeType k = start; k < end; k++) {
				SizeType jp = we.getCol(k);
				temporary(jp, alpha) += psiV.fastAccess(i0, x)*we.getValue(k);
			}
		}
	}

	void transformVector2(VectorWithOffsetType& psiDest,
	                      const VectorWithOffsetType& psiSrc,
	                      const LeftRightSuperType& lrs,
	                      const VectorSizeType& nk) const
	{
		if (wftOptions_.twoSiteDmrg && wftOptions_.accel != WftOptions::ACCEL_PATCHES)
			return transformVector2FromInfinite(psiDest,psiSrc,lrs,nk);

		typename ProgramGlobals::DirectionEnum dir2 = ProgramGlobals::EXPAND_SYSTEM;
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			transformVectorParallel(psiDest,psiSrc,lrs,i0,nk,dir2);
		}
	}

	void transformVector2FromInfinite(VectorWithOffsetType& psiDest,
	                                  const VectorWithOffsetType& psiSrc,
	                                  const LeftRightSuperType& lrs,
	                                  const VectorSizeType& nk) const
	{
		PsimagLite::OstringStream msg;
		msg<<" Destination sectors "<<psiDest.sectors();
		msg<<" Source sectors "<<psiSrc.sectors();
		progress_.printline(msg,std::cout);
		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		SparseMatrixType we;
		dmrgWaveStruct_.we.toSparse(we);
		SparseMatrixType ws;
		dmrgWaveStruct_.ws.toSparse(ws);
		SparseMatrixType wsT;
		transposeConjugate(wsT,ws);
		VectorType psiV;
		for (SizeType srcI = 0; srcI < psiSrc.sectors(); ++srcI) {
			SizeType srcII = psiSrc.sector(srcI);
			psiSrc.extract(psiV,srcII);
			SizeType offset = psiSrc.offset(srcII);
			SizeType qn = psiSrc.qn(srcII);
			for (SizeType ii=0;ii<psiDest.sectors();ii++) {
				SizeType i0 = psiDest.sector(ii);
				if (qn != psiDest.qn(i0)) continue;
				SizeType start = psiDest.offset(i0);
				SizeType final = psiDest.effectiveSize(i0)+start;
				VectorType dest(final-start,0.0);
				if (srcI > 0) psiDest.extract(dest,i0);
				if (wftOptions_.accel == WftOptions::ACCEL_TEMP) {
					transformTemp2FromInfinite(dest,
					                           start,
					                           psiV,
					                           offset,
					                           lrs,
					                           nk,
					                           ws,
					                           we);
				} else if (wftOptions_.accel == WftOptions::ACCEL_BLOCKS) {
					wftAccelBlocks_.systemFromInfinite(psiDest,
					                                   i0,
					                                   psiSrc,
					                                   srcII,
					                                   lrs,
					                                   nk);
					continue;
				} else {
					tVector2FromInfinite(dest,start,psiV,offset,lrs,nk,wsT,we);
				}

				psiDest.setDataInSector(dest,i0);
			}
		}
	}

	void tVector2FromInfinite(VectorType& dest,
	                          SizeType destOffset,
	                          const VectorType& psiV,
	                          SizeType offset,
	                          const LeftRightSuperType& lrs,
	                          const VectorSizeType& nk,
	                          const SparseMatrixType& wsT,
	                          const SparseMatrixType& we) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();

		assert(nip==dmrgWaveStruct_.ws.cols());

		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		for (SizeType x=0;x<dest.size();x++) {
			SizeType isn,jen;
			pack1.unpack(isn,jen,(SizeType)lrs.super().permutation(x+destOffset));
			SizeType is,jpl;
			pack2.unpack(is,jpl,(SizeType)lrs.left().permutation(isn));

			dest[x] += createAux2bFromInfinite(psiV,offset,is,jpl,jen,wsT,we,nk);

		}
	}

	SparseElementType createAux2bFromInfinite(const VectorType& psiV,
	                                          SizeType offset,
	                                          SizeType is,
	                                          SizeType jpl,
	                                          SizeType jen,
	                                          const SparseMatrixType& wsT,
	                                          const SparseMatrixType& we,
	                                          const VectorSizeType& nk) const
	{
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		SparseElementType sum=0;
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk;

		MatrixOrIdentityType weRef(wftOptions_.twoSiteDmrg && ni>volumeOfNk,we);
		SizeType start = wsT.getRowPtr(is);
		SizeType end = wsT.getRowPtr(is+1);
		for (SizeType k2=weRef.getRowPtr(jen);k2<weRef.getRowPtr(jen+1);k2++) {
			int jpr = weRef.getColOrExit(k2);
			// jpr < 0 could be due to an m smaller than h, the Hilbert size of one site
			// this is checked against elsewhere
			assert(jpr >= 0);
			SizeType jp = dmrgWaveStruct_.lrs.right().permutationInverse(jpl + jpr*volumeOfNk);
			SparseElementType sum2 = 0;
			for (SizeType k = start;k < end;k++) {
				SizeType ip = wsT.getCol(k);
				SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ip + jp*nalpha);
				if (y < offset || y - offset >= psiV.size()) continue;
				sum2 += wsT.getValue(k)*psiV[y-offset]*weRef.getValue(k2);
			}

			sum += sum2;
		}

		return sum;
	}

	void transformTemp2FromInfinite(VectorType& dest,
	                                SizeType destOffset,
	                                const VectorType& psiV,
	                                SizeType offset,
	                                const LeftRightSuperType& lrs,
	                                const VectorSizeType& nk,
	                                const SparseMatrixType& ws,
	                                const SparseMatrixType& we) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();
		SizeType ni = dmrgWaveStruct_.lrs.right().size()/volumeOfNk;
		MatrixOrIdentityType weRef(wftOptions_.twoSiteDmrg && ni>volumeOfNk, we);

		assert(nip==dmrgWaveStruct_.ws.cols());

		MatrixType temporary;
		buildTemporarySystem(temporary, psiV, offset, ws);

		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		for (SizeType x=0;x<dest.size();x++) {
			SizeType isn,jen;
			pack1.unpack(isn,jen,(SizeType)lrs.super().permutation(x+destOffset));
			SizeType is = 0;
			SizeType jpl = 0;
			pack2.unpack(is,jpl,(SizeType)lrs.left().permutation(isn));
			SparseElementType sum = 0;
			for (SizeType k2=weRef.getRowPtr(jen);k2<weRef.getRowPtr(jen+1);k2++) {
				int jpr = weRef.getColOrExit(k2);
				// jpr < 0 could be due to an m smaller than h, the Hilbert size of one site
				// this is checked against elsewhere
				assert(jpr >= 0);
				SizeType jp = dmrgWaveStruct_.lrs.right().permutationInverse(jpl + jpr*volumeOfNk);
				sum += temporary(is, jp)*weRef.getValue(k2);
			}

			dest[x] += sum;
		}
	}

	void buildTemporarySystem(MatrixType& temporary,
	                          const VectorType& psiV,
	                          SizeType offset,
	                          const SparseMatrixType& ws) const
	{
		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		temporary.resize(nalpha, dmrgWaveStruct_.lrs.super().size()/nalpha);
		temporary.setTo(0.0);

		PackIndicesType packSuper(nalpha);
		SizeType total = psiV.size();

		for (SizeType y = 0; y < total; ++y) {
			SizeType ip = 0;
			SizeType jp = 0;
			packSuper.unpack(ip, jp, dmrgWaveStruct_.lrs.super().permutation(y + offset));
			SizeType start = ws.getRowPtr(ip);
			SizeType end = ws.getRowPtr(ip + 1);
			for (SizeType k = start;k < end; k++) {
				SizeType is = ws.getCol(k);
				temporary(is, jp) += ws.getValue(k)*psiV[y];
			}
		}
	}

	void transformVector1bounce(VectorWithOffsetType& psiDest,
	                            const VectorWithOffsetType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            const VectorSizeType& nk) const
	{
		SparseMatrixType ws;
		dmrgWaveStruct_.ws.toSparse(ws);
		MatrixOrIdentityType wsRef(wftOptions_.twoSiteDmrg, ws);
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			tVector1bounce(psiDest,psiSrc,lrs,i0,nk,wsRef);
		}
	}

	void tVector1bounce(VectorWithOffsetType& psiDest,
	                    const VectorWithOffsetType& psiSrc,
	                    const LeftRightSuperType& lrs,
	                    SizeType i0,
	                    const VectorSizeType& nk,
	                    const MatrixOrIdentityType& wsRef) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();
		PsimagLite::OstringStream msg;
		msg<<" We're bouncing on the right, so buckle up!";
		progress_.printline(msg,std::cout);

		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);

		SizeType nalpha=dmrgWaveStruct_.lrs.left().permutationInverse().size();
		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);

		SizeType nip2 = (wftOptions_.twoSiteDmrg) ? dmrgWaveStruct_.ws.cols() : nip;

		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;
			SizeType ip,beta,kp,jp;
			pack1.unpack(ip,beta,(SizeType)lrs.super().permutation(x+start));
			for (SizeType k=wsRef.getRowPtr(ip);k<wsRef.getRowPtr(ip+1);k++) {
				int ip2 = wsRef.getColOrExit(k);
				if (ip2 < 0) continue;
				pack2.unpack(kp,jp,(SizeType)lrs.right().permutation(beta));
				SizeType ipkp = dmrgWaveStruct_.lrs.left().permutationInverse(ip2 + kp*nip2);
				SizeType y = dmrgWaveStruct_.lrs.super().permutationInverse(ipkp + jp*nalpha);
				psiDest.fastAccess(i0,x) += psiSrc.slowAccess(y)*wsRef.getValue(k);
			}
		}
	}

	void transformVector2bounce(VectorWithOffsetType& psiDest,
	                            const VectorWithOffsetType& psiSrc,
	                            const LeftRightSuperType& lrs,
	                            const VectorSizeType& nk) const
	{
		SparseMatrixType we;
		dmrgWaveStruct_.we.toSparse(we);
		MatrixOrIdentityType weRef(wftOptions_.twoSiteDmrg, we);
		for (SizeType ii=0;ii<psiDest.sectors();ii++) {
			SizeType i0 = psiDest.sector(ii);
			tVector2bounce(psiDest,psiSrc,lrs,i0,nk,weRef);
		}
	}

	void tVector2bounce(VectorWithOffsetType& psiDest,
	                    const VectorWithOffsetType& psiSrc,
	                    const LeftRightSuperType& lrs,
	                    SizeType i0,
	                    const VectorSizeType& nk,
	                    const MatrixOrIdentityType& weRef) const
	{
		SizeType volumeOfNk = DmrgWaveStructType::volumeOf(nk);
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();

		PsimagLite::OstringStream msg;
		msg<<" We're bouncing on the left, so buckle up!";
		progress_.printline(msg,std::cout);

		assert(dmrgWaveStruct_.lrs.super().permutationInverse().size()==psiSrc.size());

		SizeType start = psiDest.offset(i0);
		SizeType total = psiDest.effectiveSize(i0);
		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);

		for (SizeType x=0;x<total;x++) {
			psiDest.fastAccess(i0,x) = 0.0;

			SizeType ip,alpha,kp,jp;
			pack1.unpack(alpha,jp,(SizeType)lrs.super().permutation(x+start));
			pack2.unpack(ip,kp,(SizeType)lrs.left().permutation(alpha));

			for (SizeType k=weRef.getRowPtr(jp);k<weRef.getRowPtr(jp+1);k++) {
				int jp2 = weRef.getColOrExit(k);
				if (jp2 < 0) continue;
				SizeType kpjp = dmrgWaveStruct_.lrs.right().
				        permutationInverse(kp + jp2*volumeOfNk);

				SizeType y = dmrgWaveStruct_.lrs.super().
				        permutationInverse(ip + kpjp*nip);
				psiDest.fastAccess(i0,x) += psiSrc.slowAccess(y) * weRef.getValue(k);
			}
		}
	}

	SizeType findIold(const VectorWithOffsetType& psiSrc,
	                  const VectorWithOffsetType& psiDest,
	                  SizeType iNew) const
	{
		SizeType sectors = psiSrc.sectors();
		for (SizeType i = 0; i < sectors; ++i) {
			if (psiSrc.qn(i) == psiDest.qn(iNew))
				return psiSrc.sector(i);
		}

		err("WaveFunctionTransfLocal::findIold(): Cannot find sector in old vector\n");
		throw PsimagLite::RuntimeError("UNREACHABLE\n");
	}

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptions& wftOptions_;
	WftAccelBlocksType wftAccelBlocks_;
	PsimagLite::ProgressIndicator progress_;
}; // class WaveFunctionTransfLocals
} // namespace Dmrg

/*@}*/
#endif

