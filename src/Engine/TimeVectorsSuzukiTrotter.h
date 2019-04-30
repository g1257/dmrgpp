/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

/*! \file TimeVectorsSuzukiTrotter.h
 *
 *
 */
#ifndef TIME_VECTORS_SUZUKI_TROTTER
#define TIME_VECTORS_SUZUKI_TROTTER
#include <iostream>
#include "TimeVectorsBase.h"
#include "VectorWithOffsets.h"
#include "MatrixOrIdentity.h"
#include "Sort.h"
#include "Utils.h"

namespace Dmrg {

template<typename TargetParamsType,
         typename ModelType,
         typename WaveFunctionTransfType,
         typename LanczosSolverType,
         typename VectorWithOffsetType>
class TimeVectorsSuzukiTrotter : public  TimeVectorsBase<
        TargetParamsType,
        ModelType,
        WaveFunctionTransfType,
        LanczosSolverType,
        VectorWithOffsetType> {

	typedef TimeVectorsBase<TargetParamsType,
	ModelType,
	WaveFunctionTransfType,
	LanczosSolverType,
	VectorWithOffsetType> BaseType;
	typedef typename BaseType::PairType PairType;
	typedef typename TargetParamsType::RealType RealType;
	typedef typename TargetParamsType::SparseMatrixType SparseMatrixType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType BlockType;
	typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef VectorComplexOrRealType TargetVectorType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type
	VectorVectorWithOffsetType;
	typedef typename ModelType::HilbertBasisType HilbertBasisType;
	typedef typename ModelType::HilbertBasisType::value_type HilbertStateType;

public:

	TimeVectorsSuzukiTrotter(const SizeType& currentTimeStep,
	                         const VectorRealType& times,
	                         VectorVectorWithOffsetType& targetVectors,
	                         const ModelType& model,
	                         const WaveFunctionTransfType& wft,
	                         const LeftRightSuperType& lrs,
	                         const RealType& E0)
	    : BaseType(times),
	      progress_("TimeVectorsSuzukiTrotter"),
	      currentTimeStep_(currentTimeStep),
	      times_(times),
	      targetVectors_(targetVectors),
	      model_(model),
	      wft_(wft),
	      lrs_(lrs),
	      E0_(E0),
	      twoSiteDmrg_(wft_.options().twoSiteDmrg)
	{}

	virtual void calcTimeVectors(const PairType& startEnd,
	                             RealType Eg,
	                             const VectorWithOffsetType& phi,
	                             ProgramGlobals::DirectionEnum systemOrEnviron,
	                             bool allOperatorsApplied,
	                             const VectorSizeType& block,
	                             const TargetParamsType&)
	{
		PsimagLite::OstringStream msg;
		msg<<"EXPERIMENTAL: using SuzukiTrotter";

		RealType norma = norm(phi);
		if (norma<1e-10) return;
		msg<<" Norm of phi= "<<norma;
		progress_.printline(msg,std::cout);

		// set non-zero sectors
		targetVectors_[0] = phi;

		bool returnFlag = false;
		for (SizeType i=1;i<times_.size();i++) {
			if (targetVectors_[i].size()==0 || !allOperatorsApplied) {
				targetVectors_[i] = phi;
				returnFlag = true;
			}
		}

		if (returnFlag) return;

		// skip odd links if expanding system and
		// skip even links if expanding environ
		SizeType sitesPerBlock = model_.params().sitesPerBlock;
		SizeType lastIndexLeft = lrs_.left().block().size();
		assert(lastIndexLeft>0);
		lastIndexLeft--;
		SizeType site = static_cast<SizeType>(lrs_.left().block()[lastIndexLeft]/
		                                      sitesPerBlock);
		bool oddLink = (site & 1);
		bool b1 = (oddLink && systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		bool b2 = (!oddLink && systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON);
		if (b2 && lrs_.left().block().size() == sitesPerBlock) b2=false;

		wftAll(block);

		if (b1 || b2) return;

		bool areAllLinksSeen = allLinksSeen();
		PsimagLite::OstringStream msg2;
		msg2<<"LINKS SEEN ";
		for (SizeType i=0;i<linksSeen_.size();i++)
			msg2<<linksSeen_[i]<<" ";
		progress_.printline(msg2,std::cout);

		if (!areAllLinksSeen) {
			for (SizeType i = 0; i < block.size(); ++i)
				linksSeen_.push_back(lastIndexLeft+i);
		} else {
			PsimagLite::OstringStream msg3;
			msg3<<"ALL LINKS SEEN";
			progress_.printline(msg3,std::cout);
			return;
		}

		SparseMatrixType transformS;
		wft_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toSparse(transformS);
		SparseMatrixType transformST;
		transposeConjugate(transformST,transformS);

		SparseMatrixType transformE;
		wft_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toSparse(transformE);
		SparseMatrixType transformET;
		transposeConjugate(transformET,transformE);

		SizeType hilbertSize = model_.hilbertSize(block[0]);
		if (systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM &&
		    lrs_.right().size()==hilbertSize) {
			transformE.makeDiagonal(hilbertSize,1);
			transformET.makeDiagonal(hilbertSize,1);
		}

		if (systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON &&
		    lrs_.left().size()==hilbertSize) {
			transformS.makeDiagonal(hilbertSize,1);
			transformST.makeDiagonal(hilbertSize,1);
		}

		for (SizeType i=startEnd.first+1;i<startEnd.second;i++) {
			VectorWithOffsetType src = targetVectors_[i];
			// Only time differences here (i.e. times_[i] not times_[i]+currentTime_)
			calcTargetVector(targetVectors_[i],
			                 Eg,
			                 src,
			                 systemOrEnviron,
			                 times_[i],
			                 transformS,
			                 transformST,
			                 transformE,
			                 transformET);
			assert(targetVectors_[i].size()==targetVectors_[0].size());
		}
	}

	RealType time() const { return currentTimeStep_*BaseType::tau(); }

	virtual void timeHasAdvanced()
	{
		linksSeen_.clear();
		PsimagLite::OstringStream msg;
		msg<<"ALL LINKS CLEARED";
		progress_.printline(msg,std::cout);
	}

private:

	bool allLinksSeen() const
	{
		SizeType nsites = model_.geometry().numberOfSites();
		assert(nsites>0);
		SizeType start = model_.params().sitesPerBlock - 1;
		for (SizeType i=start;i<nsites-1;i++) {
			VectorSizeType::const_iterator it =
			        find(linksSeen_.begin(),linksSeen_.end(),i);
			if (it == linksSeen_.end()) return false;
		}
		return true;
	}

	void wftAll(const VectorSizeType& block)
	{
		for (SizeType i=1;i<times_.size();i++)
			wftOne(i,block);
	}

	void wftOne(SizeType i,const VectorSizeType& block)
	{
		VectorWithOffsetType phiNew = targetVectors_[0];

		// OK, now that we got the partition number right, let's wft:
		VectorSizeType nk;
		setNk(nk,block);
		// generalize for su(2)
		wft_.setInitialVector(phiNew,targetVectors_[i],lrs_,nk);
		phiNew.collapseSectors();
		assert(norm(phiNew)>1e-6);
		targetVectors_[i]=phiNew;
	}

	void calcTargetVector(VectorWithOffsetType& target,
	                      RealType Eg,
	                      const VectorWithOffsetType& phi,
	                      const ProgramGlobals::DirectionEnum systemOrEnviron,
	                      const RealType& time,
	                      const SparseMatrixType& S,
	                      const SparseMatrixType& ST,
	                      const SparseMatrixType& E,
	                      const SparseMatrixType& ET)
	{
		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i0 = phi.sector(ii);
			SizeType total = phi.effectiveSize(i0);
			TargetVectorType result(total,0.0);
			calcTimeVectorsSuzukiTrotter(result,
			                             Eg,
			                             phi,
			                             systemOrEnviron,
			                             i0,
			                             time,
			                             S,
			                             ST,
			                             E,
			                             ET);
			//NOTE: targetVectors_[0] = exp(iHt) |phi>
			target.setDataInSector(result,i0);
		}
	}

	void calcTimeVectorsSuzukiTrotter(TargetVectorType& result,
	                                  RealType,
	                                  const VectorWithOffsetType& phi,
	                                  const ProgramGlobals::DirectionEnum systemOrEnviron,
	                                  SizeType i0,
	                                  const RealType& time,
	                                  const SparseMatrixType& transformS,
	                                  const SparseMatrixType& transformST,
	                                  const SparseMatrixType& transformE,
	                                  const SparseMatrixType& transformET) const
	{
		SizeType offset = phi.offset(i0);
		TargetVectorType phi0(result.size());
		phi.extract(phi0,i0);

		// NOTE: result =  exp(iHt) |phi0>
		SizeType ns = lrs_.left().size();
		PackIndicesType packSuper(ns);

		VectorSizeType block;
		calcBlock(block);

		MatrixComplexOrRealType m;
		getMatrix(m,systemOrEnviron,block,time);

		VectorSizeType iperm;
		suzukiTrotterPerm(iperm,block);
		for (SizeType i=0;i<phi0.size();i++) {
			SizeType xp=0,yp=0;
			packSuper.unpack(xp,yp,lrs_.super().permutation(i+offset));
			if (systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				timeVectorSystem(result,
				                 phi0,
				                 xp,
				                 yp,
				                 packSuper,
				                 block,
				                 m,
				                 i,
				                 offset,
				                 transformE,
				                 transformET,
				                 iperm);
			} else {
				timeVectorEnviron(result,
				                  phi0,
				                  xp,
				                  yp,
				                  packSuper,
				                  block,
				                  m,
				                  i,
				                  offset,
				                  transformS,
				                  transformST,
				                  iperm);
			}
		}
	}

	void timeVectorSystem(TargetVectorType& result,
	                      const TargetVectorType& phi0,
	                      SizeType xp,
	                      SizeType yp,
	                      const PackIndicesType& packSuper,
	                      const BlockType& block,
	                      const MatrixComplexOrRealType& m,
	                      SizeType i,
	                      SizeType offset,
	                      const SparseMatrixType& transform,
	                      const SparseMatrixType& transformT,
	                      const VectorSizeType& iperm) const
	{
		const LeftRightSuperType& oldLrs = lrs_;
		SizeType hilbertSize = model_.hilbertSize(block[0]);
		SizeType ns = lrs_.left().size();
		SizeType nx = ns/hilbertSize;
		PackIndicesType packLeft(nx);
		PackIndicesType packRight(hilbertSize);

		if (!twoSiteDmrg_) {
			assert(transform.cols()==lrs_.right().size());
			assert(transform.rows()==oldLrs.right().permutationInverse().size());
		}

		MatrixOrIdentityType transformT1(!twoSiteDmrg_,transformT);
		MatrixOrIdentityType transform1(!twoSiteDmrg_,transform);
		for (SizeType k=transformT1.getRowPtr(yp);k<transformT1.getRowPtr(yp+1);k++) {
			SizeType x1=0,x2p=0;
			packLeft.unpack(x1,x2p,lrs_.left().permutation(xp));
			int yfull = transformT1.getColOrExit(k);
			if (yfull<0) yfull = yp;
			SizeType y1p=0,y2=0;
			packRight.unpack(y1p,y2,oldLrs.right().permutation(yfull));
			for (SizeType x2=0;x2<hilbertSize;x2++) {
				for (SizeType y1=0;y1<hilbertSize;y1++) {
					SizeType yfull2 = packRight.pack(y1,
					                                 y2,
					                                 oldLrs.right().permutationInverse());
					for (SizeType k2=transform1.getRowPtr(yfull2);
					     k2<transform1.getRowPtr(yfull2+1);
					     k2++) {
						int y = transform1.getColOrExit(k2);
						if (y<0) y = yfull2;
						SizeType x = packLeft.pack(x1,
						                           x2,
						                           lrs_.left().permutationInverse());
						SizeType j = packSuper.pack(x,
						                            y,
						                            lrs_.super().permutationInverse());
						ComplexOrRealType tmp = m(iperm[x2+y1*hilbertSize],
						        iperm[x2p+y1p*hilbertSize]);
						if (PsimagLite::norm(tmp)<1e-12) continue;
						if (j<offset || j >= offset+phi0.size())
							throw PsimagLite::RuntimeError("j out of bounds\n");
						result[j-offset] += tmp*phi0[i]*transformT1.getValue(k)*
						        transform1.getValue(k2);
					}
				}
			}
		}
	}

	void timeVectorEnviron(TargetVectorType& result,
	                       const TargetVectorType& phi0,
	                       SizeType xp,
	                       SizeType yp,
	                       const PackIndicesType& packSuper,
	                       const BlockType& block,
	                       const MatrixComplexOrRealType& m,
	                       SizeType i,
	                       SizeType offset,
	                       const SparseMatrixType& transform,
	                       const SparseMatrixType& transformT,
	                       const VectorSizeType& iperm) const
	{
		const LeftRightSuperType& oldLrs = lrs_;
		SizeType hilbertSize = model_.hilbertSize(block[0]);
		SizeType ns = oldLrs.left().permutationInverse().size();
		SizeType nx = ns/hilbertSize;
		PackIndicesType packLeft(nx);
		PackIndicesType packRight(hilbertSize);

		if (!twoSiteDmrg_) {
			assert(transform.cols()==lrs_.left().size());
			assert(transform.rows()==oldLrs.left().permutationInverse().size());
		}

		MatrixOrIdentityType transformT1(!twoSiteDmrg_,transformT);
		MatrixOrIdentityType transform1(!twoSiteDmrg_,transform);

		for (SizeType k=transformT1.getRowPtr(xp);k<transformT1.getRowPtr(xp+1);k++) {
			int xfull = transformT1.getColOrExit(k);
			if (xfull<0) xfull = xp;
			SizeType x1=0,x2p=0;
			packLeft.unpack(x1,x2p,oldLrs.left().permutation(xfull));
			assert(x2p<hilbertSize);
			SizeType y1p=0,y2=0;
			packRight.unpack(y1p,y2,lrs_.right().permutation(yp));
			for (SizeType x2=0;x2<hilbertSize;x2++) {
				for (SizeType y1=0;y1<hilbertSize;y1++) {
					SizeType xfull2 = packLeft.pack(x1,
					                                x2,
					                                oldLrs.left().permutationInverse());
					for (SizeType k2=transform1.getRowPtr(xfull2);
					     k2<transform1.getRowPtr(xfull2+1);
					     k2++) {
						int x = transform1.getColOrExit(k2);
						if (x<0) x = xfull2;
						SizeType y = packRight.pack(y1,
						                            y2,
						                            lrs_.right().permutationInverse());
						SizeType j = packSuper.pack(x,
						                            y,
						                            lrs_.super().permutationInverse());

						ComplexOrRealType tmp = m(iperm[x2+y1*hilbertSize],
						        iperm[x2p+y1p*hilbertSize]);
						if (PsimagLite::norm(tmp)<1e-12) continue;
						if (j < offset || j >= offset+phi0.size())
							throw PsimagLite::RuntimeError("j out of bounds (environ)\n");

						result[j-offset] += tmp*phi0[i]*transformT1.getValue(k)*
						        transform1.getValue(k2);
					}
				}
			}
		}
	}

	void suzukiTrotterPerm(VectorSizeType&,
	                       const VectorSizeType&) const
	{
		err("suzukiTrotter no longer supported (sorry!)\n");
	}

	void getMatrix(MatrixComplexOrRealType& m,
	               const ProgramGlobals::DirectionEnum systemOrEnviron,
	               const BlockType& block,
	               const RealType& time) const
	{
		SparseMatrixType hmatrix;
		RealType factorForDiagonals =
		        (systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? 1.0 : 0.0;
		if (systemOrEnviron==ProgramGlobals::DirectionEnum::EXPAND_ENVIRON && block[0] == 0)
			factorForDiagonals = 1.0;

		if (fabs(factorForDiagonals)>1e-6) {
			PsimagLite::OstringStream msg;
			msg<<"LINKS factors="<<factorForDiagonals;
			msg<<" added for diagonals on sites ";
			for (SizeType i = 0; i < block.size(); ++i)
				msg<<block[i]<<" ";
			progress_.printline(msg,std::cout);
		}

		err("ST not supported\n");
		// model_.hamiltonianOnLink(hmatrix,block,currentTime_,factorForDiagonals);
		crsMatrixToFullMatrix(m,hmatrix);
		assert(isHermitian(m));
		m *= (-time);
		exp(m);
	}

	void setNk(VectorSizeType& nk, const  VectorSizeType& block) const
	{
		for (SizeType i=0;i<block.size();i++)
			nk.push_back(model_.hilbertSize(block[i]));
	}

	void calcBlock(VectorSizeType& block) const
	{
		const VectorSizeType& blockLeft = lrs_.left().block();
		const VectorSizeType& blockRight = lrs_.right().block();
		SizeType smax = blockLeft[blockLeft.size()-1];
		SizeType emin = blockRight[0];

		for (SizeType i = 0; i < blockLeft.size(); ++i) {
			SizeType ind = blockLeft[i];
			for (SizeType j=0; j < blockRight.size(); ++j) {
				SizeType jnd = blockRight[j];
				assert(ind != jnd);
				if (!model_.geometry().connected(smax, emin, ind, jnd))
					continue;
				block.push_back(ind);
				block.push_back(jnd);
			}
		}
		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(block.size());
		sort.sort(block,iperm);
	}

	PsimagLite::ProgressIndicator progress_;
	const SizeType& currentTimeStep_;
	const VectorRealType& times_;
	VectorVectorWithOffsetType& targetVectors_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
	RealType E0_;
	bool twoSiteDmrg_;
	VectorSizeType linksSeen_;
}; //class TimeVectorsSuzukiTrotter
} // namespace Dmrg
/*@}*/
#endif

