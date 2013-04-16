/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

namespace Dmrg {

template<typename TargettingParamsType,
		 typename ModelType,
		 typename WaveFunctionTransfType,
		 typename LanczosSolverType,
		 typename VectorWithOffsetType>
class TimeVectorsSuzukiTrotter : public  TimeVectorsBase<
		TargettingParamsType,
		ModelType,
		WaveFunctionTransfType,
		LanczosSolverType,
		VectorWithOffsetType> {

	typedef TimeVectorsBase<TargettingParamsType,ModelType,WaveFunctionTransfType,LanczosSolverType,VectorWithOffsetType> BaseType;
	typedef typename BaseType::PairType PairType;
	typedef typename TargettingParamsType::RealType RealType;
	typedef typename TargettingParamsType::SparseMatrixType SparseMatrixType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType BlockType;
	typedef MatrixOrIdentity<SparseMatrixType> MatrixOrIdentityType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef std::vector<ComplexOrRealType> VectorComplexOrRealType;
	typedef std::vector<RealType> VectorRealType;
	typedef VectorComplexOrRealType TargetVectorType;

public:

	TimeVectorsSuzukiTrotter(RealType& currentTime,
							 const TargettingParamsType& tstStruct,
							 const VectorRealType& times,
							 std::vector<VectorWithOffsetType>& targetVectors,
							 const ModelType& model,
							 const WaveFunctionTransfType& wft,
							 const LeftRightSuperType& lrs,
							 const RealType& E0)
		: progress_("TimeVectorsSuzukiTrotter",0),
		  currentTime_(currentTime),
		  tstStruct_(tstStruct),
		  times_(times),
		  targetVectors_(targetVectors),
		  model_(model),
		  wft_(wft),
		  lrs_(lrs),
		  E0_(E0),
		  twoSiteDmrg_(wft_.twoSiteDmrg())
	{}

	virtual void calcTimeVectors(const PairType& startEnd,
	                             RealType Eg,
								 const VectorWithOffsetType& phi,
								 size_t systemOrEnviron)
	{
		static bool firstcall = true;
		std::ostringstream msg;
		msg<<"EXPERIMENTAL: using SuzukiTrotter";

		RealType norma = std::norm(phi);
		if (norma<1e-10) return;
		msg<<" Norm of phi= "<<norma;
		progress_.printline(msg,std::cout);

		// set non-zero sectors
		targetVectors_[0] = phi;

		for (size_t i=1;i<times_.size();i++)
			if (targetVectors_[i].size()==0)
				targetVectors_[i] = phi;

		if (firstcall) {
			firstcall=false;
			return;
		}
		// skip odd links if expanding system and
		// skip even links if expanding environ
		size_t lastIndexLeft = lrs_.left().block().size();
		assert(lastIndexLeft>0);
		lastIndexLeft--;
		bool oddLink = (lrs_.left().block()[lastIndexLeft] & 1);
		bool b1 = (oddLink && systemOrEnviron==ProgramGlobals::EXPAND_SYSTEM);
		bool b2 = (!oddLink && systemOrEnviron==ProgramGlobals::EXPAND_ENVIRON);
		if (b2 && lrs_.left().block().size()==1)
			b2=false;

		wftAll(lrs_.left().block()[lastIndexLeft]);

		if (b1 || b2) return;

		bool areAllLinksSeen = allLinksSeen();
		std::ostringstream msg2;
		msg2<<"LINKS SEEN ";
		for (size_t i=0;i<linksSeen_.size();i++)
			msg2<<linksSeen_[i]<<" ";
		progress_.printline(msg2,std::cout);

		if (!areAllLinksSeen) {
			linksSeen_.push_back(lastIndexLeft);
		} else {
			std::ostringstream msg3;
			msg3<<"ALL LINKS SEEN";
			progress_.printline(msg3,std::cout);
			return;
		}

		for (size_t i=startEnd.first+1;i<startEnd.second;i++) {
			VectorWithOffsetType src = targetVectors_[i];
			// Only time differences here (i.e. times_[i] not times_[i]+currentTime_)
			calcTargetVector(targetVectors_[i],Eg,src,systemOrEnviron,times_[i]);
			assert(targetVectors_[i].size()==targetVectors_[0].size());
		}
	}

	virtual void timeHasAdvanced()
	{
		linksSeen_.clear();
		std::ostringstream msg;
		msg<<"ALL LINKS CLEARED";
		progress_.printline(msg,std::cout);
	}

private:

	bool allLinksSeen() const
	{
		size_t nsites = model_.geometry().numberOfSites();
		assert(nsites>0);
		for (size_t i=0;i<nsites-1;i++) {
			std::vector<size_t>::const_iterator it = find(linksSeen_.begin(),linksSeen_.end(),i);
			if (it == linksSeen_.end()) return false;
		}
		return true;
	}

	void wftAll(size_t site)
	{
		size_t nk = model_.hilbertSize(site);
		for (size_t i=1;i<times_.size();i++)
			wftOne(i,nk);
	}

	void wftOne(size_t i,size_t nk)
	{
		if (targetVectors_[i].size()==0)
			targetVectors_[i] = targetVectors_[0];

		VectorWithOffsetType phiNew;
		phiNew.populateSectors(lrs_.super());

		// OK, now that we got the partition number right, let's wft:
		wft_.setInitialVector(phiNew,targetVectors_[i],lrs_,nk); // generalize for su(2)
		phiNew.collapseSectors();
		targetVectors_[i]=phiNew;
	}

	void calcTargetVector(VectorWithOffsetType& target,
						  RealType Eg,
						  const VectorWithOffsetType& phi,
						  size_t systemOrEnviron,
						  const RealType& time)
	{
		for (size_t ii=0;ii<phi.sectors();ii++) {
			size_t i0 = phi.sector(ii);
			size_t total = phi.effectiveSize(i0);
			TargetVectorType result(total,0.0);
			calcTimeVectorsSuzukiTrotter(result,Eg,phi,systemOrEnviron,i0,time);
			//NOTE: targetVectors_[0] = exp(iHt) |phi>
			target.setDataInSector(result,i0);
		}
	}

	void calcTimeVectorsSuzukiTrotter(TargetVectorType& result,
									  RealType Eg,
									  const VectorWithOffsetType& phi,
									  size_t systemOrEnviron,
									  size_t i0,
									  const RealType& time) const
	{
		size_t offset = phi.offset(i0);
		TargetVectorType phi0(result.size());
		phi.extract(phi0,i0);

		// NOTE: result =  exp(iHt) |phi0>
		size_t ns = lrs_.left().size();
		PackIndicesType packSuper(ns);
		std::vector<size_t> block(2);

		size_t lastIndexLeft = lrs_.left().block().size();
		assert(lastIndexLeft>0);
		lastIndexLeft--;
		block[0]=lrs_.left().block()[lastIndexLeft];
		block[1]=lrs_.right().block()[0];

		MatrixComplexOrRealType m;
		getMatrix(m,systemOrEnviron,block,time);

		SparseMatrixType transformS = wft_.stackTransform(ProgramGlobals::SYSTEM);
		SparseMatrixType transformST;
		transposeConjugate(transformST,transformS);

		SparseMatrixType transformE = wft_.stackTransform(ProgramGlobals::ENVIRON);
		SparseMatrixType transformET;
		transposeConjugate(transformET,transformE);

		size_t hilbertSize = model_.hilbertSize(block[0]);
		if (systemOrEnviron==ProgramGlobals::EXPAND_SYSTEM && lrs_.right().size()==hilbertSize) {
			transformE.makeDiagonal(hilbertSize,1);
			transformET.makeDiagonal(hilbertSize,1);
		}

		for (size_t i=0;i<phi0.size();i++) {
			size_t xp=0,yp=0;
			packSuper.unpack(xp,yp,lrs_.super().permutation(i+offset));
			if (systemOrEnviron==ProgramGlobals::EXPAND_SYSTEM) {
				timeVectorSystem(result,phi0,xp,yp,packSuper,block,m,i,offset,transformE,transformET);
			} else {
				timeVectorEnviron(result,phi0,xp,yp,packSuper,block,m,i,offset,transformS,transformST);
			}
		}
	}

	void timeVectorSystem(TargetVectorType& result,
									   const TargetVectorType& phi0,
									   size_t xp,
									   size_t yp,
									   const PackIndicesType& packSuper,
									   const BlockType& block,
									   const MatrixComplexOrRealType& m,
									   size_t i,
									   size_t offset,
									   const SparseMatrixType& transform,
									   const SparseMatrixType& transformT) const
	{
		std::vector<size_t> iperm;
		suzukiTrotterPerm(iperm,block);

		const LeftRightSuperType& oldLrs = lrs_;
		size_t hilbertSize = model_.hilbertSize(block[0]);
		size_t ns = lrs_.left().size();
		size_t nx = ns/hilbertSize;
		PackIndicesType packLeft(nx);
		PackIndicesType packRight(hilbertSize);

		if (!twoSiteDmrg_) {
			assert(transform.col()==lrs_.right().size());
			assert(transform.row()==oldLrs.right().permutationInverse().size());
		}

		MatrixOrIdentityType transformT1(!twoSiteDmrg_,transformT);
		MatrixOrIdentityType transform1(!twoSiteDmrg_,transform);
		for (size_t k=transformT1.getRowPtr(yp);k<transformT1.getRowPtr(yp+1);k++) {
			size_t x1=0,x2p=0;
			packLeft.unpack(x1,x2p,lrs_.left().permutation(xp));

			size_t yfull = transformT1.getCol(k);
			size_t y1p=0,y2=0;
			packRight.unpack(y1p,y2,oldLrs.right().permutation(yfull));

			for (size_t x2=0;x2<hilbertSize;x2++) {
				for (size_t y1=0;y1<hilbertSize;y1++) {
					size_t yfull2 = packRight.pack(y1,y2,oldLrs.right().permutationInverse());
					for (size_t k2=transform1.getRowPtr(yfull2);k2<transform1.getRowPtr(yfull2+1);k2++) {
						size_t y = transform1.getCol(k2);
						size_t x = packLeft.pack(x1,x2,lrs_.left().permutationInverse());
						size_t j = packSuper.pack(x,y,lrs_.super().permutationInverse());
						ComplexOrRealType tmp = m(iperm[x2+y1*hilbertSize],iperm[x2p+y1p*hilbertSize]);
						if (std::norm(tmp)==0) continue;
						assert(j>=offset && j<offset+phi0.size());
						result[j-offset] += tmp*phi0[i]*transformT1.getValue(k)*transform1.getValue(k2);
					}
				}
			}
		}
	}

	void timeVectorEnviron(TargetVectorType& result,
										const TargetVectorType& phi0,
										size_t xp,
										size_t yp,
										const PackIndicesType& packSuper,
										const BlockType& block,
										const MatrixComplexOrRealType& m,
										size_t i,
										size_t offset,
										const SparseMatrixType& transform,
										const SparseMatrixType& transformT) const
	{
		std::vector<size_t> iperm;
		suzukiTrotterPerm(iperm,block);

		const LeftRightSuperType& oldLrs = lrs_;
		size_t hilbertSize = model_.hilbertSize(block[0]);
		size_t ns = oldLrs.left().permutationInverse().size();
		size_t nx = ns/hilbertSize;
		PackIndicesType packLeft(nx);
		PackIndicesType packRight(hilbertSize);

		if (!twoSiteDmrg_) {
			assert(transform.col()==lrs_.left().size());
			assert(transform.row()==oldLrs.left().permutationInverse().size());
		}

		MatrixOrIdentityType transformT1(!twoSiteDmrg_,transformT);
		MatrixOrIdentityType transform1(!twoSiteDmrg_,transform);

		for (size_t k=transformT1.getRowPtr(xp);k<transformT1.getRowPtr(xp+1);k++) {
			size_t xfull = transformT1.getCol(k);
			size_t x1=0,x2p=0;
			packLeft.unpack(x1,x2p,oldLrs.left().permutation(xfull));
			assert(x2p<hilbertSize);

			size_t y1p=0,y2=0;
			packRight.unpack(y1p,y2,lrs_.right().permutation(yp));

			for (size_t x2=0;x2<hilbertSize;x2++) {
				for (size_t y1=0;y1<hilbertSize;y1++) {
					size_t xfull2 = packLeft.pack(x1,x2,oldLrs.left().permutationInverse());
					for (size_t k2=transform1.getRowPtr(xfull2);k2<transform1.getRowPtr(xfull2+1);k2++) {
						size_t x = transform1.getCol(k2);
						size_t y = packRight.pack(y1,y2,lrs_.right().permutationInverse());
						size_t j = packSuper.pack(x,y,lrs_.super().permutationInverse());
						ComplexOrRealType tmp = m(iperm[x2+y1*hilbertSize],iperm[x2p+y1p*hilbertSize]);
						if (std::norm(tmp)==0) continue;
						assert(j>=offset && j<offset+phi0.size());
						//								if (j<offset || j>=offset+phi0.size()) continue;
						result[j-offset] += tmp*phi0[i]*transformT1.getValue(k)*transform1.getValue(k2);
					}
				}
			}
		}
	}

	void suzukiTrotterPerm(std::vector<size_t>& iperm,const std::vector<size_t>& block) const
	{
		typename ModelType::HilbertBasisType  basis;
		std::vector<size_t> q;
		model_.setNaturalBasis(basis,q,block);
		iperm.resize(basis.size());
		for (size_t i=0;i<basis.size();i++) {
			size_t ket = basis[i];
			int jj = PsimagLite::isInVector(basis,ket);
			assert(jj>=0);
			iperm[ket] = jj;
		}
	}

	void getMatrix(MatrixComplexOrRealType& m,size_t systemOrEnviron,const BlockType& block,const RealType& time) const
	{
		assert(block.size()==2);
		SparseMatrixType hmatrix;
		RealType factorForDiagonals = (systemOrEnviron==ProgramGlobals::EXPAND_SYSTEM) ? 1.0 : 0.0;
		if (systemOrEnviron==ProgramGlobals::EXPAND_ENVIRON && block[0]==0) factorForDiagonals = 1.0;

		if (fabs(factorForDiagonals)>1e-6) {
			std::ostringstream msg;
			msg<<"LINKS factors="<<factorForDiagonals<<" added for diagonals on sites ";
			msg<<block[0]<<" and "<<block[1];
			progress_.printline(msg,std::cout);
		}
		model_.hamiltonianOnLink(hmatrix,block,currentTime_,factorForDiagonals);
		crsMatrixToFullMatrix(m,hmatrix);
		assert(isHermitian(m));
		m *= (-time);
		exp(m);
	}

	PsimagLite::ProgressIndicator progress_;
	RealType& currentTime_;
	const TargettingParamsType& tstStruct_;
	const VectorRealType& times_;
	std::vector<VectorWithOffsetType>& targetVectors_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
	const RealType& E0_;
	bool twoSiteDmrg_;
	std::vector<size_t> linksSeen_;
}; //class TimeVectorsSuzukiTrotter
} // namespace Dmrg
/*@}*/
#endif
