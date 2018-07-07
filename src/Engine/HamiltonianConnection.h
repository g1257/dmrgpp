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
/** \file HamiltonianConnection.h
*/

#ifndef HAMILTONIAN_CONNECTION_H
#define HAMILTONIAN_CONNECTION_H

#include "LinkProductStruct.h"
#include "CrsMatrix.h"
#include "Concurrency.h"
#include <cassert>
#include "ProgramGlobals.h"
#include "HamiltonianAbstract.h"

namespace Dmrg {

template<typename GeometryType,typename ModelHelperType,typename LinkProductType>
class HamiltonianConnection {

public:

	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef LinkProductStruct<SparseElementType> LinkProductStructType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Concurrency ConcurrencyType;
	typedef HamiltonianAbstract HamiltonianAbstractType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	HamiltonianConnection(const GeometryType& geometry,
	                      const ModelHelperType& modelHelper,
	                      const LinkProductStructType* lps = 0,
	                      typename PsimagLite::Vector<SparseElementType>::Type* x = 0,
	                      const typename PsimagLite::Vector<SparseElementType>::Type* y = 0)
	    : geometry_(geometry),
	      modelHelper_(modelHelper),
	      lps_(*lps),x_(*x),y_(*y),
	      systemBlock_(modelHelper.leftRightSuper().left().block()),
	      envBlock_(modelHelper.leftRightSuper().right().block()),
	      smax_(*std::max_element(systemBlock_.begin(),systemBlock_.end())),
	      emin_(*std::min_element(envBlock_.begin(),envBlock_.end())),
	      xtemp_(ConcurrencyType::storageSize(ConcurrencyType::codeSectionParams.npthreads)),
	      total_(0),
	      hamAbstract_(modelHelper_.leftRightSuper().super().block())
	{}

	~HamiltonianConnection()
	{
		SizeType n = garbage_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
	}

	bool compute(SizeType x,
	             SparseMatrixType* matrixBlock,
	             LinkProductStructType* lps,
	             SizeType& total) const
	{
		bool flag = false;

		const VectorSizeType& hItems = hamAbstract_.item(x);
		if (!geometry_.connected(smax_, emin_, hItems)) return flag;

		ProgramGlobals::ConnectionEnum type = geometry_.connectionKind(smax_, hItems);

		if (type == ProgramGlobals::SYSTEM_SYSTEM ||
		        type == ProgramGlobals::ENVIRON_ENVIRON) return flag;

		SparseMatrixType mBlock;

		AdditionalDataType additionalData;
		VectorSizeType edofs(LinkProductType::edofs());

		for (SizeType term = 0; term < geometry_.terms(); ++term) {
			geometry_.fillAdditionalData(additionalData, term, hItems);
			SizeType dofsTotal = LinkProductType::dofs(term, additionalData);
			for (SizeType dofs = 0; dofs < dofsTotal; ++dofs) {
				LinkProductType::connectorDofs(edofs,
				                               term,
				                               dofs,
				                               additionalData);
				SparseElementType tmp = geometry_(smax_,
				                                  emin_,
				                                  hItems,
				                                  edofs,
				                                  term);

				if (tmp == static_cast<RealType>(0.0)) continue;

				tmp = geometry_.vModifier(term, tmp, modelHelper_.time());

				flag = true;

				if (lps!=0) {
					assert(lps->typesaved.size() > total);
					lps->xsaved[total] = x;
					lps->typesaved[total]=type;
					lps->tmpsaved[total]=tmp;
					lps->termsaved[total]=term;
					lps->dofssaved[total]=dofs;
					total++;
				} else {
					calcBond(mBlock,
					         x,
					         type,
					         tmp,
					         term,
					         dofs,
					         additionalData);
					*matrixBlock += mBlock;
				}
			}
		}

		return flag;
	}

	void doTask(SizeType taskNumber ,SizeType threadNum)
	{
		if (xtemp_[threadNum].size() != x_.size())
			xtemp_[threadNum].resize(x_.size(),0.0);

		SparseElementType tmp = 0.0;

		if (taskNumber == 0) {
			modelHelper_.hamiltonianLeftProduct(xtemp_[threadNum],y_);
			return;
		}

		if (taskNumber == 1) {
			modelHelper_.hamiltonianRightProduct(xtemp_[threadNum],y_);
			return;
		}

		assert(taskNumber > 1);
		taskNumber -= 2;
		AdditionalDataType additionalData;
		SizeType xx = 0;
		ProgramGlobals::ConnectionEnum type;
		SizeType term = 0;
		SizeType dofs =0;
		prepare(taskNumber,xx,type,tmp,term,dofs,additionalData);

		linkProduct(xtemp_[threadNum],y_,xx,type,tmp,term,dofs,additionalData);
	}

	void tasks(SizeType total) { total_ = total; }

	SizeType tasks() const { return total_; }

	void sync()
	{
		SizeType total = 0;
		for (SizeType threadNum = 0; threadNum < xtemp_.size(); threadNum++)
			if (xtemp_[threadNum].size() == x_.size()) total++;

		typename PsimagLite::Vector<SparseElementType>::Type x(x_.size(),0);
		for (SizeType threadNum = 0; threadNum < total; threadNum++)
			for (SizeType i=0;i<x_.size();i++)
				x[i]+=xtemp_[threadNum][i];

		if (!ConcurrencyType::isMpiDisabled("HamiltonianConnection"))
			PsimagLite::MPI::allReduce(x);

		for (SizeType i=0;i<x_.size();i++)
			x_[i] += x[i];
	}

	void prepare(SizeType ix,
	             SizeType& x,
	             ProgramGlobals::ConnectionEnum& type,
	             SparseElementType& tmp,
	             SizeType& term,
	             SizeType& dofs,
	             AdditionalDataType& additionalData) const
	{
		x = lps_.xsaved[ix];
		type=lps_.typesaved[ix];
		term = lps_.termsaved[ix];
		dofs = lps_.dofssaved[ix];
		tmp=lps_.tmpsaved[ix];
		geometry_.fillAdditionalData(additionalData,
		                             term,
		                             lps_.hamiltonianAbstract().items(x));
	}

	LinkType getKron(const SparseMatrixType** A,
	                 const SparseMatrixType** B,
	                 SizeType xx,
	                 ProgramGlobals::ConnectionEnum type,
	                 const SparseElementType& valuec,
	                 SizeType term,
	                 SizeType dofs,
	                 const AdditionalDataType& additionalData) const
	{
		int offset = modelHelper_.leftRightSuper().left().block().size();
		PairType ops;
		std::pair<char,char> mods('N','C');
		ProgramGlobals::FermionOrBosonEnum fermionOrBoson=ProgramGlobals::FERMION;
		SizeType angularMomentum=0;
		SizeType category=0;
		RealType angularFactor=0;
		bool isSu2 = modelHelper_.isSu2();
		SparseElementType value = valuec;
		LinkProductType::valueModifier(value,term,dofs,isSu2,additionalData);

		LinkProductType::setLinkData(term,
		                             dofs,
		                             isSu2,
		                             fermionOrBoson,
		                             ops,
		                             mods,
		                             angularMomentum,
		                             angularFactor,
		                             category,
		                             additionalData);
		LinkType link2(i,
		               j,
		               type,
		               value,
		               dofs,
		               fermionOrBoson,
		               ops,
		               mods,
		               angularMomentum,
		               angularFactor,
		               category);
		SizeType sysOrEnv = (link2.type==ProgramGlobals::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SYSTEM : ProgramGlobals::ENVIRON;
		SizeType envOrSys = (link2.type==ProgramGlobals::SYSTEM_ENVIRON) ?
		            ProgramGlobals::ENVIRON : ProgramGlobals::SYSTEM;
		SizeType site1Corrected =(link2.type==ProgramGlobals::SYSTEM_ENVIRON) ?
		            link2.site1 : link2.site1-offset;
		SizeType site2Corrected =(link2.type==ProgramGlobals::SYSTEM_ENVIRON) ?
		            link2.site2-offset : link2.site2;

		*A = &reducedOperator(link2.mods.first,
		                      site1Corrected,
		                      link2.ops.first,
		                      sysOrEnv);
		*B = &reducedOperator(link2.mods.second,
		                      site2Corrected,
		                      link2.ops.second,
		                      envOrSys);

		assert(isNonZeroMatrix(**A));
		assert(isNonZeroMatrix(**B));

		return link2;
	}

	template<typename SomeConcurrencyType,typename SomeOtherConcurrencyType>
	void sync(SomeConcurrencyType& conc,SomeOtherConcurrencyType& conc2)
	{
		conc.reduce(x_,conc2);
	}

private:

	const SparseMatrixType& reducedOperator(char modifier,
	                                        SizeType i,
	                                        SizeType sigma,
	                                        SizeType type) const
	{
		if (!ModelHelperType::BasisType::useSu2Symmetry())
			return reducedOperatorLocal(modifier, i, sigma, type);
		else
			return reducedOperatorSu2(modifier, i, sigma, type);
	}

	const SparseMatrixType& reducedOperatorLocal(char modifier,
	                                             SizeType i,
	                                             SizeType sigma,
	                                             SizeType type) const
	{
		const SparseMatrixType* m = 0;
		if (type == ProgramGlobals::SYSTEM) {
			PairType ii = modelHelper_.leftRightSuper().left().getOperatorIndices(i,sigma);
			m = &(modelHelper_.leftRightSuper().left().getOperatorByIndex(ii.first).data);
		} else {
			assert(type == ProgramGlobals::ENVIRON);
			PairType ii =modelHelper_.leftRightSuper().right().getOperatorIndices(i,sigma);
			m =&(modelHelper_.leftRightSuper().right().getOperatorByIndex(ii.first).data);
		}

		if (modifier == 'N') return *m;

		SizeType typeIndex = (type == ProgramGlobals::SYSTEM) ? 0 : 1;
		SizeType dummySite = 0; // FIXME for Immm
		SizeType opsPerSiteMax = modelHelper_.leftRightSuper().left().operatorsPerSite(dummySite);
		SizeType packed = sigma + typeIndex*opsPerSiteMax + i*2*opsPerSiteMax;
		int seen = indexOfItem(seen_, packed);
		if (seen >= 0)
			return *(garbage_[seen]);

		SparseMatrixType* mc = new SparseMatrixType;
		transposeConjugate(*mc, *m);
		garbage_.push_back(mc);
		seen_.push_back(packed);
		return *mc;
	}

	const SparseMatrixType& reducedOperatorSu2(char modifier,
	                                           SizeType i,
	                                           SizeType sigma,
	                                           SizeType type) const
	{
		if (type == ProgramGlobals::SYSTEM) {
			PairType ii = modelHelper_.leftRightSuper().left().getOperatorIndices(i,
			                                                                      sigma);
			return modelHelper_.leftRightSuper().left().getReducedOperatorByIndex(modifier,
			                                                                      ii).data;
		}

		assert(type == ProgramGlobals::ENVIRON);
		PairType ii =modelHelper_.leftRightSuper().right().getOperatorIndices(i,
		                                                                      sigma);
		return modelHelper_.leftRightSuper().right().getReducedOperatorByIndex(modifier,
		                                                                       ii).data;
	}

	//! Adds a connector between system and environment
	SizeType calcBond(SparseMatrixType &matrixBlock,
	                  SizeType xx,
	                  ProgramGlobals::ConnectionEnum type,
	                  const SparseElementType& valuec,
	                  SizeType term,
	                  SizeType dofs,
	                  const AdditionalDataType& additionalData) const
	{
		SparseMatrixType const* A = 0;
		SparseMatrixType const* B = 0;
		LinkType link2 = getKron(&A,
		                         &B,
		                         xx,
		                         type,
		                         valuec,
		                         term,
		                         dofs,
		                         additionalData);
		modelHelper_.fastOpProdInter(*A,*B,matrixBlock,link2);

		return matrixBlock.nonZeros();
	}

	//! Computes x+=H_{ij}y where H_{ij} is a Hamiltonian that connects system and environment
	void linkProduct(typename PsimagLite::Vector<SparseElementType>::Type& x,
	                 const typename PsimagLite::Vector<SparseElementType>::Type& y,
	                 SizeType xx,
	                 ProgramGlobals::ConnectionEnum type,
	                 const SparseElementType &valuec,
	                 SizeType term,
	                 SizeType dofs,
	                 const AdditionalDataType& additionalData) const
	{
		SparseMatrixType const* A = 0;
		SparseMatrixType const* B = 0;
		LinkType link2 = getKron(&A,&B,xx,type,valuec,term,dofs,additionalData);
		modelHelper_.fastOpProdInter(x,y,*A,*B,link2);
	}

	bool isNonZeroMatrix(const SparseMatrixType& m) const
	{
		if (m.rows() > 0 && m.cols() > 0) return true;
		return false;
	}

	static SizeType indexOfItem(const VectorSizeType& v, SizeType x)
	{
		SizeType n = v.size();
		for (SizeType i = 0; i < n; ++i)
			if (v[i] == x) return i;

		throw PsimagLite::RuntimeError("indexOfItem(): item not found " + ttos(x) + "\n");
	}

	const GeometryType& geometry_;
	const ModelHelperType& modelHelper_;
	const LinkProductStructType& lps_;
	typename PsimagLite::Vector<SparseElementType>::Type& x_;
	const typename PsimagLite::Vector<SparseElementType>::Type& y_;
	const VectorSizeType& systemBlock_;
	const VectorSizeType& envBlock_;
	SizeType smax_,emin_;
	VectorVectorType xtemp_;
	SizeType total_;
	HamiltonianAbstractType hamAbstract_;
	mutable typename PsimagLite::Vector<SparseMatrixType*>::Type garbage_;
	mutable VectorSizeType seen_;
}; // class HamiltonianConnection
} // namespace Dmrg

/*@}*/
#endif // HAMILTONIAN_CONNECTION_H

