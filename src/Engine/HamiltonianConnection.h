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

#include "CrsMatrix.h"
#include "Concurrency.h"
#include <cassert>
#include "ProgramGlobals.h"
#include "HamiltonianAbstract.h"
#include "SuperGeometry.h"
#include "Vector.h"
#include "VerySparseMatrix.h"
#include "ProgressIndicator.h"
#include "OperatorStorage.h"

namespace Dmrg {

// Keep this class independent of x and y in x = H*y
// For things that depend on x and y use ParallelHamiltonianConnection.h
template<typename ModelLinksType, typename ModelHelperType_>
class HamiltonianConnection {

public:

	typedef ModelHelperType_ ModelHelperType;
	typedef typename ModelLinksType::GeometryType GeometryType;
	typedef SuperGeometry<GeometryType> SuperGeometryType;
	typedef HamiltonianAbstract<SuperGeometryType> HamiltonianAbstractType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::OperatorStorageType OperatorStorageType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef VerySparseMatrix<ComplexOrRealType> VerySparseMatrixType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Concurrency ConcurrencyType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::ParamsForKroneckerDumperType ParamsForKroneckerDumperType;
	typedef typename LeftRightSuperType::KroneckerDumperType KroneckerDumperType;
	typedef typename PsimagLite::Vector<LinkType>::Type VectorLinkType;
	typedef typename ModelLinksType::HermitianEnum HermitianEnum;

	HamiltonianConnection(SizeType m,
	                      const LeftRightSuperType& lrs,
	                      const GeometryType& geometry,
	                      const ModelLinksType& lpb,
	                      RealType targetTime,
	                      const ParamsForKroneckerDumperType* pKroneckerDumper)
	    : modelHelper_(m, lrs),
	      superGeometry_(geometry),
	      modelLinks_(lpb),
	      targetTime_(targetTime),
	      kroneckerDumper_(pKroneckerDumper,lrs,m),
	      progress_("HamiltonianConnection"),
	      systemBlock_(modelHelper_.leftRightSuper().left().block()),
	      envBlock_(modelHelper_.leftRightSuper().right().block()),
	      smax_(*std::max_element(systemBlock_.begin(),systemBlock_.end())),
	      emin_(*std::min_element(envBlock_.begin(),envBlock_.end())),
	      hamAbstract_(superGeometry_,
	                   smax_,
	                   emin_,
	                   modelHelper_.leftRightSuper().super().block()),
	      totalOnes_(hamAbstract_.items())
	{
		lps_.reserve(ProgramGlobals::MAX_LPS);
		SizeType nitems = hamAbstract_.items();
		for (SizeType x = 0; x < nitems; ++x)
			totalOnes_[x] = cacheConnections(x);

		SizeType last = lrs.super().block().size();
		assert(last > 0);
		--last;
		SizeType numberOfSites = geometry.numberOfSites();
		assert(numberOfSites > 0);
		bool superIsReallySuper = (lrs.super().block()[0] == 0 &&
		        lrs.super().block()[last] == numberOfSites - 1);

		if (!superIsReallySuper)
			return; // <-- CONDITIONAL EARLY EXIT HERE

		PsimagLite::OstringStream msg;
		msg<<"LinkProductStructSize="<<lps_.size();
		progress_.printline(msg,std::cout);

		PsimagLite::OstringStream msg2;
		// add left and right contributions
		msg2<<"PthreadsTheoreticalLimitForThisPart="<<(lps_.size() + 2);

		// The theoretical maximum number of pthreads that are useful
		// is equal to C + 2, where
		// C = number of connection = 2*G*M
		// where G = geometry factor
		// and   M = model factor
		// G = 1 for a chain
		// G = leg for a ladder with leg legs. (For instance, G=2 for a 2-leg ladder).
		// M = 2 for the Hubbard model
		//
		// In general, M = \sum_{term=0}^{terms} dof(term)
		// where terms and dof(term) is model dependent.
		// To find M for a model, go to the Model directory and see LinkProduct*.h file.
		// the return of function terms() for terms.
		// For dof(term) see function dof(SizeType term,...).
		// For example, for the HubbardModelOneBand, one must look at
		// src/Model/HubbardOneBand/LinkProductHubbardOneBand.h
		// In that file, terms() returns 1, so that terms=1
		// Therefore there is only one term: term = 0.
		// And dof(0,...) = 2, as you can see in  LinkProductHubbardOneBand.h.
		// Then M = 2.
		progress_.printline(msg2,std::cout);
	}

	void matrixBond(VerySparseMatrixType& matrix) const
	{
		SizeType matrixRank = matrix.rows();
		VerySparseMatrixType matrix2(matrixRank, matrixRank);
		SizeType nitems = totalOnes_.size();

		SizeType x = 0;
		for (SizeType xx = 0; xx < nitems; ++xx) {
			SparseMatrixType matrixBlock(matrixRank, matrixRank);
			for (SizeType i = 0; i < totalOnes_[xx]; ++i) {
				SparseMatrixType mBlock;
				OperatorStorageType const* A = 0;
				OperatorStorageType const* B = 0;
				const LinkType& link2 = getKron(&A, &B, x++);
				modelHelper_.fastOpProdInter(A->getCRS(), B->getCRS(), mBlock, link2);

				matrixBlock += mBlock;
			}

			VerySparseMatrixType vsm(matrixBlock);
			matrix2+=vsm;
		}

		matrix += matrix2;
	}

	const LinkType& getKron(const OperatorStorageType** A,
	                        const OperatorStorageType** B,
	                        SizeType xx) const
	{
		assert(xx < lps_.size());
		const LinkType& link2 = lps_[xx];

		assert(link2.type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON ||
		       link2.type == ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM);

		const ProgramGlobals::SysOrEnvEnum sysOrEnv =
		        (link2.type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::SYSTEM : ProgramGlobals::SysOrEnvEnum::ENVIRON;
		const ProgramGlobals::SysOrEnvEnum envOrSys =
		        (link2.type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::ENVIRON : ProgramGlobals::SysOrEnvEnum::SYSTEM;

		SizeType i = PsimagLite::indexOrMinusOne(modelHelper_.leftRightSuper().super().block(),
		                                         link2.site1);
		SizeType j = PsimagLite::indexOrMinusOne(modelHelper_.leftRightSuper().super().block(),
		                                         link2.site2);

		int offset = modelHelper_.leftRightSuper().left().block().size();

		SizeType site1Corrected = (link2.type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            i : i - offset;
		SizeType site2Corrected = (link2.type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            j - offset : j;

		*A = &modelHelper_.reducedOperator(link2.mods.first,
		                                   site1Corrected,
		                                   link2.ops.first,
		                                   sysOrEnv);
		*B = &modelHelper_.reducedOperator(link2.mods.second,
		                                   site2Corrected,
		                                   link2.ops.second,
		                                   envOrSys);

		assert(isNonZeroMatrix(**A));
		assert(isNonZeroMatrix(**B));

		(*A)->checkValidity();
		(*B)->checkValidity();

		return link2;
	}

	KroneckerDumperType& kroneckerDumper() const
	{
		return kroneckerDumper_;
	}

	const ModelHelperType& modelHelper() const { return modelHelper_; }

	SizeType tasks() const {return lps_.size(); }

private:

	SizeType cacheConnections(SizeType x)
	{
		const VectorSizeType& hItems = hamAbstract_.item(x);
		assert(superGeometry_.connected(smax_, emin_, hItems));

		ProgramGlobals::ConnectionEnum type = superGeometry_.connectionKind(smax_, hItems);

		assert(type != ProgramGlobals::ConnectionEnum::SYSTEM_SYSTEM &&
		        type != ProgramGlobals::ConnectionEnum::ENVIRON_ENVIRON);

		assert(hItems.size() == 2);

		VectorSizeType edofs(2, 0);
		SizeType totalOne = 0;
		SizeType geometryTerms = superGeometry_.geometry().terms();
		for (SizeType termIndex = 0; termIndex < geometryTerms; ++termIndex) {

			if (!modelLinks_.areSitesCompatibleForThisTerm(termIndex, hItems[0], hItems[1]))
				continue;

			const typename ModelLinksType::TermType& term = modelLinks_.term(termIndex);

			SizeType dofsTotal = term.size();
			for (SizeType dofs = 0; dofs < dofsTotal; ++dofs) {

				const typename ModelLinksType::TermType::OneLinkType& oneLink = term(dofs);

				edofs[0] = oneLink.orbs.first;
				edofs[1] = oneLink.orbs.second;
				ComplexOrRealType tmp = superGeometry_(smax_,
				                                       emin_,
				                                       hItems,
				                                       edofs,
				                                       termIndex);

				if (tmp == static_cast<RealType>(0.0)) continue;

				tmp = superGeometry_.geometry().vModifier(termIndex, tmp, targetTime_);

				oneLink.modifier(tmp);

				LinkType link2(hItems[0],
				        hItems[1],
				        type,
				        tmp,
				        oneLink.fermionOrBoson,
				        oneLink.indices,
				        oneLink.mods,
				        oneLink.angularMomentum,
				        oneLink.angularFactor,
				        oneLink.category);

				++totalOne;
				lps_.push_back(link2);

				// add h.c. parts if needed
				if (connectionIsHermitian(link2)) continue;

				link2.value = PsimagLite::conj(tmp);

				if (oneLink.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION)
					link2.value *= (-1.0);

				char saved = oneLink.mods.first;
				link2.mods.first = oneLink.mods.second;
				link2.mods.second = saved;

				++totalOne;
				lps_.push_back(link2);

			}
		}

		return totalOne;
	}

	bool connectionIsHermitian(const LinkType& link) const
	{
		return (link.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION) ?
		            linkIsHermitianFermion(link) : linkIsHermitianBoson(link);
	}

	bool linkIsHermitianFermion(const LinkType& link) const
	{
		assert(link.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION);

		HermitianEnum h1 = modelLinks_.getHermitianProperty(link.ops.first, link.site1);
		HermitianEnum h2 = modelLinks_.getHermitianProperty(link.ops.second, link.site2);

		bool isHermit1 = (h1 == ModelLinksType::HERMIT_PLUS);
		bool isHermit2 = (h2 == ModelLinksType::HERMIT_PLUS);
		bool isAnti1 = (h1 == ModelLinksType::HERMIT_MINUS);
		bool isAnti2 = (h2 == ModelLinksType::HERMIT_MINUS);
		bool b1 = (isHermit1 && isAnti2);
		bool b2 = (isAnti1 && isHermit2);
		return (b1 || b2);
	}

	bool linkIsHermitianBoson(const LinkType& link) const
	{
		assert(link.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::BOSON);

		HermitianEnum h1 = modelLinks_.getHermitianProperty(link.ops.first, link.site1);
		HermitianEnum h2 = modelLinks_.getHermitianProperty(link.ops.second, link.site2);

		bool isHermit1 = (h1 == ModelLinksType::HERMIT_PLUS);
		bool isHermit2 = (h2 == ModelLinksType::HERMIT_PLUS);

		return (isHermit1 && isHermit2);
	}

	const ModelHelperType modelHelper_;
	SuperGeometryType superGeometry_;
	const ModelLinksType& modelLinks_;
	RealType targetTime_;
	mutable KroneckerDumperType kroneckerDumper_;
	PsimagLite::ProgressIndicator progress_;
	VectorLinkType lps_;
	const VectorSizeType& systemBlock_;
	const VectorSizeType& envBlock_;
	SizeType smax_;
	SizeType emin_;
	HamiltonianAbstractType hamAbstract_;
	VectorSizeType totalOnes_;
}; // class HamiltonianConnection
} // namespace Dmrg

/*@}*/
#endif // HAMILTONIAN_CONNECTION_H

