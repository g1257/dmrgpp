/*
Copyright (c) 2008-2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 1.0.0]
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

/*! \file ObserverHelper.h
 *
 *  A class to read and serve precomputed data to the observer
 *
 */
#ifndef PRECOMPUTED_H
#define PRECOMPUTED_H
#include "SparseVector.h"
#include "ProgramGlobals.h"
#include "TimeSerializer.h"
#include "DmrgSerializer.h"
#include "VectorWithOffsets.h" // to include norm
#include "VectorWithOffset.h" // to include norm
#include "GetBraOrKet.h"

namespace Dmrg {

template<typename IoInputType_,
         typename MatrixType_,
         typename VectorType_,
         typename VectorWithOffsetType_,
         typename LeftRightSuperType>
class ObserverHelper {

public:

	typedef IoInputType_ IoInputType;
	typedef MatrixType_ MatrixType;
	typedef VectorType_ VectorType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef SizeType IndexType;
	typedef typename VectorType::value_type FieldType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef DmrgSerializer<LeftRightSuperType,VectorWithOffsetType> DmrgSerializerType;
	typedef typename DmrgSerializerType::FermionSignType FermionSignType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<short int>::Type VectorShortIntType;
	typedef PsimagLite::GetBraOrKet GetBraOrKetType;

	enum class SaveEnum {YES, NO};

	ObserverHelper(IoInputType& io,
	               SizeType start,
	               SizeType nf,
	               SizeType trail,
	               bool withLegacyBugs)
	    : io_(io),
	      withLegacyBugs_(withLegacyBugs),
	      noMoreData_(false),
	      numberOfSites_(0)
	{
		typename BasisWithOperatorsType::VectorBoolType odds;
		io_.read(odds, "OddElectronsOneSite");
		SizeType n = odds.size();
		signsOneSite_.resize(n);
		for (SizeType i = 0; i < n; ++i)
			signsOneSite_[i] = (odds[i]) ? -1 : 1;

		if (nf > 0)
			if (!init(start, start + nf, SaveEnum::YES))
				return;

		if (trail > 0)
			if (!init(start, start + trail, SaveEnum::NO))
				return;
	}

	~ObserverHelper()
	{
		for (SizeType i = 0; i < dSerializerV_.size(); ++i) {
			delete dSerializerV_[i];
			dSerializerV_[i] = 0;
		}

		for (SizeType i = 0; i < timeSerializerV_.size(); ++i) {
			delete timeSerializerV_[i];
			timeSerializerV_[i] = 0;
		}
	}

	const SizeType& numberOfSites() const { return numberOfSites_; }

	bool endOfData() const { return noMoreData_; }

	void transform(SparseMatrixType& ret,
	               const SparseMatrixType& O2,
	               SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->transform(ret, O2);
	}

	SizeType cols(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->cols();
	}

	SizeType rows(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->rows();
	}

	short int signsOneSite(SizeType site) const
	{
		assert(site < signsOneSite_.size());
		return signsOneSite_[site];
	}

	const FermionSignType& fermionicSignLeft(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->fermionicSignLeft();
	}

	const FermionSignType& fermionicSignRight(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->fermionicSignRight();
	}

	const LeftRightSuperType& leftRightSuper(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->leftRightSuper();
	}

	ProgramGlobals::DirectionEnum direction(SizeType ind) const
	{
		checkIndex(ind);
		return dSerializerV_[ind]->direction();
	}

	const VectorWithOffsetType& wavefunction(SizeType ind,
	                                         SizeType levelIndex,
	                                         SizeType sectorIndex) const
	{
		checkIndex(ind);

		return dSerializerV_[ind]->wavefunction(levelIndex, sectorIndex);
	}

	RealType time(SizeType ind) const
	{
		if (timeSerializerV_.size() == 0) return 0.0;
		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->time();
	}

	SizeType site(SizeType ind) const
	{
		if (timeSerializerV_.size() == 0) {
			checkIndex(ind);
			return dSerializerV_[ind]->site();
		}

		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->site();
	}

	SizeType size() const { return dSerializerV_.size(); }

	const VectorWithOffsetType& getVectorFromBracketId(const PsimagLite::GetBraOrKet& braOrKet,
	                                                   SizeType index) const
	{
		if (braOrKet.isPvector()) {
			const SizeType pIndex = braOrKet.pIndex();
			return timeVector(index, pIndex);
		}

		const SizeType levelIndex = braOrKet.levelIndex();
		return wavefunction(index, levelIndex, braOrKet.sectorIndex());
	}

	const VectorWithOffsetType& timeVector(SizeType braketId,
	                                       SizeType ind) const
	{
		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->vector(braketId);
	}

	bool withLegacyBugs() const
	{
		return withLegacyBugs_;
	}

	friend std::ostream& operator<<(std::ostream& os, ObserverHelper& p)
	{
		for (SizeType i=0;i<p.SpermutationInverse_.size();i++) {
			os<<"i="<<i<<"\n";
			os<<"\tS.size="<<p.SpermutationInverse_[i].size();
			os<<" "<<p.Spermutation_[i].size()<<"\n";
			os<<"\tSE.size="<<p.SEpermutationInverse_[i].size();
			os<<" "<<p.SEpermutation_[i].size()<<"\n";
			os<<"\tElectrons.size="<<p.electrons_[i].size()<<"\n";
			os<<"\tTransform="<<p.transform_[i].n_row()<<"x";
			os<<p.transform_[i].n_col()<<"\n";
			os<<"\tWF.size="<<p.wavefunction_[i].size()<<"\n";
		}

		return os;
	}

private:

	bool init(SizeType start, SizeType end, SaveEnum saveOrNot)
	{
		PsimagLite::String prefix = "Serializer";
		SizeType total = 0;
		io_.read(total, prefix + "/Size");
		if (start >= end || start >= total || end > total) return false;

		for (SizeType i = start; i < end; ++i) {

			DmrgSerializerType* dSerializer = new DmrgSerializerType(io_,
			                                                         prefix + "/" + ttos(i),
			                                                         false,
			                                                         true);


			SizeType tmp = dSerializer->leftRightSuper().sites();
			if (tmp > 0 && numberOfSites_ == 0) numberOfSites_ = tmp;

			if (saveOrNot == SaveEnum::YES)
				dSerializerV_.push_back(dSerializer);
			else
				delete dSerializer;

			try {
				PsimagLite::String prefix("/TargetingCommon/" + ttos(i));
				TimeSerializerType* ts = new TimeSerializerType(io_, prefix);
				if (saveOrNot == SaveEnum::YES)
					timeSerializerV_.push_back(ts);
				else
					delete ts;
			} catch(...) {}

			std::cerr<<__FILE__<<" read "<<i<<" out of "<<total<<"\n";
		}

		noMoreData_ = (end == total);
		return (dSerializerV_.size() > 0);
	}

	static SizeType braketStringToNumber(const PsimagLite::String& str)
	{
		GetBraOrKetType ketOrBra(str);
		if (!ketOrBra.isPvector())
			return 0;

		return ketOrBra.levelIndex();
	}

	void checkIndex(SizeType ind) const
	{
		if (ind >= dSerializerV_.size())
			err("Index " + ttos(ind) + " bigger than " + ttos(dSerializerV_.size()));

		if (dSerializerV_[ind]) return;

		err("dSerializerV_ at index " + ttos(ind) + " point to 0x0\n");
	}

	IoInputType& io_;
	typename PsimagLite::Vector<DmrgSerializerType*>::Type dSerializerV_;
	typename PsimagLite::Vector<TimeSerializerType*>::Type timeSerializerV_;
	const bool withLegacyBugs_;
	bool noMoreData_;
	VectorShortIntType signsOneSite_;
	SizeType numberOfSites_;
};  // ObserverHelper
} // namespace Dmrg

/*@}*/
#endif
