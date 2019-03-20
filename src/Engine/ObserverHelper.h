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

// temporary class only <<--------- FIXME DELETE
class PointerForSerializer {

public:

	PointerForSerializer(SizeType n)
	    : pos_(n)
	{}

	void setPointer(SizeType pos)
	{
		pos_ = pos;
	}

	SizeType get() const { return  pos_; }

private:

	SizeType pos_;
};

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
	typedef PointerForSerializer PointerForSerializerType;
	typedef GetBraOrKet GetBraOrKetType;

	enum class SaveEnum {YES, NO};

	ObserverHelper(IoInputType& io,
	               SizeType start,
	               SizeType nf,
	               SizeType trail,
	               bool withLegacyBugs)
	    : io_(io),
	      withLegacyBugs_(withLegacyBugs),
	      noMoreData_(false),
	      dSsize_(0),
	      timeSsize_(0),
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
		for (SizeType i = 0; i < dSsize_; ++i)
			delete dSerializerV_[i];

		for (SizeType i = 0; i < timeSerializerV_.size(); ++i)
			delete timeSerializerV_[i];

		dSerializerV_.clear();
		timeSerializerV_.clear();
		dSsize_ = timeSsize_ = 0;
	}

	const SizeType& numberOfSites() const { return numberOfSites_; }

	bool endOfData() const { return noMoreData_; }

	void transform(SparseMatrixType& ret,
	               const SparseMatrixType& O2,
	               const PointerForSerializerType& ind) const
	{
		assert(ind.get() < dSerializerV_.size());
		return dSerializerV_[ind.get()]->transform(ret,O2);
	}

	SizeType cols(const PointerForSerializerType& ind) const
	{
		assert(ind.get() < dSerializerV_.size());
		return dSerializerV_[ind.get()]->cols();
	}

	SizeType rows(const PointerForSerializerType& ind) const
	{
		assert(ind.get() < dSerializerV_.size());
		return dSerializerV_[ind.get()]->rows();
	}

	short int signsOneSite(SizeType site) const
	{
		assert(site < signsOneSite_.size());
		return signsOneSite_[site];
	}

	const FermionSignType& fermionicSignLeft(const PointerForSerializerType& ind) const
	{
		assert(ind.get() < dSerializerV_.size());
		return dSerializerV_[ind.get()]->fermionicSignLeft();
	}

	const FermionSignType& fermionicSignRight(const PointerForSerializerType& ind) const
	{
		assert(ind.get() < dSerializerV_.size());
		return dSerializerV_[ind.get()]->fermionicSignRight();
	}

	const LeftRightSuperType& leftRightSuper(const PointerForSerializerType& ind) const
	{
		return dSerializerV_[ind.get()]->leftRightSuper();
	}

	ProgramGlobals::DirectionEnum direction(const PointerForSerializerType& ind) const
	{
		assert(ind.get() < dSerializerV_.size());
		return dSerializerV_[ind.get()]->direction();
	}

	const VectorWithOffsetType& wavefunction(const PointerForSerializerType& ind) const
	{
		assert(ind.get() < dSerializerV_.size());
		return dSerializerV_[ind.get()]->wavefunction();
	}

	RealType time(const PointerForSerializerType& index) const
	{
		if (timeSsize_ == 0) return 0.0;
		assert(index.get() < dSerializerV_.size());
		SizeType ind = index.get();
		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->time();
	}

	SizeType site(const PointerForSerializerType& index) const
	{
		assert(index.get() < dSerializerV_.size());
		SizeType ind = index.get();

		if (timeSsize_ == 0) {
			assert(ind < dSerializerV_.size());
			assert(dSerializerV_[ind]);
			return dSerializerV_[ind]->site();
		}

		assert(ind < timeSerializerV_.size());
		assert(timeSerializerV_[ind]);
		return timeSerializerV_[ind]->site();
	}

	SizeType size() const
	{
		return dSsize_; //-1;
	}

	const VectorWithOffsetType& getVectorFromBracketId(PsimagLite::String braOrKet,
	                                                   const PointerForSerializerType& index) const
	{
		SizeType braketId = braketStringToNumber(braOrKet);
		// braketId == 0 means GS
		if (braketId == 0)
			return wavefunction(index);

		// braketId > 0 then it means the "time vector" number braketId - 1
		assert(braketId > 0);
		return timeVector(braketId - 1, index);
	}

	const VectorWithOffsetType& timeVector(SizeType braketId,
	                                       const PointerForSerializerType& index) const
	{
		assert(index.get() < dSerializerV_.size());
		SizeType ind = index.get();
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

		dSsize_ = dSerializerV_.size();
		timeSsize_ = timeSerializerV_.size();
		noMoreData_ = (end == total);
		return (dSsize_ > 0);
	}

	static SizeType braketStringToNumber(const PsimagLite::String& str)
	{
		if (str == "gs") return 0;
		if (str == "time") return 1; // == "P0", "time" is legacy notation
		int x = GetBraOrKetType::getPtype(str);
		if (x >= 0) return x;

		PsimagLite::String msg("ObserverHelper::braketStringToNumber:");
		throw PsimagLite::RuntimeError(msg + " must be gs or time or P\\d+\n");
	}

	IoInputType& io_;
	typename PsimagLite::Vector<DmrgSerializerType*>::Type dSerializerV_;
	typename PsimagLite::Vector<TimeSerializerType*>::Type timeSerializerV_;
	bool withLegacyBugs_;
	bool noMoreData_;
	SizeType dSsize_;
	SizeType timeSsize_;
	VectorShortIntType signsOneSite_;
	SizeType numberOfSites_;
};  // ObserverHelper
} // namespace Dmrg

/*@}*/
#endif
