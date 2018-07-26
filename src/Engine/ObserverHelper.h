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

namespace Dmrg {
template<
        typename IoInputType_,
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
	typedef typename LeftRightSuperType::BasisWithOperatorsType
	BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef DmrgSerializer<LeftRightSuperType,VectorWithOffsetType> DmrgSerializerType;
	typedef typename DmrgSerializerType::FermionSignType FermionSignType;

	enum {LEFT_BRAKET=0,RIGHT_BRAKET=1};
	enum SaveEnum {SAVE_YES, SAVE_NO};

	ObserverHelper(IoInputType& io,
	               SizeType nf,
	               SizeType trail,
	               SizeType numberOfPthreads,
	               bool hasTimeEvolution,
	               bool verbose,
	               bool withLegacyBugs)
	    : io_(io),
	      currentPos_(numberOfPthreads),
	      verbose_(verbose),
	      withLegacyBugs_(withLegacyBugs),
	      bracket_(2,0),
	      noMoreData_(false)
	{
		if (nf > 0)
			if (!init(hasTimeEvolution,nf,SAVE_YES))
				return;

		if (trail > 0)
			if (!init(hasTimeEvolution,trail,SAVE_NO))
				return;
	}

	~ObserverHelper()
	{
		for (SizeType i=0;i<dSerializerV_.size();i++) {
			DmrgSerializerType* p = dSerializerV_[i];
			delete p;
		}
	}

	bool endOfData() const { return noMoreData_; }

	void setPointer(SizeType threadId,SizeType pos)
	{
		assert(threadId<currentPos_.size());
		currentPos_[threadId]=pos;
	}

	SizeType getPointer(SizeType threadId) const
	{
		assert(threadId<currentPos_.size());
		return currentPos_[threadId];
	}

	void setBrakets(SizeType left,SizeType right)
	{
		bracket_[LEFT_BRAKET]=left;
		bracket_[RIGHT_BRAKET]=right;
	}

	void transform(SparseMatrixType& ret,const SparseMatrixType& O2,size_t threadId) const
	{
		assert(checkPos(threadId));
		return dSerializerV_[currentPos_[threadId]]->transform(ret,O2);
	}

	SizeType cols(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return dSerializerV_[currentPos_[threadId]]->cols();
	}

	SizeType rows(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return dSerializerV_[currentPos_[threadId]]->rows();
	}

	const FermionSignType& fermionicSignLeft(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return dSerializerV_[currentPos_[threadId]]->fermionicSignLeft();
	}

	const FermionSignType& fermionicSignRight(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return dSerializerV_[currentPos_[threadId]]->fermionicSignRight();
	}

	const LeftRightSuperType& leftRightSuper(SizeType threadId) const
	{
		if (!checkPos(threadId))
			err("No more data\n");
		return dSerializerV_[currentPos_[threadId]]->leftRightSuper();
	}

	ProgramGlobals::DirectionEnum direction(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return dSerializerV_[currentPos_[threadId]]->direction();
	}

	const VectorWithOffsetType& wavefunction(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return dSerializerV_[currentPos_[threadId]]->wavefunction();
	}

	RealType time(SizeType threadId) const
	{
		if (timeSerializerV_.size() == 0) return 0.0;
		assert(checkPos(threadId));
		return timeSerializerV_[currentPos_[threadId]].time();
	}

	SizeType site(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return  (timeSerializerV_.size()==0) ?
		            dSerializerV_[currentPos_[threadId]]->site()
		        : timeSerializerV_[currentPos_[threadId]].site();
		}

		SizeType size() const
		{
		return dSerializerV_.size(); //-1;
	}

	SizeType marker(SizeType threadId) const
	{
		assert(checkPos(threadId));
		return timeSerializerV_[currentPos_[threadId]].marker();
	}

	const VectorWithOffsetType& getVectorFromBracketId(SizeType leftOrRight,
	                                                   SizeType threadId) const
	{
		SizeType braketId = bracket_[leftOrRight];
		// braketId == 0 means GS
		if (braketId == 0) {
			return wavefunction(threadId);
		}

		// braketId > 0 then it means the "time vector" number braketId - 1
		assert(braketId > 0);
		return timeVector(braketId - 1, threadId);
	}

	const VectorWithOffsetType& timeVector(SizeType braketId,
	                                       SizeType threadId) const
	{
		assert(checkPos(threadId));
		return timeSerializerV_[currentPos_[threadId]].vector(braketId);
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

	bool init(bool hasTimeEvolution,
	          SizeType nf,
	          SaveEnum saveOrNot)
	{
		PsimagLite::String prefix = "Serializer";
		SizeType total = 0;
		io_.read(total, prefix + "/Size");
		for (SizeType i = 0; i < total; ++i) {
			if (nf > 0 && i == nf) break;

			DmrgSerializerType* dSerializer = new DmrgSerializerType(io_,
			                                                         prefix + "/" + ttos(i),
			                                                         false,
			                                                         true);
			if (saveOrNot == SAVE_YES)
				dSerializerV_.push_back(dSerializer);
			else
				delete dSerializer;
			if (hasTimeEvolution) {
				TimeSerializerType ts(io_, ""); // FIXME
				if (saveOrNot == SAVE_YES)
					timeSerializerV_.push_back(ts);
			}

			std::cerr<<__FILE__<<" read "<<i<<" out of "<<total<<"\n";
		}

		noMoreData_ = true;
		return (dSerializerV_.size() > 0);
	}

	void integrityChecks()
	{
		if (dSerializerV_.size()!=timeSerializerV_.size())
			throw PsimagLite::RuntimeError("Error 1\n");
		if (dSerializerV_.size()==0) return;
		for (SizeType x=0;x<dSerializerV_.size()-1;x++) {
			SizeType n = dSerializerV_[x]->leftRightSuper().super().size();
			if (n==0) continue;
			if (n!=timeSerializerV_[x].size())
				throw PsimagLite::RuntimeError("Error 2\n");
		}

	}

	void getTransform(MatrixType& transform,int ns)
	{
		io_.read(transform, "TRANSFORM_sites",ns);
		io_.rewind();
	}

	void getWaveFunction(VectorType& wavefunction,SizeType ns)
	{
		VectorWithOffsetType tmpV;
		tmpV.read(io_,"WAVEFUNCTION_sites=",ns);
		tmpV.toSparse(wavefunction);
		io_.rewind();
	}

	void getDirection(int& x,int ns)
	{
		io_.readline(x,"DIRECTION=",ns);
		io_.rewind();
	}

	bool checkPos(SizeType threadId) const
	{
		if (threadId>=currentPos_.size())
			return checkFailedThread(threadId);

		SizeType pos = currentPos_[threadId];

		if (pos>=dSerializerV_.size())
			return checkFailed1(threadId,pos);

		bool hasTimeE = (timeSerializerV_.size()>0);

		if (!hasTimeE) return true;
		if (pos>=timeSerializerV_.size())
			return checkFailed2(threadId,pos);
		return true;
	}

	bool checkFailedThread(SizeType threadId) const
	{
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += " thread=" + ttos(threadId);
		str += " >= currentPos.size=" + ttos(currentPos_.size());
		str += "\n";
		std::cerr<<str;
		return false;
	}

	bool checkFailed1(SizeType threadId,SizeType pos) const
	{
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += " thread=" + ttos(threadId) + " currentPos=" + ttos(pos);
		str += " >= serializer.size=" + ttos(dSerializerV_.size());
		str += "\n";
		std::cerr<<str;
		return false;
	}

	bool checkFailed2(SizeType threadId,SizeType pos) const
	{
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += " thread=" + ttos(threadId) + " currentPos=" + ttos(pos);
		str += " >= time serializer.size=" + ttos(timeSerializerV_.size());
		str += "\n";
		std::cerr<<str;
		return false;
	}

	IoInputType& io_;
	typename PsimagLite::Vector<DmrgSerializerType*>::Type dSerializerV_;
	typename PsimagLite::Vector<TimeSerializerType>::Type timeSerializerV_;
	typename PsimagLite::Vector<SizeType>::Type currentPos_; // it's a vector: one per pthread
	bool verbose_;
	bool withLegacyBugs_;
	typename PsimagLite::Vector<SizeType>::Type bracket_;
	bool noMoreData_;
};  //ObserverHelper
} // namespace Dmrg

/*@}*/
#endif
