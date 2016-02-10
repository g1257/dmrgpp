/*
Copyright (c) 2008-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

/*! \file Observer.h
 *
 *  A class to perform post-processing calculation of observables
 *
 */
#ifndef DMRG_OBSERVE_H
#define DMRG_OBSERVE_H
#include "ObserverHelper.h"
#include "OnePointCorrelations.h"
#include "CorrelationsSkeleton.h"
#include "TwoPointCorrelations.h"
#include "FourPointCorrelations.h"
#include "VectorWithOffsets.h" // for operator*
#include "VectorWithOffset.h" // for operator*
#include "Profiling.h"
#include "Parallel4PointDs.h"
#include "MultiPointCorrelations.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "Utils.h"

namespace Dmrg {

template<typename VectorWithOffsetType_,typename ModelType_,typename IoInputType>
class Observer {

	typedef typename VectorWithOffsetType_::value_type FieldType;
	typedef PsimagLite::SparseVector<FieldType> VectorType;
	typedef typename ModelType_::RealType RealType;
	typedef PsimagLite::Matrix<FieldType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename ModelType_::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename ModelType_::ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef ObserverHelper<IoInputType,
	MatrixType,
	VectorType,
	VectorWithOffsetType_,
	LeftRightSuperType> ObserverHelperType;
	typedef CorrelationsSkeleton<ObserverHelperType,ModelType_> CorrelationsSkeletonType;
	typedef OnePointCorrelations<ObserverHelperType> OnePointCorrelationsType;
	typedef TwoPointCorrelations<CorrelationsSkeletonType> TwoPointCorrelationsType;
	typedef FourPointCorrelations<CorrelationsSkeletonType> FourPointCorrelationsType;
	typedef MultiPointCorrelations<CorrelationsSkeletonType> MultiPointCorrelationsType;
	typedef PsimagLite::Profiling ProfilingType;

	static SizeType const GROW_RIGHT = CorrelationsSkeletonType::GROW_RIGHT;
	static SizeType const GROW_LEFT = CorrelationsSkeletonType::GROW_LEFT;
	static SizeType const DIAGONAL = CorrelationsSkeletonType::DIAGONAL;
	static SizeType const NON_DIAGONAL = CorrelationsSkeletonType::NON_DIAGONAL;

	enum {GS_VECTOR=ObserverHelperType::GS_VECTOR,
		  TIME_VECTOR=ObserverHelperType::TIME_VECTOR};
	enum {LEFT_BRACKET=ObserverHelperType::LEFT_BRACKET,
		  RIGHT_BRACKET=ObserverHelperType::RIGHT_BRACKET};

public:

	typedef ModelType_ ModelType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;

	Observer(IoInputType& io,
	         SizeType nf,
	         bool hasTimeEvolution,
	         const ModelType& model,
	         bool verbose=false)
	    : helper_(io,nf,model.params().nthreads,hasTimeEvolution,verbose),
	      verbose_(verbose),
	      onepoint_(helper_),
	      skeleton_(helper_,model,verbose),
	      twopoint_(helper_,skeleton_),
	      fourpoint_(helper_,skeleton_)
	{}

	SizeType size() const { return helper_.size(); }

	RealType time(SizeType threadId) const { return helper_.time(threadId); }

	SizeType site(SizeType threadId) const { return helper_.site(threadId); }

	SizeType marker(SizeType threadId) const { return helper_.marker(threadId); }

	void setPointer(SizeType threadId,SizeType x) { helper_.setPointer(threadId,x); }

	bool endOfData() const { return helper_.endOfData(); }

	// return true if
	// we're at site 1 or n-2
	bool isAtCorner(SizeType numberOfSites,SizeType threadId) const
	{
		bool es = (helper_.direction(threadId) == ProgramGlobals::EXPAND_SYSTEM);
		if (es && helper_.site(threadId) ==  numberOfSites-2) return true;
		if (!es && helper_.site(threadId) == 1) return true;
		return false;
	}

	void setBrackets(const PsimagLite::String& left,const PsimagLite::String& right)
	{
		helper_.setBrackets(bracketStringToNumber(left),
		                    bracketStringToNumber(right));
	}

	MatrixType correlations(const MatrixType& O1,
	                        const MatrixType& O2,
	                        int fermionicSign,
	                        SizeType rows,
	                        SizeType cols)
	{
		return twopoint_(O1,O2,fermionicSign,rows,cols);
	}

	template<typename SomeBracketType>
	void threePoint(const SomeBracketType& bracket,
	                SizeType rows,
	                SizeType cols)
	{
		SizeType site0 = bracket.site(0);
		std::cout<<"#site0="<<site0<<"\n";
		SizeType threadId = 0;
		MatrixType m0;
		MatrixType m1;
		MatrixType m2;
		crsMatrixToFullMatrix(m0,bracket.op(0).data);
		crsMatrixToFullMatrix(m1,bracket.op(1).data);
		crsMatrixToFullMatrix(m2,bracket.op(2).data);
		int fermionSign = bracket.op(0).fermionSign;
		for (SizeType site1 = 0; site1 < rows; ++site1) {
			for (SizeType site2 = 0; site2 < cols; ++site2) {
				try {
					typename MatrixType::value_type tmp = fourpoint_.threePoint('N',site0,m0,
					                                                            'N',site1,m1,
					                                                            'N',site2,m2,
					                                                            fermionSign,
					                                                            threadId);
					std::cout<<site1<<" "<<site2<<" "<<tmp<<"\n";
				} catch(std::exception&) {}
			}
		}
	}

	template<typename SomeBracketType>
	void fourPoint(const SomeBracketType& bracket,
	               SizeType rows,
	               SizeType cols)
	{
		SizeType site0 = bracket.site(0);
		SizeType site1 = bracket.site(1);
		std::cout<<"#site0="<<site0<<"\n#site1="<<site1<<"\n";
		if (site0 >= site1)
			throw PsimagLite::RuntimeError("FourPoint: needs ordered distinct points\n");

		SizeType threadId = 0;
		MatrixType m0;
		MatrixType m1;
		MatrixType m2;
		MatrixType m3;
		crsMatrixToFullMatrix(m0,bracket.op(0).data);
		crsMatrixToFullMatrix(m1,bracket.op(1).data);
		crsMatrixToFullMatrix(m2,bracket.op(2).data);
		crsMatrixToFullMatrix(m3,bracket.op(3).data);

		RealType fermionSign = bracket.op(0).fermionSign;
		typename FourPointCorrelationsType::SparseMatrixType O2gt;
		fourpoint_.firstStage(O2gt,'N',site0,m0,'N',site1,m1,fermionSign,threadId);

		for (SizeType site2 = site1+1; site2 < rows; ++site2) {
			for (SizeType site3 = site2+1; site3 < cols; ++site3) {
				typename MatrixType::value_type tmp = fourpoint_.
				        secondStage(O2gt,site1,'N',site2,m2,'N',site3,m3,fermionSign,threadId);
					std::cout<<site2<<" "<<site3<<" "<<tmp<<"\n";
			}
		}
	}

	template<typename SomeBracketType>
	VectorMatrixType ladder(const SomeBracketType& bracket,
	                        SizeType rows,
	                        SizeType,
	                        SizeType threadId)
	{
		SizeType rowsOver2 = utils::exactDivision(rows,2);

		MatrixType m0,m1,m2,m3;
		crsMatrixToFullMatrix(m0,bracket.op(0).data);
		crsMatrixToFullMatrix(m1,bracket.op(1).data);
		crsMatrixToFullMatrix(m2,bracket.op(2).data);
		crsMatrixToFullMatrix(m3,bracket.op(3).data);
		int f = bracket.op(0).fermionSign;
		VectorMatrixType v(4);
		for (SizeType i = 0; i < v.size(); ++i) v[i].resize(rowsOver2,rowsOver2);

		for (SizeType i = 0; i < rowsOver2; i += 2) {
			for (SizeType j = rowsOver2; j < rows; j += 2) {
				SizeType jj = j - rowsOver2;
				v[0](i,jj) = ladder_(m0,i,m1,i+1,m2,j,m3,j+1,f,threadId);
				v[1](i,jj) = f*ladder_(m0,i,m1,i+1,m3,j,m2,j+1,f,threadId);
				v[2](i,jj) = f*ladder_(m1,i,m0,i+1,m2,j,m3,j+1,f,threadId);
				v[3](i,jj) = ladder_(m1,i,m0,i+1,m3,j,m2,j+1,f,threadId);
			}
		}

		return v;
	}

	template<typename SomeModelType>
	void fourPointDeltas(MatrixType& fpd,
	                     const typename PsimagLite::Vector<SizeType>::Type& gammas,
	                     const SomeModelType& model)
	{
		if (gammas.size()!=4) {
			std::cerr<<"Observer: fourPointDeltas(...):  wrong number of gammas ";
			std::cerr<<" expected "<<4<<" got "<<gammas.size()<<"\n";
			throw PsimagLite::RuntimeError("Observer::fourPointDeltas(...)\n");
		}

		SizeType nsites = 2*fpd.n_row();
		assert(fpd.n_row()==fpd.n_col());

		SizeType hs = model.hilbertSize(0);
		SizeType nx = 0;
		while (hs) {
			hs>>=1;
			nx++;
		}

		nx /= 2;

		assert(fpd.n_row()>1);
		typedef std::pair<SizeType,SizeType> PairType;

		typename PsimagLite::Vector<PairType>::Type pairs;
		for (SizeType i=0;i<fpd.n_row();i++) {
			if (2*i+1>=nsites) continue;
			for (SizeType j=i+1;j<fpd.n_col();j++) {
				if (2*j+1>=nsites) continue;
				pairs.push_back(PairType(i,j));
			}
		}

		typedef Parallel4PointDs<ModelType,FourPointCorrelationsType>
		        Parallel4PointDsType;
		typedef PsimagLite::Parallelizer<Parallel4PointDsType> ParallelizerType;
		ParallelizerType threaded4PointDs(PsimagLite::Concurrency::npthreads,
		                                  PsimagLite::MPI::COMM_WORLD);

		Parallel4PointDsType helper4PointDs(fpd,fourpoint_,model,gammas,pairs);

		threaded4PointDs.loopCreate(pairs.size(),helper4PointDs);
	}

	template<typename ApplyOperatorType>
	FieldType onePoint(SizeType site,
	                   const typename ApplyOperatorType::OperatorType& A,
	                   typename ApplyOperatorType::BorderEnum corner)
	{
		return onepoint_.template operator()<ApplyOperatorType>(site,A,corner);
	}

	template<typename ApplyOperatorType>
	FieldType onePointHookForZero(SizeType site,
	                              const typename ApplyOperatorType::OperatorType& A,
	                              bool corner = false)
	{
		return onepoint_.template hookForZero<ApplyOperatorType>(site,A,corner);
	}

	template<typename VectorLikeType>
	typename PsimagLite::EnableIf
	<PsimagLite::IsVectorLike<VectorLikeType>::True,void>::Type
	multiCorrelations(VectorLikeType& result,
	                  const MatrixType& O,
	                  SizeType rows,
	                  SizeType cols)
	{
		SizeType nthreads = 1;
		MultiPointCorrelationsType multi(nthreads,helper_,skeleton_);
		multi(result,O,rows,cols);
	}

private:

	SizeType bracketStringToNumber(const PsimagLite::String& str) const
	{
		if (str == "gs") return GS_VECTOR;
		if (str == "time") return TIME_VECTOR;
		PsimagLite::String msg("Observer::bracketStringToNumber:");
		throw PsimagLite::RuntimeError(msg + " must be gs or time\n");
	}

	FieldType ladder_(const MatrixType& O1,
	                  SizeType i1,
	                  const MatrixType& O2,
	                  SizeType i2,
	                  const MatrixType& O3,
	                  SizeType j1,
	                  const MatrixType& O4,
	                  SizeType j2,
	                  int fermionicSign,
	                  SizeType threadId)
	{
		char mod = 'N';
		return fourpoint_(mod,
		                  i1,
		                  O1,
		                  mod,
		                  i2,
		                  O2,
		                  mod,
		                  j1,
		                  O3,
		                  mod,
		                  j2,
		                  O4,
		                  fermionicSign,
		                  threadId);
	}

	ObserverHelperType helper_;
	bool verbose_;
	OnePointCorrelationsType onepoint_;
	CorrelationsSkeletonType skeleton_;
	TwoPointCorrelationsType twopoint_;
	FourPointCorrelationsType fourpoint_;
};  //class Observer
} // namespace Dmrg

/*@}*/
#endif

