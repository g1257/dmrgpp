/*
Copyright (c) 2008-2015, UT-Battelle, LLC
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
#include "Parallel4PointDs.h"
#include "MultiPointCorrelations.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "Utils.h"
#include "ExpressionCalculator.h"
#include "PredicateAwesome.h"

namespace Dmrg {

template<typename VectorWithOffsetType_,typename ModelType_,typename IoInputType>
class Observer {

public:

	typedef typename VectorWithOffsetType_::value_type FieldType;
	typedef PsimagLite::SparseVector<FieldType> VectorType;
	typedef typename ModelType_::RealType RealType;
	typedef PsimagLite::Matrix<FieldType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename ModelType_::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename ModelType_::ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelType_::ParametersType ParametersType;
	typedef ObserverHelper<IoInputType, MatrixType, VectorType, VectorWithOffsetType_,
	LeftRightSuperType> ObserverHelperType;
	typedef CorrelationsSkeleton<ObserverHelperType,ModelType_> CorrelationsSkeletonType;
	typedef OnePointCorrelations<ObserverHelperType> OnePointCorrelationsType;
	typedef TwoPointCorrelations<CorrelationsSkeletonType> TwoPointCorrelationsType;
	typedef FourPointCorrelations<CorrelationsSkeletonType> FourPointCorrelationsType;
	typedef MultiPointCorrelations<CorrelationsSkeletonType> MultiPointCorrelationsType;
	typedef typename CorrelationsSkeletonType::BraketType BraketType;
	typedef ModelType_ ModelType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef Parallel4PointDs<ModelType,FourPointCorrelationsType> Parallel4PointDsType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::PredicateAwesome<> PredicateAwesomeType;

	class ManyPointAction {

	public:

		typedef PsimagLite::ExpressionCalculator<int> ExpressionCalculatorType;
		typedef PsimagLite::PrepassData<int> PrepassDataType;
		typedef PsimagLite::ExpressionPrepass<PrepassDataType> ExpressionPrepassType;

		ManyPointAction(bool hasNonTrivialAction,
		                PsimagLite::String actionString)
		    : nonTrivial_(hasNonTrivialAction),
		      actionString_(actionString)
		{}

		bool operator()(SizeType s0, SizeType s1, SizeType s2, SizeType s3) const
		{
			if (!nonTrivial_) return true;

			PredicateAwesomeType pred(actionString_, '~');
			return pred.isTrue("%0", s0, "%1", s1, "%2", s2, "%3", s3);
		}

	private:

		bool nonTrivial_;
		PsimagLite::String actionString_;
	};

	Observer(IoInputType& io,
	         SizeType start,
	         SizeType nf,
	         SizeType trail,
	         const ParametersType& params)
	    : helper_(io,
	              start,
	              nf,
	              trail,
	              !params.options.isSet("fixLegacyBugs")),
	      onepoint_(helper_),
	      skeleton_(helper_, true),
	      twopoint_(skeleton_),
	      fourpoint_(skeleton_)
	{}

	const ObserverHelperType& helper() const { return helper_; }

	// return true if
	// we're at site 1 or n-2
	bool isAtCorner(SizeType numberOfSites, SizeType ptr) const
	{
		const bool es = (helper_.direction(ptr) == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		if (es && helper_.site(ptr) ==  numberOfSites - 2) return true;
		return (!es && helper_.site(ptr) == 1);
	}

	void twoPoint(MatrixType& storage, const BraketType& braket) const
	{
		assert(braket.points() == 2);

		SizeType flag = 0;

		try {
			braket.site(0);
			flag |= 1;
		} catch (std::exception&) {}

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) {}

		SparseMatrixType m0 = braket.op(0).getCRS();
		SparseMatrixType m1 = braket.op(1).getCRS();
		ProgramGlobals::FermionOrBosonEnum fermionSign = braket.op(0).fermionOrBoson();

		SizeType site1 = 0;
		SizeType sites= storage.n_col();
		PsimagLite::String str("twopoint: Give no site, first site, or all sites\n");

		switch (flag) {

		case 0: // no sites given
			return twopoint_(storage, m0, m1, fermionSign, braket.bra(), braket.ket());

		case 1: //first site given
			for (site1 = braket.site(0); site1 < sites; ++site1)
				storage(braket.site(0),site1) = twopoint_.calcCorrelation(braket.site(0),
				                                                          site1,
				                                                          braket.op(0).getCRS(),
				                                                          braket.op(1).getCRS(),
				                                                          fermionSign,
				                                                          braket.bra(),
				                                                          braket.ket());
			return;

		case 3:
			storage(braket.site(0),braket.site(1)) = twopoint_.calcCorrelation(braket.site(0),
			                                                                   braket.site(1),
			                                                                   braket.op(0).getCRS(),
			                                                                   braket.op(1).getCRS(),
			                                                                   fermionSign,
			                                                                   braket.bra(),
			                                                                   braket.ket());
			return;

		default:
			err(str);
		}
	}

	void twoPoint(MatrixType& m,
	              const SparseMatrixType& O1,
	              const SparseMatrixType& O2,
	              ProgramGlobals::FermionOrBosonEnum fermionicSign,
	              PsimagLite::String bra,
	              PsimagLite::String ket) const
	{
		twopoint_(m, O1, O2, fermionicSign, bra, ket);
	}

	FieldType threePoint(const BraketType& braket,
	                     SizeType rows,
	                     SizeType cols,
	                     bool needsPrinting) const
	{
		assert(braket.points() == 3);

		SizeType flag = 0;
		try {
			braket.site(0);
			flag |= 1;
		} catch (std::exception&) {}

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) {}

		try {
			braket.site(2);
			flag |= 4;
		} catch (std::exception&) {}

		if (flag != 0 && flag != 7 && flag != 1) {
			PsimagLite::String str("threePoint: ");
			str += "Give no site, first site, or all sites\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (flag == 7) {
			typename MatrixType::value_type tmp = fourpoint_.threePoint(braket.site(0),
			                                                            braket.site(1),
			                                                            braket.site(2),
			                                                            braket);
			if (!needsPrinting) return tmp;
			std::cout<<"Fixed all sites\n";
			std::cout<<braket.site(0)<<" ";
			std::cout<<braket.site(1)<<" "<<braket.site(2)<<"  "<<tmp<<"\n";
			return tmp;
		}

		if (!needsPrinting)
			err("Observer::threePoint with !needsPrinting only for all sites fixed\n");

		if (flag == 1) {
			SizeType site0 = braket.site(0);
			std::cout<<"Fixed site0= "<<site0<<"\n";
			for (SizeType site1 = site0+1; site1 < rows; ++site1) {
				for (SizeType site2 = site1+1; site2 < cols; ++site2) {
					typename MatrixType::value_type tmp = fourpoint_.threePoint(site0,
					                                                            site1,
					                                                            site2,
					                                                            braket);
					std::cout<<site1<<" "<<site2<<"  "<<tmp<<"\n";

				}
			}
		}

		assert(flag == 0);
		for (SizeType site0 = 0; site0 < rows; ++site0) {
			for (SizeType site1 = site0+1; site1 < rows; ++site1) {
				for (SizeType site2 = site1+1; site2 < cols; ++site2) {

					typename MatrixType::value_type tmp = fourpoint_.threePoint(site0,
					                                                            site1,
					                                                            site2,
					                                                            braket);
					std::cout<<site0<<" "<<site1<<" "<<site2<<"  "<<tmp<<"\n";
				}
			}
		}

		return 0;
	}

	const FourPointCorrelationsType& fourpoint() const {return fourpoint_; }

	void fourPoint(const BraketType& braket,
	               SizeType rows,
	               SizeType cols,
	               const ManyPointAction& myaction)
	{
		assert(braket.points() == 4);

		SizeType flag = 0;
		try {
			braket.site(0);
			flag |= 1;
		} catch (std::exception&) {}

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) {}

		try {
			braket.site(2);
			flag |= 4;
		} catch (std::exception&) {}

		try {
			braket.site(3);
			flag |= 8;
		} catch (std::exception&) {}

		if (flag != 0 && flag != 3 && flag != 15) {
			PsimagLite::String str("fourPoint: ");
			str += "Give no site, first two sites, or all sites\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (flag == 15) {
			std::cout<<"Fixed all sites\n";
			SizeType site0 = braket.site(0);
			SizeType site1 = braket.site(1);
			typename FourPointCorrelationsType::SparseMatrixType O2gt;
			fourpoint_.firstStage(O2gt,'N',site0,'N',site1,braket,0,1);
			typename MatrixType::value_type tmp = fourpoint_.secondStage(O2gt,
			                                                             site1,
			                                                             'N',braket.site(2),
			                                                             'N',braket.site(3),
			                                                             braket,
			                                                             2,
			                                                             3);
			std::cout<<site0<<" "<<site1<<" ";
			std::cout<<braket.site(2)<<" "<<braket.site(3)<<"  "<<tmp<<"\n";
			return;
		}

		if (flag == 3) {
			SizeType site0 = braket.site(0);
			SizeType site1 = braket.site(1);
			std::cout<<"Fixed site0= "<<site0<<"\n";
			std::cout<<"Fixed site1= "<<site1<<"\n";
			typename FourPointCorrelationsType::SparseMatrixType O2gt;
			fourpoint_.firstStage(O2gt,'N',site0,'N',site1,braket,0,1);

			for (SizeType site2 = site1+1; site2 < rows; ++site2) {
				for (SizeType site3 = site2+1; site3 < cols; ++site3) {
					typename MatrixType::value_type tmp = fourpoint_.secondStage(O2gt,
					                                                             site1,
					                                                             'N',
					                                                             site2,
					                                                             'N',
					                                                             site3,
					                                                             braket,
					                                                             2,
					                                                             3);
					std::cout<<site2<<" "<<site3<<" "<<tmp<<"\n";
				}
			}

			return;
		}

		assert(flag == 0);
		for (SizeType site0 = 0; site0 < rows; ++site0) {
			for (SizeType site1 = site0+1; site1 < cols; ++site1) {
				for (SizeType site2 = site1+1; site2 < rows; ++site2) {
					for (SizeType site3 = site2+1; site3 < cols; ++site3) {
						if (!myaction(site0, site1, site2, site3)) continue;
						typename FourPointCorrelationsType::SparseMatrixType O2gt;
						fourpoint_.firstStage(O2gt,'N',site0,'N',site1,braket,0,1);
						typename MatrixType::value_type tmp = fourpoint_.secondStage(O2gt,
						                                                             site1,
						                                                             'N',
						                                                             site2,
						                                                             'N',
						                                                             site3,
						                                                             braket,
						                                                             2,
						                                                             3);
						std::cout<<site0<<" "<<site1<<" ";
						std::cout<<site2<<" "<<site3<<" "<<tmp<<"\n";
					}
				}
			}
		}
	}

	FieldType anyPoint(const BraketType& braket, bool needsPrinting) const
	{
		assert(braket.points() >= 4);

		SizeType flag = 0;
		try {
			braket.site(0);
			flag |= 1;
		} catch (std::exception&) {}

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) {}

		try {
			braket.site(2);
			flag |= 4;
		} catch (std::exception&) {}

		try {
			braket.site(3);
			flag |= 8;
		} catch (std::exception&) {}

		if (flag != 15) {
			PsimagLite::String str("fourPoint: ");
			str += "Give all sites\n";
			throw PsimagLite::RuntimeError(str);
		}

		FieldType tmp = fourpoint_.anyPoint(braket);

		if (!needsPrinting) return tmp; // <<-- early exit

		std::cout<<"Fixed all sites\n";
		for (SizeType i = 0; i < braket.points(); ++i)
			std::cout<<braket.site(i)<<" ";

		std::cout<<tmp<<"\n";

		return tmp;
	}

	void fourPointDeltas(MatrixType& fpd,
	                     const typename PsimagLite::Vector<SizeType>::Type& gammas,
	                     const ModelType& model)
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


		typedef PsimagLite::Parallelizer<Parallel4PointDsType> ParallelizerType;
		ParallelizerType threaded4PointDs(PsimagLite::Concurrency::codeSectionParams);

		Parallel4PointDsType helper4PointDs(fpd,
		                                    fourpoint_,
		                                    model,
		                                    gammas,
		                                    pairs,
		                                    Parallel4PointDsType::MODE_NORMAL);

		threaded4PointDs.loopCreate(helper4PointDs);
	}

	template<typename ApplyOperatorType>
	FieldType onePoint(SizeType site,
	                   const typename ApplyOperatorType::OperatorType& A,
	                   typename ApplyOperatorType::BorderEnum corner,
	                   const PsimagLite::GetBraOrKet& bra,
	                   const PsimagLite::GetBraOrKet& ket) const
	{
		return onepoint_.template operator()<ApplyOperatorType>(site, A, corner, bra, ket);
	}

	template<typename ApplyOperatorType>
	FieldType onePointHookForZero(SizeType site,
	                              const typename ApplyOperatorType::OperatorType& A,
	                              const PsimagLite::GetBraOrKet& bra,
	                              const PsimagLite::GetBraOrKet& ket) const
	{
		return onepoint_.template hookForZero<ApplyOperatorType>(site, A, bra, ket);
	}

	template<typename VectorLikeType>
	typename PsimagLite::EnableIf
	<PsimagLite::IsVectorLike<VectorLikeType>::True,void>::Type
	multiCorrelations(VectorLikeType& result,
	                  const SparseMatrixType& O,
	                  SizeType rows,
	                  SizeType cols,
	                  PsimagLite::String bra,
	                  PsimagLite::String ket)
	{
		MultiPointCorrelationsType multi(skeleton_);
		multi(result, O, rows, cols, bra, ket);
	}

private:

	const ObserverHelperType helper_;
	const OnePointCorrelationsType onepoint_;
	const CorrelationsSkeletonType skeleton_;
	const TwoPointCorrelationsType twopoint_;
	const FourPointCorrelationsType fourpoint_;
};  //class Observer
} // namespace Dmrg

/*@}*/
#endif

