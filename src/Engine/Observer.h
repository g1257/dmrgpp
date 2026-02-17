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
#include "Concurrency.h"
#include "CorrelationsSkeleton.h"
#include "FourPointCorrelations.h"
#include "ManyPointAction.h"
#include "MultiPointCorrelations.h"
#include "ObserverHelper.h"
#include "OnePointCorrelations.h"
#include "Parallel4PointDs.h"
#include "Parallelizer.h"
#include "PredicateAwesome.h"
#include "TwoPointCorrelations.h"
#include "Utils.h"
#include "VectorWithOffset.h" // for operator*
#include "VectorWithOffsets.h" // for operator*

namespace Dmrg {

template <typename ObserverHelperType, typename ModelType_> class Observer {

public:

	using ModelType                  = ModelType_;
	using VectorWithOffsetType       = typename ObserverHelperType::VectorWithOffsetType;
	using IoInputType                = typename ObserverHelperType::IoInputType;
	using FieldType                  = typename VectorWithOffsetType::value_type;
	using VectorType                 = PsimagLite::SparseVector<FieldType>;
	using RealType                   = typename ModelType_::RealType;
	using MatrixType                 = PsimagLite::Matrix<FieldType>;
	using VectorMatrixType           = typename PsimagLite::Vector<MatrixType>::Type;
	using BasisWithOperatorsType     = typename ModelType_::BasisWithOperatorsType;
	using SparseMatrixType           = typename BasisWithOperatorsType::SparseMatrixType;
	using LeftRightSuperType         = typename ModelType_::ModelHelperType::LeftRightSuperType;
	using ParametersType             = typename ModelType_::ParametersType;
	using CorrelationsSkeletonType   = CorrelationsSkeleton<ObserverHelperType, ModelType_>;
	using OnePointCorrelationsType   = OnePointCorrelations<ObserverHelperType, ModelType_>;
	using TwoPointCorrelationsType   = TwoPointCorrelations<CorrelationsSkeletonType>;
	using FourPointCorrelationsType  = FourPointCorrelations<CorrelationsSkeletonType>;
	using MultiPointCorrelationsType = MultiPointCorrelations<CorrelationsSkeletonType>;
	using BraketType                 = typename CorrelationsSkeletonType::BraketType;
	using Parallel4PointDsType       = Parallel4PointDs<ModelType, FourPointCorrelationsType>;
	using VectorStringType           = PsimagLite::Vector<PsimagLite::String>::Type;
	using ManyPointActionType        = ManyPointAction;

	Observer(IoInputType&      io,
	         SizeType          start,
	         SizeType          nf,
	         SizeType          trail,
	         const ModelType_& model)
	    : helper_(io,
	              start,
	              nf,
	              trail,
	              model.params().options.isSet("keepLegacyBugs"),
	              model.params().options.isSet("observeReadOnDemand"))
	    , onepoint_(helper_, model)
	    , skeleton_(helper_, model, true)
	    , twopoint_(skeleton_)
	    , fourpoint_(skeleton_)
	{ }

	const ObserverHelperType& helper() const { return helper_; }

	// return true if
	// we're at site 1 or n-2
	bool isAtCorner(SizeType numberOfSites, SizeType ptr) const
	{
		const bool es
		    = (helper_.direction(ptr) == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		if (es && helper_.site(ptr) == numberOfSites - 2)
			return true;
		return (!es && helper_.site(ptr) == 1);
	}

	void twoPoint(MatrixType&                storage,
	              const BraketType&          braket,
	              bool                       needsPrinting,
	              const ManyPointActionType& action) const
	{
		assert(braket.points() == 2);

		SizeType flag = 0;

		try {
			braket.site(0);
			flag |= 1;
		} catch (std::exception&) { }

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) { }

		ProgramGlobals::FermionOrBosonEnum fermionSign = braket.op(0).fermionOrBoson();

		SizeType site1 = 0;
		SizeType sites = storage.n_col();

		switch (flag) {

		case 0: // no sites given
			twopoint_(storage, braket, fermionSign, braket.bra(), braket.ket(), action);
			if (needsPrinting)
				std::cout << (storage);

			break;

		case 1: // first site given
			if (action.nonTrivial())
				err("For non trivial action, give no sites\n");

			for (site1 = braket.site(0); site1 < sites; ++site1) {
				storage(braket.site(0), site1)
				    = twopoint_.calcCorrelation(braket.site(0),
				                                site1,
				                                braket,
				                                fermionSign,
				                                braket.bra(),
				                                braket.ket());
				if (needsPrinting)
					std::cout << storage(braket.site(0), site1) << " ";
			}

			if (needsPrinting)
				std::cout << "\n";

			break;

		case 3:
			if (action.nonTrivial())
				err("For non trivial action, give no sites\n");

			storage(braket.site(0), braket.site(1))
			    = twopoint_.calcCorrelation(braket.site(0),
			                                braket.site(1),
			                                braket,
			                                fermionSign,
			                                braket.bra(),
			                                braket.ket());
			if (needsPrinting)
				std::cout << storage(braket.site(0), braket.site(1)) << "\n";
			;

			break;

		default:
			err("twopoint: Give no site, first site, or all sites\n");
		}
	}

	void twoPoint(MatrixType&                        m,
	              const BraketType                   braket,
	              ProgramGlobals::FermionOrBosonEnum fermionicSign,
	              PsimagLite::String                 bra,
	              PsimagLite::String                 ket,
	              const ManyPointActionType&         action) const
	{
		twopoint_(m, braket, fermionicSign, bra, ket, action);
	}

	FieldType
	threePoint(const BraketType& braket, SizeType rows, SizeType cols, bool needsPrinting) const
	{
		assert(braket.points() == 3);

		SizeType flag = 0;
		try {
			braket.site(0);
			flag |= 1;
		} catch (std::exception&) { }

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) { }

		try {
			braket.site(2);
			flag |= 4;
		} catch (std::exception&) { }

		if (flag != 0 && flag != 7 && flag != 1) {
			PsimagLite::String str("threePoint: ");
			str += "Give no site, first site, or all sites\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (flag == 7) {
			typename MatrixType::value_type tmp = fourpoint_.threePoint(
			    braket.site(0), braket.site(1), braket.site(2), braket);
			if (!needsPrinting)
				return tmp;
			std::cout << "Fixed all sites\n";
			std::cout << braket.site(0) << " ";
			std::cout << braket.site(1) << " " << braket.site(2) << "  " << tmp << "\n";
			return tmp;
		}

		if (!needsPrinting)
			err("Observer::threePoint with !needsPrinting only for all sites fixed\n");

		if (flag == 1) {
			SizeType site0 = braket.site(0);
			std::cout << "Fixed site0= " << site0 << "\n";
			for (SizeType site1 = site0 + 1; site1 < rows; ++site1) {
				for (SizeType site2 = site1 + 1; site2 < cols; ++site2) {
					typename MatrixType::value_type tmp
					    = fourpoint_.threePoint(site0, site1, site2, braket);
					std::cout << site1 << " " << site2 << "  " << tmp << "\n";
				}
			}
		}

		assert(flag == 0);
		for (SizeType site0 = 0; site0 < rows; ++site0) {
			for (SizeType site1 = site0 + 1; site1 < rows; ++site1) {
				for (SizeType site2 = site1 + 1; site2 < cols; ++site2) {

					typename MatrixType::value_type tmp
					    = fourpoint_.threePoint(site0, site1, site2, braket);
					std::cout << site0 << " " << site1 << " " << site2 << "  "
					          << tmp << "\n";
				}
			}
		}

		return 0;
	}

	const FourPointCorrelationsType& fourpoint() const { return fourpoint_; }

	void fourPoint(const BraketType&      braket,
	               SizeType               rows,
	               SizeType               cols,
	               const ManyPointAction& myaction)
	{
		assert(braket.points() == 4);

		SizeType flag = 0;
		try {
			braket.site(0);
			flag |= 1;
		} catch (std::exception&) { }

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) { }

		try {
			braket.site(2);
			flag |= 4;
		} catch (std::exception&) { }

		try {
			braket.site(3);
			flag |= 8;
		} catch (std::exception&) { }

		if (flag != 0 && flag != 3 && flag != 15) {
			PsimagLite::String str("fourPoint: ");
			str += "Give no site, first two sites, or all sites\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (flag == 15) {
			std::cout << "Fixed all sites\n";
			BraketType braketOrdered = braket;
			int        sign          = orderBraketIfNeeded(braketOrdered);
			SizeType   site0         = braketOrdered.site(0);
			SizeType   site1         = braketOrdered.site(1);
			typename FourPointCorrelationsType::SparseMatrixType O2gt;
			const bool                                           finalTransform = true;
			fourpoint_.firstStage(
			    O2gt, 'N', site0, 'N', site1, braketOrdered, 0, 1, finalTransform);
			typename MatrixType::value_type tmp
			    = fourpoint_.secondStage(O2gt,
			                             site1,
			                             'N',
			                             braketOrdered.site(2),
			                             'N',
			                             braketOrdered.site(3),
			                             braketOrdered,
			                             2,
			                             3);

			// sign due to reordering
			tmp *= static_cast<RealType>(sign);
			std::cout << braket.site(0) << " " << braket.site(1) << " ";
			std::cout << braket.site(2) << " " << braket.site(3) << "  " << tmp << "\n";
			return;
		}

		if (flag == 3) {
			SizeType site0 = braket.site(0);
			SizeType site1 = braket.site(1);
			std::cout << "Fixed site0= " << site0 << "\n";
			std::cout << "Fixed site1= " << site1 << "\n";
			typename FourPointCorrelationsType::SparseMatrixType O2gt;
			const bool                                           finalTransform = true;
			fourpoint_.firstStage(
			    O2gt, 'N', site0, 'N', site1, braket, 0, 1, finalTransform);

			for (SizeType site2 = site1 + 1; site2 < rows; ++site2) {
				for (SizeType site3 = site2 + 1; site3 < cols; ++site3) {
					typename MatrixType::value_type tmp
					    = fourpoint_.secondStage(
					        O2gt, site1, 'N', site2, 'N', site3, braket, 2, 3);
					std::cout << site2 << " " << site3 << " " << tmp << "\n";
				}
			}

			return;
		}

		assert(flag == 0);
		const bool finalTransform = true;
		for (SizeType site0 = 0; site0 < rows; ++site0) {
			for (SizeType site1 = site0 + 1; site1 < cols; ++site1) {
				for (SizeType site2 = site1 + 1; site2 < rows; ++site2) {
					for (SizeType site3 = site2 + 1; site3 < cols; ++site3) {
						if (!myaction(site0, site1, site2, site3))
							continue;
						typename FourPointCorrelationsType::SparseMatrixType
						    O2gt;
						fourpoint_.firstStage(O2gt,
						                      'N',
						                      site0,
						                      'N',
						                      site1,
						                      braket,
						                      0,
						                      1,
						                      finalTransform);
						typename MatrixType::value_type tmp
						    = fourpoint_.secondStage(O2gt,
						                             site1,
						                             'N',
						                             site2,
						                             'N',
						                             site3,
						                             braket,
						                             2,
						                             3);
						std::cout << site0 << " " << site1 << " ";
						std::cout << site2 << " " << site3 << " " << tmp
						          << "\n";
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
		} catch (std::exception&) { }

		try {
			braket.site(1);
			flag |= 2;
		} catch (std::exception&) { }

		try {
			braket.site(2);
			flag |= 4;
		} catch (std::exception&) { }

		try {
			braket.site(3);
			flag |= 8;
		} catch (std::exception&) { }

		if (flag != 15) {
			PsimagLite::String str("fourPoint: ");
			str += "Give all sites\n";
			throw PsimagLite::RuntimeError(str);
		}

		FieldType tmp = fourpoint_.anyPoint(braket);

		if (!needsPrinting)
			return tmp; // <<-- early exit

		std::cout << "Fixed all sites\n";
		for (SizeType i = 0; i < braket.points(); ++i)
			std::cout << braket.site(i) << " ";

		std::cout << tmp << "\n";

		return tmp;
	}

	void fourPointDeltas(MatrixType&                                        fpd,
	                     const typename PsimagLite::Vector<SizeType>::Type& gammas,
	                     const ModelType&                                   model)
	{
		if (gammas.size() != 4) {
			std::cerr << "Observer: fourPointDeltas(...):  wrong number of gammas ";
			std::cerr << " expected " << 4 << " got " << gammas.size() << "\n";
			throw PsimagLite::RuntimeError("Observer::fourPointDeltas(...)\n");
		}

		SizeType nsites = 2 * fpd.n_row();
		assert(fpd.n_row() == fpd.n_col());

		SizeType hs = model.hilbertSize(0);
		while (hs) {
			hs >>= 1;
		}

		assert(fpd.n_row() > 1);
		using PairType = std::pair<SizeType, SizeType>;

		typename PsimagLite::Vector<PairType>::Type pairs;
		for (SizeType i = 0; i < fpd.n_row(); i++) {
			if (2 * i + 1 >= nsites)
				continue;
			for (SizeType j = i + 1; j < fpd.n_col(); j++) {
				if (2 * j + 1 >= nsites)
					continue;
				pairs.push_back(PairType(i, j));
			}
		}

		using ParallelizerType = PsimagLite::Parallelizer<Parallel4PointDsType>;
		ParallelizerType threaded4PointDs(PsimagLite::Concurrency::codeSectionParams);

		Parallel4PointDsType helper4PointDs(
		    fpd, fourpoint_, model, gammas, pairs, Parallel4PointDsType::MODE_NORMAL);

		threaded4PointDs.loopCreate(helper4PointDs);
	}

	template <typename ApplyOperatorType>
	FieldType onePoint(SizeType                                        ptr,
	                   const typename ApplyOperatorType::OperatorType& A,
	                   SizeType                                        site,
	                   typename ApplyOperatorType::BorderEnum          corner,
	                   const PsimagLite::GetBraOrKet&                  bra,
	                   const PsimagLite::GetBraOrKet&                  ket) const
	{
		return onepoint_.template operator()<ApplyOperatorType>(
		    ptr, A, site, corner, bra, ket);
	}

	template <typename ApplyOperatorType>
	FieldType onePointHookForZero(SizeType                                        ptr,
	                              const typename ApplyOperatorType::OperatorType& A,
	                              SizeType                                        splitSize,
	                              const PsimagLite::GetBraOrKet&                  bra,
	                              const PsimagLite::GetBraOrKet&                  ket) const
	{
		return onepoint_.template hookForZero<ApplyOperatorType>(
		    ptr, A, splitSize, bra, ket);
	}

	template <typename VectorLikeType>
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<VectorLikeType>::True, void>::Type
	multiCorrelations(VectorLikeType&         result,
	                  const SparseMatrixType& O,
	                  SizeType                rows,
	                  SizeType                cols,
	                  PsimagLite::String      bra,
	                  PsimagLite::String      ket)
	{
		MultiPointCorrelationsType multi(skeleton_);
		multi(result, O, rows, cols, bra, ket);
	}

private:

	static int orderBraketIfNeeded(BraketType& braket)
	{
		SizeType nsites = braket.points();
		if (nsites == 1)
			return 1;

		std::vector<SizeType> sites(nsites);
		for (SizeType i = 0; i < nsites; ++i) {
			sites[i] = braket.site(i);
		}

		PsimagLite::Sort<std::vector<SizeType>> sort;
		std::vector<SizeType>                   sitesSorted = sites;
		std::vector<SizeType>                   iperm(sites.size());
		sort.sort(sitesSorted, iperm);

		if (sites == sitesSorted)
			return 1;

		// throws if mixed fermion or boson
		ProgramGlobals::FermionOrBosonEnum statistics = getStatistics(braket);

		braket.reorder(iperm);

		return (statistics == ProgramGlobals::FermionOrBosonEnum::BOSON)
		    ? 1
		    : parityOfPermutation(iperm);
	}

	static ProgramGlobals::FermionOrBosonEnum getStatistics(const BraketType& braket)
	{
		SizeType nsites = braket.points();
		if (nsites == 0)
			return ProgramGlobals::FermionOrBosonEnum::BOSON;

		ProgramGlobals::FermionOrBosonEnum statistics = braket.op(0).fermionOrBoson();
		for (SizeType i = 1; i < nsites; ++i) {
			ProgramGlobals::FermionOrBosonEnum tmp = braket.op(i).fermionOrBoson();
			if (tmp != statistics)
				err("Inconsistent statistics for braket\n");
		}

		return statistics;
	}

	static int parityOfPermutation(const std::vector<SizeType>& perm)
	{
		return (swapCountSmall(perm) & 1) ? -1 : 1;
	}

	// Inspired by
	// https://stackoverflow.com/questions/20702782/efficiently-determine-the-parity-of-a-permutation
	static SizeType swapCountSmall(const std::vector<SizeType>& perm)
	{
		SizeType          n     = perm.size();
		SizeType          swaps = 0;
		unsigned long int seen  = 0;
		for (SizeType i = 0; i < n; ++i) {
			unsigned long int mask = (1L << i);
			if ((seen & mask) != 0)
				continue;
			seen |= mask;
			for (SizeType j = perm[i]; (seen & (1L << j)) == 0; j = perm[j]) {
				seen |= (1L << j);
				++swaps;
			}
		}

		return swaps;
	}

	const ObserverHelperType        helper_;
	const OnePointCorrelationsType  onepoint_;
	const CorrelationsSkeletonType  skeleton_;
	const TwoPointCorrelationsType  twopoint_;
	const FourPointCorrelationsType fourpoint_;
}; // class Observer
} // namespace Dmrg

/*@}*/
#endif
