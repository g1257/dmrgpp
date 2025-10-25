/*
Copyright (c) 2009-2013-2019, UT-Battelle, LLC
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

/*! \file TargetingExpression.h
 * TBW FIXME TODO
 */

#ifndef TARGETING_EXPRESSION_H
#define TARGETING_EXPRESSION_H
#include "CanonicalExpression.h"
#include "GroupOfOneTimeEvolutions.h"
#include "Pvector.h"
#include "SpecForTargetingExpression.h"
#include "TargetingBase.h"
#include <iostream>
#include <stdexcept>
#include "KetForTargetingExpression.hh"

namespace Dmrg
{

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingExpression : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_>
{

	typedef TargetingBase<LanczosSolverType_, VectorWithOffsetType_> BaseType;
	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename BaseType::ModelType ModelType;
	typedef typename BaseType::CheckpointType CheckpointType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BasisType::QnType QnType;
	typedef typename TargetingCommonType::StageEnumType StageEnumType;
	typedef Pvector<typename VectorWithOffsetType_::value_type> PvectorType;
	typedef SpecForTargetingExpression<BaseType> SpecForTargetingExpressionType;
	typedef typename SpecForTargetingExpressionType::AlgebraType AlgebraType;
	typedef typename SpecForTargetingExpressionType::AssignAndDestroy AssignAndDestroyType;
	typedef PsimagLite::CanonicalExpression<SpecForTargetingExpressionType, AssignAndDestroyType>
	    CanonicalExpressionType;
	typedef AuxForTargetingExpression<BaseType> AuxForTargetingExpressionType;
	typedef typename TargetingCommonType::VectorRealType VectorRealType;
	typedef typename AuxForTargetingExpressionType::VectorStringType VectorStringType;
	typedef typename AuxForTargetingExpressionType::VectorVectorWithOffsetType
	    VectorVectorWithOffsetType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename AlgebraType::VectorSizeType VectorSizeType;
	typedef typename BaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef GroupOfOneTimeEvolutions<Pvectors<BaseType>> GroupOfOneTimeEvolutionsType;
	using ComplexOrRealType = typename ModelType::ComplexOrRealType;
	using KetType = KetForTargetingExpression<ComplexOrRealType>;
	using TermType = typename AlgebraType::TermType;

public:

	TargetingExpression(const LeftRightSuperType& lrs,
	    const CheckpointType& checkPoint,
	    const WaveFunctionTransfType& wft,
	    const QnType&,
	    InputValidatorType& io)
	    : BaseType(lrs, checkPoint, wft, 0)
	    , progress_("TargetingExpression")
	    , gsWeight_(0.3)
	    , gsWeightActual_(gsWeight_)
	    , pvectors_(io, this->common().aoe(), lrs)
	{
		io.readline(gsWeight_, "GsWeight=");
	}

	SizeType sites() const { return 0; }

	SizeType targets() const { return pvectors_.targets(); }

	RealType normSquared(SizeType i) const
	{
		return PsimagLite::real(this->tv(i) * this->tv(i));
	}

	RealType weight(SizeType i) const
	{
		assert(this->common().aoe().noStageIs(StageEnumType::DISABLED));
		assert(i < weights_.size());
		return weights_[i];
	}

	RealType gsWeight() const
	{
		return gsWeightActual_;
	}

	void evolve(const VectorRealType& energies,
	    ProgramGlobals::DirectionEnum direction,
	    const BlockType& block1,
	    const BlockType&,
	    SizeType loopNumber)
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE)
			return;

		this->common().setAllStagesTo(StageEnumType::WFT_NOADVANCE);
		assert(block1.size() == 1);
		const SizeType site = block1[0];
		const SizeType numberOfSites = this->common().aoe().model().superGeometry().numberOfSites();
		const SizeType total = this->common().aoe().tvs();
		assert(total <= pvectors_.targets());
		if (site != 0 && site + 1 != numberOfSites)
			this->common().aoeNonConst().wftSome(site, 0, total);

		assert(energies.size() > 0);
		computePvectors(direction, energies[0], site); // may alter the number of tvs

		computeAllWeights();
		this->common().printNormsAndWeights(gsWeightActual_, weights_);

		const bool doBorderIfBorder = true;
		auto testLambda = [this](const PsimagLite::GetBraOrKet& bra,
				      const PsimagLite::GetBraOrKet& ket) {
			if (!hasTimeEvolution(bra) || !hasTimeEvolution(ket))
				return;

			RealType braTime = getTimeForKet(bra);
			RealType ketTime = getTimeForKet(ket);
			const RealType globalTime = this->common().time();

			if (braTime != ketTime)
				err("BraTime= " + ttos(braTime) + " but ketTime= " + ttos(ketTime) + "\n");

			if (braTime != globalTime)
				err("Bra and kets have some time= " + ttos(braTime) + " but globalTime= " + ttos(globalTime) + "\n");
		};

		if (loopNumber >= this->model().params().finiteLoop.size() - 1) {
			if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				if (site >= numberOfSites - 1)
					return;
			} else {
				if (site < 1)
					return;
			}
		}

		this->common().cocoon(block1, direction, doBorderIfBorder, &testLambda); // in-situ

		// border trigger below

		// avoid self-triggers:
		if (site == 0 || site == numberOfSites - 1)
			return;

		if (site > 1 && site < numberOfSites - 2)
			return;

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			if (site == 1)
				return;
		} else {
			if (site == numberOfSites - 2)
				return;
		}

		SizeType x = (site == 1) ? 0 : numberOfSites - 1;
		BlockType block(1, x);
		evolve(energies, direction, block, block, loopNumber);
	}

	void read(typename TargetingCommonType::IoInputType& io,
	    PsimagLite::String prefix)
	{
		this->common().readGSandNGSTs(io, prefix, "Expression");
	}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	    PsimagLite::IoSelector::Out& io,
	    PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, prefix, block, "Expression");
	}

	bool hasTimeEvolution(const PsimagLite::GetBraOrKet& ket) const
	{
		return (ket.isPvector()) ? pvectors_.hasTimeEvolution(ket.pIndex()) : false;
	}

	RealType getTimeForKet(const PsimagLite::GetBraOrKet& ket) const
	{
		if (ket.isPvector())
			return timeEvolve_.getTimeForKet(ket.pIndex());

		if (ket.isRvector())
			throw PsimagLite::RuntimeError("R vectors cannot be tested\n");

		return 0;
	}

private:

	void computePvectors(ProgramGlobals::DirectionEnum dir, RealType Eg, SizeType site)
	{
		if (allOrigPvectorsDone()) {
			const SizeType tvs = this->common().aoe().tvs();
			if (tvs == pvectors_.origPvectors())
				return;
			if (tvs < pvectors_.origPvectors())
				err("TVS could not have decreased ?!\n");

			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg << "All user-provided P vectors finished";
			progress_.printline(msgg, std::cout);
			std::cerr << PsimagLite::AnsiColor::green;
			std::cerr << "All user-provided P vectors finished\n";
			std::cerr << PsimagLite::AnsiColor::reset;

			this->common().aoeNonConst().targetVectorsResize(pvectors_.origPvectors());

			return;
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "P0=" << pvectors_(0).lastName();
		progress_.printline(msgg, std::cout);

		CanonicalExpressionType canonicalExpression(opSpec_);
		SizeType total = pvectors_.targets();

		AuxForTargetingExpressionType aux(pvectors_, timeEvolve_, dir, Eg, site);
		AlgebraType opEmpty(aux);
		bool needsTrimming = false;
		PsimagLite::String allpvectors;
		for (SizeType i = 0; i < total; ++i) {

			if (pvectors_(i).isDone())
				continue;

			aux.setPindexOutput(i);
			AlgebraType tmp(aux);

			canonicalExpression(tmp, pvectors_(i).lastName(), opEmpty, aux);
			tmp.finalize();

			AlgebraType thispBefore(tmp);
			finalize(aux.tempVectors(), aux.tempNames(), i, thispBefore);
			PsimagLite::String thispAfter = pvectors_(i).lastName();

			int pIndexOrMinusOne = (tmp.size() == 1) ?  tmp.term(0).pIndex(): -1;
			if (pIndexOrMinusOne >= 0) {
				SizeType x = pIndexOrMinusOne;
				if (x != i) {
					VectorWithOffsetType_& v0 = this->tvNonConst(i);
					v0 = this->tv(x);
					v0 *= tmp.term(0).ket().factor();
				} else {
					std::cerr << "Ignoring self assignment P";
					std::cerr << i << "=P" << x << "\n";
				}

				pvectors_.setAsDone(i);
				continue;
			}

			if (!tmp.hasSummationKetAndNoMult()) {
				allpvectors += thispAfter;
				continue;
			}

			simplifyTerms(thispBefore);
			if (thispBefore.toString() != pvectors_(i).lastName()) {
				needsTrimming = true;
                                compressExpression(thispBefore);
				// checkNoUncompressedExists(compr);
				pvectors_.pushString(i, thispBefore.toString());
			}

			allpvectors += thispBefore.toString();
		}

		if (needsTrimming)
			pvectors_.trimPvectors(allpvectors);
	}

	void finalize(const VectorVectorWithOffsetType& tempVectors,
	    const VectorStringType& tempNames,
	    SizeType pVectorIndex,
	    const AlgebraType& tempExpr)
	{
		const SizeType ntemps = tempNames.size();

		if (ntemps == 0)
			return;

		// find tempNames_ in pVectors_ and trim tempNames_ accordingly
		VectorBoolType removed_(ntemps);
		VectorSizeType tempToP(ntemps, 10000);
		for (SizeType i = 0; i < ntemps; ++i) {
			int x = pvectors_.findInOrigNames(tempNames[i]);
			if (x < 0)
				continue;
			this->tvNonConst(x) = tempVectors[i];
			pvectors_.setAsDone(x);
			removed_[i] = true;
			tempToP[i] = x;
		}

		// these surviving tempNames_ need storage, add them
		for (SizeType i = 0; i < ntemps; ++i) {
			if (removed_[i])
				continue;
			int x = pvectors_.findInAnyNames(tempNames[i]);
			if (x >= 0)
				continue;
			auto lambda = [this, i, &tempToP, &tempNames](SizeType ind) {
				tempToP[i] = ind;
				return this->expandExpression(tempNames[i], tempToP);
			};

			pvectors_.createNew(tempVectors[i], lambda);
		}

		PsimagLite::String newpstring = compressExpression(tempExpr);
		const PsimagLite::String selfName = "|P" + ttos(pVectorIndex) + ">";
		if (newpstring == selfName)
			pvectors_.setAsDone(pVectorIndex);
		else
			pvectors_.pushString(pVectorIndex, newpstring);
	}

	bool allOrigPvectorsDone() const
	{
		for (SizeType i = 0; i < pvectors_.origPvectors(); ++i)
			if (!pvectors_(i).isDone())
				return false;
		return true;
	}

	// replace "|!m" + something ==> "|P" + number
	void compressExpression(AlgebraType& expr) const
	{
		SizeType n = expr.size();
		for (SizeType i = 0; i < n; ++i) {
			const KetType& ket = expr.term(i).ket();
			if (ket.kind() == KetType::Kind::M) {
				const SizeType ind = pvectors_.findPforThisExpression(ket.name());
				expr.setKet(i, "|P" + ttos(ind) + ">");
			}
		}
	}

	// replace "R" + i ==> "P" + tempToP[i]
	PsimagLite::String expandExpression(PsimagLite::String str,
	    const VectorSizeType& tempToP) const
	{
		SizeType i = 0;
		const SizeType len = str.length();
		if (len < 4)
			return str;
		PsimagLite::String expanded;
		for (; i < len; ++i) {
			if (i + 4 < len && str[i] == '|' && str[i + 1] == 'R') {
				SizeType j = i + 2;
				PsimagLite::String buffer;
				for (; j < len; ++j) {
					if (str[j] == '>')
						break;
					buffer += str[j];
				}

				const SizeType ind = PsimagLite::atoi(buffer);
				if (ind >= tempToP.size())
					err("tempToP.size() >= index\n");
				buffer = "|P" + ttos(tempToP[ind]);
				i = j + 1;
				expanded += buffer;
				continue;
			}

			expanded += str[i];
		}

		return expanded;
	}

	// replace |+PXpPY> ==> |PX>, update |PX>,
	// FIMXE: ask aoe destroying PY if no longer referencable
	void simplifyTerms(AlgebraType& expr)
	{
		SizeType n = expr.size();
		for (SizeType i = 0; i < n; ++i) {
			const KetType& ket = expr.term(i).ket();
			if (ket.kind() != KetType::Kind::S) {
				continue;
			}

			auto opaque = ket.fillSumStruct();

			pvectors_.sumPvectors(opaque[0].first,
			    opaque[0].second,
			    opaque[1].first,
			    opaque[1].second,
			    ket.name());

			expr.setKet(i,  "|P" + ttos(opaque[0].first) + ">");
		}
	}

	static void checkNoUncompressedExists(PsimagLite::String str)
	{
		SizeType i = 0;
		const SizeType len = str.length();
		if (len < 4)
			return;
		for (; i < len; ++i) {
			if (i + 4 < len && str.substr(i, 3) == "|!m")
				err("Uncompressed exists in " + str + "\n");
		}
	}

	void computeAllWeights()
	{
		const SizeType n = this->common().aoe().tvs();
		assert(n <= pvectors_.targets());
		weights_.resize(n);
		std::fill(weights_.begin(), weights_.end(), 0);
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			RealType norma = norm(this->tv(i));
			if (norma < 1e-6)
				continue;
			weights_[i] = pvectors_(i).weight() / norma;
			sum += weights_[i];
		}

		gsWeightActual_ = 1 - sum;

		if (gsWeightActual_ >= gsWeight_)
			return; // <--- EARLY EXIT HERE

		assert(sum > 1e-6);
		RealType factor = (1 - gsWeight_) / sum;
		for (SizeType i = 0; i < n; ++i)
			weights_[i] *= factor;

		gsWeightActual_ = gsWeight_;
	}

	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	RealType gsWeightActual_;
	VectorRealType weights_;
	SpecForTargetingExpressionType opSpec_;
	Pvectors<BaseType> pvectors_;
	GroupOfOneTimeEvolutionsType timeEvolve_;
}; // class TargetingExpression
} // namespace Dmrg
/*@}*/
#endif
