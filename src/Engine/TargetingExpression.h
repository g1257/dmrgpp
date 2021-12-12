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
#include <iostream>
#include "TargetingBase.h"
#include <stdexcept>
#include "Pvector.h"
#include "SpecForTargetingExpression.h"
#include "CanonicalExpression.h"
#include "GroupOfOneTimeEvolutions.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingExpression : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

	typedef TargetingBase<LanczosSolverType_,VectorWithOffsetType_> BaseType;
	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename BaseType::ModelType ModelType;
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
	typedef GroupOfOneTimeEvolutions<Pvectors<BaseType> > GroupOfOneTimeEvolutionsType;

public:

	TargetingExpression(const LeftRightSuperType& lrs,
	                    const ModelType& model,
	                    const WaveFunctionTransfType& wft,
	                    const QnType&,
	                    InputValidatorType& io)
	    : BaseType(lrs, model, wft, 0),
	      progress_("TargetingExpression"),
	      gsWeight_(0.3),
	      pvectors_(io, this->common().aoe(), lrs)
	{
		io.readline(gsWeight_, "GsWeight=");
	}

	SizeType sites() const { return 0; }

	SizeType targets() const { return pvectors_.targets(); }

	RealType normSquared(SizeType i) const
	{
		return PsimagLite::real(this->tv(i)*
		                        this->tv(i));
	}

	RealType weight(SizeType i) const
	{
		assert(this->common().aoe().noStageIs(StageEnumType::DISABLED));
		VectorRealType weights;
		RealType gsWeight = 0;
		computeAllWeights(gsWeight, weights);
		assert(i < weights.size());
		return weights[i];
	}

	RealType gsWeight() const
	{
		VectorRealType weights;
		RealType gsWeight = 0;
		computeAllWeights(gsWeight, weights);
		return gsWeight;
	}

	void evolve(const VectorRealType& energies,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType&,
	            SizeType)
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE) return;

		this->common().setAllStagesTo(StageEnumType::WFT_NOADVANCE);
		assert(block1.size() == 1);
		const SizeType site = block1[0];
		const SizeType total = this->common().aoe().tvs();
		assert(total <= pvectors_.targets());
		this->common().aoeNonConst().wftSome(site, 0, total);

		assert(energies.size() > 0);
		computePvectors(direction, energies[0]); // may alter the number of tvs

		VectorRealType weight;
		RealType gsWeight = 0;
		computeAllWeights(gsWeight, weight);
		this->common().printNormsAndWeights(gsWeight, weight);

		const bool doBorderIfBorder = true;
		auto testLambda = [this](const PsimagLite::GetBraOrKet& bra,
		        const PsimagLite::GetBraOrKet& ket)
		{
			if (!hasTimeEvolution(bra) || !hasTimeEvolution(ket)) return;

			RealType braTime = getTimeForKet(bra);
			RealType ketTime = getTimeForKet(ket);
			if (braTime == ketTime) return;
			err("BraTime= " + ttos(braTime) + " but ketTime= " + ttos(ketTime) + "\n");
		};

		this->common().cocoon(block1, direction, doBorderIfBorder, &testLambda); // in-situ
	}

	void read(typename TargetingCommonType::IoInputType& io,
	          PsimagLite::String prefix)
	{
		this->common().read(io, prefix);
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

	void computePvectors(ProgramGlobals::DirectionEnum dir, RealType Eg)
	{
		if (allOrigPvectorsDone()) {
			const SizeType tvs = this->common().aoe().tvs();
			if (tvs == pvectors_.origPvectors()) return;
			if (tvs < pvectors_.origPvectors())
				err("TVS could not have decreased ?!\n");

			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"All user-provided P vectors finished";
			progress_.printline(msgg, std::cout);
			std::cerr<<PsimagLite::AnsiColor::green;
			std::cerr<<"All user-provided P vectors finished\n";
			std::cerr<<PsimagLite::AnsiColor::reset;

			this->common().aoeNonConst().targetVectorsResize(pvectors_.origPvectors());

			return;
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"P0="<<pvectors_(0).lastName();
		progress_.printline(msgg, std::cout);

		CanonicalExpressionType canonicalExpression(opSpec_);
		SizeType total = pvectors_.targets();

		AuxForTargetingExpressionType aux(pvectors_, timeEvolve_, dir, Eg);
		AlgebraType opEmpty(aux);
		bool needsTrimming = false;
		PsimagLite::String allpvectors;
		for (SizeType i = 0; i < total; ++i) {

			if (pvectors_(i).isDone()) continue;

			aux.setPindexOutput(i);
			AlgebraType tmp(aux);

			canonicalExpression(tmp, pvectors_(i).lastName(), opEmpty, aux);
			tmp.finalize();

			PsimagLite::String thispBefore = tmp.toString();
			finalize(aux.tempVectors(), aux.tempNames(), i, thispBefore);
			PsimagLite::String thispAfter = pvectors_(i).lastName();

			int x = tmp.pIndex();
			if (x >= 0) {
				if (static_cast<SizeType>(x) == i) err("Self assigment\n");
				VectorWithOffsetType_& v0 = this->tvNonConst(i);
				v0 =  this->tv(x);
				pvectors_.setAsDone(i);
				continue;
			}

			if (!tmp.hasSummationKetAndNoMult()) {
				allpvectors += thispAfter;
				continue;
			}

			PsimagLite::String newpstring = simplifyTerms(thispBefore);
			if (newpstring != pvectors_(i).lastName()) {
				needsTrimming = true;
				PsimagLite::String compr = compressExpression(newpstring);
				//checkNoUncompressedExists(compr);
				pvectors_.pushString(i, compr);
				newpstring += compr;
			}

			allpvectors += newpstring;
		}

		if (needsTrimming) pvectors_.trimPvectors(allpvectors);
	}

	void finalize(const VectorVectorWithOffsetType& tempVectors,
	              const VectorStringType& tempNames,
	              SizeType pVectorIndex,
	              PsimagLite::String tempExpr)
	{
		const SizeType ntemps = tempNames.size();

		if (ntemps == 0) return;

		// find tempNames_ in pVectors_ and trim tempNames_ accordingly
		VectorBoolType removed_(ntemps);
		VectorSizeType tempToP(ntemps, 10000);
		for (SizeType i = 0; i < ntemps; ++i) {
			int x = pvectors_.findInOrigNames(tempNames[i]);
			if (x < 0) continue;
			this->tvNonConst(x) = tempVectors[i];
			pvectors_.setAsDone(x);
			removed_[i] = true;
			tempToP[i] = x;
		}

		// these surviving tempNames_ need storage, add them
		for (SizeType i = 0; i < ntemps; ++i) {
			if (removed_[i]) continue;
			int x = pvectors_.findInAnyNames(tempNames[i]);
			if (x >= 0) continue;
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
			if (!pvectors_(i).isDone()) return false;
		return true;
	}

	// replace "|!m" + something ==> "|P" + number
	PsimagLite::String compressExpression(PsimagLite::String str) const
	{
		SizeType i = 0;
		const SizeType len = str.length();
		if (len < 4) return str;
		PsimagLite::String result;
		for (; i < len; ++i) {
			if (i + 4 < len && str.substr(i,3) == "|!m") {
				SizeType j = i + 3;
				PsimagLite::String buffer;
				for (;j < len; ++j) {
					buffer += str[j];
					if (str[j] == '>') break;
				}

				const SizeType ind = pvectors_.findPforThisExpression(buffer);
				buffer = "|P" + ttos(ind) + ">";
				i = j;
				result += buffer;
				continue;
			}

			result += str[i];
		}

		return result;
	}

	// replace "R" + i ==> "P" + tempToP[i]
	PsimagLite::String expandExpression(PsimagLite::String str,
	                                    const VectorSizeType& tempToP) const
	{
		SizeType i = 0;
		const SizeType len = str.length();
		if (len < 4) return str;
		PsimagLite::String expanded;
		for (; i < len; ++i) {
			if (i + 4 < len && str[i] == '|' && str[i + 1] == 'R') {
				SizeType j = i + 2;
				PsimagLite::String buffer;
				for (;j < len; ++j) {
					if (str[j] == '>') break;
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
	PsimagLite::String simplifyTerms(PsimagLite::String str)
	{
		SizeType i = 0;
		const SizeType len = str.length();
		if (len < 4) return str;
		PsimagLite::String simplified;
		for (; i < len; ++i) {
			if (i + 4 < len && str.substr(i, 4) == "|!aP") {
				SizeType j = i + 4;
				PsimagLite::String buffer;
				for (;j < len; ++j) {
					if (str[j] == 'p') break;
					buffer += str[j];
				}

				SizeType ind0 = PsimagLite::atoi(buffer);
				buffer = "";
				++j;
				assert(str[j] == 'P');
				++j;
				for (;j < len; ++j) {
					if (str[j] == '>') break;
					buffer += str[j];
				}

				SizeType ind1 = PsimagLite::atoi(buffer);

				// before reordering ind0 and ind1
				PsimagLite::String p0PlusP1 = "|P" + ttos(ind0) + ">+|P" + ttos(ind1) + ">";

				// ask aoe to sum ind0 and ind1 and put it into ind0
				assert(ind0 != ind1);
				if (ind0 > ind1) {
					const SizeType tmp = ind0;
					ind0 = ind1;
					ind1 = tmp;
				}

				pvectors_.sumPvectors(ind0, ind1, p0PlusP1);

				buffer = "|P" + ttos(ind0) + ">";
				i = j;
				simplified += buffer;
				continue;
			}

			simplified += str[i];
		}

		return simplified;
	}

	static void checkNoUncompressedExists(PsimagLite::String str)
	{
		SizeType i = 0;
		const SizeType len = str.length();
		if (len < 4) return;
		for (; i < len; ++i) {
			if (i + 4 < len && str.substr(i, 3) == "|!m")
				err("Uncompressed exists in " + str + "\n");
		}
	}

	void computeAllWeights(RealType& gsWeight, VectorRealType& weight) const
	{
		const SizeType n = this->common().aoe().tvs();
		assert(n <= pvectors_.targets());
		weight.resize(n);
		std::fill(weight.begin(), weight.end(), 0);
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			RealType norma = norm(this->tv(i));
			if (norma < 1e-6) continue;
			weight[i] = pvectors_(i).weight()/norma;
			sum += weight[i];
		}

		gsWeight = 1 - sum;

		if (gsWeight >= gsWeight_) return; // <--- EARLY EXIT HERE

		assert(sum > 1e-6);
		RealType factor = (1 - gsWeight_)/sum;
		for (SizeType i = 0; i < n; ++i)
			weight[i] *= factor;

		gsWeight = gsWeight_;
	}

	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	SpecForTargetingExpressionType opSpec_;
	Pvectors<BaseType> pvectors_;
	GroupOfOneTimeEvolutionsType timeEvolve_;
};     //class TargetingExpression
} // namespace Dmrg
/*@}*/
#endif

