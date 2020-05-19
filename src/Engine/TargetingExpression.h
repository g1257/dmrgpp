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
	typedef typename PsimagLite::Vector<PvectorType*>::Type VectorPvectorType;
	typedef SpecForTargetingExpression<BaseType> SpecForTargetingExpressionType;
	typedef typename SpecForTargetingExpressionType::AlgebraType AlgebraType;
	typedef PsimagLite::CanonicalExpression<SpecForTargetingExpressionType>
	CanonicalExpressionType;
	typedef AuxForTargetingExpression<BaseType> AuxForTargetingExpressionType;
	typedef typename TargetingCommonType::VectorRealType VectorRealType;
	typedef typename AuxForTargetingExpressionType::VectorStringType VectorStringType;
	typedef typename AuxForTargetingExpressionType::VectorVectorWithOffsetType
	VectorVectorWithOffsetType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename AlgebraType::VectorSizeType VectorSizeType;

public:

	TargetingExpression(const LeftRightSuperType& lrs,
	                    const ModelType& model,
	                    const WaveFunctionTransfType& wft,
	                    const QnType&,
	                    InputValidatorType& io)
	    : BaseType(lrs,model,wft,0),
	      progress_("TargetingExpression"),
	      gsWeight_(0.3),
	      origPvectors_(0)
	{
		io.readline(gsWeight_, "GsWeight=");
		pvectorsFromInput(io);
	}

	~TargetingExpression()
	{
		for (SizeType i = 0; i < pVectors_.size(); ++i) {
			delete pVectors_[i];
			pVectors_[i] = 0;
		}
	}

	SizeType sites() const { return 0; }

	SizeType targets() const { return pVectors_.size(); }

	RealType normSquared(SizeType i) const
	{
		return PsimagLite::real(this->common().aoe().targetVectors()[i]*
		                        this->common().aoe().targetVectors()[i]);
	}

	RealType weight(SizeType i) const
	{
		assert(this->common().aoe().noStageIs(StageEnumType::DISABLED));
		assert(i < pVectors_.size());
		return pVectors_[i]->weight();
	}

	RealType gsWeight() const
	{
		return (this->common().aoe().noStageIs(StageEnumType::DISABLED)) ? gsWeight_ : 1.0;
	}

	void evolve(const VectorRealType&,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType&,
	            SizeType)
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE) return;

		this->common().setAllStagesTo(StageEnumType::WFT_NOADVANCE);
		assert(block1.size() == 1);
		const SizeType site = block1[0];
		const SizeType total = this->common().aoe().targetVectors().size();
		assert(total <= pVectors_.size());
		this->common().aoe().wftSome(site, 0, total);

		computePvectors(direction); // may alter the number of tvs

		const SizeType n = this->common().aoe().targetVectors().size();
		assert(n <= pVectors_.size());
		VectorRealType weight(n);
		for (SizeType i = 0; i < n; ++i)
			weight[i] = pVectors_[i]->weight();

		this->common().printNormsAndWeights(gsWeight_, weight);

		const bool doBorderIfBorder = true;
		this->common().cocoon(block1, direction, doBorderIfBorder); // in-situ
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

private:

	void pvectorsFromInput(InputValidatorType& io)
	{
		SizeType total = 0;
		io.readline(total, "Pvectors=");
		pVectors_.resize(total);
		PsimagLite::String tmp;
		RealType sum = 0.0;
		for (SizeType i = 0; i < total; ++i) {
			io.readline(tmp, "P" + ttos(i) + "=");
			pVectors_[i] = new PvectorType(tmp);
			sum += pVectors_[i]->weight();
		}

		if (sum == 0.0) return;

		RealType factor = (1.0 - gsWeight_)/sum;
		for (SizeType i = 0; i < total; ++i)
			pVectors_[i]->multiplyWeight(factor);

		origPvectors_ = pVectors_.size();
	}

	void computePvectors(ProgramGlobals::DirectionEnum dir)
	{
		if (allOrigPvectorsDone()) {
			const SizeType tvs = this->common().aoe().targetVectors().size();
			if (tvs == origPvectors_) return;
			if (tvs < origPvectors_)
				err("TVS could not have decreased ?!\n");

			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"All user-provided P vectors finished";
			progress_.printline(msgg, std::cout);
			std::cerr<<PsimagLite::AnsiColor::green;
			std::cerr<<"All user-provided P vectors finished\n";
			std::cerr<<PsimagLite::AnsiColor::reset;

			this->common().aoe().targetVectorsResize(origPvectors_);

			return;
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"P0="<<pVectors_[0]->lastName();
		progress_.printline(msgg, std::cout);

		CanonicalExpressionType canonicalExpression(opSpec_);
		SizeType total = pVectors_.size();

		AuxForTargetingExpressionType aux(this->common().aoe(),
		                                  this->model(),
		                                  this->lrs(),
		                                  dir);
		const AlgebraType opEmpty(aux);
		bool needsTrimming = false;
		PsimagLite::String allpvectors;
		for (SizeType i = 0; i < total; ++i) {

			if (pVectors_[i]->lastName() == "DONE") continue;

			AlgebraType tmp(aux);

			canonicalExpression(tmp, pVectors_[i]->lastName(), opEmpty, aux);
			tmp.finalize();

			PsimagLite::String thispBefore = tmp.toString();
			finalize(aux.tempVectors(), aux.tempNames(), i, thispBefore);
			PsimagLite::String thispAfter = pVectors_[i]->lastName();

			int x = tmp.pIndex();
			if (x >= 0) {
				if (static_cast<SizeType>(x) == i) err("Self assigment\n");
				VectorWithOffsetType_& v0 = this->common().aoe().targetVectors(i);
				v0 =  this->common().aoe().targetVectors(x);
				pVectors_[i]->pushString("DONE");
				continue;
			}

			if (!tmp.hasSummationKetAndNoMult()) {
				allpvectors += thispAfter;
				continue;
			}

			PsimagLite::String newpstring = simplifyTerms(thispBefore);
			if (newpstring != pVectors_[i]->lastName()) {
				needsTrimming = true;
				PsimagLite::String compr = compressExpression(newpstring);
				//checkNoUncompressedExists(compr);
				pVectors_[i]->pushString(compr);
				newpstring += compr;
			}

			allpvectors += newpstring;
		}

		if (needsTrimming) trimPvectors(allpvectors);
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
			int x = findInOrigNames(tempNames[i]);
			if (x < 0) continue;
			this->common().aoe().targetVectors(x) = tempVectors[i];
			pVectors_[x]->pushString("DONE");
			removed_[i] = true;
			tempToP[i] = x;
		}

		// these surviving tempNames_ need storage, add them
		for (SizeType i = 0; i < ntemps; ++i) {
			if (removed_[i]) continue;
			int x = findInAnyNames(tempNames[i]);
			if (x >= 0) continue;
			const SizeType ind = this->common().aoe().createPvector(tempVectors[i]);
			tempToP[i] = ind;
			const PsimagLite::String ename = expandExpression(tempNames[i], tempToP);
			PvectorType* pnew = new PvectorType(ename);
			pnew->pushString("DONE");
			pVectors_.push_back(pnew);

			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"P["<<ind<<"]="<<ename<<" created";
			progress_.printline(msgg, std::cout);
		}

		PsimagLite::String newpstring = compressExpression(tempExpr);
		const PsimagLite::String selfName = "|P" + ttos(pVectorIndex) + ">";
		if (newpstring == selfName) newpstring = "DONE";
		pVectors_[pVectorIndex]->pushString(newpstring);
	}

	bool allOrigPvectorsDone() const
	{
		for (SizeType i = 0; i < origPvectors_; ++i)
			if (pVectors_[i]->lastName() != "DONE") return false;
		return true;
	}

	int findInOrigNames(PsimagLite::String str) const
	{
		return findInNames(str, origPvectors_);
	}

	int findInAnyNames(PsimagLite::String str) const
	{
		return findInNames(str, pVectors_.size());
	}

	int findInNames(PsimagLite::String str, SizeType end) const
	{
		assert(end <= pVectors_.size());
		for (SizeType i = 0; i < end; ++i)
			if (pVectors_[i]->hasAnyName(str)) return i;

		return -1;
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

				const SizeType ind = findPforThisExpression(buffer);
				buffer = "|P" + ttos(ind) + ">";
				i = j;
				result += buffer;
				continue;
			}

			result += str[i];
		}

		return result;
	}

	SizeType findPforThisExpression(PsimagLite::String str) const
	{
		const SizeType pvectors = pVectors_.size();
		for (SizeType i = 0; i < pvectors; ++i) {
			if (!pVectors_[i]->hasAnyName(str)) continue;
			return i;
		}

		throw PsimagLite::RuntimeError("findPforThisExpression: P not found\n");
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

				sumPvectors(ind0, ind1, p0PlusP1);

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

	void sumPvectors(SizeType ind0, SizeType ind1, PsimagLite::String p0PlusP1)
	{
		assert(ind0 < ind1);
		VectorWithOffsetType_& v0 = this->common().aoe().targetVectors(ind0);
		const VectorWithOffsetType_& v1 =  this->common().aoe().targetVectors(ind1);
		v0 += v1;
		pVectors_[ind0]->sum(*(pVectors_[ind1]), p0PlusP1);
	}

	void trimPvectors(PsimagLite::String str)
	{
		const SizeType tvs = this->common().aoe().targetVectors().size();
		VectorBoolType used(tvs, false);
		assert(origPvectors_ <= tvs);
		for (SizeType i = 0; i < origPvectors_; ++i)
			used[i] = true;

		findUsedPvectors(used, str);

		for (SizeType i = 0; i < tvs; ++i) {
			if (used[i]) continue;
			this->common().aoe().destroyPvector(i);
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"P["<<i<<"] destroyed";
			progress_.printline(msgg, std::cout);
		}

		// this line is commented out because then the index of aoe.targetVectors
		// will be different than the index of pVectors
		// aoe.targetVectors is resized after all vectors have been computed
		//this->common().aoe().trimVectors();
		const SizeType tvsFinal = this->common().aoe().targetVectors().size();
		if (tvs != tvsFinal) {
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Number of target vectors is "<<tvsFinal<<" now";
			progress_.printline(msgg, std::cout);
		}
	}

	void findUsedPvectors(VectorBoolType& used, PsimagLite::String str) const
	{
		const SizeType len = str.size();
		for (SizeType i = 0; i < len; ++i) {

			if (str.substr(i, 2) != "|P") continue;

			PsimagLite::String buffer;
			SizeType j = i + 2;
			for (; j < len; ++j) {
				if (str[j] >= 48 && str[j] <= 57) {
					buffer += str[j];
				} else if (str[j] == '>') {
					break;
				} else {
					buffer = "";
					break;
				}
			}

			i = j;
			if (buffer == "") continue;
			const SizeType ind = PsimagLite::atoi(buffer);
			used[ind] = true;
		}
	}

	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	SizeType origPvectors_;
	VectorPvectorType pVectors_;
	SpecForTargetingExpressionType opSpec_;
};     //class TargetingExpression
} // namespace Dmrg
/*@}*/
#endif

