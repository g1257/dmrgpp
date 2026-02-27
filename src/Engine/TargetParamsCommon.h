/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file TargetParamsCommon.h
 *
 *  FIXME
 */
#ifndef TARGET_PARAMS_COMMON_H
#define TARGET_PARAMS_COMMON_H
#include "AlgebraicStringToNumber.h"
#include "Io/IoSerializerStub.h"
#include "ProgramGlobals.h"
#include "TargetParamsBase.h"
#include <iostream>
#include <stdexcept>
#include <vector>

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType> class TargetParamsCommon : public TargetParamsBase<ModelType> {

public:

	using RealType           = typename ModelType::RealType;
	using BaseType           = TargetParamsBase<ModelType>;
	using OperatorType       = typename ModelType::OperatorType;
	using PairType           = typename OperatorType::PairType;
	using SparseMatrixType   = typename OperatorType::StorageType;
	using ComplexOrRealType  = typename SparseMatrixType::value_type;
	using MatrixType         = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorType         = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using VectorSizeType     = PsimagLite::Vector<SizeType>::Type;
	using VectorMatrixType   = typename PsimagLite::Vector<MatrixType>::Type;
	using VectorOperatorType = typename PsimagLite::Vector<OperatorType>::Type;
	using InputValidatorType = typename ModelType::InputValidatorType;
	using PairSizeType       = std::pair<SizeType, SizeType>;
	using VectorStringType   = PsimagLite::Vector<PsimagLite::String>::Type;

	TargetParamsCommon(InputValidatorType& io,
	                   PsimagLite::String  targeting,
	                   const ModelType&    model)
	    : BaseType(targeting)
	    , sites_(0)
	    , startingLoops_(0)
	    , concatenation_(BaseType::ConcatEnum::PRODUCT)
	    , noOperator_(false)
	    , skipTimeZero_(false)
	    , isEnergyForExp_(false)
	    , gsWeight_(0.0)
	    , energyForExp_(0.0)
	    , io_(io)
	    , model_(model)
	    , sectorLevel_(PairSizeType(0, 0))
	{
		/*PSIDOC TargetParamsCommon
		\item[TSPSites] [VectorInteger] The first number is the number of numbers
		to follow. The following numbers are the
		sites $\pi'(0),\pi'(1),\cdots$ where the operators
		to build state $|\phi\rangle$ should be applied; see Eq.~(4) of
		\cite{re:alvarez11}. These sites must be
		ordered by appearance in the DMRG sweeping.
		\item[TSPLoops] [VectorInteger]
		The first number is the number of numbers
		to follow. The following numbers are
		the delay (in units of finite loops) before evolving in time. Delaying the
		time evolution helps converge the state $|\phi\rangle$.
		\item[TSPProductOrSum] [String] Either \verb!product! or \verb!sum! indicating
		whether the operators $B$ in Eq.~(4) of \cite{re:alvarez11} should be multiplied
		or summed. Note that Eq.~(4) shows only multiplication.
		\item[TSPSkipTimeZero] [Integer] Either 0 or 1 to indicate whether to skip the
		application of the time evolution at $t=0$. It is ignored unless
		$|\phi\rangle$ is the ground state.
		\item[TSPEnergyForExp] [RealType] Energy to use as origin for the exponential
		in the time evolution.
		*/

		if (targeting == "TargetingExpression")
			return;

		io.readline(gsWeight_, "GsWeight=");

		readOperatorsToApplyIfAny();

		try {
			int tmp = 0;
			io.readline(tmp, "TSPSkipTimeZero=");
			skipTimeZero_ = (tmp > 0);
		} catch (std::exception&) { }

		try {
			io.readline(energyForExp_, "TSPEnergyForExp=");
			isEnergyForExp_ = true;
		} catch (std::exception&) { }

		try {
			VectorType tmpVector;
			io.read(tmpVector, "TSPOperatorMultiplier");
			multiplyOperators(tmpVector);
			std::cout << "TSPOperatorMultiplier found with " << tmpVector.size()
			          << " entries.\n";
		} catch (std::exception&) { }

		bool               hasApplyTo = false;
		PsimagLite::String tmp;
		try {
			io.readline(tmp, "TSPApplyTo=");
			hasApplyTo = true;
		} catch (std::exception&) { }

		if (hasApplyTo) {
			PsimagLite::GetBraOrKet gbok(tmp);
			sectorLevel_ = PairSizeType(gbok.sectorIndex(), gbok.levelIndex());
		}

		noOperator_ = isNoOperator();
		checkBorderOperators();
		checkSizesOfOperators();
	}

	SizeType memResolv(PsimagLite::MemResolv&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	virtual SizeType sites() const { return sites_.size(); }

	virtual SizeType sites(SizeType i) const
	{
		assert(i < sites_.size());
		return sites_[i];
	}

	virtual int findIndexOfSite(SizeType site, SizeType start) const
	{
		if (start >= sites_.size())
			return -1;
		VectorSizeType::const_iterator it
		    = std::find(sites_.begin() + start, sites_.end(), site);
		if (it == sites_.end())
			return -1;
		return it - sites_.begin();
	}

	virtual void setOperator(SizeType i, SizeType j, const OperatorType& op)
	{
		assert(i < sites_.size());
		sites_[i] = j;
		assert(i < aOperators_.size());
		aOperators_[i] = op;
	}

	virtual const VectorSizeType& startingLoops() const { return startingLoops_; }

	virtual typename BaseType::ConcatEnum concatenation() const { return concatenation_; }

	virtual const VectorOperatorType& aOperators() const { return aOperators_; }

	virtual bool noOperator() const { return noOperator_; }

	virtual void noOperator(bool x) { noOperator_ = x; }

	virtual bool skipTimeZero() const { return skipTimeZero_; }

	virtual bool isEnergyForExp() const { return isEnergyForExp_; }

	virtual RealType energyForExp() const { return energyForExp_; }

	virtual RealType gsWeight() const { return gsWeight_; }

	virtual SizeType sectorIndex() const { return sectorLevel_.first; }

	virtual SizeType levelIndex() const { return sectorLevel_.second; }

	void write(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);

		ioSerializer.write(label + "/sites_", sites_);
		ioSerializer.write(label + "/startingLoops_", startingLoops_);
		ioSerializer.write(label + "/concatenation_", concatenation_);
		ioSerializer.write(label + "/noOperator_", noOperator_);
		ioSerializer.write(label + "/skipTimeZero_", skipTimeZero_);
		ioSerializer.write(label + "/isEnergyForExp_", isEnergyForExp_);
		ioSerializer.write(label + "/gsWeight_", gsWeight_);
		ioSerializer.write(label + "/energyForExp_", energyForExp_);
		ioSerializer.write(label + "/aOperators_", aOperators_);
	}

	friend std::ostream& operator<<(std::ostream& os, const TargetParamsCommon& t)
	{
		os << "TargetParams.operators=" << t.aOperators_.size() << "\n";
		for (SizeType i = 0; i < t.aOperators_.size(); i++) {
			os << "TargetParams.operator " << i << "\n";
			os << t.aOperators_[i];
		}

		os << "TargetParams.site=" << t.sites_;
		os << "TargetParams.startingLoop=" << t.startingLoops_ << "\n";

		return os;
	}

private:

	void multiplyOperators(const VectorType& tmpVector)
	{
		if (aOperators_.size() != tmpVector.size()) {
			PsimagLite::String str(__FILE__);
			str += "\nFATAL: ";
			str += " TSPOperatorMultiplier of size " + ttos(tmpVector.size());
			str += " but " + ttos(aOperators_.size()) + " expected.\n";
		}

		for (SizeType i = 0; i < aOperators_.size(); i++)
			aOperators_[i] *= tmpVector[i];
	}

	bool isNoOperator() const
	{
		if (aOperators_.size() != 1)
			return false;
		return (isTheIdentity(aOperators_[0].getStorage())
		        && aOperators_[0].fermionOrBoson()
		            == ProgramGlobals::FermionOrBosonEnum::BOSON);
	}

	void checkSizesOfOperators() const
	{
		if (sites_.size() != aOperators_.size()) {
			PsimagLite::String str(__FILE__);
			str += " Listed sites is " + ttos(sites_.size());
			str += " but operators found is " + ttos(aOperators_.size()) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		for (SizeType i = 0; i < aOperators_.size(); i++) {
			SizeType n       = aOperators_[i].getStorage().rows();
			SizeType hilbert = model_.hilbertSize(sites_[i]);
			if (n != hilbert) {
				PsimagLite::String str(__FILE__);
				str += " Operator rank is " + ttos(n) + " but ";
				str += ttos(hilbert) + " was expected\n";
				throw PsimagLite::RuntimeError(str);
			}
		}
	}

	void checkBorderOperators() const
	{
		if (sites_.size() == 0)
			return;

		SizeType linSize = model_.superGeometry().numberOfSites();

		if (hasOperatorAt(0) && !hasOperatorAt(1)) {
			errorBorderOperators(0);
		}

		if (hasOperatorAt(linSize - 1) && !hasOperatorAt(linSize - 2)) {
			errorBorderOperators(linSize - 1);
		}
	}

	void checkSites() const
	{
		SizeType linSize = model_.superGeometry().numberOfSites();
		for (SizeType i = 0; i < sites_.size(); ++i) {
			if (sites_[i] >= linSize) {
				PsimagLite::String str(__FILE__);
				str += " TSPSites: The " + ttos(i) + "-th site is ";
				str += ttos(sites_[i]) + " is larger than the total ";
				str += "number of sites " + ttos(linSize) + "\n";
				throw PsimagLite::RuntimeError(str);
			}
		}
	}

	bool hasOperatorAt(SizeType site) const
	{
		for (SizeType i = 0; i < sites_.size(); ++i) {
			if (sites_[i] == site)
				return true;
		}
		return false;
	}

	void errorBorderOperators(SizeType site) const
	{
		SizeType linSize = model_.superGeometry().numberOfSites();
		SizeType site2   = (site == 0) ? 1 : linSize - 2;

		PsimagLite::String str("ERROR: Operators at border site: Please ");
		str += "add the identity operator at site " + ttos(site2) + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	void vecstringToVecnumbers(VectorSizeType& nums, const VectorStringType& strs)
	{
		using AlgebraicStringToNumberType = AlgebraicStringToNumber<RealType>;
		const SizeType numberOfSites      = model_.superGeometry().numberOfSites();
		const SizeType n                  = strs.size();
		if (n == 0) {
			return;
		}

		nums.resize(n);

		AlgebraicStringToNumberType algebraicStringToNumber("TSPSites", numberOfSites);
		std::cout << "TSPSites=[";
		for (SizeType i = 0; i < n; ++i) {
			nums[i] = algebraicStringToNumber.procLength(strs[i]);
			std::cout << nums[i];
			if (i + 1 < n)
				std::cout << ", ";
		}

		std::cout << "];\n";
	}

	void readOperatorsToApplyIfAny()
	{
		VectorStringType sitesStr;
		try {
			io_.read(sitesStr, "TSPSites");
		} catch (std::exception&) { }

		vecstringToVecnumbers(sites_, sitesStr);
		checkSites();

		bool has_label = false;
		try {
			io_.read(startingLoops_, "TSPLoops");
			has_label = true;
		} catch (std::exception&) { }

		checkError(sitesStr.size() > 0, has_label, "TSPLoops");

		if (sites_.size() != startingLoops_.size()) {
			PsimagLite::String str(__FILE__);
			str += " Listed sites is " + ttos(sites_.size());
			str += " but delay loops found is " + ttos(startingLoops_.size()) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String productOrSum = "product";
		has_label                       = false;
		try {
			io_.readline(productOrSum, "TSPProductOrSum=");
			has_label = true;
		} catch (std::exception&) { }

		checkError(sitesStr.size() > 1, has_label, "TSPProductOrSum");

		//! Concatenation specifies what to do with
		//! operators at different sites, add them or multiply them
		if (productOrSum == "product") {
			this->concatenation_ = BaseType::ConcatEnum::PRODUCT;
		} else if (productOrSum == "sum") {
			this->concatenation_ = BaseType::ConcatEnum::SUM;
		} else {
			PsimagLite::String s(__FILE__);
			s += " : Unknown concatentation " + productOrSum + "\n";
			s += "Possible values: product sum\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		if (sites_.size() == 0) {
			// Here we add a bogus operator: the identity in the middle of the lattice
			const SizeType numberOfSites = model_.superGeometry().numberOfSites();
			// it's OK if numberOfSites is odd also
			sites_.push_back(numberOfSites / 2);
			startingLoops_.push_back(0);
			OperatorType identity = model_.naturalOperator("identity", 0, 0);
			aOperators_.push_back(identity);
		} else {
			readOperators();
		}
	}

	void readOperators()
	{
		for (SizeType i = 0; i < sites_.size(); ++i) {
			PsimagLite::String prefix2 = (io_.isAinur()) ? "TSPOp" + ttos(i) + ":" : "";
			OperatorType       myOp(
                            io_, model_, OperatorType::MUST_BE_NONZERO, prefix2, sites_[i]);
			aOperators_.push_back(myOp);
		}
	}

	static void checkError(bool label_is_needed, bool has_label, const std::string& label)
	{
		if (label_is_needed) {
			if (!has_label) {
				err("Label " + label + " is needed in input file\n");
			}
		} else {
			if (has_label) {
				std::string msg = "Label " + label + " will be ignored\n";
				std::cerr << msg;
				std::cout << msg;
			}
		}
	}

	VectorSizeType                sites_;
	VectorSizeType                startingLoops_;
	typename BaseType::ConcatEnum concatenation_;
	bool                          noOperator_;
	bool                          skipTimeZero_;
	bool                          isEnergyForExp_;
	RealType                      gsWeight_;
	RealType                      energyForExp_;
	VectorOperatorType            aOperators_;
	InputValidatorType&           io_;
	const ModelType&              model_;
	PairSizeType                  sectorLevel_;
}; // class TargetParamsCommon
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_COMMON_H
