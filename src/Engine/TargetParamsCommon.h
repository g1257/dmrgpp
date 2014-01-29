/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file TargetParamsCommon.h
 *
 *  FIXME
 */
#ifndef TARGET_PARAMS_COMMON_H
#define TARGET_PARAMS_COMMON_H
#include <vector>
#include <stdexcept>
#include <iostream>
#include "CookedOperator.h"
#include "TargetParamsBase.h"

namespace Dmrg {
//! Coordinates reading of TargetSTructure from input file
template<typename ModelType>
class TargetParamsCommon : public TargetParamsBase<ModelType> {

public:

	typedef typename ModelType::RealType RealType;

	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;

	enum {PRODUCT,SUM};

	template<typename IoInputter>
	TargetParamsCommon(IoInputter& io,const ModelType& model)
	    : sites_(0),
	      startingLoops_(0),
	      concatenation_(PRODUCT),
	      noOperator_(false),
	      model_(model)
	{
		io.read(sites_,"TSPSites");
		io.read(startingLoops_,"TSPLoops");

		PsimagLite::String productOrSum_ = "product";
		try {
			io.readline(productOrSum_,"TSPProductOrSum=");
		} catch (std::exception& e) {
			PsimagLite::String s(__FILE__);
			s += "\n FATAL: Must provide TSPProductOrSum=.\n";
			s += "Please add TSPProductOrSum=product or TSPProductOrSum=sum  ";
			s += "immediately below the TSPLoops= line in the input file\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		//! Concatenation specifies what to do with
		//! operators at different sites, add them or multiply them
		if (productOrSum_ == "product") {
			this->concatenation_ = PRODUCT;
		} else if (productOrSum_ == "sum") {
			this->concatenation_ = SUM;
		} else {
			PsimagLite::String s(__FILE__);
			s += " : Unknown concatentation " + productOrSum_ + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		aOperators_.resize(sites_.size());

		CookedOperator<ModelType> cookedOperator(model_);

		for (SizeType i=0;i<sites_.size();i++) {
			OperatorType myOp(io,cookedOperator,OperatorType::MUST_BE_NONZERO);
			aOperators_[i] = myOp;
		}

		noOperator_ = isNoOperator();
		checkBorderOperators();
		checkSizesOfOperators();
	}

	virtual SizeType sites() const
	{
		return sites_.size();
	}

	virtual SizeType sites(SizeType i) const
	{
		assert(i < sites_.size());
		return sites_[i];
	}

	virtual const VectorSizeType& startingLoops() const
	{
		return startingLoops_;
	}

	virtual SizeType concatenation() const
	{
		return concatenation_;
	}

	virtual const VectorOperatorType& aOperators() const
	{
		return aOperators_;
	}

	virtual void setConcatenation(SizeType x)
	{
		concatenation_ = x;
	}

	virtual bool noOperator() const
	{
		return noOperator_;
	}

	virtual void noOperator(bool x)
	{
		noOperator_ = x;
	}

	template<typename ModelType_>
	friend std::ostream& operator<<(std::ostream& os,
	                                const TargetParamsCommon<ModelType_>& t);

private:

	bool isNoOperator() const
	{
		if (aOperators_.size()!=1) return false;
		return (isTheIdentity(aOperators_[0].data) && aOperators_[0].fermionSign);
	}

	void checkSizesOfOperators() const
	{
		if (sites_.size() != aOperators_.size() ||
		    sites_.size() != startingLoops_.size())
			throw PsimagLite::RuntimeError("CommonTargetting\n");

		for (SizeType i=0;i<aOperators_.size();i++) {
			SizeType n = aOperators_[i].data.row();
			if (n!=model_.hilbertSize(sites_[i]))
				throw PsimagLite::RuntimeError("CommonTargetting\n");
		}
	}

	void checkBorderOperators()
	{
		if (sites_.size() == 0) return;

		SizeType linSize = model_.geometry().numberOfSites();

		if (hasOperatorAt(0) && !hasOperatorAt(1)) {
			errorBorderOperators(0);
		}

		if (hasOperatorAt(linSize-1) && !hasOperatorAt(linSize-2)) {
			errorBorderOperators(linSize-1);
		}
	}

	bool hasOperatorAt(SizeType site) const
	{
		for (SizeType i = 0; i < sites_.size(); ++i) {
			if (sites_[i] == site) return true;
		}
		return false;
	}

	void errorBorderOperators(SizeType site)
	{
		SizeType linSize = model_.geometry().numberOfSites();
		SizeType site2 = (site == 0) ? 1 : linSize - 2;

		PsimagLite::String str("ERROR: Operators at border site: Please ");
		str += "add the identity operator at site " + ttos(site2) + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	VectorSizeType sites_;
	VectorSizeType startingLoops_;
	SizeType concatenation_;
	VectorOperatorType aOperators_;
	bool noOperator_;
	const ModelType& model_;
}; // class TargetParamsCommon

template<typename ModelType>
inline std::ostream&
operator<<(std::ostream& os,const TargetParamsCommon<ModelType>& t)
{
	os<<"#TargetParams.operators="<<t.aOperators_.size()<<"\n";
	for (SizeType i=0;i<t.aOperators_.size();i++) {
		os<<"#TargetParams.operator "<<i<<"\n";
		os<<t.aOperators_[i];
	}

	os<<"#TargetParams.site="<<t.sites_;
	os<<"#TargetParams.startingLoop="<<t.startingLoops_<<"\n";

	return os;
}
} // namespace Dmrg 

/*@}*/
#endif // TARGET_PARAMS_COMMON_H

