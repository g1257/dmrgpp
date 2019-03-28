/*
Copyright (c) 2014, UT-Battelle, LLC
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

/*! \file TargetParamsBase.h
 *
 *  FIXME
 */
#ifndef TARGET_PARAMS_BASE_H
#define TARGET_PARAMS_BASE_H
#include <vector>
#include <stdexcept>
#include "MemResolv.h"
#include "FreqEnum.h"

namespace Dmrg {

template<typename ModelType>
class TargetParamsBase {
public:

	typedef typename ModelType::RealType RealType;
	typedef std::pair<PsimagLite::FreqEnum, RealType> PairFreqType;
	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	enum class ConcatEnum {PRODUCT, SUM};

	enum class AlgorithmEnum {KRYLOV, CONJUGATE_GRADIENT, CHEBYSHEV, RUNGE_KUTTA, SUZUKI_TROTTER};

	virtual ~TargetParamsBase() {}

	virtual SizeType sites() const = 0;

	virtual SizeType sites(SizeType) const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented sites\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual void setOperator(SizeType, SizeType, const OperatorType&)
	{
		unimplemented("type(setOperator)");
	}

	virtual const VectorSizeType& startingLoops() const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented startingLoops\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual const VectorRealType& chebyTransform() const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented chebyTransform\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual ConcatEnum concatenation() const
	{
		throw PsimagLite::RuntimeError("concatenation");
	}

	virtual const VectorOperatorType& aOperators() const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented aOperators\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual RealType correctionA() const
	{
		return 0;
	}

	virtual SizeType type() const
	{
		return unimplementedInt("type");
	}

	virtual void type(SizeType)
	{
		unimplemented("type(SizeType)");
	}

	virtual SizeType advanceEach() const
	{
		return 0;
	}

	virtual SizeType cgSteps() const
	{
		return unimplementedInt("cgSteps");
	}

	virtual PairFreqType omega() const
	{
		PsimagLite::String s("TargetParamsBase: unimplemented omega \n");
		throw PsimagLite::RuntimeError(s);
	}

	virtual void omega(PsimagLite::FreqEnum, RealType)
	{
		unimplemented("omega(RealType)");
	}

	virtual RealType eta() const
	{
		return unimplemented("eta");
	}

	virtual RealType cgEps() const
	{
		return unimplemented("cgEps");
	}

	virtual AlgorithmEnum algorithm() const
	{
		throw PsimagLite::RuntimeError("algorithm");
	}

	virtual RealType tau() const
	{
		return unimplemented("tau");
	}

	virtual RealType maxTime() const
	{
		return unimplemented("maxTime");
	}

	virtual SizeType timeSteps() const
	{
		return unimplementedInt("timeSteps");
	}

	virtual bool noOperator() const
	{
		return static_cast<bool>(unimplemented("noOperator"));
	}

	virtual void noOperator(bool)
	{
		unimplemented("noOperator");
	}

	virtual bool skipTimeZero() const
	{
		return static_cast<bool>(unimplemented("skipTimeZero"));
	}

	virtual bool isEnergyForExp() const
	{
		return static_cast<bool>(unimplemented("isEnergyForExp"));
	}

	virtual RealType energyForExp() const
	{
		return unimplemented("energyForExp");
	}

	virtual RealType gsWeight() const
	{
		return unimplemented("gsWeight");
	}

	virtual RealType timeDirection() const
	{
		return unimplemented("timeDirection");
	}

private:

	RealType unimplemented(PsimagLite::String s) const
	{
		s = "TargetParamsBase: unimplemented " + s + "\n";
		throw PsimagLite::RuntimeError(s);
	}

	SizeType unimplementedInt(PsimagLite::String s) const
	{
		s = "TargetParamsBase: unimplemented " + s + "\n";
		throw PsimagLite::RuntimeError(s);
	}
}; // class TargetParamsBase

} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_BASE_H

