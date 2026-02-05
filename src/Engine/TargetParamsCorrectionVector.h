/*
Copyright (c) 2009-2016-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file TargetParamsCorrectionVector.h
 *
 *  This is a structure to represent the parameters of the
 *  Correction vector DMRG algorithm.
 *  Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef TARGET_PARAMS_CORRECTION_V_H
#define TARGET_PARAMS_CORRECTION_V_H

#include "FreqEnum.h"
#include "TargetParamsCommon.h"

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType>
class TargetParamsCorrectionVector : public TargetParamsCommon<ModelType> {

public:

	using BaseType         = TargetParamsCommon<ModelType>;
	using RealType         = typename ModelType::RealType;
	using PairFreqType     = typename BaseType::BaseType::PairFreqType;
	using OperatorType     = typename ModelType::OperatorType;
	using PairType         = typename OperatorType::PairType;
	using SparseMatrixType = typename OperatorType::StorageType;
	using ComplexOrReal    = typename SparseMatrixType::value_type;
	using MatrixType       = PsimagLite::Matrix<ComplexOrReal>;

	template <typename IoInputter>
	TargetParamsCorrectionVector(IoInputter&        io,
	                             PsimagLite::String targeting,
	                             const ModelType&   model)
	    : BaseType(io, targeting, model)
	    , cgSteps_(1000)
	    , firstRitz_(0)
	    , nForFraction_(1)
	    , advanceEach_(0)
	    , cgEps_(1e-6)
	{
		io.readline(correctionA_, "CorrectionA=");
		io.readline(type_, "DynamicDmrgType=");
		PsimagLite::String tmp;
		io.readline(tmp, "CorrectionVectorFreqType=");
		PsimagLite::FreqEnum freqEnum = PsimagLite::FREQ_REAL;

		if (tmp == "Matsubara") {
			freqEnum = PsimagLite::FREQ_MATSUBARA;
		} else if (tmp != "Real") {
			PsimagLite::String msg("CorrectionVectorFreqType");
			throw PsimagLite::RuntimeError(msg += "must be either Real or Matsubara\n");
		}

		RealType omega;
		io.readline(omega, "CorrectionVectorOmega=");
		omega_ = PairFreqType(freqEnum, omega);
		io.readline(eta_, "CorrectionVectorEta=");

		io.readline(tmp, "CorrectionVectorAlgorithm=");
		if (tmp == "Krylov") {
			algorithm_ = BaseType::AlgorithmEnum::KRYLOV;
		} else if (tmp == "ConjugateGradient") {
			algorithm_ = BaseType::AlgorithmEnum::CONJUGATE_GRADIENT;
		} else if (tmp == "Chebyshev") {
			algorithm_ = BaseType::AlgorithmEnum::CHEBYSHEV;
		} else if (tmp == "KrylovTime") {
			algorithm_ = BaseType::AlgorithmEnum::KRYLOVTIME;
		} else {
			PsimagLite::String str("TargetParamsCorrectionVector ");
			str += "Unknown algorithm " + tmp + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		try {
			io.readline(cgSteps_, "ConjugateGradientSteps=");
		} catch (std::exception&) { }

		try {
			io.readline(cgEps_, "ConjugateGradientEps=");
		} catch (std::exception&) { }

		try {
			int x = 0;
			io.readline(x, "TSPUseQns=");
			err("TSPUseQns= is no longer needed, please delete it from the input "
			    "file\n");
		} catch (std::exception&) { }

		try {
			io.readline(firstRitz_, "FirstRitz=");
		} catch (std::exception&) { }

		try {
			io.readline(nForFraction_, "CVnForFraction=");
		} catch (std::exception&) { }

		if (nForFraction_ > 1)
			io.readline(advanceEach_, "TSPAdvanceEach=");

		if (freqEnum == PsimagLite::FREQ_MATSUBARA && firstRitz_ != 0)
			err("FirstRitz must be 0 for Matsubara\n");
	}

	virtual RealType correctionA() const { return correctionA_; }

	virtual SizeType type() const { return type_; }

	virtual void type(SizeType x) { type_ = x; }

	virtual SizeType cgSteps() const { return cgSteps_; }

	virtual PairFreqType omega() const { return omega_; }

	virtual void omega(PsimagLite::FreqEnum freqEnum, RealType x)
	{
		omega_ = PairFreqType(freqEnum, x);
	}

	virtual RealType eta() const { return eta_; }

	virtual RealType cgEps() const { return cgEps_; }

	virtual typename BaseType::AlgorithmEnum algorithm() const { return algorithm_; }

	virtual SizeType firstRitz() const { return firstRitz_; }

	virtual SizeType nForFraction() const { return nForFraction_; }

	virtual SizeType advanceEach() const { return advanceEach_; }

private:

	SizeType                         type_;
	typename BaseType::AlgorithmEnum algorithm_;
	SizeType                         cgSteps_;
	SizeType                         firstRitz_;
	SizeType                         nForFraction_;
	SizeType                         advanceEach_;
	RealType                         correctionA_;
	PairFreqType                     omega_;
	RealType                         eta_;
	RealType                         cgEps_;
}; // class TargetParamsCorrectionVector

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const TargetParamsCorrectionVector<ModelType>& t)
{
	os << "TargetParams.type=AdaptiveDynamic\n";
	const TargetParamsCommon<ModelType>& tp = t;
	os << tp;
	os << "DynamicDmrgType=" << t.type() << "\n";
	os << "CorrectionVectorOmega=" << t.omega() << "\n";
	os << "CorrectionVectorEta=" << t.eta() << "\n";
	os << "ConjugateGradientSteps" << t.cgSteps() << "\n";
	os << "ConjugateGradientEps" << t.cgEps() << "\n";
	os << "firstRitz" << t.firstRitz() << "\n";
	os << "nForFraction" << t.nForFraction() << "\n";

	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_CORRECTION_V_H
