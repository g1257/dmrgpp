/*
Copyright (c) 2009-2011, 2013, UT-Battelle, LLC
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

/*! \file CorrectionVectorParams.h
 *
 *  This is a structure to represent the parameters of the
 *  Correction vector DMRG algorithm.
 *  Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef CORRECTION_V_PARAMS_H
#define CORRECTION_V_PARAMS_H

#include "TargetParamsCommon.h"

namespace Dmrg {
//! Coordinates reading of TargetSTructure from input file
template<typename ModelType>
class CorrectionVectorParams : public TargetParamsCommon<ModelType> {

public:

	typedef TargetParamsCommon<ModelType> BaseType;
	typedef typename ModelType::RealType RealType;

	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrReal;
	typedef PsimagLite::Matrix<ComplexOrReal> MatrixType;

	static SizeType const PRODUCT = BaseType::PRODUCT;
	static SizeType const SUM = BaseType::SUM;

	template<typename IoInputter>
	CorrectionVectorParams(IoInputter& io,const ModelType& model)
	    : BaseType(io,model),
	      cgSteps_(1000),
	      cgEps_(1e-6)
	{
		this->setConcatenation(SUM);
		io.readline(correctionA_,"CorrectionA=");
		io.readline(type_,"DynamicDmrgType=");
		io.readline(steps_,"DynamicDmrgSteps=");
		io.readline(eps_,"DynamicDmrgEps=");
		io.readline(omega_,"CorrectionVectorOmega=");
		io.readline(eta_,"CorrectionVectorEta=");

		try {
			io.readline(cgSteps_,"ConjugateGradientSteps=");
		} catch (std::exception& e) {}

		try {
			io.readline(cgEps_,"ConjugateGradientEps=");
		} catch (std::exception& e) {}
	}

	virtual RealType correctionA() const
	{
		return correctionA_;
	}

	virtual SizeType type() const
	{
		return type_;
	}

	virtual SizeType steps() const
	{
		return steps_;
	}

	virtual RealType eps() const
	{
		return eps_;
	}

	virtual SizeType cgSteps() const
	{
		return cgSteps_;
	}

	virtual RealType omega() const
	{
		return omega_;
	}

	virtual RealType eta() const
	{
		return eta_;
	}

	virtual SizeType cgEps() const
	{
		return cgEps_;
	}

private:

	SizeType type_;
	SizeType steps_;
	RealType eps_;
	RealType correctionA_;
	SizeType cgSteps_;
	RealType omega_;
	RealType eta_;
	RealType cgEps_;
}; // class CorrectionVectorParams

template<typename ModelType>
inline std::ostream& operator<<(std::ostream& os,
                                const CorrectionVectorParams<ModelType>& t)
{
	os<<"#TargetParams.type=AdaptiveDynamic\n";
	const TargetParamsCommon<ModelType>& tp = t;
	os<<tp;
	os<<"DynamicDmrgType="<<t.type()<<"\n";
	os<<"DynamicDmrgSteps="<<t.steps()<<"\n";
	os<<"DynamicDmrgEps="<<t.eps()<<"\n";
	os<<"CorrectionVectorOmega="<<t.omega()<<"\n";
	os<<"CorrectionVectorEta="<<t.eta()<<"\n";
	os<<"ConjugateGradientSteps"<<t.cgSteps()<<"\n";
	os<<"ConjugateGradientEps"<<t.cgEps()<<"\n";

	return os;
}
} // namespace Dmrg 

/*@}*/
#endif //CORRECTION_V_PARAMS_H

