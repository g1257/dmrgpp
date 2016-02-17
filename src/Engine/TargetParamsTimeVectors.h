/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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

/*! \file TargetParamsTimeVectors.h
 *
 *  This is a structure to represent the parameters of the TimeStep Evolution
 *  algorithm. Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef TARGET_PARAMS_TIME_VECTORS_H
#define TARGET_PARAMS_TIME_VECTORS_H
#include "TargetParamsCommon.h"

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template<typename ModelType>
class TargetParamsTimeVectors : public TargetParamsCommon<ModelType> {

	typedef TargetParamsCommon<ModelType> BaseType;

public:

	typedef typename ModelType::RealType RealType;

	template<typename IoInputter>
	TargetParamsTimeVectors(IoInputter& io,const ModelType& model)
	    : BaseType(io,model),
	      timeSteps_(0),
	      advanceEach_(0),
	      algorithm_(BaseType::KRYLOV),
	      tau_(0),
	      timeDirection_(1.0)
	{
		/*PSIDOC TargetParamsTimeVectors
		\item[TSPTau] [RealType], $\tau$ for the Krylov,
		see \cite{re:alvarez11} Section II.B and II.C.
		\item[TSPTimeSteps] [Integer]  $n_v$ as defined in
		\cite{re:alvarez11} Section II.B
		\item[TSPAdvanceEach] [Integer] Number of sites to sweep before
		advancing to the next time.
		\item[TSPAlgorithm] [String] Either
		\verb!Krylov! or \verb!RungeKutta! or \verb!SuzukiTrotter!\\
		Note that SuzukiTrotter is currently very experimental and unsupported.
		*/
		io.readline(tau_,"TSPTau=");
		io.readline(timeSteps_,"TSPTimeSteps=");
		io.readline(advanceEach_,"TSPAdvanceEach=");
		PsimagLite::String s="";

		try {
			io.readline(s,"TSPAlgorithm=");
			if (s=="RungeKutta" || s=="rungeKutta" || s=="rungekutta")
				algorithm_ = BaseType::RUNGE_KUTTA;
			if (s=="SuzukiTrotter" || s=="suzukiTrotter" || s=="suzukitrotter")
				algorithm_ = BaseType::SUZUKI_TROTTER;
		} catch (std::exception&) {
			PsimagLite::String s(__FILE__);
			s += "\n FATAL: TSPAlgorithm not found in input file.\n";
			s += "Please add either TSPAlgorithm=Krylov or TSPAlgorithm=RungeKutta";
			s += " or TSPAlgorithm=SuzukiTrotter just below the TSPAdvanceEach= ";
			s += " line in the input file.\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		try {
			io.readline(timeDirection_,"TSPTimeFactor=");
		} catch (std::exception&) {}
	}

	virtual SizeType timeSteps() const
	{
		return timeSteps_;
	}

	virtual SizeType advanceEach() const
	{
		return advanceEach_;
	}

	virtual SizeType algorithm() const
	{
		return algorithm_;
	}

	virtual RealType tau() const
	{
		return tau_;
	}

	virtual RealType timeDirection() const
	{
		return timeDirection_;
	}

private:

	SizeType timeSteps_;
	SizeType advanceEach_;
	SizeType algorithm_;
	RealType tau_;
	RealType timeDirection_;

}; // class TargetParamsTimeVectors

template<typename ModelType>
inline std::ostream&
operator<<(std::ostream& os,const TargetParamsTimeVectors<ModelType>& t)
{
	os<<"#TargetParams.type=TimeVectors";
	os<<"#TargetParams.tau="<<t.tau()<<"\n";
	os<<"#TargetParams.timeSteps="<<t.timeSteps()<<"\n";
	os<<"#TargetParams.advanceEach="<<t.advanceEach()<<"\n";
	os<<"#TargetParams.algorithm="<<t.algorithm()<<"\n";
	os<<"#TargetParams.timeDirection="<<t.timeDirection()<<"\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_TIME_VECTORS_H

