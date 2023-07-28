/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file ParametersModelTjMultiOrb.h
 *
 *  Contains the parameters for the Hubbard model and function to
 *  read them from a file
 *
 */
#ifndef DMRG_PARAMS_TJMULTIORB_H
#define DMRG_PARAMS_TJMULTIORB_H
#include "ParametersModelBase.h"

namespace Dmrg
{
//! ModelTjMultiOrb Parameters
template <typename RealType, typename QnType>
struct ParametersModelTjMultiOrb : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;

	template <typename IoInputType>
	ParametersModelTjMultiOrb(IoInputType& io)
	    : BaseType(io, false)
	    , reinterpretAndTruncate(0)
	{
		io.read(potentialV, "potentialV");
		io.readline(orbitals, "Orbitals=");

		try {
			io.readline(reinterpretAndTruncate, "JHundInfinity=");
		} catch (std::exception&) {
		}

		if (orbitals != 2 && reinterpretAndTruncate > 0)
			throw PsimagLite::RuntimeError("JHundInfinity>0 only possible for orbitals==2\n");

		if (reinterpretAndTruncate > 3)
			throw PsimagLite::RuntimeError("JHundInfinity must be less or equal to 3\n");
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	    SizeType,
	    PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1,
	    PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelTjMultiOrb";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/orbitals", orbitals);
		io.write(label + "/reinterpretAndTruncate", reinterpretAndTruncate);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os,
	    const ParametersModelTjMultiOrb& parameters)
	{
		os << "potentialV\n";
		os << parameters.potentialV;
		os << "orbitals=" << parameters.orbitals << "\n";
		os << "JHundInfinity=" << parameters.reinterpretAndTruncate << "\n";
		return os;
	}

	// Do not include here connection parameters
	typename PsimagLite::Vector<RealType>::Type potentialV;
	SizeType orbitals;
	SizeType reinterpretAndTruncate;
};
} // namespace Dmrg

/*@}*/
#endif
