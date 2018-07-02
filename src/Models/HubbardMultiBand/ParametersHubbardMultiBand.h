/*
Copyright (c) 2009-2018, UT-Battelle, LLC
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

/*! \file ParametersHubbardMultiBand.h
 *
 *  Contains the parameters for the HubbardMultiBand model
 *
 */
#ifndef PARAMS_HUBBARD_MULTI_BAND_H
#define PARAMS_HUBBARD_MULTI_BAND_H
#include <stdexcept>
#include <vector>
#include "Matrix.h"
#include "TargetQuantumElectrons.h"

namespace Dmrg {
//! FeAs Model Parameters
template<typename ComplexOrRealType>
struct ParametersHubbardMultiBand {
	// no connections here please!!
	// connections are handled by the geometry

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<PsimagLite::Matrix<ComplexOrRealType> >::Type VectorType;

	template<typename IoInputType>
	ParametersHubbardMultiBand(IoInputType& io)
	    : targetQuantum(io)
	{
		io.readline(orbitals,"Orbitals=");
		io.read(hubbardU,"hubbardU");
		io.read(potentialV,"potentialV");
		SizeType h = 1;
		try {
			io.readline(h, "NumberOfHoppingOrbitalMatrices=");
		} catch (std::exception&) {}

		SizeType sites = 0;
		io.readline(sites, "TotalNumberOfSites=");
		if (h != 1 && h != sites)
			err("NumberOfHoppingOrbitalMatrices=1 or =numberOfSites\n");

		for (SizeType i = 0; i < h; ++i) {
			PsimagLite::Matrix<ComplexOrRealType> m;
			io.read(m, "hopOnSite");
			hopOnSite.push_back(m);
		}
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersHubbardMultiBand";
		io.createGroup(label);
		targetQuantum.write(label, io);
		io.write(label + "/orbitals", orbitals);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/hopOnSite", hopOnSite);
	}

	TargetQuantumElectrons<RealType> targetQuantum;
	SizeType orbitals;
	VectorRealType hubbardU;
	// Onsite potential values, one for each site
	VectorRealType potentialV;
	VectorType hopOnSite;
};

//! Function that prints model parameters to stream os
template<typename RealType>
std::ostream& operator<<(std::ostream &os,
                         const ParametersHubbardMultiBand<RealType>& parameters)
{
	os<<parameters.targetQuantum;
	os<<"Orbitals="<<parameters.orbitals<<"\n";
	os<<"hubbardU\n";
	os<<parameters.hubbardU;
	os<<"potentialV\n";
	os<<parameters.potentialV;

	return os;
}
} // namespace Dmrg

/*@}*/
#endif

