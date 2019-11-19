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

/*! \file ParametersModelHeisenberg.h
 *
 *  Contains the parameters for the Heisenberg model and function
 *  to read them from a file
 *
 */
#ifndef PARAMETERSMODELHEISENBERG_H
#define PARAMETERSMODELHEISENBERG_H
#include "Vector.h"
#include "ParametersModelBase.h"

namespace Dmrg {
//! Heisenberg Model Parameters
template<typename RealType, typename QnType>
struct ParametersModelHeisenberg : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	// no connectors here, connectors are handled by the geometry
	template<typename IoInputType>
	ParametersModelHeisenberg(IoInputType& io)
	    : BaseType(io, false)
	{
		io.readline(twiceTheSpin,"HeisenbergTwiceS=");
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");

		bool hasField = false;
		try {
			magneticFieldV.resize(nsites);
			io.read(magneticFieldV, "MagneticField");
			hasField = true;
		} catch (std::exception&) {}

		if (hasField) {
			magneticFieldDirection = "z";
			try {
				io.readline(magneticFieldDirection, "MagneticFieldDirection=");
			} catch (std::exception&) {
				std::cerr<<"WARNING: MagneticFieldDirection= not given, assuming z\n";
			}

			if (magneticFieldDirection != "x" && magneticFieldDirection != "z")
				err("magneticFieldDirection must be in {x, z}\n");
		} else {
			magneticFieldV.clear();
		}

		try {
			io.read(anisotropyD,"AnisotropyD");
		} catch (std::exception&) {}

		try {
			io.read(anisotropyE,"AnisotropyE");
		} catch (std::exception&) {}
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
		PsimagLite::String label = label1 + "/ParametersModelHeisenberg";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/twiceTheSpin", twiceTheSpin);
		io.write(label + "/magneticFieldV", magneticFieldV);
		io.write(label + "/magneticFieldDirection", magneticFieldDirection);
		io.write(label + "/anisotropyD", anisotropyD);
		io.write(label + "/anisotropyE", anisotropyE);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream &os,
	                                const ParametersModelHeisenberg& parameters)
	{
		os<<"MagneticField="<<parameters.magneticFieldV<<"\n";
		os<<"magneticFieldDirection="<<parameters.magneticFieldDirection<<"\n";
		os<<"AnisotropyD="<<parameters.anisotropy<<"\n";
		os<<"HeisenbergTwiceS="<<parameters.twiceTheSpin<<"\n";
		os<<parameters.targetQuantum;
		return os;
	}

	SizeType twiceTheSpin;
	PsimagLite::String magneticFieldDirection;
	VectorRealType magneticFieldV;
	VectorRealType anisotropyD;
	VectorRealType anisotropyE;
};
} // namespace Dmrg

/*@}*/
#endif

