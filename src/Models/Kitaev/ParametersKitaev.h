/*
Copyright (c) 2009-2012-2018, UT-Battelle, LLC
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

/*! \file ParametersKitaev.h
 *
 *  Contains the parameters for the Kitaev model and function
 *  to read them from a file (started March 2018)
 *
 */
#ifndef PARAMETERS_KITAEV_H
#define PARAMETERS_KITAEV_H
#include "Vector.h"
#include "ParametersModelBase.h"

namespace Dmrg {
//! Kitaev Model Parameters
// no connectors here, connectors are handled by the geometry
template<typename RealType, typename QnType>
struct ParametersKitaev : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	template<typename IoInputType>
	ParametersKitaev(IoInputType& io)
	    : BaseType(io,false)
	{
		SizeType nsites = BaseType::targetQuantum().totalNumberOfSites;
		try {
			magneticFieldX.resize(nsites, 0.0);
			io.read(magneticFieldX,"MagneticFieldX");
			std::cerr<<"Has MagneticFieldX \n";
		} catch (std::exception&) {
			magneticFieldX.clear();
		}

		try {
			magneticFieldY.resize(nsites, 0.0);
			io.read(magneticFieldY,"MagneticFieldY");
			std::cerr<<"Has MagneticFieldY \n";
		} catch (std::exception&) {
			magneticFieldY.clear();
		}

		try {
			magneticFieldZ.resize(nsites, 0.0);
			io.read(magneticFieldZ,"MagneticFieldZ");
			std::cerr<<"Has MagneticFieldZ \n";
		} catch (std::exception&) {
			magneticFieldZ.clear();
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
		PsimagLite::String label = label1 + "/ParametersKitaev";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/magneticFieldX", magneticFieldX);
		io.write(label + "/magneticFieldY", magneticFieldY);
		io.write(label + "/magneticFieldZ", magneticFieldZ);
	}


	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream &os,
	                                const ParametersKitaev& parameters)
	{
		os<<"MagneticFieldX="<<parameters.magneticFieldX<<"\n";
		os<<"MagneticFieldY="<<parameters.magneticFieldY<<"\n";
		os<<"MagneticFieldZ="<<parameters.magneticFieldZ<<"\n";
		os<<parameters.targetQuantum;
		return os;
	}

	VectorRealType magneticFieldX;
	VectorRealType magneticFieldY;
	VectorRealType magneticFieldZ;
};
} // namespace Dmrg

/*@}*/
#endif

