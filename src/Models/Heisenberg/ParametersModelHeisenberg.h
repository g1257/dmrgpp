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
#include "../../Engine/ParametersModelBase.h"
#include "Vector.h"

namespace Dmrg {
//! Heisenberg Model Parameters
template <typename RealType, typename QnType>
struct ParametersModelHeisenberg : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	// no connectors here, connectors are handled by the geometry
	template <typename IoInputType>
	ParametersModelHeisenberg(IoInputType& io)
	    : BaseType(io, false)
	    , twiceTheSpinBorder(0)
	{
		PsimagLite::String model;
		io.readline(model, "Model=");

		io.readline(twiceTheSpin, "HeisenbergTwiceS=");

		if (model == "HeisenbergMix")
			io.readline(twiceTheSpinBorder, "HeisenbergTwiceSborder=");

		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");

		try {
			magneticFieldX.resize(nsites);
			io.read(magneticFieldX, "MagneticFieldX");
		} catch (std::exception&) {
			magneticFieldX.clear();
		}

		try {
			magneticFieldZ.resize(nsites);
			io.read(magneticFieldZ, "MagneticFieldZ");
		} catch (std::exception&) {
			magneticFieldZ.clear();
		}

		// throw if supplying MagneticField label
		bool invalidLabel = false;
		try {
			VectorRealType tmpVector;
			io.read(tmpVector, "MagneticField=");
			invalidLabel = true;
		} catch (std::exception&) { }

		if (invalidLabel) {
			throw PsimagLite::RuntimeError(
			    "MagneticField label is no longer supported.\n"
			    + PsimagLite::String("Please use MagneticField[XZ] instead\n"));
		}

		// throw if supplying MagneticFieldDirection label
		try {
			PsimagLite::String tmpStr;
			io.readline(tmpStr, "MagneticFieldDirection=");
			invalidLabel = true;
		} catch (std::exception&) { }

		if (invalidLabel) {
			throw PsimagLite::RuntimeError(
			    "MagneticFieldDirection label is no longer supported.\n"
			    + PsimagLite::String("Please use MagneticField[XZ] instead\n"));
		}

		try {
			io.read(anisotropyD, "AnisotropyD");
		} catch (std::exception&) { }

		try {
			io.read(anisotropyE, "AnisotropyE");
		} catch (std::exception&) { }
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelHeisenberg";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/twiceTheSpin", twiceTheSpin);
		io.write(label + "/magneticFieldX", magneticFieldX);
		io.write(label + "/magneticFieldZ", magneticFieldZ);
		io.write(label + "/anisotropyD", anisotropyD);
		io.write(label + "/anisotropyE", anisotropyE);
	}

	static void checkMagneticField(SizeType s, unsigned char c, SizeType n)
	{
		if (s == 0 || s == n)
			return;

		PsimagLite::String msg("ModelHeisenberg: If provided, ");
		msg += " MagneticField" + ttos(c) + " must be a vector of " + ttos(n)
		    + " entries.\n";
		err(msg);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os,
	                                const ParametersModelHeisenberg& parameters)
	{
		if (!parameters.magneticFieldX.empty())
			os << "MagneticFieldX=" << parameters.magneticFieldX << "\n";
		if (!parameters.magneticFieldZ.empty())
			os << "MagneticFieldZ=" << parameters.magneticFieldZ << "\n";

		os << "AnisotropyD=" << parameters.anisotropy << "\n";
		os << "HeisenbergTwiceS=" << parameters.twiceTheSpin << "\n";
		os << parameters.targetQuantum;
		return os;
	}

	SizeType twiceTheSpin;
	SizeType twiceTheSpinBorder;
	VectorRealType magneticFieldX;
	VectorRealType magneticFieldZ;
	VectorRealType anisotropyD;
	VectorRealType anisotropyE;
};
} // namespace Dmrg

/*@}*/
#endif
