/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file ParametersModelFeAs.h
 *
 *  Contains the parameters for the FeAs model and function to read them from a JSON file
 *
 */
#ifndef LANCZOS_PARAMS_MODELFEAS_H
#define LANCZOS_PARAMS_MODELFEAS_H
#include "PsimagLite.h"

namespace LanczosPlusPlus {
//! FeAs Model Parameters
template <typename ComplexOrRealType> struct ParametersModelFeAs {
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	// no connections here please!!
	// connections are handled by the geometry

	enum class IntEnum
	{
		INT_PAPER33,
		INT_V,
		INT_CODE2,
		INT_IMPURITY,
		INT_KSPACE,
		INT_ORBITAL0
	};

	template <typename IoInputType>
	ParametersModelFeAs(IoInputType& io)
	    : feAsMode(IntEnum::INT_PAPER33)
	    , coulombV(0)
	    , anisotropyD(0)
	{
		io.readline(orbitals, "Orbitals=");
		io.read(hubbardU, "hubbardU");
		io.read(potentialV, "potentialV");

		bool decayInInputFile = false;
		try {
			PsimagLite::String tmp;
			io.readline(tmp, "Decay=");
			decayInInputFile = true;
		} catch (std::exception& e) { }

		if (decayInInputFile) {
			PsimagLite::String str("Please use FeAsMode= instead of Decay=");
			str += " in input file\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String tmp;
		io.readline(tmp, "FeAsMode=");
		feAsMode = convertToEnum(tmp);

		if (feAsMode == IntEnum::INT_V || feAsMode == IntEnum::INT_CODE2) {
			SizeType tmp = orbitals * orbitals;
			if (feAsMode == IntEnum::INT_CODE2)
				tmp *= 2;
			if (hubbardU.size() != tmp) {
				PsimagLite::String str("FeAsMode: expecting ");
				str += ttos(tmp) + " U values\n";
				throw PsimagLite::RuntimeError(str);
			}
		}

		if (feAsMode == IntEnum::INT_V) {
			if (orbitals != 3)
				throw PsimagLite::RuntimeError("FeAsMode: expecting 3 orbitals\n");
			io.readline(coulombV, "CoulombV=");
		}

		if (feAsMode == IntEnum::INT_PAPER33 || feAsMode == IntEnum::INT_IMPURITY) {
			if (hubbardU.size() < 4 || hubbardU.size() > 6) {
				PsimagLite::String str("FeAsMode: expecting");
				str += " 4 or 5 or 6 U values\n";
				throw PsimagLite::RuntimeError(str);
			}

			if (hubbardU.size() == 4 || hubbardU.size() == 5) {
				hubbardU.resize(6);
				hubbardU[4] = hubbardU[2];
				hubbardU[5] = 0.0;
			}

			try {
				io.read(spinOrbit, "SpinOrbit");
			} catch (std::exception&) { }

			std::cout << "U[0]=" << hubbardU[0] << " =U\n";
			std::cout << "U[1]=" << hubbardU[1] << " =U'-J/2\n";
			std::cout << "U[2]=" << hubbardU[2];
			std::cout << " = factor for 1/2(S+_aS-_b + S-_aS+_b) term\n";
			std::cout << "U[3]=" << hubbardU[3] << " =-J\n";
			std::cout << "U[4]=" << hubbardU[4] << " = factor for Sz_aSz_b term\n";
			std::cout << "U[5]=" << hubbardU[5] << " = factor for \\sum_\\sigma ";
			std::cout << "n_{a\\sigma}*n_{b\\sigma} term\n";
		}

		if (feAsMode == IntEnum::INT_KSPACE) {
			if (hubbardU.size() != 1) {
				PsimagLite::String str("FeAsMode: expecting");
				str += " just 1 U values\n";
				throw PsimagLite::RuntimeError(str);
			}
		}

		try {
			io.readline(anisotropyD, "AnisotropyD=");
		} catch (std::exception&) { }
	}

	static IntEnum convertToEnum(PsimagLite::String x)
	{
		if (x == "INT_PAPER33")
			return IntEnum::INT_PAPER33;

		if (x == "INT_V")
			return IntEnum::INT_V;

		if (x == "INT_CODE2")
			return IntEnum::INT_CODE2;

		if (x == "INT_IMPURITY")
			return IntEnum::INT_IMPURITY;

		if (x == "INT_KSPACE")
			return IntEnum::INT_KSPACE;

		// if (x == "INT_ORBITAL0")
		//	return IntEnum::INT_ORBITAL0;

		PsimagLite::String all = "INT_PAPER33 INT_V INT_CODE2 INT_IMPURITY";
		all += PsimagLite::String(" INT_KSPACE") + " INT_ORBITAL0";
		throw PsimagLite::RuntimeError("FeAsMode= can only be one of " + all + "\n");
	}

	static PsimagLite::String modeString(IntEnum x)
	{
		switch (x) {
		case IntEnum::INT_PAPER33:
			return "INT_PAPER33";
		case IntEnum::INT_V:
			return "INT_V";
		case IntEnum::INT_CODE2:
			return "INT_CODE2";
		case IntEnum::INT_IMPURITY:
			return "INT_IMPURITY";
		case IntEnum::INT_KSPACE:
			return "INT_KSPACE";
		case IntEnum::INT_ORBITAL0:
			return "INT_ORBITAL0";
		}

		return "UNKNOWN";
	}

	SizeType orbitals;
	// Hubbard U values (one for each site)
	typename PsimagLite::Vector<RealType>::Type hubbardU;
	// Onsite potential values, one for each site
	typename PsimagLite::Vector<RealType>::Type potentialV;
	IntEnum                                     feAsMode;
	RealType                                    coulombV;
	RealType                                    anisotropyD;
	PsimagLite::Matrix<ComplexOrRealType>       spinOrbit;
	// target number of electrons  in the system
	int nOfElectrons;
}; // struct ParametersModelFeAs

//! Function that prints model parameters to stream os
template <typename RealTypeType>
std::ostream& operator<<(std::ostream& os, const ParametersModelFeAs<RealTypeType>& parameters)
{
	os << "orbitals=" << parameters.orbitals << "\n";
	os << "hubbardU\n";
	os << parameters.hubbardU;
	os << "potentialV\n";
	os << parameters.potentialV;
	os << "SpinOrbit\n";
	os << parameters.spinOrbit;
	os << "FeAsMode=" << parameters.modeString(parameters.feAsMode) << "\n";
	os << "CoulombV=" << parameters.coulombV << "\n";
	os << "AnisotropyD=" << parameters.anisotropyD << "\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif
