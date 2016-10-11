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

/*! \file ParametersModelFeAs.h
 *
 *  Contains the parameters for the FeAs model and function to read them from a JSON file
 *
 */
#ifndef PARAMETERSMODELFEAS_H
#define PARAMETERSMODELFEAS_H
#include <stdexcept>
#include <vector>
#include "Matrix.h"
#include "TargetQuantumElectrons.h"

namespace Dmrg {
//! FeAs Model Parameters
template<typename ComplexOrRealType>
struct ParametersModelFeAs {
	// no connections here please!!
	// connections are handled by the geometry

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	enum IntEnum {INT_PAPER33,
		          INT_V,
		          INT_CODE2,
		          INT_IMPURITY,
		          INT_KSPACE,
		          INT_ORBITAL0};

	static PsimagLite::String modeString(IntEnum x)
	{
		switch (x) {
		case INT_PAPER33:
			return "INT_PAPER33";
		case INT_V:
			return "INT_V";
		case INT_CODE2:
			return "INT_CODE2";
		case INT_IMPURITY:
			return "INT_IMPURITY";
		case INT_KSPACE:
			return "INT_KSPACE";
		case INT_ORBITAL0:
			return "INT_ORBITAL0";
		}

		return "UNKNOWN";
	}

	static IntEnum convertToEnum(PsimagLite::String x)
	{
		if (x == "INT_PAPER33")
			return INT_PAPER33;

		if (x == "INT_V")
			return INT_V;

		if (x == "INT_CODE2")
			return INT_CODE2;

		if (x == "INT_IMPURITY")
			return INT_IMPURITY;

		if (x == "INT_KSPACE")
			return INT_KSPACE;

		if (x == "INT_ORBITAL0")
			return INT_ORBITAL0;

		PsimagLite::String all = "INT_PAPER33 INT_V INT_CODE2 INT_IMPURITY";
		all += PsimagLite::String(" INT_KSPACE") + " INT_ORBITAL0";
		throw PsimagLite::RuntimeError("FeAsMode= can only be one of " + all + "\n");
	}

	template<typename IoInputType>
	ParametersModelFeAs(IoInputType& io)
	    : targetQuantum(io),
	      minElectronsPerSite(0),
	      potentialT(0),
	      feAsMode(INT_PAPER33),
	      coulombV(0),
	      magneticField(0,0),
	      jzSymmetry(false)
	{
		io.readline(orbitals,"Orbitals=");
		io.read(hubbardU,"hubbardU");
		io.read(potentialV,"potentialV");

		bool decayInInputFile = false;
		try {
			PsimagLite::String tmp;
			io.readline(tmp,"Decay=");
			feAsMode = convertToEnum(tmp);
			decayInInputFile = true;
		} catch (std::exception&) {}

		if (decayInInputFile) {
			PsimagLite::String str("Please use FeAsMode= instead of Decay=");
			str += " in input file\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String tmp;
		io.readline(tmp,"FeAsMode=");
		feAsMode = convertToEnum(tmp);

		if (feAsMode == INT_V || feAsMode == INT_CODE2) {
			SizeType tmp = orbitals * orbitals;
			if (feAsMode == INT_CODE2) tmp *= 2;
			if (hubbardU.size() != tmp) {
				PsimagLite::String str("FeAsMode: expecting ");
				str += ttos(tmp) + " U values\n";
				throw PsimagLite::RuntimeError(str);
			}
		}

		if (feAsMode == INT_V) {
			if (orbitals != 3)
				throw PsimagLite::RuntimeError("FeAsMode: expecting 3 orbitals\n");
			io.readline(coulombV,"CoulombV=");
		}

		if (feAsMode == INT_PAPER33 || feAsMode == INT_IMPURITY) {
			if (hubbardU.size() != 4 && hubbardU.size() != 5) {
				PsimagLite::String str("FeAsMode: expecting");
				str +=  " 4 or 5 U values\n";
				throw PsimagLite::RuntimeError(str);
			}

			if (hubbardU.size() == 4) {
				hubbardU.resize(5);
				hubbardU[4] = hubbardU[2];
			}

			std::cout<<"U[0]="<<hubbardU[0]<<" U\n";
			std::cout<<"U[1]="<<hubbardU[1]<<" U'-J/2\n";
			std::cout<<"U[2]="<<hubbardU[2]<<" -2J for S+S- + S-S+ term\n";
			std::cout<<"U[3]="<<hubbardU[3]<<" -J\n";
			std::cout<<"U[4]="<<hubbardU[4]<<" -2J for SzSz term\n";
		}

		if (feAsMode == INT_KSPACE) {
			if (hubbardU.size() != 1) {
				PsimagLite::String str("FeAsMode: expecting");
				str +=  " just 1 U value\n";
				throw PsimagLite::RuntimeError(str);
			}
		}
		try {
			io.readMatrix(magneticField,"magneticField");
		} catch (std::exception& e) {}

		try {
			io.read(potentialT,"potentialT");
		} catch (std::exception& e) {}

		if (magneticField.n_row()!=0 && magneticField.n_row()!=3)
			throw PsimagLite::RuntimeError("MagneticField: if present must have 3 rows\n");
		if (magneticField.n_row()!=0 && magneticField.n_col()!=potentialV.size())
			throw PsimagLite::RuntimeError("MagneticField: Expecting columns == sites\n");

		try {
			io.readline(minElectronsPerSite,"MinElectronsPerSite=");
		} catch (std::exception& e) {}

		try {
			io.readMatrix(spinOrbit,"SpinOrbit");
		} catch (std::exception&) {}

		try {
			io.readline(jzSymmetry,"JzSymmetry=");
		} catch (std::exception&) {}
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	//serializr start class ParametersModelFeAs

	TargetQuantumElectrons<RealType> targetQuantum;

	//serializr normal orbitals
	SizeType orbitals;
	//serializr normal minElectronsPerSite
	SizeType minElectronsPerSite;
	// Hubbard U values (one for each site)
	//serializr normal hubbardU
	typename PsimagLite::Vector<RealType>::Type hubbardU;
	// Onsite potential values, one for each site
	//serializr normal potentialV
	typename PsimagLite::Vector<RealType>::Type potentialV;
	//serializr normal potentialT
	typename PsimagLite::Vector<RealType>::Type potentialT;
	//serializr normal feAsMode
	IntEnum feAsMode;
	//serializr normal coulombV
	RealType coulombV;
	//serializr normal magneticField
	PsimagLite::Matrix<RealType> magneticField;
	PsimagLite::Matrix<ComplexOrRealType> spinOrbit;
	SizeType jzSymmetry;
};

//! Function that prints model parameters to stream os
template<typename RealType>
std::ostream& operator<<(std::ostream &os,
                         const ParametersModelFeAs<RealType>& parameters)
{
	os<<parameters.targetQuantum;
	os<<"Orbitals="<<parameters.orbitals<<"\n";
	os<<"hubbardU\n";
	os<<parameters.hubbardU;
	os<<"potentialV\n";
	os<<parameters.potentialV;
	if (parameters.magneticField.n_row()>0) {
		os<<"magneticField\n";
		os<<parameters.magneticField;
	}
	if (parameters.spinOrbit.n_row()>0) {
		os<<"SpinOrbit\n";
		os<<parameters.spinOrbit;
	}

	if (parameters.jzSymmetry>0) {
		os<<"using jzSymmetry, works only for 3 orbitals \n";
		os<<parameters.jzSymmetry;
	}
	os<<"FeAsMode=";
	os<<ParametersModelFeAs<RealType>::modeString(parameters.feAsMode)<<"\n";
	if (parameters.feAsMode == ParametersModelFeAs<RealType>::INT_V)
		os<<"CoulombV="<<parameters.coulombV<<"\n";

	if (parameters.potentialT.size()>0) {
		os<<"potentialT\n";
		os<<parameters.potentialT;
	}

	os<<"MinElectronsPerSite="<<parameters.minElectronsPerSite<<"\n";

	return os;
}
} // namespace Dmrg

/*@}*/
#endif

