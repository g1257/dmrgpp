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

namespace Dmrg {
//! FeAs Model Parameters
template<typename Field>
struct ParametersModelFeAs {
	// no connections here please!!
	// connections are handled by the geometry

	template<typename IoInputType>
	ParametersModelFeAs(IoInputType& io)
	    : potentialT(0),feAsMode(0),coulombV(0),magneticField(0,0),minElectronsPerSite(0)
	{
		io.readline(orbitals,"Orbitals=");
		io.read(hubbardU,"hubbardU");
		io.read(potentialV,"potentialV");


		bool decayInInputFile = false;
		try {
			io.readline(feAsMode,"Decay=");
			decayInInputFile = true;
		} catch (std::exception& e) {}

		if (decayInInputFile) {
			PsimagLite::String str("Please use FeAsMode= instead of Decay=");
			str += " in input file\n";
			throw PsimagLite::RuntimeError(str);
		}

		try {
			io.readline(feAsMode,"FeAsMode=");
		} catch (std::exception& e) {}

		if (feAsMode > 2)
			throw PsimagLite::RuntimeError("FeAsMode: expecting 0 or 1 or 2\n");

		if (feAsMode > 0) {
			SizeType tmp = orbitals * orbitals;
			if (feAsMode == 2) tmp *= 2;
			if (hubbardU.size() != tmp) {
				PsimagLite::String str("FeAsMode: expecting");
				str += ttos(tmp) + " U values\n";
				throw PsimagLite::RuntimeError(str);
			}
		}

		if (feAsMode == 1) {
			if (orbitals != 3)
				throw PsimagLite::RuntimeError("FeAsMode: expecting 3 orbitals\n");
			io.readline(coulombV,"CoulombV=");
		}

		try {
			io.readMatrix(magneticField,"MagneticField");
		} catch (std::exception& e) {}

		try {
			io.read(potentialT,"potentialT");
		} catch (std::exception& e) {}

		if (magneticField.n_row()!=0 && magneticField.n_row()!=3)
			throw PsimagLite::RuntimeError("Magnetic Field: if present must have 3 rows\n");
		if (magneticField.n_row()!=0 && magneticField.n_col()!=potentialV.size())
			throw PsimagLite::RuntimeError("Magnetic Field: Expecting columns equal sites\n");

		try {
			io.readline(minElectronsPerSite,"MinElectronsPerSite=");
		} catch (std::exception& e) {}
	}

	SizeType orbitals;
	// Hubbard U values (one for each site)
	typename PsimagLite::Vector<Field>::Type hubbardU;
	// Onsite potential values, one for each site
	typename PsimagLite::Vector<Field>::Type potentialV;
	typename PsimagLite::Vector<Field>::Type potentialT;
	SizeType feAsMode;
	Field coulombV;
	PsimagLite::Matrix<Field> magneticField;
	SizeType minElectronsPerSite;
};

//! Function that prints model parameters to stream os
template<typename FieldType>
std::ostream& operator<<(std::ostream &os,const ParametersModelFeAs<FieldType>& parameters)
{
	os<<"hubbardU\n";
	os<<parameters.hubbardU;
	os<<"potentialV\n";
	os<<parameters.potentialV;
	if (parameters.magneticField.n_row()>0) {
		os<<"MagneticField\n";
		os<<parameters.magneticField;
	}

	os<<"FeAsMode="<<parameters.feAsMode<<"\n";
	if (parameters.feAsMode == 1)
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
