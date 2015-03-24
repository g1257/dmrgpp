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
template<typename RealType>
struct ParametersModelFeAs {
	// no connections here please!!
	// connections are handled by the geometry

	template<typename IoInputType>
	ParametersModelFeAs(IoInputType& io)
	    : targetQuantum(io),
	      minElectronsPerSite(0),
	      potentialT(0),
	      feAsMode(0),
	      coulombV(0),
	      magneticField(0,0)
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

		io.readline(feAsMode,"FeAsMode=");

		if (feAsMode > 4)
			throw PsimagLite::RuntimeError("FeAsMode: expecting 0 to 4\n");

		if (feAsMode == 1 || feAsMode == 2) {
			SizeType tmp = orbitals * orbitals;
			if (feAsMode == 2) tmp *= 2;
			if (hubbardU.size() != tmp) {
				PsimagLite::String str("FeAsMode: expecting ");
				str += ttos(tmp) + " U values\n";
				throw PsimagLite::RuntimeError(str);
			}
		}

		if (feAsMode == 1) {
			if (orbitals != 3)
				throw PsimagLite::RuntimeError("FeAsMode: expecting 3 orbitals\n");
			io.readline(coulombV,"CoulombV=");
		}

		if (feAsMode == 0 || feAsMode == 3) {
			if (hubbardU.size() != 4) {
				PsimagLite::String str("FeAsMode: expecting");
				str +=  " 4 U values\n";
				throw PsimagLite::RuntimeError(str);
			}
		}

		if (feAsMode == 4) {
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
			throw PsimagLite::RuntimeError("Magnetic RealType: if present must have 3 rows\n");
		if (magneticField.n_row()!=0 && magneticField.n_col()!=potentialV.size())
			throw PsimagLite::RuntimeError("Magnetic RealType: Expecting columns equal sites\n");

		try {
			io.readline(minElectronsPerSite,"MinElectronsPerSite=");
		} catch (std::exception& e) {}
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "ParametersModelFeAs";

		const char* start = reinterpret_cast<const char *>(this);
		const char* end = reinterpret_cast<const char *>(&minElectronsPerSite);
		SizeType total = mres.memResolv(&orbitals, end-start, str + " orbitals");

		start = end;
		end = reinterpret_cast<const char *>(&hubbardU);
		total += mres.memResolv(&minElectronsPerSite, end-start, str + " minElectronsPerSite");

		start = end;
		end = reinterpret_cast<const char *>(&potentialV);
		total += mres.memResolv(&hubbardU, end-start, str + " hubbardU");

		start = end;
		end = reinterpret_cast<const char *>(&potentialT);
		total += mres.memResolv(&potentialV, end-start, str + " potentialV");

		start = end;
		end = reinterpret_cast<const char *>(&feAsMode);
		total += mres.memResolv(&potentialT, end-start, str + " potentialT");

		start = end;
		end = reinterpret_cast<const char *>(&coulombV);
		total += mres.memResolv(&feAsMode, end-start, str + " feAsMode");

		start = end;
		end = reinterpret_cast<const char *>(&magneticField);
		total += mres.memResolv(&coulombV, end-start, str + " coulombV");

		total += mres.memResolv(&magneticField,
		                        sizeof(*this) - total,
		                        str + " magneticField");

		return total;
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
	SizeType feAsMode;
	//serializr normal coulombV
	RealType coulombV;
	//serializr normal magneticField
	PsimagLite::Matrix<RealType> magneticField;
};

//! Function that prints model parameters to stream os
template<typename RealType>
std::ostream& operator<<(std::ostream &os,
                         const ParametersModelFeAs<RealType>& parameters)
{
	os<<parameters.targetQuantum;
	os<<"hubbardU\n";
	os<<parameters.hubbardU;
	os<<"potentialV\n";
	os<<parameters.potentialV;
	if (parameters.magneticField.n_row()>0) {
		os<<"magneticField\n";
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

