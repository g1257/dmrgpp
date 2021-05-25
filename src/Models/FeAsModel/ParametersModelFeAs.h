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
#include "ParametersModelBase.h"

namespace Dmrg {

//! FeAs Model Parameters
template<typename ComplexOrRealType, typename QnType>
struct ParametersModelFeAs : public ParametersModelBase<ComplexOrRealType, QnType> {
	// no connections here please!!
	// connections are handled by the geometry

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef ParametersModelBase<ComplexOrRealType, QnType> BaseType;

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
	    : BaseType(io, false),
	      potentialT(0),
	      feAsMode(INT_PAPER33),
	      coulombV(0),
	      magneticField(0,0),
	      jzSymmetry(false),
	      orbDependence(false)
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

		try {
			io.readline(orbDependence,"OrbDependence=");
		} catch (std::exception&) {}

		if (feAsMode == INT_PAPER33 || feAsMode == INT_IMPURITY) {
			if (!orbDependence){
				if (hubbardU.size() != 4 && hubbardU.size() != 5) {
					PsimagLite::String str("FeAsMode: expecting");
					str +=  " 4 or 5 U values\n";
					throw PsimagLite::RuntimeError(str);
				}
			} else {
				if (orbitals==2 && hubbardU.size() != 5) {
					PsimagLite::String str("FeAsMode: expecting");
					str +=  " 5 U values with OrbDependence and 2 orbitals\n";
					throw PsimagLite::RuntimeError(str);
				}
				if (orbitals==3 && hubbardU.size() != 12) {
					PsimagLite::String str("FeAsMode: expecting");
					str +=  " 12 U values with OrbDependence and 3 orbitals\n";
					throw PsimagLite::RuntimeError(str);
				}
			}

			if (hubbardU.size() == 4) {
				hubbardU.resize(5);
				hubbardU[4] = hubbardU[2];
			}

			if (!orbDependence){
				std::cout<<"U[0]="<<hubbardU[0]<<" U\n";
				std::cout<<"U[1]="<<hubbardU[1]<<" U'-J/2\n";
				std::cout<<"U[2]="<<hubbardU[2]<<" -2J for S+S- + S-S+ term\n";
				std::cout<<"U[3]="<<hubbardU[3]<<" -J\n";
				std::cout<<"U[4]="<<hubbardU[4]<<" -2J for SzSz term\n";
			} else {
				if (orbitals==2) {
					std::cout<<"U[0,0]="<<hubbardU[0]<<" U_0\n";
					std::cout<<"U[0,1]="<<hubbardU[1]<<" U_1\n";
					std::cout<<"U[1,0]="<<hubbardU[2]<<" U00\n";
					std::cout<<"U[2,0]="<<hubbardU[3]<<" S.S\n";
					std::cout<<"U[3,0]="<<hubbardU[4]<<" PairHop\n";
				} else {
					std::cout<<"U[0,0]="<<hubbardU[0]<<" U_0\n";
					std::cout<<"U[0,1]="<<hubbardU[1]<<" U_1\n";
					std::cout<<"U[0,2]="<<hubbardU[2]<<" U_2\n";
					std::cout<<"U[1,0]="<<hubbardU[3]<<" Un0n1\n";
					std::cout<<"U[1,1]="<<hubbardU[4]<<" Un0n2\n";
					std::cout<<"U[1,2]="<<hubbardU[5]<<" Un1n2\n";
					std::cout<<"U[2,0]="<<hubbardU[6]<<" S0S1\n";
					std::cout<<"U[2,1]="<<hubbardU[7]<<" S0S2\n";
					std::cout<<"U[2,2]="<<hubbardU[8]<<" S1S2\n";
					std::cout<<"U[3,0]="<<hubbardU[9]<<" P0P1\n";
					std::cout<<"U[3,1]="<<hubbardU[10]<<" P0P2\n";
					std::cout<<"U[3,2]="<<hubbardU[11]<<" P1P2\n";
				}
			}
		}

		if (feAsMode == INT_KSPACE) {
			if (hubbardU.size() != 1) {
				PsimagLite::String str("FeAsMode: expecting");
				str +=  " just 1 U value\n";
				throw PsimagLite::RuntimeError(str);
			}
		}
		try {
			io.read(magneticField, "magneticField");
		} catch (std::exception&) {}

		try {
			io.read(potentialT,"potentialT");
		} catch (std::exception&) {}

		if (magneticField.rows()!=0 && magneticField.rows()!=3)
			throw PsimagLite::RuntimeError("MagneticField: if present must have 3 rows\n");
		if (magneticField.rows()!=0 && magneticField.cols()!=potentialV.size())
			throw PsimagLite::RuntimeError("MagneticField: Expecting columns == sites\n");

		try {
			io.readline(jzSymmetry,"JzSymmetry=");
		} catch (std::exception&) {}

		bool hasSpinOrbitMatrix = false;
		try {
			io.read(spinOrbit, "SpinOrbit");
			hasSpinOrbitMatrix = true;
		} catch (std::exception&) {}

		if (!hasSpinOrbitMatrix && jzSymmetry > 0)
			err("jzSymmetry > 0 needs SpinOrbit matrix in input file\n");

		if (hasSpinOrbitMatrix && jzSymmetry == 0)
			err("SpinOrbit matrix found in input but jzSymmetry set to 0\n");

		std::cout<<"JzSymmetry="<<jzSymmetry<<std::endl;
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
		PsimagLite::String label = label1 + "/ParametersModelFeAs";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/orbitals", orbitals);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/potentialT", potentialT);
		io.write(label + "/feAsMode", feAsMode);
		io.write(label + "/coulombV", coulombV);
		magneticField.write(label + "/magneticField", io);
		spinOrbit.write(label + "/spinOrbit", io);
		io.write(label + "/jzSymmetry", jzSymmetry);
		io.write(label + "/orbDependence", orbDependence);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream &os,
	                                const ParametersModelFeAs& parameters)
	{
		os<<"Orbitals="<<parameters.orbitals<<"\n";
		os<<"hubbardU\n";
		os<<parameters.hubbardU;
		os<<"potentialV\n";
		os<<parameters.potentialV;
		if (parameters.magneticField.rows()>0) {
			os<<"magneticField\n";
			os<<parameters.magneticField;
		}
		if (parameters.spinOrbit.rows()>0) {
			os<<"SpinOrbit\n";
			os<<parameters.spinOrbit;
		}

		if (parameters.jzSymmetry>0) {
			os<<"using jzSymmetry, works only for 3 orbitals \n";
			os<<parameters.jzSymmetry;
		}
		if (parameters.orbDependence>0) {
			os<<"using OrbDependence with Orbitals="<<parameters.orbitals<<" \n";
			os<<parameters.orbDependence;
		}
		os<<"FeAsMode=";
		os<<ParametersModelFeAs<RealType, QnType>::modeString(parameters.feAsMode)<<"\n";
		if (parameters.feAsMode == ParametersModelFeAs<RealType, QnType>::INT_V)
			os<<"CoulombV="<<parameters.coulombV<<"\n";

		if (parameters.potentialT.size()>0) {
			os<<"potentialT\n";
			os<<parameters.potentialT;
		}

		return os;
	}

	SizeType orbitals;
	typename PsimagLite::Vector<RealType>::Type hubbardU;
	typename PsimagLite::Vector<RealType>::Type potentialV;
	typename PsimagLite::Vector<RealType>::Type potentialT;
	IntEnum feAsMode;
	RealType coulombV;
	PsimagLite::Matrix<RealType> magneticField;
	PsimagLite::Matrix<ComplexOrRealType> spinOrbit;
	SizeType jzSymmetry;
	SizeType orbDependence;
};
} // namespace Dmrg

/*@}*/
#endif

