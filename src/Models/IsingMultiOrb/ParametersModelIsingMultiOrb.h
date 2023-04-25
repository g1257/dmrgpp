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

/*! \file ParametersModelIsingMultiOrb.h
 *
 *  Contains the parameters for the Ising model with many orbitals
 *
 */
#ifndef PARAMETERSMODELISINGMULTIORB_H
#define PARAMETERSMODELISINGMULTIORB_H
#include "Vector.h"
#include "../../Engine/ParametersModelBase.h"

namespace Dmrg {
//! IsingMultiOrb Model Parameters
template<typename RealType, typename QnType>
struct ParametersModelIsingMultiOrb : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Matrix<RealType> MatrixRealType;
	// no connectors here, connectors are handled by the geometry
	template<typename IoInputType>
	ParametersModelIsingMultiOrb(IoInputType& io)
	    : BaseType(io, false),
	      magneticFieldX(0, 0),
	      magneticFieldZ(0, 0),
	      onsitelinksSzSz(0, 0),
	      TimeSchedule(0, 0),
	      hasTimeSchedule_(false),
	      ta(0.0)
	{
		PsimagLite::String model;
		io.readline(model, "Model=");
		io.readline(orbitals,"Orbitals=");
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");

		try {
			io.read(magneticFieldX, "MagneticFieldX");
		} catch (std::exception&) {}

		if (magneticFieldX.rows()>0 && magneticFieldX.rows()!=orbitals)
			throw PsimagLite::RuntimeError("MagneticFieldX: Expecting rows == orbitals\n");
		if (magneticFieldX.cols()>0 && magneticFieldX.cols()!=nsites)
			throw PsimagLite::RuntimeError("MagneticFieldX: Expecting columns == sites\n");

		try {
			io.read(magneticFieldZ, "MagneticFieldZ");
		} catch (std::exception&) {}

		if (magneticFieldZ.rows()>0 && magneticFieldZ.rows()!=orbitals)
			throw PsimagLite::RuntimeError("MagneticFieldZ: Expecting rows == orbitals\n");
		if (magneticFieldZ.cols()>0 && magneticFieldZ.cols()!=nsites)
			throw PsimagLite::RuntimeError("MagneticFieldZ: Expecting columns == sites\n");

		// throw if supplying MagneticField label
		bool invalidLabel = false;
		try {
			VectorRealType tmpVector;
			io.read(tmpVector, "MagneticField=");
			invalidLabel = true;
		} catch (std::exception&) {}

		if (invalidLabel) {
			throw PsimagLite::RuntimeError("MagneticField label is no longer supported.\n" +
			                               PsimagLite::String("Please use MagneticField[XZ] instead\n"));
		}

		// throw if supplying MagneticFieldDirection label
		try {
			PsimagLite::String tmpStr;
			io.readline(tmpStr, "MagneticFieldDirection=");
			invalidLabel = true;
		} catch (std::exception&) {}

		if (invalidLabel) {
			throw PsimagLite::RuntimeError("MagneticFieldDirection label is no longer supported.\n" +
			                               PsimagLite::String("Please use MagneticField[XZ] instead\n"));
		}

		try {
			io.read(onsitelinksSzSz,"OnSiteLinksSzSz");
		} catch (std::exception&) {}


		const SizeType orbs1 = combinations(orbitals, 2);
		if (onsitelinksSzSz.cols()>0 && onsitelinksSzSz.cols()!=nsites)
			throw PsimagLite::RuntimeError("OnSiteLinksSzSz: Expecting cols == sites\n");
		if (onsitelinksSzSz.rows()>0 && onsitelinksSzSz.rows()!=orbs1)
			throw PsimagLite::RuntimeError("OnSiteLinksSzSz: Expecting rows == combinations(orbitals,2)\n");

		bool hasTimeSchedule = false;
		try {
			io.read(TimeSchedule,"TimeSchedule");
			hasTimeSchedule = true;
		} catch (std::exception&) {}

		if (hasTimeSchedule ==true) { // Check Input to see if everything is correct

			if (TimeSchedule.cols()>0 && TimeSchedule.cols()!=3)
				throw PsimagLite::RuntimeError("TimeSchedule: Expecting cols == 3, [s,Gamma(s),J(s)]\n");

			bool hasta = false;
			try {
				io.readline(ta,"ta=");
				hasta = true;
			} catch (std::exception&) {}

			RealType tau =0.0;
			bool hastau = false;
			try {
				io.readline(tau, "TSPTau=");
				hastau = true;
			} catch (std::exception&) {}

			if (hastau && !hasta)
				throw PsimagLite::RuntimeError("TimeSchedule: TSPTau is set, "
			                                   "so you must have ta=something>0 in the input!\n");
			if (hasta && hastau)
				if(ta<0 || fabs(ta)<1e-5)
					throw PsimagLite::RuntimeError("TimeSchedule: ta is negative or too small!\n");

			hasTimeSchedule_ = true;
		}

	}

	static SizeType combinations(SizeType n, SizeType r)
	{
		if (r > n)
			return 0;
		if (r == 0 || r == n)
			return 1;
		return combinations(n - 1, r - 1) + combinations(n - 1, r);
	}

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelIsingMultiOrb";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/orbitals", orbitals);
		magneticFieldX.write(label + "/magneticFieldX", io);
		magneticFieldZ.write(label + "/magneticFieldZ", io);
		onsitelinksSzSz.write(label + "/onsitelinksSzSz", io);
	}

	static void checkMagneticField(unsigned char c, SizeType s, SizeType n, SizeType ss, SizeType orbs)
	{
		if (s == 0 || s == n) return;
		if (ss == 0 || ss == orbs) return;

		PsimagLite::String msg("ModelIsingMultiOrb: If provided, ");
		msg += " MagneticField" + ttos(c) + " must be a matrix of (rows) " +
		        ttos(orbs) + "orbitals times (cols) " + ttos(n) +" entries.\n";
		err(msg);
	}

	static void checkOnSiteLinksSzSz(SizeType s, SizeType n, SizeType ss, SizeType orbs1)
	{
		if (s == 0 || s == n) return;
		if (ss == 0 || ss == orbs1) return;
		PsimagLite::String msg("ModelIsingMultiOrb: If provided, ");
		msg += " OnsiteLinksSzSz must be a matrix of (rows) " +
		        ttos(orbs1) + "terms times (cols) " + ttos(n) +" entries.\n";
		err(msg);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream &os,
	                                const ParametersModelIsingMultiOrb& parameters)
	{
		if (parameters.magneticFieldX.cols()>0)
			os<<"MagneticFieldX="<<parameters.magneticFieldX<<"\n";
		if (parameters.magneticFieldZ.cols()>0)
			os<<"MagneticFieldZ="<<parameters.magneticFieldZ<<"\n";
		if (parameters.onsitelinksSzSz.cols()>0)
			os<<"OnSiteLinksSzSz="<<parameters.onsitelinksSzSz<<"\n";
		
		os<<parameters.targetQuantum;
		return os;
	}

	SizeType orbitals;
	MatrixRealType magneticFieldX;
	MatrixRealType magneticFieldZ;
	MatrixRealType onsitelinksSzSz;
	MatrixRealType TimeSchedule;
	bool hasTimeSchedule_;
	RealType ta;
};
} // namespace Dmrg

/*@}*/
#endif

