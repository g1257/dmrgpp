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

/*! \file ParametersModelHubbard.h
 *
 *  Contains the parameters for the Hubbard model and function to read them from a file
 *
 */
#ifndef PARAMETERSMODELHUBBARD_H
#define PARAMETERSMODELHUBBARD_H
#include "InputCheck.h"
#include "InputNg.h"
#include "ParametersModelBase.h"

namespace Dmrg
{
//! Hubbard Model Parameters
template <typename RealType, typename QnType>
struct ParametersModelHubbard : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef PsimagLite::InputNg<InputCheck>::Readable IoInputType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	ParametersModelHubbard(IoInputType& io)
	    : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		hubbardU.resize(nsites, 0.0);
		potentialV.resize(2 * nsites, 0.0);
		io.read(hubbardU, "hubbardU");
		io.read(potentialV, "potentialV");
		try {
			anisotropy.resize(nsites, 0.0);
			io.read(anisotropy, "AnisotropyD");
			std::cerr << "Has AnisotropyD\n";
		} catch (std::exception&) {
			anisotropy.clear();
		}

		if (anisotropy.size() > 0 && anisotropy.size() != nsites) {
			throw PsimagLite::RuntimeError("AnisotropyD size must be " + ttos(nsites) +
			    " entries long, if provided.\n");
		}

		try {
			magneticX.resize(nsites, 0.0);
			io.read(magneticX, "MagneticFieldX");
			std::cerr << "Has MagneticFieldX\n";
		} catch (std::exception&) {
			magneticX.clear();
		}

		try {
			io.readline(potentialA, "PotentialA=");
			std::cout<<"PotentialA="<<potentialA<<"\n";
		} catch (std::exception&) {}

		if (magneticX.size() > 0 && magneticX.size() != nsites) {
			throw PsimagLite::RuntimeError("MagneticFieldX size must be " +
			    ttos(nsites) + " entries long, if provided.\n");
		}

		onSiteHaddLegacy = readOldT(io, nsites);
	}

	void write(PsimagLite::String label1,
	    PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelHubbard";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/anisotropy", anisotropy);
		io.write(label + "/magneticX", magneticX);
	}

	VectorStringType readOldT(IoInputType& io, SizeType nsites)
	{
		typename PsimagLite::Vector<RealType>::Type potentialTlegacy;
		try {
			io.read(potentialTlegacy, "PotentialT");
			std::cerr << "Has PotentialT\n";
		} catch (std::exception&) {
			return VectorStringType();
		}

		RealType omega = 0;
		try {
			io.readline(omega, "omega=");
		} catch (std::exception&) {
		}

		RealType phase = 0;
		try {
			io.readline(phase, "phase=");
		} catch (std::exception&) {
		}

		// c means cosine below
		const PsimagLite::String function = "*(c:+:*:%t:" + ttos(omega) + ":" + ttos(phase) + ")*";
		const PsimagLite::String nup = function + "nup";
		const PsimagLite::String ndown = function + "ndown";

		VectorStringType potentialTv(nsites);
		for (SizeType site = 0; site < nsites; ++site) {
			const RealType val = potentialTlegacy[site];
			const PsimagLite::String plusSignOrNot = ((val < 0) && (site > 0)) ? "+" : "";
			PsimagLite::String expression = ttos(val) + nup + " + ";
			expression += plusSignOrNot + ttos(val) + ndown;
			potentialTv[site] = ProgramGlobals::killSpaces(expression);
		}

		return potentialTv;
	}

	typename PsimagLite::Vector<RealType>::Type hubbardU;
	typename PsimagLite::Vector<RealType>::Type potentialV;
	typename PsimagLite::Vector<RealType>::Type anisotropy;
	typename PsimagLite::Vector<RealType>::Type magneticX;
	RealType potentialA = 0.;

	// for time-dependent H:
	VectorStringType onSiteHaddLegacy;
};
} // namespace Dmrg

/*@}*/
#endif
