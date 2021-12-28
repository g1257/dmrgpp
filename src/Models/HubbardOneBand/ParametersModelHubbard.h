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
#include "ParametersModelBase.h"
#include "InputNg.h"
#include "InputCheck.h"

namespace Dmrg {
//! Hubbard Model Parameters
template<typename RealType, typename QnType>
struct ParametersModelHubbard : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef PsimagLite::InputNg<InputCheck>::Readable IoInputType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	ParametersModelHubbard(IoInputType& io) : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		hubbardU.resize(nsites, 0.0);
		potentialV.resize(2*nsites, 0.0);
		io.read(hubbardU,"hubbardU");
		io.read(potentialV,"potentialV");
		try {
			anisotropy.resize(nsites, 0.0);
			io.read(anisotropy,"AnisotropyD");
			std::cerr<<"Has AnisotropyD\n";
		} catch (std::exception&) {
			anisotropy.clear();
		}

		try {
			magneticX.resize(nsites, 0.0);
			io.read(magneticX,"MagneticFieldX");
			std::cerr<<"Has MagneticFieldX\n";
		} catch (std::exception&) {
			magneticX.clear();
		}

		VectorStringType tmp = readOldT(io, nsites);
		PsimagLite::String tmp2;

		try {
			io.readline(tmp2, "AddOnSiteHamiltonian=");
			onSiteHadd.resize(nsites);
			stringToVectorOfStrings(onSiteHadd, tmp2);
		} catch (std::exception&) {}

		if (tmp2 != "") {
			if (tmp.size() > 0)
				err("AddOnSiteHamiltonian: You cannot give both legacy and standard entries\n");
		} else {
			onSiteHadd.swap(tmp);
		}
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

	PsimagLite::String killSpaces(PsimagLite::String str)
	{
		PsimagLite::String buffer;
		const SizeType n = str.length();
		for (SizeType i = 0; i < n; ++i)
			if (str[i] != ' ') buffer += str[i];
		return buffer;
	}

	VectorStringType readOldT(IoInputType& io, SizeType nsites)
	{
		typename PsimagLite::Vector<RealType>::Type potentialTlegacy;
		try {
			io.read(potentialTlegacy, "PotentialT");
			std::cerr<<"Has PotentialT\n";
		} catch (std::exception&) {
			return VectorStringType();
		}

		RealType omega = 0;
		try {
			io.readline(omega,"omega=");
		} catch (std::exception&) {}

		RealType phase = 0;
		try {
			io.readline(phase,"phase=");
		} catch (std::exception&) {}

		// c means cosine below
		const PsimagLite::String function = "*(c:+:*:%t:" +
		        ttos(omega) + ":" + ttos(phase) + ")*";
		const PsimagLite::String nup = function + "nup";
		const PsimagLite::String ndown = function + "ndown";

		VectorStringType potentialTv(nsites);
		for (SizeType site = 0; site < nsites; ++site) {
			const RealType val = potentialTlegacy[site];
			const PsimagLite::String plusSignOrNot = ((val < 0) && (site > 0)) ? "+" : "";
			PsimagLite::String expression = ttos(val) + nup + " + ";
			expression += plusSignOrNot + ttos(val) + ndown;
			potentialTv[site] = killSpaces(expression);
		}

		return potentialTv;
	}

	void stringToVectorOfStrings(VectorStringType& vec, PsimagLite::String str)
	{
		if (str[0] == '[') {
			stringToVectorOfStringsCommaMode(vec, str);
		} else {
			stringToVectorOfStringsPlusMode(vec, str);
		}
	}

	void stringToVectorOfStringsPlusMode(VectorStringType& vec, PsimagLite::String str)
	{
		str = killSpaces(str);

		// break on plus
		const SizeType nsites = vec.size();
		VectorStringType tokens;
		PsimagLite::split(tokens, str, "+");
		const SizeType n = tokens.size();
		for (SizeType i = 0; i < n; ++i) {
			std::pair<PsimagLite::String, SizeType> oneSummand = getSiteAndContent(tokens[i]);
			const SizeType site = oneSummand.second;
			if (site >= nsites)
				err("You provided a site " + ttos(site) + " > " + ttos(nsites) + "\n");
			vec[site] = oneSummand.first;
		}
	}

	std::pair<PsimagLite::String, SizeType> getSiteAndContent(PsimagLite::String str)
	{
		const SizeType n = str.length();
		SizeType status = 0; // 0 = closed, 1 = open
		SizeType site = 0;
		bool foundSite = false;
		PsimagLite::String buffer;
		PsimagLite::String content;
		for (SizeType i = 0; i < n; ++i) {
			const char c = str[i];
			if (c == '[') {
				if (status == 1)
					err("Nested brakets found\n");
				status = 1; // open
				continue;
			} else if (c == ']') {
				if (status != 1)
					err("Closing braket without opening one\n");
				site = PsimagLite::atoi(buffer);
				buffer = "";
				foundSite = true;
				status = 0; // closing
				continue;
			}

			if (status == 1) buffer += c;
			else content += c;
		}

		if (!foundSite)
			err("A term for AddOnSiteHamiltonian was given without a site\n");

		return std::pair<PsimagLite::String, SizeType>(content, site);
	}

	void stringToVectorOfStringsCommaMode(VectorStringType& vec, PsimagLite::String str)
	{
		const SizeType last = str.length() - 1;
		if (str.length() < 3 || str[0] != '[' || str[last] != ']')
			err("Expected [...] in comma mode\n");

		str = str.substr(1, str.length() - 2); // remove [ and ]

		// break on ,
		const SizeType nsites = vec.size();
		VectorStringType tokens;
		PsimagLite::split(tokens, str, ",");
		const SizeType n = tokens.size();
		if (n != nsites)
			err("Expected " + ttos(nsites) + " entries but got " + ttos(n) + "\n");
		vec.swap(tokens);
	}

	typename PsimagLite::Vector<RealType>::Type hubbardU;
	typename PsimagLite::Vector<RealType>::Type potentialV;
	typename PsimagLite::Vector<RealType>::Type anisotropy;
	typename PsimagLite::Vector<RealType>::Type magneticX;

	// for time-dependent H:
	VectorStringType onSiteHadd;
};
} // namespace Dmrg

/*@}*/
#endif

