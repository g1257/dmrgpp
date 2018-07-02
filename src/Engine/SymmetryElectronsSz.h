/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

/*! \file SymmetryElectronsSz.h
 *
 *
 */
#ifndef DMRG_SYMM_ELECTRONS_SZ_H
#define  DMRG_SYMM_ELECTRONS_SZ_H
#include "Vector.h"
#include "TypeToString.h"
#include "ProgramGlobals.h"
#include "Io/IoSelector.h"

namespace Dmrg {

template<typename RealType>
class SymmetryElectronsSz {

	typedef PsimagLite::IoSelector::Out IoOutType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairSizeType;

public:

	void set(const VectorPairSizeType& j,
	         const VectorSizeType& f,
	         const VectorSizeType& e,
	         const VectorSizeType& o)
	{
		jmValues = j;
		flavors = f;
		electrons = e;
		other = o;
	}

	SizeType electronsMax() const
	{
		return *(std::max_element(electrons.begin(),electrons.end()));
	}

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/SymmetryElectronsSz";
		io.createGroup(label);
		io.write(label + "/electrons", electrons);
		io.write(label + "/other", other);
		io.write(label + "/jmValues", jmValues);
		io.write(label + "/flavors", flavors);
	}

	VectorSizeType electrons;
	VectorSizeType other;
	VectorPairSizeType jmValues;
	VectorSizeType flavors;
}; // struct SymmetryElectronsSz

template<typename IoOutputter, typename RealType>
void write(IoOutputter& io,
          const SymmetryElectronsSz<RealType>& bd,
          typename PsimagLite::EnableIf<
          PsimagLite::IsOutputLike<IoOutputter>::True, int>::Type = 0)
{
	io.write(bd.electronsUp,"bdElectronsUp");
	io.write(bd.electronsDown,"bdElectronsDown");
	io.write(bd.jmValues,"bdJmValues");
	io.write(bd.flavors,"bdFlavors=");
	PsimagLite::String msg("Symmetry=ElectronsSz\n");
	io.print(msg);
}

} // namespace Dmrg

/*@}*/
#endif

