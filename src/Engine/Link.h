/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file Link.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LINK_H
#define LINK_H

#include "ProgramGlobals.h"

namespace Dmrg {

template<typename FieldType>
struct Link {

	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef std::pair<char,char> PairCharType;
	typedef typename PsimagLite::Real<FieldType>::Type RealType;

	Link(SizeType i,
	     SizeType j,
	     ProgramGlobals::ConnectionEnum type1,
	     const FieldType& value1,
	     SizeType dofs1,
	     ProgramGlobals::FermionOrBosonEnum fOb,
	     const PairSizeType& ops1,
	     const PairCharType& mods1,
	     SizeType aM,RealType aF,
	     SizeType cat)
	    : site1(i),site2(j),
	      type(type1),value(value1),dofs(dofs1),
	      fermionOrBoson(fOb),ops(ops1),mods(mods1),
	      angularMomentum(aM),angularFactor(aF),
	      category(cat)
	{}

	SizeType site1,site2;
	ProgramGlobals::ConnectionEnum type;
	FieldType value;
	SizeType dofs;
	ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
	std::pair<SizeType,SizeType> ops; // operator indices
	PairCharType mods; // operator modifiers
	SizeType angularMomentum;
	RealType angularFactor;
	SizeType category;

}; // struct Link
} // namespace Dmrg 

/*@}*/
#endif // LINK_H
