/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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
/** \ingroup PsimagLite */
/*@{*/

/*! \file ChainEx.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef CHAIN_EX_H
#define CHAIN_EX_H
#include "GeometryBase.h"
#include "String.h"

namespace PsimagLite {

template<typename InputType>
class ChainEx : public GeometryBase<InputType> {

public:

	enum { DIRECTION_X, DIRECTION_NNN };

	ChainEx(SizeType linSize,InputType& io) : linSize_(linSize)
	{}

	virtual SizeType maxConnections() const { return 2; }

	virtual SizeType dirs() const { return 2; }

	SizeType handle(SizeType i,SizeType j) const
	{
		this->unimplemented("handle");
		return (i<j) ? i : j;
	}

	SizeType getVectorSize(SizeType dirId) const
	{
		this->unimplemented("getVectorSize");
		return linSize_-1;
	}

	bool connected(SizeType i1,SizeType i2) const
	{
		if (i1==i2) return false;
		bool b1 = this->neighbors(i1,i2);
		bool b2 = (fabs(i1-i2) == 2 && (i1 & 1) == 0);
		return (b1 | b2);
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1,SizeType i2) const
	{
		bool b1 = this->neighbors(i1,i2);
		bool b2 = (fabs(i1-i2) == 2 && (i1 & 1) == 0);
		assert(b1 ^ b2);

		if (b1) return DIRECTION_X;
		return DIRECTION_NNN;
	}

	bool fringe(SizeType i,SizeType smax,SizeType emin) const
	{
		bool b1 = (i==smax || i==emin);
		if (smax & 1) {
			return (b1 || (i == smax - 1));
		}

		return (b1 || (i == emin + 1));
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax,SizeType emin,SizeType siteNew2) const
	{
		return smax+1;
	}

	String label() const
	{
		return "chainEx";
	}

	SizeType findReflection(SizeType site) const
	{
		return linSize_ - site -1;
	}

	SizeType length(SizeType i) const
	{
		this->unimplemented("length");
		return linSize_;
	}

	SizeType translate(SizeType site,SizeType dir,SizeType amount) const
	{
		this->unimplemented("getVectorSize");

		site+=amount;
		while (site>=linSize_) site -= linSize_;
		return site;
	}

private:

	SizeType linSize_;
}; // class ChainEx
} // namespace PsimagLite 

/*@}*/
#endif // CHAIN_EX_H

