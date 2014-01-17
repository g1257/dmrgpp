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

/*! \file Chain.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef CHAIN_H
#define CHAIN_H
#include "GeometryBase.h"
#include "String.h"

namespace PsimagLite {

template<typename InputType>
class Chain : public GeometryBase<InputType> {

public:

	enum { DIRECTION_X };

	Chain(SizeType linSize,InputType& io)
		: linSize_(linSize),isPeriodic_(false)
	{
		try {
			int x = 0;
			io.readline(x,"IsPeriodicX=");
			isPeriodic_ = (x > 0) ? true : false;
			std::cout<<"periodic\n";
		} catch (std::exception& e) {}
	}

	virtual SizeType maxConnections() const { return 1; }

	virtual SizeType dirs() const { return 1; }

	SizeType handle(SizeType i,SizeType j) const
	{
		return (i<j) ? i : j;
	}

	SizeType getVectorSize(SizeType dirId) const
	{
		if (dirId!=DIRECTION_X) throw RuntimeError("Chain must have direction 0\n");
		return linSize_-1;
	}

	bool connected(SizeType i1,SizeType i2) const
	{
		if (i1==i2) return false;
		return this->neighbors(i1,i2,isPeriodic_,linSize_-1);
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1,SizeType i2) const
	{
		return DIRECTION_X;
	}

	bool fringe(SizeType i,SizeType smax,SizeType emin) const
	{
		return (i==smax || i==emin);
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax,SizeType emin,SizeType siteNew2) const
	{
		return smax+1;
	}

	String label() const
	{
		return "chain";
	}

	SizeType findReflection(SizeType site) const
	{
		return linSize_ - site -1;
	}

	SizeType length(SizeType i) const
	{
		assert(i==0);
		return linSize_;
	}

	SizeType translate(SizeType site,SizeType dir,SizeType amount) const
	{
		assert(dir==0);

		site+=amount;
		while (site>=linSize_) site -= linSize_;
		return site;
	}

private:

	SizeType linSize_;
	bool isPeriodic_;
}; // class Ladder
} // namespace PsimagLite 

/*@}*/
#endif // LADDER_H

