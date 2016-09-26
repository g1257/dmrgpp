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

/*! \file GeometryBase.h
 *
 *  Well, I need to read that chapter in
 *  Alexandrescu's "Modern C++ design" again to have
 *  a decent factory here, but this will have to do for now
 *
 */
#ifndef GEOMETRY_BASE_H
#define GEOMETRY_BASE_H

#include "InputNg.h"
#include "MemResolv.h"

namespace PsimagLite {

template<typename ComplexOrRealType, typename InputType>
class GeometryBase {

	typedef std::pair<SizeType,SizeType> PairType;
	typedef Matrix<ComplexOrRealType> MatrixType;

	struct AdditionalData {
		AdditionalData() : type1(0),type2(0),TYPE_C(GeometryBase::TYPE_C) {}

		SizeType type1;
		SizeType type2;
		SizeType TYPE_C;
	};

public:

	enum {TYPE_O,TYPE_C};

	typedef AdditionalData AdditionalDataType;

	virtual ~GeometryBase()
	{}

	template<class Archive>
	void serialize(Archive&, const unsigned int)
	{}

	virtual SizeType memResolv(MemResolv& mres,
	                           SizeType x,
	                           String msg) const = 0;

	virtual SizeType dirs() const = 0;

	virtual SizeType handle(SizeType i,SizeType j) const = 0;

	virtual SizeType getVectorSize(SizeType dirId) const = 0;

	virtual bool connected(SizeType i1,SizeType i2) const = 0;

	virtual SizeType calcDir(SizeType i1,SizeType i2) const = 0;

	virtual bool fringe(SizeType i,SizeType smax,SizeType emin) const = 0;

	virtual SizeType getSubstituteSite(SizeType smax,
	                                   SizeType emin,
	                                   SizeType siteNew2) const = 0;

	virtual String label() const = 0;

	virtual SizeType length(SizeType i) const = 0;

	virtual SizeType translate(SizeType site,SizeType dir,SizeType amount) const = 0;

	virtual SizeType maxConnections() const = 0;

	virtual SizeType findReflection(SizeType site) const = 0;

	virtual void set(MatrixType&, SizeType&) const
	{
		throw RuntimeError("GeometryBase::set() unimplemented for derived class\n");
	}

	virtual void fillAdditionalData(AdditionalDataType&,
	                                SizeType,
	                                SizeType) const
	{}

	virtual int index(SizeType i1,SizeType edof1,SizeType edofTotal) const
	{
		return edof1+i1*edofTotal;
	}

	virtual SizeType matrixRank(SizeType linSize,SizeType maxEdof) const
	{
		return linSize*maxEdof;
	}

	virtual int signChange(SizeType, SizeType) const
	{
		return 1;
	}

protected:

	SizeType unimplemented(const String& str) const
	{
		String str2 = "unimplemented " + str + "\n";
		throw RuntimeError(str2);
	}

	bool neighbors(SizeType i1,SizeType i2,bool periodic = false,SizeType period = 1) const
	{
		SizeType imin = (i1<i2) ? i1 : i2;
		SizeType imax = (i1>i2) ? i1 : i2;
		bool b = (imax-imin==1);
		if (!periodic) return b;
		bool b2 = (imax-imin == period);
		return (b || b2);
	}
}; // class GeometryBase
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_BASE_H

