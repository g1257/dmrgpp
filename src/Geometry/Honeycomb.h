/*
Copyright (c) 2009-2013-2018, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.]
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

/*! \file Honeycomb.h
 *
 *  HoneyComb lattice implementation (started March 2018)
 */
#ifndef PSI_GEOMETRY_HONEYCOMB_H
#define PSI_GEOMETRY_HONEYCOMB_H

#include <stdexcept>
#include "GeometryBase.h"

namespace PsimagLite {

template<typename ComplexOrRealType, typename InputType>
class Honeycomb : public GeometryBase<ComplexOrRealType, InputType> {

public:

	enum Dir {DIR_X = 0, DIR_Y = 1, DIR_Z = 2};
	enum GeLe {GREATER_OR_EQUAL, LESS_OR_EQUAL};

	Honeycomb(SizeType linSize, InputType& io)
	    : linSize_(linSize), periodicY_(false), periodicX_(false)
	{
		io.readline(ly_,"HoneycombLy=");
		int x = 0;
		io.readline(x, "IsPeriodicY=");
		periodicY_ = (x > 0);

		try {
			io.readline(x, "IsPeriodicX=");
			periodicX_ = (x > 0);
		} catch (std::exception&) {}
	}

	virtual SizeType maxConnections() const
	{
		throw RuntimeError("Honeycomb::maxConnections() unimplemented\n");
	}

	virtual SizeType dirs() const { return 3; }

	virtual SizeType length(SizeType) const
	{
		return this->unimplemented("length");
	}

	virtual SizeType translate(SizeType,SizeType,SizeType) const
	{
		return this->unimplemented("translate");
	}

	SizeType getVectorSize(SizeType dirId) const
	{
		throw RuntimeError("Honeycomb::getVectorSize() unimplemented\n");
	}

	bool connected(SizeType i1, SizeType i2) const
	{
		return connectedInternal(i1, i2).first;
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1, SizeType i2) const
	{
		std::pair<bool, Dir> bdir = connectedInternal(i1, i2);
		assert(bdir.first);
		return bdir.second;
	}

	bool fringe(SizeType i,SizeType smax,SizeType emin) const
	{
		if (i <= smax)
			return isThereAneighbor(i, emin, linSize_);

		if (i >= emin)
			return isThereAneighbor(i, 0, smax);

		return false;
	}

	// assumes i1 and i2 are connected
	SizeType handle(SizeType i1,SizeType i2) const
	{
		throw RuntimeError("Honeycomb::handle() unimplemented\n");
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax,SizeType emin,SizeType siteNew2) const
	{
		throw RuntimeError("Honeycomb::getSubstituteSite() unimplemented\n");
	}

	String label() const
	{
		return "Honeycomb";
	}

	SizeType findReflection(SizeType) const
	{
		throw RuntimeError("findReflection: unimplemented (sorry)\n");
	}

	template<class Archive>
	void serialize(Archive & ar, const unsigned int)
	{
		ar & boost::serialization::base_object<GeometryBase<ComplexOrRealType, InputType> >(*this);
		ar & linSize_;
		ar & ly_;
	}

	SizeType memResolv(MemResolv&,
	                   SizeType,
	                   String) const
	{
		return 0;
	}

private:

	std::pair<bool, Dir> connectedInternal(SizeType ii1, SizeType ii2) const
	{
		std::pair<bool, Dir> falseDir(false, DIR_X);
		if (ii1 == ii2) return falseDir;

		bool normal = (ii1 < ii2);
		SizeType i1 = (normal) ? ii1 : ii2;
		SizeType i2 = (normal) ? ii2 : ii1;

		SizeType x1 = 0;
		SizeType y1 = 0;
		getCoordinates(x1, y1, i1);

		SizeType x2 = 0;
		SizeType y2 = 0;
		getCoordinates(x2, y2, i2);

		if (isDirectionX(x1, y1, x2, y2))
			return std::pair<bool, Dir>(true, DIR_X);

		if (isDirectionY(x1, y1, x2, y2))
			return std::pair<bool, Dir>(true, DIR_Y);

		if (isDirectionZ(x1, y1, x2, y2))
			return std::pair<bool, Dir>(true, DIR_Z);

		return falseDir;
	}

	void getCoordinates(SizeType& x, SizeType& y, SizeType i) const
	{
		throw RuntimeError("Honeycomb::getCoordinates() unimplemented\n");
	}

	// assumes i1 < i2
	bool isDirectionZ(SizeType x1, SizeType y1, SizeType x2, SizeType y2) const
	{
		if (x1 != x2) return false;
		assert(y1 < y2);
		SizeType d = y2 - y1;
		if ((d == 1) && (y1 & 1)) return true;

		if (!periodicY_) return false;

		// periodic in y direction
		return (d + 1 == ly_);
	}

	// assumes i1 < i2
	bool isDirectionY(SizeType x1, SizeType y1, SizeType x2, SizeType y2) const
	{
		if (x1 == x2)
			return isDirectionYsameCol(y1, y2);

		SizeType dy = abs(y2 - y1);
		if (dy != 1)
			return false;

		assert(x1 < x2);
		SizeType dx = x2 - x1;
		if ((dx == 1) && isMultipleOf3(y1) && (y2 < y1))
			return true;

		if (!periodicX_) return false;

		// periodic in x direction
		SizeType lx = linSize_/ly_;
		if (dx + 1 != lx)
			return false;

		return (isMultipleOf3(y2) && (y1 < y2));
	}

	// assumes i1 < i2
	bool isDirectionX(SizeType x1, SizeType y1, SizeType x2, SizeType y2) const
	{
		if (x1 == x2)
			return isDirectionXsameCol(y1, y2);

		SizeType dy = abs(y2 - y1);
		if (dy != 1)
			return false;

		assert(x1 < x2);
		SizeType dx = x2 - x1;
		if ((dx == 1) && isMultipleOf4(y1) && (y1 < y2))
			return true;

		if (!periodicX_) return false;

		// periodic in x direction
		SizeType lx = linSize_/ly_;
		if (dx + 1 != lx)
			return false;

		return (isMultipleOf4(y2) && (y2 < y1));
	}

	bool isDirectionYsameCol(SizeType y1, SizeType y2) const
	{
		throw RuntimeError("Honeycomb::isDirectionYsameCol() unimplemented\n");
	}

	bool isDirectionXsameCol(SizeType y1, SizeType y2) const
	{
		throw RuntimeError("Honeycomb::isDirectionXsameCol() unimplemented\n");
	}

	bool isThereAneighbor(SizeType ind, SizeType start, SizeType end) const
	{
		for (SizeType i = start; i < end; ++i) {
			if (connectedInternal(ind, i).first) return true;
		}

		return false;
	}

	static bool isMultipleOf3(SizeType y)
	{
		return ((y % 3) == 0);
	}

	static bool isMultipleOf4(SizeType y)
	{
		return ((y % 4) == 0);
	}

	SizeType linSize_;
	SizeType ly_;
	bool periodicY_;
	bool periodicX_;
};
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_H

