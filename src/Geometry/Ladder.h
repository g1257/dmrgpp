/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.4.0]
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

/*! \file Ladder.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LADDER_H
#define LADDER_H

#include "GeometryBase.h"
#include "String.h"

namespace PsimagLite {

template<typename InputType>
class Ladder : public GeometryBase<InputType> {
public:
	enum {DIRECTION_X,DIRECTION_Y};

	Ladder(SizeType linSize,InputType& io)
	    : linSize_(linSize)
	{

		io.readline(leg_,"LadderLeg=");
		if (leg_!=2) {
			std::cerr<<"WARNING: LadderLeg!=2 is experimental!\n";
		}
		try {
			int x = 0;
			io.readline(x,"PeriodicY=");
			isPeriodicY_ = (x > 0) ? true : false;
			if (leg_==2) throw RuntimeError("LadderLeg==2 cannot have PeriodicY set\n");
			std::cerr<<"INFO: PeriodicY="<<isPeriodicY_<<"\n";
		} catch (std::exception& e) {
			if (leg_>2) throw RuntimeError("LadderLeg>2 must have PeriodicY= line\n");
		}


		if (leg_ & 1)
			throw RuntimeError("Ladder: leg must be even\n");

		if (leg_ == 2)
			isPeriodicY_ = false;

		if (linSize % leg_ !=0)
			throw RuntimeError("Ladder: leg must divide number of sites\n");
	}


	virtual SizeType maxConnections() const { return 4; }

	virtual SizeType dirs() const { return 2; }

	SizeType getVectorSize(SizeType dirId) const
	{
		if (dirId==DIRECTION_X)
			return linSize_-leg_;
		else if (dirId==DIRECTION_Y)
			return (isPeriodicY_) ? linSize_ : linSize_ - linSize_/leg_;

		throw RuntimeError("Unknown direction\n");
	}

	bool connected(SizeType i1,SizeType i2) const
	{
		if (i1==i2) return false;
		SizeType c1 = i1/leg_;
		SizeType c2 = i2/leg_;
		SizeType r1 = i1%leg_;
		SizeType r2 = i2%leg_;
		if (c1==c2) return this->neighbors(r1,r2,isPeriodicY_,leg_-1);
		if (r1==r2) return this->neighbors(c1,c2);
		return false;
	}

	SizeType calcDir(SizeType i1,SizeType i2) const
	{
		assert(connected(i1,i2));
		if (sameColumn(i1,i2)) return DIRECTION_Y;
		return DIRECTION_X;
	}

	bool fringe(SizeType i,SizeType smax,SizeType emin) const
	{
		SizeType c = smax % leg_;
		SizeType r = 2 + 2*c;
		if (c>=leg_/2) r = r - leg_;

		if (smax+1 == emin) r = leg_; // finite loops

		bool a = (i<emin && i>=smax - r + 1);
		bool b = (i>smax && i<=emin + r - 1);
		return (a || b);
	}

	// siteEnv is fringe in the environment
	SizeType getSubstituteSite(SizeType smax,SizeType emin,SizeType siteEnv) const
	{
		if (smax+1 == emin) return siteEnv; // finite loops

		SizeType c = smax % leg_;
		SizeType s = int(emin/leg_) - int(smax/leg_);
		assert(s>0);
		if (c>=leg_/2) s--;
		return  siteEnv - s*leg_;
	}

	SizeType handle(SizeType i1,SizeType i2) const
	{
		SizeType dir = calcDir(i1,i2);
		SizeType imin = (i1<i2) ? i1 : i2;
		SizeType y = imin/leg_;
		switch(dir) {
		case DIRECTION_X:
			return imin;
		case DIRECTION_Y:
			if (!isPeriodicY_) return imin-imin/leg_;
			if (imin ==0 || imin % leg_ == 0) imin = (i1>i2) ? i1 : i2;
			return imin-imin/leg_ + y;
		}
		throw RuntimeError("hanlde: Unknown direction\n");
	}

	bool sameColumn(SizeType i1,SizeType i2) const
	{
		SizeType c1 = i1/leg_;
		SizeType c2 = i2/leg_;
		if (c1 == c2) return true;
		return false;
	}

	bool sameRow(SizeType i1,SizeType i2) const
	{
		SizeType c1 = i1%leg_;
		SizeType c2 = i2%leg_;
		if (c1 == c2) return true;
		return false;
	}

	String label() const
	{
		return "ladder";
	}

	SizeType length(SizeType i) const
	{
		assert(i<2);
		return (i==1) ? leg_ : SizeType(linSize_/leg_);
	}

	SizeType translate(SizeType site,SizeType dir,SizeType amount) const
	{
		assert(dir<2);
		SizeType x = SizeType(site/leg_);
		SizeType y = site % leg_;
		SizeType lx = SizeType(linSize_/leg_);
		if (dir==DIRECTION_X) x = translateInternal(x,lx,amount);
		else y = translateInternal(y,leg_,amount);
		SizeType ind = y + x*leg_;
		assert(ind<linSize_);
		return ind;
	}

	SizeType findReflection(SizeType site) const
	{
		SizeType x = SizeType(site/leg_);
		SizeType y = site % leg_;
		SizeType lx = SizeType(linSize_/leg_);
		SizeType ind = y + (lx-x-1)*leg_;
		assert(ind<linSize_);
		return ind;
	}

private:

	SizeType translateInternal(SizeType c,SizeType l,SizeType amount) const
	{
		c += amount;
		while (c>=l) c-=l;
		return c;
	}

	SizeType linSize_;
	SizeType leg_;
	bool isPeriodicY_;
}; // class Ladder
} // namespace PsimagLite 

/*@}*/
#endif // LADDER_H

