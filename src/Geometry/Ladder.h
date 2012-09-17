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

#include "GeometryUtils.h"

namespace PsimagLite {
	
	class Ladder  {
	public:
		enum {DIRECTION_X,DIRECTION_Y};

		Ladder(size_t linSize,size_t leg,bool isPeriodicY)
		: linSize_(linSize),leg_(leg),isPeriodicY_(isPeriodicY)
		{
			if (leg & 1) throw std::runtime_error("Ladder: leg must be even\n");
			if (leg == 2)  isPeriodicY_ = false;
			//if (leg>2) std::cerr<<"isPeriodicY="<<isPeriodicY_<<"\n";
			if (linSize % leg !=0) throw std::runtime_error("Ladder: leg must divide number of sites\n");
		}

		size_t getVectorSize(size_t dirId) const
		{
			if (dirId==DIRECTION_X) return linSize_-leg_;
			else if (dirId==DIRECTION_Y) return (isPeriodicY_) ? linSize_ : linSize_ - linSize_/leg_;

			throw std::runtime_error("Unknown direction\n");
		}

		bool connected(size_t i1,size_t i2) const
		{
			if (i1==i2) return false;
			size_t c1 = i1/leg_;
			size_t c2 = i2/leg_;
			size_t r1 = i1%leg_;
			size_t r2 = i2%leg_;
			if (c1==c2) return GeometryUtils::neighbors(r1,r2,isPeriodicY_,leg_-1);
			if (r1==r2) return GeometryUtils::neighbors(c1,c2);
			return false;
		}

		size_t calcDir(size_t i1,size_t i2) const
		{
			assert(connected(i1,i2));
			if (sameColumn(i1,i2)) return DIRECTION_Y;
			return DIRECTION_X;
		}

		bool fringe(size_t i,size_t smax,size_t emin) const
		{
			size_t c = smax % leg_;
			size_t r = 2 + 2*c;
			if (c>=leg_/2) r = r - leg_;

			if (smax+1 == emin) r = leg_; // finite loops

			bool a = (i<emin && i>=smax - r + 1);
			bool b = (i>smax && i<=emin + r - 1);
			return (a || b);
		}

		// siteEnv is fringe in the environment
		size_t getSubstituteSite(size_t smax,size_t emin,size_t siteEnv) const
		{
			if (smax+1 == emin) return siteEnv; // finite loops

			size_t c = smax % leg_;
			size_t s = int(emin/leg_) - int(smax/leg_);
			assert(s>0);
			if (c>=leg_/2) s--;
			return  siteEnv - s*leg_;
		}

		size_t handle(size_t i1,size_t i2) const
		{
			size_t dir = calcDir(i1,i2);
			size_t imin = (i1<i2) ? i1 : i2;
			size_t y = imin/leg_;
			switch(dir) {
			case DIRECTION_X:
				return imin;
			case DIRECTION_Y:
				if (!isPeriodicY_) return imin-imin/leg_;
				if (imin ==0 || imin % leg_ == 0) imin = (i1>i2) ? i1 : i2;
				return imin-imin/leg_ + y;
			}
			throw std::runtime_error("hanlde: Unknown direction\n");
		}

		bool sameColumn(size_t i1,size_t i2) const
		{
			size_t c1 = i1/leg_;
			size_t c2 = i2/leg_;
			if (c1 == c2) return true;
			return false;
		}

		bool sameRow(size_t i1,size_t i2) const
		{
			size_t c1 = i1%leg_;
			size_t c2 = i2%leg_;
			if (c1 == c2) return true;
			return false;
		}

		std::string label() const
		{
			return "ladder";
		}
		
		size_t findReflection(size_t site) const
		{
			throw std::runtime_error("findReflection: unimplemented (sorry)\n");
		}

	private:

		size_t linSize_;
		size_t leg_;
		bool isPeriodicY_;
	}; // class Ladder
} // namespace PsimagLite 

/*@}*/
#endif // LADDER_H

