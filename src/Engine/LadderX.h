// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.4.0]
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file LadderX.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LADDERX_H
#define LADDERX_H

#include "GeometryBase.h"
#include "Ladder.h"

namespace Dmrg {
	
	class LadderX : public GeometryBase {
			typedef Ladder LadderType;
		public:
			enum {DIRECTION_X=LadderType::DIRECTION_X,DIRECTION_Y=LadderType::DIRECTION_Y,DIRECTION_XPY,DIRECTION_XMY};

			LadderX(size_t linSize,size_t leg) : ladder_(linSize,leg),linSize_(linSize),leg_(leg)
			{
			}

			virtual size_t getVectorSize(size_t dirId) const
			{
				switch (dirId) {
				case DIRECTION_XPY:
						return linSize_ - linSize_/leg_;
				case DIRECTION_XMY:
					return linSize_ - linSize_/leg_;
				}
				return ladder_.getVectorSize(dirId);
			}

			virtual bool connected(size_t i1,size_t i2) const
			{
				if (i1==i2) return false;
				if (ladder_.connected(i1,i2)) return true;
				size_t c1 = i1/leg_;
				size_t c2 = i2/leg_;
				size_t r1 = i1%leg_;
				size_t r2 = i2%leg_;
				if (c1==c2) return this->neighbors(r1,r2);
				if (r1==r2) return this->neighbors(c1,c2);
				return (this->neighbors(r1,r2) && this->neighbors(c1,c2));
			}

			// assumes i1 and i2 are connected
			virtual size_t calcDir(size_t i1,size_t i2) const
			{
				if (ladder_.sameColumn(i1,i2)) return DIRECTION_Y;
				if (ladder_.sameRow(i1,i2)) return DIRECTION_X;
				size_t imin = (i1<i2) ? i1 : i2;
				if (imin&1) return DIRECTION_XPY;
				return DIRECTION_XMY;
			}

			virtual bool fringe(size_t i,size_t smax,size_t emin) const
			{
				bool a = (i<emin && i>=smax-1);
				bool b = (i>smax && i<=emin+1);
				if (smax & 1) return (a || b);
				a = (i<emin && i>=smax-2);
				b = (i>smax && i<=emin+2);
				return (a || b);
			}

			// assumes i1 and i2 are connected
			virtual size_t handle(size_t i1,size_t i2) const
			{
				size_t dir = calcDir(i1,i2);
				size_t imin = (i1<i2) ? i1 : i2;
				switch(dir) {
					case DIRECTION_X:
						return imin;
					case DIRECTION_Y:
						return imin-imin/leg_;
					case DIRECTION_XPY: // only checked for leg_=2
						return (imin-1)/leg_;
					case DIRECTION_XMY:// only checked for leg_=2
						return imin/leg_;
				}
				throw std::runtime_error("handle: Unknown direction\n");
			}

			// siteNew2 is fringe in the environment
			virtual size_t getSubstituteSite(size_t smax,size_t emin,size_t siteNew2) const
			{
				return smax+siteNew2-emin+1;
			}

			virtual std::string label() const
			{
				return "ladderx";
			}

		private:

			LadderType ladder_; // owner
			size_t linSize_;
			size_t leg_;
	}; // class LadderBath
} // namespace Dmrg 

/*@}*/
#endif // GEOMETRY_H
