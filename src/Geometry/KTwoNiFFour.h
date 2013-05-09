/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file KTwoNiFFour.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef KTWONIFFOUR_H
#define KTWONIFFOUR_H
#include <stdexcept>
#include <cassert>
#include "String.h"

namespace PsimagLite {

	class KTwoNiFFour  {

		typedef std::pair<int,int> PairType;

		enum {TYPE_O,TYPE_C};
		enum {SUBTYPE_X,SUBTYPE_Y};
		enum {DIR_X,DIR_Y,DIR_XPY,DIR_XMY};

	public:

		struct AdditionalData {
			AdditionalData() : type1(0),type2(0),TYPE_C(KTwoNiFFour::TYPE_C) {}

			size_t type1;
			size_t type2;
			size_t TYPE_C;
		};

		KTwoNiFFour(size_t linSize,int signChange) 
		: linSize_(linSize),signChange_(signChange)
		{
			std::cerr<<"SIGN CHANGE="<<signChange<<"\n";
		}

		size_t getVectorSize(size_t dirId) const
		{
			assert(false);
			throw RuntimeError("getVectorSize: unimplemented\n");
		}

		bool connected(size_t i1,size_t i2) const
		{
			if (i1==i2) return false;
			PairType type1 = findTypeOfSite(i1);
			PairType type2 = findTypeOfSite(i2);
			//! 4 possibilities
			//! c-c
			if (type1.first == type2.first && type1.first == TYPE_C) return false;
			//! o-o
			if (type1.first == type2.first && type1.first == TYPE_O) {
				if (type1.second==type2.second) return false;
				size_t newi1 = (type1.second==SUBTYPE_X) ? i1 : i2;
				size_t newi2 = (type1.second==SUBTYPE_X) ? i2 : i1;
				if (newi1>newi2) {
					assert(newi1>=4);
					if (newi1-2==newi2 || newi1-3==newi2) return true;
					return false;
				}
				if (newi2-1==newi1) return true;
				if (newi2>=2 && newi2-2==newi1) return true;
				return false;
			}
			//! o-c or c-o
			size_t newi1 = (type1.first==TYPE_O) ? i1 : i2;
			size_t newi2 = (type1.first==TYPE_O) ? i2 : i1;
			assert(newi2>=3);
			if (newi2-1==newi1 || newi2-2==newi1 || newi2-3==newi1 || newi2+1==newi1)
				return true;
			return false;
		}

		// assumes i1 and i2 are connected
		size_t calcDir(size_t i1,size_t i2) const
		{
			PairType type1 = findTypeOfSite(i1);
			PairType type2 = findTypeOfSite(i2);
			//! o-o
			if (type1.first == type2.first && type1.first == TYPE_O) {
				assert(type1.second!=type2.second);
				size_t newi1 = (type1.second==SUBTYPE_X) ? i1 : i2;
				size_t newi2 = (type1.second==SUBTYPE_X) ? i2 : i1;
				if (newi1>newi2) {
					assert(newi1>=4);
					size_t distance = newi1-newi2;
					if (distance==2) return DIR_XPY; 
					assert(distance==3);
					return DIR_XMY;
				}
				size_t distance = newi2-newi1;
				if (distance==1)  return DIR_XPY; 
				assert(distance==2);
				return DIR_XMY;
			}
			//! o-c or c-o
			size_t newi1 = (type1.first==TYPE_O) ? i1 : i2;
			type1 = findTypeOfSite(newi1);
			return (type1.second==SUBTYPE_X) ? DIR_X : DIR_Y;
		}

		bool fringe(size_t i,size_t smax,size_t emin) const
		{
			size_t r = smax%4;
			switch (r) {
			case 0:
				return fringe0(i,smax,emin);
			case 1:
				return fringe1(i,smax,emin);
			case 2:
				return fringe2(i,smax,emin);
			case 3:
				return fringe3(i,smax,emin);
			}
			assert(false);
			return false;
		}

		// assumes i1 and i2 are connected
		size_t handle(size_t i1,size_t i2) const
		{
			assert(false);
			return 0;
		}

		// siteNew2 is fringe in the environment
		size_t getSubstituteSite(size_t smax,size_t emin,size_t siteNew2) const
		{
			size_t r = smax%4;
			switch (r) {
			case 0:
				return subs0(smax,emin,siteNew2);
			case 1:
				return subs1(smax,emin,siteNew2);
			case 2:
				return subs2(smax,emin,siteNew2);
			case 3:
				return subs3(smax,emin,siteNew2);
			}
			assert(false);
			return 0;
		}

		String label() const
		{
			return "KTwoNiFFour";
		}

		size_t maxConnections() const
		{
			return 6;
		}
		
		size_t findReflection(size_t site) const
		{
			throw RuntimeError("findReflection: unimplemented (sorry)\n");
		}

		void fillAdditionalData(AdditionalData& additionalData,size_t ind,size_t jnd) const
		{
			additionalData.type1 = findTypeOfSite(ind).first;
			additionalData.type2 = findTypeOfSite(jnd).first;
			additionalData.TYPE_C = TYPE_C;
		}

		size_t matrixRank() const
		{
			size_t sites = linSize_;
			size_t no = 0;
			size_t nc = 0;
			for (size_t i=0;i<sites;i++) {
				size_t type1 = findTypeOfSite(i).first;
				if (type1==TYPE_C) nc++;
				else no++;
			}
			return 2*no+nc;
		}
		
		int index(size_t i,size_t orb) const
		{
			size_t type1 = findTypeOfSite(i).first;
			if (type1==TYPE_C && orb>0) return -1;
			if (type1==TYPE_C || orb==0) return i;
			size_t sites = linSize_;
			size_t tmp = (i+1)/4;
			assert(sites+i>=tmp);
			return sites+i-tmp;
		}
		
		int signChange(size_t i1,size_t i2) const
		{
			// change the sign of Cu-O
			size_t newi1 = std::min(i1,i2);
			size_t newi2 = std::max(i1,i2);
			PairType type1 = findTypeOfSite(newi1);
			PairType type2 = findTypeOfSite(newi2);
			int sign1 = 1;
			if (type1.first!=type2.first) {
				
				int diff =  newi2-newi1;
				assert(diff==1 || diff==2 || diff==3);
				if (diff<2) sign1 = -1;
			}
			
			if (isInverted(i1) || isInverted(i2)) return signChange_*sign1;
			return sign1;
		}

	private:

		//! Given as a pair, first number is the type,
		//! If first number == TYPE_C then second number is bogus
		//! If first number == TYPE_O then second number is the subtype
		PairType findTypeOfSite(size_t site) const
		{
			size_t sitePlusOne = site + 1;
			size_t r = sitePlusOne%4;
			if (r==0) return PairType(TYPE_C,0);
			
			if (r==1) return PairType(TYPE_O,SUBTYPE_X);
			return PairType(TYPE_O,SUBTYPE_Y);
		}
		
		bool isInverted(size_t i) const
		{
			size_t j = i+4;
			return ((j%8)==0);
		}
		
		bool fringe0(size_t i,size_t smax,size_t emin) const
		{
			return (i==smax || i==emin);
		}
		
		size_t subs0(size_t smax,size_t emin,size_t i) const
		{
			return smax+3;
		}

		bool fringe1(size_t i,size_t smax,size_t emin) const
		{
			bool b1 = (i==smax || i==smax-1);
			bool b2 = (i==emin || i==emin+1 || i==emin+2);
			return (b1 || b2);
		}
		
		size_t subs1(size_t smax,size_t emin,size_t i) const
		{
			if (i==emin) return smax+3;
			if (i==emin+1) return smax+2;
			if (i==emin+2) return smax+4;
			assert(false);
			return 0;
		}
		
		bool fringe2(size_t i,size_t smax,size_t emin) const
		{
			return (i==smax-2 || i==emin+2);
		}
		
		size_t subs2(size_t smax,size_t emin,size_t i) const
		{
			assert(i==emin+2);
			return smax+1;
		}

		bool fringe3(size_t i,size_t smax,size_t emin) const
		{
			return (i==emin || i==smax || i==smax-1 || i==smax-2);
		}
		
		size_t subs3(size_t smax,size_t emin,size_t i) const
		{
			assert(i==emin);
			return smax+1;
		}

		size_t linSize_;
		int signChange_;
	}; // class KTwoNiFFour
} // namespace PsimagLite 

/*@}*/
#endif // KTWONIFFOUR_H
