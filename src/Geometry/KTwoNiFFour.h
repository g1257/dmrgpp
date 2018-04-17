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
#include "GeometryBase.h"

namespace PsimagLite {

template<typename ComplexOrRealType, typename InputType>
class KTwoNiFFour : public GeometryBase<ComplexOrRealType, InputType> {

	typedef std::pair<int,int> PairType;
	typedef GeometryBase<ComplexOrRealType, InputType> GeometryBaseType;
	typedef typename GeometryBaseType::AdditionalDataType AdditionalDataType;

	enum {TYPE_O = GeometryBaseType::TYPE_O, TYPE_C = GeometryBaseType::TYPE_C};

	enum {SUBTYPE_X,SUBTYPE_Y};

	enum {DIR_X,DIR_Y,DIR_XPY,DIR_XMY};

public:

	KTwoNiFFour() {}

	KTwoNiFFour(SizeType linSize,InputType& io)
	    : linSize_(linSize)
	{
		io.readline(signChange_,"SignChange=");
		std::cout<<"KTwoNiFFour: SIGN CHANGE="<<signChange_<<"\n";
	}

	template<class Archive>
	void write(Archive & ar, const unsigned int)
	{
		ar & boost::serialization::base_object<GeometryBase<ComplexOrRealType, InputType> >(*this);
		ar & linSize_;
		ar & signChange_;
	}

	SizeType memResolv(MemResolv& mres,
	                   SizeType,
	                   String msg) const
	{
		String str = msg;
		str += "KTwoNiFFour";
		const char* start = (const char *)this;
		const char* end = (const char*)&linSize_;
		SizeType total = end - start;
		mres.push(MemResolv::MEMORY_TEXTPTR, total, start,str+" vptr");

		start = end;
		end = (const char*)&signChange_;
		total += mres.memResolv(&linSize_,end-start,str + " linSize");

		mres.memResolv(&signChange_,sizeof(*this) - total,str + " signChange");

		return sizeof(*this);
	}

	SizeType getVectorSize(SizeType) const
	{
		assert(false);
		throw RuntimeError("getVectorSize: unimplemented\n");
	}

	virtual SizeType dirs() const { return 4; }

	virtual SizeType length(SizeType) const
	{
		return this->unimplemented("length");
	}

	virtual SizeType translate(SizeType,SizeType,SizeType) const
	{
		return this->unimplemented("translate");
	}

	bool connected(SizeType i1,SizeType i2) const
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
			SizeType newi1 = (type1.second==SUBTYPE_X) ? i1 : i2;
			SizeType newi2 = (type1.second==SUBTYPE_X) ? i2 : i1;
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
		SizeType newi1 = (type1.first==TYPE_O) ? i1 : i2;
		SizeType newi2 = (type1.first==TYPE_O) ? i2 : i1;
		assert(newi2>=3);
		if (newi2-1==newi1 || newi2-2==newi1 || newi2-3==newi1 || newi2+1==newi1)
			return true;
		return false;
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1,SizeType i2) const
	{
		PairType type1 = findTypeOfSite(i1);
		PairType type2 = findTypeOfSite(i2);
		//! o-o
		if (type1.first == type2.first && type1.first == TYPE_O) {
			assert(type1.second!=type2.second);
			SizeType newi1 = (type1.second==SUBTYPE_X) ? i1 : i2;
			SizeType newi2 = (type1.second==SUBTYPE_X) ? i2 : i1;
			if (newi1>newi2) {
				assert(newi1>=4);
				SizeType distance = newi1-newi2;
				if (distance==2) return DIR_XPY;
				assert(distance==3);
				return DIR_XMY;
			}
			SizeType distance = newi2-newi1;
			if (distance==1)  return DIR_XPY;
			assert(distance==2);
			return DIR_XMY;
		}
		//! o-c or c-o
		SizeType newi1 = (type1.first==TYPE_O) ? i1 : i2;
		type1 = findTypeOfSite(newi1);
		return (type1.second==SUBTYPE_X) ? DIR_X : DIR_Y;
	}

	bool fringe(SizeType i,SizeType smax,SizeType emin) const
	{
		SizeType r = smax%4;
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
	SizeType handle(SizeType,SizeType) const
	{
		assert(false);
		return 0;
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax,SizeType emin,SizeType siteNew2) const
	{
		SizeType r = smax%4;
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

	SizeType maxConnections() const
	{
		return 6;
	}

	SizeType findReflection(SizeType) const
	{
		throw RuntimeError("findReflection: unimplemented (sorry)\n");
	}

	void fillAdditionalData(AdditionalDataType& additionalData,
	                        SizeType ind,
	                        SizeType jnd) const
	{
		additionalData.type1 = findTypeOfSite(ind).first;
		additionalData.type2 = findTypeOfSite(jnd).first;
		additionalData.TYPE_C = TYPE_C;
	}

	SizeType matrixRank(SizeType,SizeType) const
	{
		SizeType sites = linSize_;
		SizeType no = 0;
		SizeType nc = 0;
		for (SizeType i=0;i<sites;i++) {
			SizeType type1 = findTypeOfSite(i).first;
			if (type1==TYPE_C) nc++;
			else no++;
		}
		return 2*no+nc;
	}

	int index(SizeType i,SizeType orb,SizeType) const
	{
		SizeType type1 = findTypeOfSite(i).first;
		if (type1==TYPE_C && orb>0) return -1;
		if (type1==TYPE_C || orb==0) return i;
		SizeType sites = linSize_;
		SizeType tmp = (i+1)/4;
		assert(sites+i>=tmp);
		return sites+i-tmp;
	}

	int signChange(SizeType i1,SizeType i2) const
	{
		// change the sign of Cu-O
		SizeType newi1 = std::min(i1,i2);
		SizeType newi2 = std::max(i1,i2);
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
	PairType findTypeOfSite(SizeType site) const
	{
		SizeType sitePlusOne = site + 1;
		SizeType r = sitePlusOne%4;
		if (r==0) return PairType(TYPE_C,0);

		if (r==1) return PairType(TYPE_O,SUBTYPE_X);
		return PairType(TYPE_O,SUBTYPE_Y);
	}

	bool isInverted(SizeType i) const
	{
		SizeType j = i+4;
		return ((j%8)==0);
	}

	bool fringe0(SizeType i,SizeType smax,SizeType emin) const
	{
		return (i==smax || i==emin);
	}

	SizeType subs0(SizeType smax,SizeType,SizeType) const
	{
		return smax+3;
	}

	bool fringe1(SizeType i,SizeType smax,SizeType emin) const
	{
		bool b1 = (i==smax || i==smax-1);
		bool b2 = (i==emin || i==emin+1 || i==emin+2);
		return (b1 || b2);
	}

	SizeType subs1(SizeType smax,SizeType emin,SizeType i) const
	{
		if (i==emin) return smax+3;
		if (i==emin+1) return smax+2;
		if (i==emin+2) return smax+4;
		assert(false);
		return 0;
	}

	bool fringe2(SizeType i,SizeType smax,SizeType emin) const
	{
		return (i==smax-2 || i==emin+2);
	}

	SizeType subs2(SizeType smax,SizeType emin,SizeType i) const
	{
		assert(i==emin+2);
		return smax+1;
	}

	bool fringe3(SizeType i,SizeType smax,SizeType emin) const
	{
		return (i==emin || i==smax || i==smax-1 || i==smax-2);
	}

	SizeType subs3(SizeType smax,SizeType emin,SizeType i) const
	{
		assert(i==emin);
		return smax+1;
	}

	SizeType linSize_;
	int signChange_;
}; // class KTwoNiFFour
} // namespace PsimagLite

/*@}*/
#endif // KTWONIFFOUR_H

