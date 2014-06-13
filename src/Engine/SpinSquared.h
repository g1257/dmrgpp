/*
Copyright (c) 2009, UT-Battelle, LLC
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

/*! \file SpinSquared.h
 *
 *  Encapsulates the calculation of S^2 from creation
 *  and anihilation operators. Works in conjunction with SpinSquaredHelper
 *
 */
#ifndef SPIN_SQUARED_H
#define SPIN_SQUARED_H

namespace Dmrg {
template<typename CallbackType>
class SpinSquared {

	typedef typename CallbackType::Word Word;
	typedef typename CallbackType::FieldType FieldType;

public:

	enum {SPIN_UP=0,SPIN_DOWN=1};

	enum {ORBITAL_A=0,ORBITAL_B=1};

	SpinSquared(CallbackType& callback,
	            int NUMBER_OF_ORBITALS1,
	            int DEGREES_OF_FREEDOM1)
	    : callback_(callback),
	      NUMBER_OF_ORBITALS(NUMBER_OF_ORBITALS1),
	      DEGREES_OF_FREEDOM(DEGREES_OF_FREEDOM1)
	{}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType x,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "SpinSquared";

		const char* start = reinterpret_cast<const char *>(this);
		const char* end = reinterpret_cast<const char *>(&NUMBER_OF_ORBITALS);
		SizeType total = end - start;
		mres.push(SomeMemResolvType::MEMORY_HEAPPTR,
		          SomeMemResolvType::SIZEOF_HEAPREF,
		          this,
		          str + " ref to callback");

		mres.memResolv(&callback_, 0, str + " callback");

		start = end;
		end = reinterpret_cast<const char *>(&DEGREES_OF_FREEDOM);
		total += mres.memResolv(&NUMBER_OF_ORBITALS,
		                        end-start,
		                        str + " NUMBER_OF_ORBITALS");

		total += mres.memResolv(&DEGREES_OF_FREEDOM,
		                        sizeof(*this) - total,
		                        str + " DEGREES_OF_FREEDOM");

		return total;
	}

	void doOnePairOfSitesA(const Word& ket,SizeType i,SizeType j) const
	{
		FieldType value = 0.5;
		Word bra=ket;
		int ret = c(bra,j,ORBITAL_A,SPIN_UP);
		ret = cDagger(bra,j,ORBITAL_A,SPIN_DOWN,ret);
		ret = c(bra,i,ORBITAL_A,SPIN_DOWN,ret);
		ret = cDagger(bra,i,ORBITAL_A,SPIN_UP,ret);
		if (ret==0) callback_(ket,bra,value);
		if (NUMBER_OF_ORBITALS==1) return;

		bra=ket;
		ret = c(bra,j,ORBITAL_B,SPIN_UP);
		ret = cDagger(bra,j,ORBITAL_B,SPIN_DOWN,ret);
		ret = c(bra,i,ORBITAL_A,SPIN_DOWN,ret);
		ret = cDagger(bra,i,ORBITAL_A,SPIN_UP,ret);
		if (ret==0) callback_(ket,bra,value);

		bra=ket;
		ret = c(bra,j,ORBITAL_A,SPIN_UP);
		ret = cDagger(bra,j,ORBITAL_A,SPIN_DOWN,ret);
		ret = c(bra,i,ORBITAL_B,SPIN_DOWN,ret);
		ret = cDagger(bra,i,ORBITAL_B,SPIN_UP,ret);
		if (ret==0) callback_(ket,bra,value);

		bra=ket;
		ret = c(bra,j,ORBITAL_B,SPIN_UP);
		ret = cDagger(bra,j,ORBITAL_B,SPIN_DOWN,ret);
		ret = c(bra,i,ORBITAL_B,SPIN_DOWN,ret);
		ret = cDagger(bra,i,ORBITAL_B,SPIN_UP,ret);
		if (ret==0) callback_(ket,bra,value);
	}

	void doOnePairOfSitesB(const Word& ket,SizeType i,SizeType j) const
	{
		FieldType value = 0.5;
		Word bra=ket;
		int ret = c(bra,j,ORBITAL_A,SPIN_DOWN);
		ret = cDagger(bra,j,ORBITAL_A,SPIN_UP,ret);
		ret = c(bra,i,ORBITAL_A,SPIN_UP,ret);
		ret = cDagger(bra,i,ORBITAL_A,SPIN_DOWN,ret);
		if (ret==0) callback_(ket,bra,value);
		if (NUMBER_OF_ORBITALS==1) return;

		bra=ket;
		ret = c(bra,j,ORBITAL_B,SPIN_DOWN);
		ret = cDagger(bra,j,ORBITAL_B,SPIN_UP,ret);
		ret = c(bra,i,ORBITAL_A,SPIN_UP,ret);
		ret = cDagger(bra,i,ORBITAL_A,SPIN_DOWN,ret);
		if (ret==0) callback_(ket,bra,value);

		bra=ket;
		ret = c(bra,j,ORBITAL_A,SPIN_DOWN);
		ret = cDagger(bra,j,ORBITAL_A,SPIN_UP,ret);
		ret = c(bra,i,ORBITAL_B,SPIN_UP,ret);
		ret = cDagger(bra,i,ORBITAL_B,SPIN_DOWN,ret);
		if (ret==0) callback_(ket,bra,value);

		bra=ket;
		ret = c(bra,j,ORBITAL_B,SPIN_DOWN);
		ret = cDagger(bra,j,ORBITAL_B,SPIN_UP,ret);
		ret = c(bra,i,ORBITAL_B,SPIN_UP,ret);
		ret = cDagger(bra,i,ORBITAL_B,SPIN_DOWN,ret);
		if (ret==0) callback_(ket,bra,value);
	}

	void doDiagonal(const Word& ket,SizeType i,SizeType j) const
	{
		FieldType value = spinZ(ket,i)*spinZ(ket,j);
		callback_(ket,ket,value);
	}

	FieldType spinZ(const Word& ket,SizeType i) const
	{
		int sum=0;
		for (SizeType gamma=0;gamma<SizeType(NUMBER_OF_ORBITALS);gamma++) {
			sum += n(ket,i,gamma,SPIN_UP);
			sum -= n(ket,i,gamma,SPIN_DOWN);
		}
		return sum*0.5;
	}

private:

	int c(Word& ket,SizeType i,SizeType gamma,SizeType spin,int ret=0) const
	{
		if (ret<0) return ret;
		SizeType shift_ = (DEGREES_OF_FREEDOM*i);
		Word mask = 1<<(shift_ + gamma+spin*NUMBER_OF_ORBITALS);
		if ((ket & mask)==0) return -1;
		ket ^= mask;
		return 0;
	}

	int cDagger(Word& ket,SizeType i,SizeType gamma,SizeType spin,int ret=0) const
	{
		if (ret<0) return ret;
		SizeType shift_ = (DEGREES_OF_FREEDOM*i);
		Word mask = 1<<(shift_+gamma+spin*NUMBER_OF_ORBITALS);
		if ((ket & mask)>0) return -1;
		ket |= mask;
		return 0;
	}

	int n(const Word& ket,SizeType i,SizeType gamma,SizeType spin) const
	{
		SizeType shift_ = (DEGREES_OF_FREEDOM*i);
		Word mask = 1<<(shift_+gamma+spin*NUMBER_OF_ORBITALS);
		if ((ket & mask)>0) return 1;
		return 0;
	}

	//serializr start class SpinSquared
	//serializr ref callback_ this
	CallbackType& callback_;
	//serializr normal NUMBER_OF_ORBITALS
	int const NUMBER_OF_ORBITALS;
	//serializr normal DEGREES_OF_FREEDOM
	int const DEGREES_OF_FREEDOM;
}; // SpinSquared
} // namespace Dmrg

/*@}*/
#endif

