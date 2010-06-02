// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file ClebschGordan.h
 *
 *  
 *
 */
#ifndef CLEBSCH_GORDAN_H
#define CLEBSCH_GORDAN_H

#include "Utils.h"
#include "ProgramLimits.h"

namespace Dmrg {
	//! This is a class to compute ClebschGordan Coefficients
	//! Don't use this class directly, use ClebschGordanCached instead, it'll improve performance
	// Parts taken from Ref.S = http://caps.gsfc.nasa.gov/simpson/software/cg_f90.txt
	//( David G. Simpson NASA, Goddard Space Flight Center, Greenbelt, Maryland  20771)

	template<typename FieldType>
	class ClebschGordan {
			typedef FieldType LongType;
		public:
			typedef std::pair<size_t,size_t> PairType;
			ClebschGordan(size_t dummy=0) 
			{
				static int firstcall=1;
				if (firstcall) createFactorials();
				firstcall=0;
			}

			// receiving format is (2*j,j+m)
			FieldType operator()(const PairType& jm,const PairType& jm1,const PairType& jm2) const
			{
				FieldType j,m,j1,m1,j2,m2;
				convert_(j,m,jm);
				convert_(j1,m1,jm1);
				convert_(j2,m2,jm2);
				FieldType ret=cg(j,m,j1,m1,j2,m2);
				return ret;
			}

		private:
			static std::vector<LongType> factorial_;
			
			bool passesHurdles(FieldType j,FieldType m,FieldType j1,FieldType m1,FieldType j2,FieldType m2) const
			{
				//From Ref.S
				if (isfrac(j1+j2+j) || isfrac(j1+m1) || isfrac(j2+m2) ||
					isfrac(j+m) || isfrac(-j1+j-m2) || isfrac(-j2+j+m1)) return false;
						
				if (m != (m1+m2)) return false; 
				return true;
			}

			//From Ref.S
			bool isNonZero(FieldType j3,FieldType m3,FieldType j1,FieldType m1,FieldType j2,FieldType m2) const
			{
				if ( (j3 < fabs(j1-j2)) || (j3 > (j1+j2)) || (fabs(m1) > j1)  ||
					(fabs(m2) > j2) || (fabs(m3) > j3)) return false;
				return true;
			}

			void convert_(FieldType& j,FieldType& m,const PairType& jm) const
			{
				j=0.5*jm.first;
				m=jm.second-j;
			}

			// From Ref.S
			bool isfrac(FieldType x) const
			{
				FieldType eps=1e-8;
				if ((fabs(x)-int(fabs(x)))>eps) return true;
				return false;	
			}

// 			\delta_{m,m_1+m_2} \sqrt{\frac{(2j+1)(j+j_1-j_2)!(j-j_1+j_2)!(j_1+j_2-j)! }{(j_1+j_2+j+1)!}} \ \times
// 
// 			\sqrt{(j+m)!(j-m)!(j_1-m_1)!(j_1+m_1)!(j_2-m_2)!(j_2+m_2)!}\ \times
// 
// 			\sum_k \frac{(-1)^k}{k!(j_1+j_2-j-k)!(j_1-m_1-k)!(j_2+m_2-k)!(j-j_2+m_1+k)!(j-j_1-m_2+k)!}.
			FieldType cg(FieldType j,FieldType m,FieldType j1,FieldType m1,FieldType j2,FieldType m2) const 
			{
				if (!passesHurdles(j,m,j1,m1,j2,m2)) return 0;
				if (!isNonZero(j,m,j1,m1,j2,m2)) return 0;
				
				// now we consider  m>=0 and assume all hurdles have passed
				return cg_f1(j,m,j1,m1,j2,m2)*cg_f2(j,m,j1,m1,j2,m2);
				
			}
			
			// From Ref. S
			FieldType cg_f1(FieldType j3,FieldType m3,FieldType j1,FieldType m1,FieldType j2,FieldType m2) const 
			{
				FieldType c = sqrt((j3+j3+1)/fact_(nint(j1+j2+j3+1)));
				c *= sqrt(fact_(nint(j1+j2-j3))*fact_(nint(j2+j3-j1))*fact_(nint(j3+j1-j2)));
				c *= sqrt(fact_(nint(j1+m1))*fact_(nint(j1-m1))*fact_(nint(j2+m2))*fact_(nint(j2-m2))*fact_(nint(j3+m3))
						*fact_(nint(j3-m3)));
				return c;
			}

			// From Ref. S
			FieldType cg_f2(FieldType j3,FieldType m3,FieldType j1,FieldType m1,FieldType j2,FieldType m2)  const
			{
				FieldType sumk = 0;
				for (size_t k=0;k<factorial_.size();k++) {
					if (j1+j2-j3-k <0) continue;
					if (j3-j1-m2+k <0) continue;
					if (j3-j2+m1+k <0) continue;
					if (j1-m1-k    <0) continue;
					if (j2+m2-k    <0) continue;
					FieldType term = fact_(nint(j1+j2-j3-k))*fact_(nint(j3-j1-m2+k))*
							 fact_(nint(j3-j2+m1+k))*
							 fact_(nint(j1-m1-k))*fact_(nint(j2+m2-k))*fact_(k);
					if (k%2==1) term = -term;
					sumk += 1.0/term;
				}
				return sumk;
			}
			
			LongType nint(FieldType t) const
			{
				return (LongType)rint(t); // see man rint 
						// CONFORMING TO C99, POSIX.1-2001.
			}

			void createFactorials() 
			{
				factorial_[0]=1;
				for (size_t i=1;i<factorial_.size();i++) factorial_[i]=i*factorial_[i-1];
			}

			LongType fact_(LongType x) const { return factorial_[(int)x]; }

			int parityOf(const FieldType& f)
			{
				int x = int(f);
				if (x%2==0) return 1;
				return -1;
			}
	}; // ClebschGordan
	
	
	template<typename FieldType>
	std::vector<FieldType> ClebschGordan<FieldType>::factorial_(ProgramLimits::NumberOfFactorials);
} // namespace Dmrg

/*@}*/
#endif
