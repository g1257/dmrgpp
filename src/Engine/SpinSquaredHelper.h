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

/*! \file SpinSquaredHelper.h
 *
 *  Helper class for SpinSquared
 *
 */
#ifndef SPIN_SQ_HELPER_H
#define SPIN_SQ_HELPER_H

namespace Dmrg {	
	template<typename FieldType_,typename Word_>
	class SpinSquaredHelper {
		public:
			typedef FieldType_ FieldType;
			typedef Word_ Word;
		
			SpinSquaredHelper() : data_(0),ketSaved_(0) { }
		
			void operator()(const Word& ket,const Word& bra,const FieldType& value)
			{
				if (ket!=bra) {
					throw std::runtime_error("SpinSquaredHelper::operator(): ket!=bra\n");
				}
				if (ket!=ketSaved_) data_=0;
				ketSaved_=ket;
				data_ += value;
			}

			//! receives m, returns (2*j,m+j)
			std::pair<size_t,size_t> getJmPair(const FieldType& m) const
			{
				size_t j = getJvalue();
				size_t mtilde = getMvalue(m,j);
				return std::pair<size_t,size_t>(j,mtilde);
			}

			void clear() { data_=0; }

		private:
			FieldType data_;
			Word ketSaved_;

			int perfectSquareOrCrash(const FieldType& t) const
			{
				FieldType r = sqrt(t);
				int ri=int(r);
				if (ri!=r) std::runtime_error("SpinSquaredHelper:: sqrt(1+4d) not an integer\n");
				return ri;
			}

			size_t getJvalue() const
			{
				if (data_<0) std::runtime_error("SpinSquaredHelper::getJvalue(): d<0\n");
				int tmp = perfectSquareOrCrash(1.0+4.0*data_);
				if (tmp<1) std::runtime_error("SpinSquaredHelper::getJvalue(): sqrt(1+4d)<1\n");
				size_t ret = tmp-1;
				return ret;	
			}

			size_t getMvalue(const FieldType& m,size_t j) const
			{
				FieldType tmp = m+j*0.5;
				if (tmp<0)  std::runtime_error("SpinSquaredHelper::getMvalue(): j+m <0\n");
				size_t ret = size_t(tmp);
				if (ret!=tmp) std::runtime_error("SpinSquaredHelper::getMvalue(): j+m not size_t\n");
				return ret;
			}
		}; // class SpinSquaredHelper
} // namespace Dmrg

/*@}*/
#endif
