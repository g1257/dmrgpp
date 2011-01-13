
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

#ifndef BASISHUBBARDLANCZOS_H
#define BASISHUBBARDLANCZOS_H

#include "Utils.h"

namespace Dmrg {
	
	class BasisHubbardLanczos {
	public:
		
		typedef unsigned int long long WordType;
		static size_t nsite_;
		static psimag::Matrix<size_t> comb_;
		static std::vector<WordType> bitmask_; 
		
		BasisHubbardLanczos(size_t nsite, size_t npart) : npart_(npart)
		{
			if (nsite_>0 && nsite!=nsite_)
				throw std::runtime_error("BasisHubbardLanczos: All basis must have same number of sites\n");
			nsite_ = nsite;
			doCombinatorial();
			doBitmask();

			/* compute size of basis */
			size_t hilbert=1;
			int n=nsite;
			size_t m=1;
			for (;m<=npart;n--,m++)
				hilbert=hilbert*n/m;

			if (data_.size()!=hilbert) {
				data_.clear();
				data_.resize(hilbert);
			}

			if (npart==0) {
				data_[0]=0;
				return;
			}
			
			/* define basis states */
			WordType ket = (1ul<<npart)-1;
			for (size_t i=0;i<hilbert;i++) {
				data_[i] = ket;
				n=m=0;
				for (;(ket&3)!=1;n++,ket>>=1) {
					m += ket&1;
				}
				ket = ((ket+1)<<n) ^ ((1<<m)-1);
			}
			size_ = hilbert;
		} 
		

		size_t size() const { return size_; } 

		const WordType& operator[](size_t i) const
		{
			return data_[i];
		} 

		size_t perfectIndex(WordType state) const
		{
			size_t n=0;
			for (size_t b=0,c=1;state>0;b++,state>>=1)
				if (state&1) n += comb_(b,c++);

			return n;
		} 

		static const WordType& bitmask(size_t i)
		{
			return bitmask_[i];
		} 

		static int bitcnt (WordType b)
		{
		#if (ULONG_MAX == 0xfffffffful)
			b = (0x55555555 & b) + (0x55555555 & (b >> 1));
			b = (0x33333333 & b) + (0x33333333 & (b >> 2));
			b = (0x0f0f0f0f & b) + (0x0f0f0f0f & (b >> 4));
			b = (0x00ff00ff & b) + (0x00ff00ff & (b >> 8));
			b = (0x0000ffff & b) + (0x0000ffff & (b >> 16));

			return (int) b;
		#else
			b = (0x5555555555555555 & b) + (0x5555555555555555 & (b >> 1));
			b = (0x3333333333333333 & b) + (0x3333333333333333 & (b >> 2));
			b = (0x0f0f0f0f0f0f0f0f & b) + (0x0f0f0f0f0f0f0f0f & (b >> 4));
			b = (0x00ff00ff00ff00ff & b) + (0x00ff00ff00ff00ff & (b >> 8));
			b = (0x0000ffff0000ffff & b) + (0x0000ffff0000ffff & (b >> 16));
			b = (0x00000000ffffffff & b) + (0x00000000ffffffff & (b >> 32));

			return (int) b;
		#endif
		}  

	private:
		

		void doCombinatorial()
		{
			/* look-up table for binomial coefficients */
			comb_.resize(nsite_,nsite_);

			for (size_t n=0;n<nsite_;n++)
				for (size_t i=0;i<nsite_;i++)
					comb_(n,i)=0;

			for (size_t n=0;n<nsite_;n++) {
				size_t m = 0;
				int j = n;
				size_t i = 1;
				size_t cnm  = 1;
				for (;m<=n/2;m++,cnm=cnm*j/i,i++,j--)
					comb_(n,m) = comb_(n,n-m) = cnm;
			}
		} 

		void doBitmask()
		{
			bitmask_.resize(nsite_);
			bitmask_[0]=1ul;
			for (size_t i=1;i<nsite_;i++)
				bitmask_[i] = bitmask_[i-1]<<1;
		} 
		
		
		size_t size_;
		size_t npart_;
		std::vector<WordType> data_;
		
	}; // class BasisHubbardLanczos

	size_t BasisHubbardLanczos::nsite_=0;
	psimag::Matrix<size_t> BasisHubbardLanczos::comb_;
	std::vector<BasisHubbardLanczos::WordType> BasisHubbardLanczos::bitmask_; 
	
} // namespace
#endif

