// BEGIN LICENSE BLOCK
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file FermionSign.h
 *
 *  FIXME documentation
 */
#ifndef FEMION_SIGN_H
#define FEMION_SIGN_H

#include "Utils.h"

namespace Dmrg {
	
	class FermionSign {
		public:
			FermionSign(const std::vector<size_t>& electrons) 
				: signs_(electrons.size())
			{
				init(electrons);
			}
			
			template<typename SomeBasisType>
			FermionSign(const SomeBasisType& basis,const std::vector<size_t>& electrons)
				: signs_(basis.size())
			{
				
				const std::vector<size_t>& basisElectrons = basis.electronsVector(SomeBasisType::BEFORE_TRANSFORM);
				if (basisElectrons.size()!=basis.permutationInverse().size()) throw std::runtime_error("Problem\n");
				size_t nx = basisElectrons.size()/electrons.size();
				std::vector<size_t> el(basisElectrons.size());
				for (size_t x=0;x<basisElectrons.size();x++) {
					size_t x0,x1;
					utils::getCoordinates(x0,x1,basis.permutation(x),nx);
					int nx0 = basisElectrons[x]-electrons[x1];
					if (nx0<0) throw std::runtime_error("FermionSign::ctor(...)\n");
					el[x] = nx0;
				}
				init(el);
			}
			
			template<typename IoInputter>
			FermionSign(IoInputter& io,bool bogus = false)
			{
				if (bogus) return;
				io.read(signs_,"#FERMIONICSIGN");
			}
			
			int operator()(size_t i,size_t f) const
			{
				return (signs_[i]) ? f : 1;
			}
			
			template<typename IoOutputter>
			void save(IoOutputter& io) const
			{
				io.printVector(signs_,"#FERMIONICSIGN");
			}
			
		private:
			
			void init(const std::vector<size_t>& electrons)
			{
				for (size_t i=0;i<signs_.size();i++)
					signs_[i] = (electrons[i] & 1);
			}
			
			std::vector<bool> signs_;
	}; // class FermionSign
} // namespace Dmrg 

/*@}*/
#endif // FEMION_SIGN_H
