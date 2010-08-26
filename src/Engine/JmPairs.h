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

/*! \file JmPairs.h
 *
 *  This is a "vector" of (2j,m+j) pairs.
 *  Repeated entries are stored only once
 *  Provides a transparent access as if it were a normal vector
 *
 */
#ifndef JMPAIRS_HEADER_H
#define JMPAIRS_HEADER_H

#include "Utils.h"
#include <algorithm>

namespace Dmrg {

	template<typename PairType_>
	class JmPairs {
		public:
			typedef PairType_ PairType;
			typedef PairType value_type;
			
			PairType operator[](size_t alpha) const 
			{
				return jmPairs_[indices_[alpha]];
			}

			JmPairs<PairType>& operator=(const std::vector<PairType>& jmpairs)
			{
				jmPairs_.clear();
				indices_.clear();
				for (size_t i=0;i<jmpairs.size();i++) {
					int x = utils::isInVector(jmPairs_,jmpairs[i]);
					if (x<0) {
						jmPairs_.push_back(jmpairs[i]);
						x=jmPairs_.size()-1;
					}
					indices_.push_back(x);
				}
				return *this;
			}
			
			//! indices_[alpha] = jm
			void push(const PairType& jm,size_t alpha)
			{
				int x = utils::isInVector(jmPairs_,jm);

				if (x<0) {
					jmPairs_.push_back(jm);
					x=jmPairs_.size()-1;
				}

				indices_.push_back(x);
			}

			void clear()
			{
				jmPairs_.clear();
				indices_.clear();
			}

			void reorder(const std::vector<size_t>& permutation)
			{
				utils::reorder(indices_,permutation);
			}

			void truncate(const std::vector<size_t>& removedIndices)
			{
				utils::truncateVector(indices_,removedIndices);
				std::vector<size_t> unusedPairs;
				findUnusedJmPairs(unusedPairs);
				removeUnusedPairs(unusedPairs);
			}

			template<typename Op>
			void maxFirst(size_t& maxvalue)
			{
				Op f;
				for (size_t i=0;i<jmPairs_.size();i++) {
					if (f(jmPairs_[i].first,maxvalue)) {
						maxvalue=jmPairs_[i].first;
					}
				}
			}

			size_t size() const { return indices_.size(); }

			void resize(size_t dummy) { }; // does nothing, safely

			template<typename IoOutputter>
			void save(IoOutputter& io) const
			{
				io.printVector(jmPairs_,"#su2JmPairs");		
				io.printVector(indices_,"#su2JmIndices");		
			}

			template<typename IoInputter>
			void load(IoInputter& io) 
			{
				io.read(jmPairs_,"#su2JmPairs");
				io.read(indices_,"#su2JmIndices");	
			}

		private:
			std::vector<PairType> jmPairs_;
			std::vector<size_t> indices_;

			void findUnusedJmPairs(std::vector<size_t>& unusedPairs)
			{
				for (size_t i=0;i<jmPairs_.size();i++) 
					if (isUnusedPair(i)) unusedPairs.push_back(i);
				
			}

			void removeUnusedPairs(const std::vector<size_t>& unusedPairs)
			{
				size_t counter=0;
				std::vector<size_t> neworder(jmPairs_.size());
				std::vector<PairType> tmpVector(jmPairs_.size()-unusedPairs.size());
				
				for (size_t i=0;i<jmPairs_.size();i++) {
					if (utils::isInVector(unusedPairs,i)>=0) continue;
					tmpVector[counter]=jmPairs_[i];
					neworder[i]=counter;
					counter++;
				}
				jmPairs_=tmpVector;
				std::vector<size_t> tmpVector2(indices_.size());
				for (size_t i=0;i<indices_.size();i++) 
					tmpVector2[i]=neworder[indices_[i]];
				indices_=tmpVector2;
			}

			bool isUnusedPair(size_t ind)
			{
				for (size_t i=0;i<indices_.size();i++) 
					if (indices_[i]==ind) return false;
				return true;
			}
	}; // JmPairs

	template<typename PairType>
	std::ostream& operator<<(std::ostream& os,JmPairs<PairType> jmPairs)
	{
		for (size_t i=0;i<jmPairs.size();i++) 
			os<<"jmPair["<<i<<"]="<<jmPairs[i]<<"\n";
		return os;
	}

	std::istream& operator>>(std::istream& is,std::pair<size_t,size_t>& pair)
	{
		
		is>>pair.first;
		is>>pair.second;
		return is;
	}
} // namespace Dmrg

/*@}*/
#endif

