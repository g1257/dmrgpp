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

/*! \file SparseVector.h
 *
 *  A class to represent sparse vectors
 *
 */
#ifndef DMRG_SPARSEVECTOR_H
#define DMRG_SPARSEVECTOR_H

#include <utility>
#include <iostream>
#include <algorithm>
#include "Vector.h" // in PsimagLite

namespace Dmrg {
	template<typename FieldType>
	struct SparseVector {
		public:
			typedef FieldType value_type;
			typedef std::pair<size_t,size_t> PairType;
			
			SparseVector(const std::vector<FieldType>& v) : size_(v.size())
			{
				FieldType zerovalue=static_cast<FieldType>(0);
				for (size_t i=0;i<v.size();i++) {
					if (v[i]!=zerovalue) {
						values_.push_back(v[i]);
						indices_.push_back(i);
					}
				}
			}
			
			void fromChunk(const std::vector<FieldType>& v,size_t offset,size_t total)
			{
				resize(total);
				for (size_t i=0;i<v.size();i++) {
					indices_.push_back(i+offset);
					values_.push_back(v[i]);
				}
			}

			SparseVector() {}
			
			SparseVector(size_t n) : size_(n) { }
			
			/*void normalize()
			{
				FieldType norma = std::norm(values_);
				for (size_t i=0;i<values_.size();i++) values_[i] /= norma;
			}*/

			void resize(size_t x)
			{
				values_.clear();
				indices_.clear();
				size_=x;
			}

			FieldType& operator[](size_t index)
			{
				int i=PsimagLite::isInVector(indices_,index);
				if (i<0) i=add(index,0);
				return values_[i];	
			}
			
			//FIXME: disable due to performance reasons
			FieldType operator[](size_t index) const
			{
				int i=PsimagLite::isInVector(indices_,index);
				if (i<0) {
					//std::cerr<<"index="<<index<<"\n";
					//utils::vectorPrint(indices_,"indices",std::cerr);
					//throw std::runtime_error("SparseVector::operator[](): index out of range.\n");
					return 0;
				}
				return values_[i];	
			}

			//! adds an index (maybe the indices should be sorted at some point)
			size_t add(int index,const FieldType& value)
			{
				indices_.push_back(index);
				values_.push_back(value);
				return values_.size()-1;
			}

			size_t size() const { return size_; }

			void print(std::ostream& os,const std::string& label) const
			{
				os<<label<<"\n";
				os<<size_<<"\n";
				os<<indices_.size()<<"\n";
				for (size_t i=0;i<indices_.size();i++) {
					os<<indices_[i]<<" "<<values_[i]<<"\n";
				}
			}

			size_t indices() const { return indices_.size(); }

			size_t index(size_t x) const { return indices_[x]; }
			
			FieldType value(size_t x) const { return values_[x]; }

			/*FieldType norm2() const
			{
				FieldType sum=0;
				for (size_t i=0;i<indices_.size();i++) 
					sum += std::conj(values_[i])*values_[i];
				return sum;
			}*/
			
			void toChunk(std::vector<FieldType>& dest,size_t i0, size_t total, bool test=false) const
			{
				if (test) {
					PairType firstLast = findFirstLast();
					if (i0>firstLast.first || i0+total<firstLast.second) 
						throw std::runtime_error("SparseVector::toChunk(...)"
							" check failed\n");
				}
				dest.resize(total);
				for (size_t i=0;i<indices_.size();i++) 
					dest[indices_[i]-i0]=values_[i];
			}
			
			template<typename SomeBasisType>
			size_t toChunk(std::vector<FieldType>& dest,const SomeBasisType& parts) const
			{
				size_t part = findPartition(parts);
				size_t offset = parts.partition(part);
				size_t total = parts.partition(part+1)-offset;
				dest.resize(total);
				for (size_t i=0;i<total;i++) dest[i]=0;
				for (size_t i=0;i<indices_.size();i++) 
					dest[indices_[i]-offset]=values_[i];
				return part;
			}
			
			template<typename SomeBasisType>
			size_t findPartition(const SomeBasisType& parts) const
			{
				PairType firstLast = findFirstLast();
				size_t ret = 0;
				for (size_t i=0;i<parts.partition();i++) {
					if (firstLast.first>=parts.partition(i)) {
						ret = i;
					} else {
						break;
					}
				}
				size_t ret2 = 1;
				for (size_t i=0;i<parts.partition();i++) {
					if (firstLast.second>parts.partition(i)) {
						ret2 = i;
					} else {
						break;
					}
				}
				if (ret!=ret2) throw std::runtime_error("SparseVector::findPartition(...)"
						"vector extends more than one partition\n");
				return ret;
			}
			
		private:
			
			PairType findFirstLast() const
			{
				return PairType(*(std::min_element(indices_.begin(),indices_.end() ) ),
						*( std::max_element(indices_.begin(), indices_.end() ) ));
			}
			
			std::vector<FieldType> values_;
			std::vector<size_t> indices_;
			size_t size_;
	}; // class SparseVector
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,const SparseVector<FieldType>& s)
	{
		s.print(os,"SparseVector");
		return os;
	}
	
} // namespace Dmrg

namespace std {
	template<typename FieldType>
	inline FieldType norm(const Dmrg::SparseVector<FieldType>& v)
	{
		FieldType sum=0;
		for (size_t i=0;i<v.indices();i++) sum += conj(v.value(i))*v.value(i);
		return sqrt(sum);
	}
	
	template<typename FieldType>
	inline FieldType norm(const Dmrg::SparseVector<std::complex<FieldType> >& v)
	{
		std::complex<FieldType> sum=0;
		for (size_t i=0;i<v.indices();i++) sum += conj(v.value(i))*v.value(i);
		return real(sqrt(sum));
	}
}
/*@}*/
#endif

