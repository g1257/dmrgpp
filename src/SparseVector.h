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

/*! \file SparseVector.h
 *
 *  A class to represent sparse vectors
 *
 */
#ifndef PsimagLite_SPARSEVECTOR_H
#define PsimagLite_SPARSEVECTOR_H

#include <utility>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "Vector.h" // in PsimagLite
#include "Sort.h"
#include "String.h"

namespace PsimagLite {
	template<typename FieldType>
	struct SparseVector {
		public:
			typedef FieldType value_type;
			typedef std::pair<SizeType,SizeType> PairType;
			
			SparseVector(const typename Vector<FieldType>::Type& v)
			: size_(v.size()),isSorted_(false)
			{
				FieldType zerovalue=static_cast<FieldType>(0);
				for (SizeType i=0;i<v.size();i++) {
					if (v[i]!=zerovalue) {
						values_.push_back(v[i]);
						indices_.push_back(i);
					}
				}
			}
			
			void fromChunk(const typename Vector<FieldType>::Type& v,SizeType offset,SizeType total)
			{
				resize(total);
				for (SizeType i=0;i<v.size();i++) {
					indices_.push_back(i+offset);
					values_.push_back(v[i]);
				}
			}

			SparseVector(SizeType n) : size_(n),isSorted_(false) { }
			
			void resize(SizeType x)
			{
				values_.clear();
				indices_.clear();
				size_=x;
			}

			//! adds an index (maybe the indices should be sorted at some point)
			SizeType add(int index,const FieldType& value)
			{
				indices_.push_back(index);
				values_.push_back(value);
				isSorted_ = false;
				return values_.size()-1;
			}

			SizeType size() const { return size_; }

			void print(std::ostream& os,const String& label) const
			{
				os<<label<<", ";
				os<<size_<<", ";
				os<<indices_.size()<<", ";
				for (SizeType i=0;i<indices_.size();i++) {
					os<<indices_[i]<<" "<<values_[i]<<",";
				}
				os<<"\n";
			}

			SizeType indices() const { return indices_.size(); }

			SizeType index(SizeType x) const { return indices_[x]; }
			
			FieldType value(SizeType x) const { return values_[x]; }

			void toChunk(typename Vector<FieldType>::Type& dest,SizeType i0, SizeType total, bool test=false) const
			{
				if (test) {
					PairType firstLast = findFirstLast();
					if (i0>firstLast.first || i0+total<firstLast.second) 
						throw RuntimeError("SparseVector::toChunk(...)"
							" check failed\n");
				}
				dest.resize(total);
				for (SizeType i=0;i<indices_.size();i++) 
					dest[indices_[i]-i0]=values_[i];
			}
			
			template<typename SomeBasisType>
			SizeType toChunk(typename Vector<FieldType>::Type& dest,const SomeBasisType& parts) const
			{
				SizeType part = findPartition(parts);
				SizeType offset = parts.partition(part);
				SizeType total = parts.partition(part+1)-offset;
				dest.resize(total);
				for (SizeType i=0;i<total;i++) dest[i]=0;
				for (SizeType i=0;i<indices_.size();i++) 
					dest[indices_[i]-offset]=values_[i];
				return part;
			}
			
			template<typename SomeBasisType>
			SizeType findPartition(const SomeBasisType& parts) const
			{
				PairType firstLast = findFirstLast();
				SizeType ret = 0;
				for (SizeType i=0;i<parts.partition();i++) {
					if (firstLast.first>=parts.partition(i)) {
						ret = i;
					} else {
						break;
					}
				}
				SizeType ret2 = 1;
				for (SizeType i=0;i<parts.partition();i++) {
					if (firstLast.second>parts.partition(i)) {
						ret2 = i;
					} else {
						break;
					}
				}
				if (ret!=ret2) throw RuntimeError("SparseVector::findPartition(...)"
						"vector extends more than one partition\n");
				return ret;
			}

			template<typename T>
			SparseVector<FieldType> operator*=(const T& val)
			{
				 for (SizeType i=0;i<values_.size();i++) values_[i] *= val;
				 return *this;
			}

			SparseVector<FieldType> operator-=(const SparseVector<FieldType>& v)
			{
				 for (SizeType i=0;i<v.values_.size();i++) {
					 values_.push_back(-v.values_[i]);
					 indices_.push_back(v.indices_[i]);
				 }
				 isSorted_ = false;
				 return *this;
			}

			bool operator==(const SparseVector<FieldType>& v) const
			{
				assert(isSorted_);
				assert(v.isSorted_);
				for (SizeType i=0;i<v.values_.size();i++) {
					if (indices_[i]!=v.indices_[i]) return false;
					FieldType val = values_[i] - v.values_[i];
					if (!isAlmostZero(val,1e-8)) return false;
				}
				return true;
			}

			// FIXME : needs performance
			FieldType scalarProduct(const SparseVector<FieldType>& v) const
			{
				assert(indices_.size()==1 || isSorted_);
				assert(v.indices_.size()==1 || v.isSorted_);
				FieldType sum = 0;
				SizeType i = 0, j = 0, index = 0;

				for (;i<indices_.size();i++) {
					index = indices_[i];
					for (;j<v.indices_.size();j++) {
						if (v.indices_[j]<index) continue;
						if (v.indices_[j]>index) break;
						if (v.indices_[j]==index) sum += values_[i] * std::conj(v.values_[j]);
					}
					if (j>0) j--;
				}

				return sum;
			}

			bool isOne()
			{
				SizeType x = 0;
				for (SizeType i=0;i<indices_.size();i++) {
					if (isAlmostZero(values_[i],1e-4)) continue;
					x++;
					if (x>1) return false;
				}
				return true;
			}

			void sort()
			{
				if (indices_.size()<2 || isSorted_) return;

				Sort<typename Vector<SizeType>::Type > sort;
				typename Vector<SizeType>::Type iperm(indices_.size());
				sort.sort(indices_,iperm);
				typename Vector<FieldType>::Type values(iperm.size());
				for (SizeType i=0;i<values_.size();i++) values[i] = values_[iperm[i]];
				values_.clear();
				FieldType sum = values[0];
				SizeType prevIndex = indices_[0];
				typename Vector<SizeType>::Type indices;

				for (SizeType i=1;i<indices_.size();i++) {

					if (indices_[i]!=prevIndex) {
						if (std::norm(sum)>1e-16) {
							values_.push_back(sum);
							indices.push_back(prevIndex);
						}
						sum = values[i];
						prevIndex = indices_[i];
					} else {
						sum += values[i];
					}
				}
				if (std::norm(sum)>1e-16) {
					values_.push_back(sum);
					indices.push_back(prevIndex);
				}
				indices_=indices;
				isSorted_ = true;

			}

			void clear()
			{
				values_.clear();
				indices_.clear();
				isSorted_=false;
			}

			template<typename T,typename T2>
			friend SparseVector<T2> operator*(const T& val,const SparseVector<T2>& sv);
			
		private:
			


			PairType findFirstLast() const
			{
				return PairType(*(std::min_element(indices_.begin(),indices_.end() ) ),
						*( std::max_element(indices_.begin(), indices_.end() ) ));
			}
			
			typename Vector<FieldType>::Type values_;
			typename Vector<SizeType>::Type indices_;
			SizeType size_;
			bool isSorted_;
	}; // class SparseVector
	
	template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,const SparseVector<FieldType>& s)
	{
		s.print(os,"SparseVector");
		return os;
	}

	template<typename T,typename T2>
	SparseVector<T2> operator*(const T& val,const SparseVector<T2>& sv)
	{
		 SparseVector<T2> res = sv;
		 for (SizeType i=0;i<res.values_.size();i++) res.values_[i] *= val;
		 return res;
	}

	template<typename T>
	T operator*(const SparseVector<T>& v1,const SparseVector<T>& v2)
	{
		return v1.scalarProduct(v2);
	}
} // namespace PsimagLite

namespace PsimagLite {
	template<typename FieldType>
	inline FieldType norm(const PsimagLite::SparseVector<FieldType>& v)
	{
		FieldType sum=0;
		for (SizeType i=0;i<v.indices();i++) sum += std::conj(v.value(i))*v.value(i);
		return sqrt(sum);
	}
	
	template<typename FieldType>
	inline FieldType norm(const PsimagLite::SparseVector<std::complex<FieldType> >& v)
	{
		std::complex<FieldType> sum=0;
		for (SizeType i=0;i<v.indices();i++) sum += std::conj(v.value(i))*v.value(i);
		return real(sqrt(sum));
	}

}
/*@}*/
#endif

