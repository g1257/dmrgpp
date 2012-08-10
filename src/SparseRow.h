/*
Copyright (c) 2009, UT-Battelle, LLC
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
/** \ingroup DMRG */
/*@{*/

/*! \file SparseRow.h
 *
 *  Sparse row allows one to accumulate non-zeros to a row
 *  of a CrsMatrix and add the row at the end
 */
#ifndef SPARSE_ROW_H
#define SPARSE_ROW_H

#include "Sort.h"
#include <cassert>

namespace PsimagLite {
	
	template<typename CrsMatrixType>
	class SparseRow {
	public:
		typedef typename CrsMatrixType::value_type ValueType;
		typedef std::vector<size_t> ColumnsType;
		typedef std::vector<ValueType> VectorType;

		void add(size_t col,ValueType value)
		{
			cols_.push_back(col);
			values_.push_back(value);
		}

		ValueType matrixVectorProduct(const VectorType& y) const
		{
			ValueType sum = 0;
			for (size_t i=0;i<cols_.size();i++)
				sum += values_[i]*y[cols_[i]];
			return sum;
		}
//		void clear()
//		{
//			cols_.clear();
//			values_.clear();
//		}

		size_t finalize(CrsMatrixType& matrix)
		{
			assert(cols_.size()==values_.size());
			if (cols_.size()==0) return 0;

			Sort<ColumnsType> s;
			ColumnsType iperm(cols_.size());
			s.sort(cols_,iperm);
			size_t prevCol = cols_[0];
			size_t counter = 0;
			ValueType value = 0;
			for (size_t i=0;i<cols_.size();i++) {
				if (cols_[i]==prevCol) {
					value += values_[iperm[i]];
					continue;
				}
				matrix.pushCol(prevCol);
				matrix.pushValue(value);
				counter++;
				value = values_[iperm[i]];
				prevCol = cols_[i];
			}
			matrix.pushCol(prevCol);
			matrix.pushValue(value);
			counter++;
			return counter;
		}

		ValueType finalize(const VectorType& y)
		{
			assert(cols_.size()==values_.size());
			if (cols_.size()==0) return 0;

			ValueType sum = 0.0;
			for (size_t i=0;i<cols_.size();i++)
				sum += values_[i] * y[cols_[i]];
			return sum;
		}

	private:
		ColumnsType cols_;
		VectorType values_;
			
	}; // class SparseRow

} // namespace Dmrg 

/*@}*/
#endif // SPARSE_ROW_H
