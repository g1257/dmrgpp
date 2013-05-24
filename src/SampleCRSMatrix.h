
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 0.0.1]
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

#ifndef SAMPLE_CRSMATRIX_HEADER_H
#define SAMPLE_CRSMATRIX_HEADER_H

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdlib> // <-- just for drand48 and srand48(seed)
#include "Sort.h"
#include "loki/TypeTraits.h"

namespace PsimagLite {
	
	template<typename T>
	class SampleCRSMatrix {
	public:
		typedef T value_type;
		

		SampleCRSMatrix(SizeType rank) : rank_(rank),rowptr_(rank+1)
		{
		}
		

		SampleCRSMatrix(SizeType rank,T seed,SizeType nonZeros,T maxValue) : rank_(rank),rowptr_(rank+1)
		{
			srand48(seed);
			typename Vector<SizeType>::Type rows,cols;
			typename Vector<T>::Type vals;
			for (SizeType i=0;i<nonZeros;i++) {
				// pick a row
				SizeType row = SizeType(drand48()*rank);
				// and a column
				SizeType col = SizeType(drand48()*rank);
				// and a value
				T val = drand48()*maxValue;
				rows.push_back(row);
				cols.push_back(col);
				vals.push_back(val);
			}
			// fill the matrix with this data:
			fillMatrix(rows,cols,vals);
		}
		

		template<typename SomeIoInputType>
		SampleCRSMatrix(SomeIoInputType& io)
		{
			io>>rank_;
			readVector(io,rowptr_);
			readVector(io,colind_);
			readVector(io,values_);
		}
		
		
		

		void setRow(int n,int v)
		{
			rowptr_[n]=v;
		}
		

		void pushCol(int i) { colind_.push_back(i); }
		

		void pushValue(const T& value) { values_.push_back(value); }
		

		void matrixVectorProduct(typename Vector<T>::Type& x, const typename Vector<T>::Type& y) const
		{
			for (SizeType i = 0; i < y.size(); i++)
				for (SizeType j = rowptr_[i]; j < rowptr_[i + 1]; j++)
					x[i] += values_[j] * y[colind_[j]];
		}
		

		SizeType rank() const { return rank_; }
		

		template<typename SomeIoOutputType>
		void save(SomeIoOutputType& io) const
		{
			io<<rank_<<"\n";
			saveVector(io,rowptr_);
			saveVector(io,colind_);
			saveVector(io,values_);
		}
		
		
	private:
		

		template<typename SomeIoOutputType,typename SomeVectorType>
		typename EnableIf<IsVectorLike<SomeVectorType>::True, void>::Type
		saveVector(SomeIoOutputType& io,const SomeVectorType& v) const
		{
			io<<v.size()<<"\n";
			for (SizeType i=0;i<v.size();i++) {
				io<<v[i]<<" ";
			}
			io<<"\n";
		}
		

		template<typename SomeIoInputType,typename SomeVectorType>
		typename EnableIf<IsVectorLike<SomeVectorType>::True, void>::Type
		readVector(SomeIoInputType& io,SomeVectorType& v) const
		{
			int size=0;
			io>>size;
			if (size<0) throw RuntimeError("readVector: size is zero\n");
			v.resize(size);
			for (SizeType i=0;i<v.size();i++) {
				io>>v[i];
			}
		}
		

		void fillMatrix(typename Vector<SizeType>::Type& rows,typename Vector<SizeType>::Type& cols,
				typename Vector<T>::Type& vals)
		{
			Sort<typename Vector<SizeType>::Type > s;
		    typename Vector<SizeType>::Type iperm(rows.size());
			s.sort(rows,iperm);
			SizeType counter = 0;
			SizeType prevRow = rows[0]+1;
			for (SizeType i=0;i<rows.size();i++) {
				SizeType row = rows[i];
				if (prevRow!=row) {
					// add new row
					rowptr_[row] = counter++;
					prevRow = row;
				}
				colind_.push_back(cols[iperm[i]]);
				values_.push_back(vals[iperm[i]]);
			}
			SizeType lastNonZeroRow = rows[rows.size()-1];
			for (SizeType i=lastNonZeroRow+1;i<=rank_;i++)
				rowptr_[i] = counter;
		}
		
		
		
		SizeType rank_;
		typename Vector<SizeType>::Type rowptr_;
		typename Vector<SizeType>::Type colind_;
		typename Vector<T>::Type values_;
		
	}; // class SampleCRSMatrix
	
	

	template<typename T>
	std::ostream& operator<<(std::ostream& os,const SampleCRSMatrix<T>& m)
	{
		m.save(os);
		return os;
	}
	
	
} // namespace Dmrg
#endif

