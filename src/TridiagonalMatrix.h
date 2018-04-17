/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file TridiagonalMatrix.h
 *
 *  A struct represent a tridiagonal matrix
 */
#ifndef TRIDIAGONAL_MATRIX_H
#define TRIDIAGONAL_MATRIX_H
#include "Matrix.h"

namespace PsimagLite {

template<typename FieldType>
class TridiagonalMatrix {

	typedef typename Vector<FieldType>::Type VectorType;
	typedef typename Real<FieldType>::Type RealType;

public:

	typedef FieldType value_type;

	TridiagonalMatrix()  { }

	template<typename IoInputType>
	TridiagonalMatrix(IoInputType& io)
	{
		io.read(a_,"#Avector");
		io.read(b_,"#Bvector");
	}

	template<typename IoOutputType>
	void write(IoOutputType& io) const
	{
		io.write(a_,"#Avector");
		io.write(b_,"#Bvector");
	}

	void resize(SizeType n,FieldType value)
	{
		resize(n);
		for (SizeType i=0;i<n;i++) a_[i]=b_[i]=value;
	}

	void resize(SizeType n)
	{
		a_.resize(n);
		b_.resize(n);
	}

	FieldType& a(SizeType i) { return a_[i]; }

	FieldType& b(SizeType i) { return b_[i]; }

	const FieldType& a(SizeType i) const { return a_[i]; }

	const FieldType& b(SizeType i) const { return b_[i]; }

	template<typename SomeMatrixType>
	void buildDenseMatrix(SomeMatrixType& m) const
	{
		m.resize(a_.size(),a_.size());
		for (SizeType i=0;i<m.rows();i++)
			for (SizeType j=0;j<m.cols();j++)
				m(i,j)=0;

		for (SizeType i=0;i<a_.size();i++) {
			m(i,i) = a_[i];
			if (i>0) m(i,i-1) = m(i-1,i) = b_[i-1];
		}
	}

	void push(const FieldType& a,const FieldType& b)
	{
		a_.push_back(a);
		b_.push_back(b);
	}

	SizeType size() const { return a_.size(); }

	template<typename SomeVectorType>
	FieldType excited(SomeVectorType &z, SizeType level) const
	{
		if (a_.size() > 4900)
			throw RuntimeError("TridiagonalMatrix::excited: too big\n");

		typedef typename SomeVectorType::value_type ElementType;
		Matrix<ElementType> m;
		buildDenseMatrix(m);
		SizeType n = m.rows();
		assert(m.cols() == n);
		typename Vector<RealType>::Type eigs(n);
		diag(m,eigs,'V');
		assert(level < n);
		for (SizeType i = 0; i < n; ++i)
			z[i] = m(i,level);

		return eigs[level];
	}

private:
	VectorType a_,b_;
}; // class TridiagonalMatrix
} // namespace PsimagLite

/*@}*/
#endif

