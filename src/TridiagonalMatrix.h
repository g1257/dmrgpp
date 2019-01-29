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

	typedef typename Real<FieldType>::Type RealType;

	static const bool diagWithLapack_ = true;

public:

	typedef typename Vector<FieldType>::Type VectorType;
	typedef typename Vector<RealType>::Type VectorRealType;
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
	void buildDenseMatrix(SomeMatrixType& m, SizeType n = 0) const
	{
		if (n == 0) n = a_.size();
		m.resize(n, n, 0);
		for (SizeType i = 0; i < n - 1; ++i) {
			m(i, i) = a_[i];
			m(i, i + 1) = b_[i + 1];
			m(i + 1, i) = PsimagLite::conj(b_[i + 1]);
		}

		m(n - 1, n - 1) = a_[n-1];
	}

	void diag(VectorRealType& eigs, SizeType nn) const
	{
		if (diagWithLapack_)
			diag2(eigs, nn);
		else
			ground(eigs, nn);
	}

	void push(const FieldType& a,const FieldType& b)
	{
		a_.push_back(a);
		b_.push_back(b);
	}

	SizeType size() const { return a_.size(); }

private:

	void diag2(VectorRealType&, SizeType) const;

	void ground(VectorRealType& groundD, SizeType nn) const
	{
		const int n = nn;

		groundD.resize(n);
		groundAllocations(n);

		const long int maxCounter = 10000;

		assert(a_.size() >= nn && b_.size() >= nn);
		for (SizeType i = 0; i < nn; ++i) {
			groundD[i] = a_[i];
			groundE_[i] = b_[i];
		}

		RealType s = 0;
		long int intCounter=0;
		int m = 0;
		int l = 0;
		for (; l < n; l++) {
			do {
				intCounter++;
				if (intCounter > maxCounter) {
					std::cerr<<"lanczos: ground: premature exit ";
					std::cerr<<"(may indicate an internal error)\n";
					break;
				}

				for (m = l; m < n - 1; m++) {
					RealType dd = fabs(groundD[m]) + fabs(groundD[m + 1]);
					if ((fabs(groundE_[m]) + dd) == dd) break;
				}

				if (m != l) {
					RealType g = (groundD[l + 1] - groundD[l])/(2.0*groundE_[l]);
					RealType r = sqrt(g*g + 1.0);
					g = groundD[m] - groundD[l] + groundE_[l]/
					        (g + (g >= 0 ? fabs(r) : -fabs(r)));
					RealType p = 0.0;
					RealType c = 1.0;
					int i = m -1;
					for (s = 1.0; i >= l; i--) {
						RealType f = s * groundE_[i];
						RealType h = c * groundE_[i];
						groundE_[i + 1] = (r = sqrt(f * f + g * g));
						if (r == 0.0) {
							groundD[i + 1] -= p;
							groundE_[m] = 0.0;
							break;
						}

						s = f / r;
						c = g / r;
						g = groundD[i + 1] - p;
						r = (groundD[i] - g) * s + 2.0 * c * h;
						groundD[i + 1] = g + (p = s * r);
						g = c * r - h;
					}

					if (r == 0.0 && i >= l) continue;
					groundD[l] -= p;
					groundE_[l] = g;
					groundE_[m] = 0.0;
				}
			} while (m != l);
		}

		std::sort(groundD.begin(), groundD.end());

		if (intCounter > maxCounter)
			throw RuntimeError(String(__FILE__) + "::ground(): internal error\n");
	}

	void groundAllocations(SizeType n) const
	{
		if (groundE_.size() != n) {
			groundE_.clear();
			groundE_.resize(n);
		}
	}

	VectorRealType a_;
	VectorRealType b_;
	mutable VectorRealType groundE_;
}; // class TridiagonalMatrix
} // namespace PsimagLite

/*@}*/
#endif

