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

/*! \file ReflectionColor
 *
 *  Critical problems:
 *
 *  - getting a l.i. set in an efficient way
 *  - support for fermionic reflections
 *
 *  Low priority work that needs to be done:
 *
 *  - support for bases that change from site to site
 *
 */
#ifndef REFLECTION_COLOR_H
#define REFLECTION_COLOR_H

#include "LAPACK.h"
#include "Matrix.h"
#include "PackIndices.h" // in PsimagLite
#include "ProgressIndicator.h"
#include "Sort.h"
#include "SparseVector.h"

namespace Dmrg
{

// FIXME: MOVE ELSEWHERE:
template <typename RealType>
bool isAlmostZero(const RealType& x, RealType eps = 1e-20)
{
	return (fabs(x) < eps);
}

// FIXME: MOVE ELSEWHERE:
template <typename RealType>
bool isAlmostZero(const std::complex<RealType>& x, RealType eps = 1e-20)
{
	return (fabs(real(x) * real(x) + imag(x) * imag(x)) < eps);
}

template <typename RealType, typename SparseMatrixType>
class ReflectionColor
{

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef SparseVector<typename VectorType::value_type> SparseVectorType;

	enum { AVAILABLE,
		NOT_AVAILABLE,
		COLOR };

public:

	ReflectionColor(const SparseMatrixType& reflection, bool idebug)
	    : reflection_(reflection)
	    , idebug_(idebug)
	{
		//		printFullMatrix(reflection,"reflection");
		findIsolated();

		printInfo();
		//		std::cerr<<"Isolated nodes="<<ipIsolated_.size()<<"\n";
		//		std::cerr<<"Connected nodes="<<ipConnected_.size()<<"\n";

		SizeType firstColor = 2;
		// SparseMatrixType A;
		// generateAmatrix(A,ipConnected_);
		SizeType ncolor = gencolorLabel(reflection, firstColor);

		//		for (SizeType i=0ilabel(ip_isolated) = 0;

		//		% print statistics
		//		% ----------------
		printInfo2(ncolor, firstColor);

		// gencolorCheck(ilabel_,ncolor);

		// gencolorPerm(ilabel_,firstColor,ncolor);

		for (SizeType i = 0; i < ipIsolated_.size(); i++)
			ilabel_[ipIsolated_[i]] = AVAILABLE;

		for (SizeType i = 0; i < ilabel_.size(); i++)
			if (ilabel_[i] == firstColor)
				ipcolor_.push_back(i);
	}

	const typename PsimagLite::Vector<SizeType>::Type& ipcolor() const { return ipcolor_; }

	const typename PsimagLite::Vector<SizeType>::Type& isolated() const { return ipIsolated_; }

private:

	void printInfo() const
	{
		if (!idebug_)
			return;
		std::cout << "ipIsolated ";
		for (SizeType i = 0; i < ipIsolated_.size(); i++)
			std::cout << ipIsolated_[i] << " ";
		std::cout << "\n";
		for (SizeType i = 0; i < ipConnected_.size(); i++)
			std::cout << ipConnected_[i] << " ";
		std::cout << "\n";
	}

	void printInfo2(SizeType ncolor, SizeType firstColor) const
	{
		if (!idebug_)
			return;
		std::cout << "ncolor=" << ncolor << "\n";
		for (SizeType icolor = firstColor; icolor <= ncolor; icolor++) {
			typename PsimagLite::Vector<SizeType>::Type ilist;
			findWithLabel(ilist, ilabel_, icolor);
			std::cout << "color=" << icolor << " size=" << ilist.size() << "\n";
		}
	}

	void generateAmatrix(SparseMatrixType& A, const typename PsimagLite::Vector<SizeType>::Type& ipConnected) const
	{
		A.resize(ipConnected.size());
		SizeType counter = 0;
		typename PsimagLite::Vector<int>::Type ipConnectedInverse(reflection_.rank(), -1);
		for (SizeType i = 0; i < ipConnected.size(); i++)
			ipConnectedInverse[ipConnected[i]] = i;

		for (SizeType i = 0; i < ipConnected.size(); i++) {
			A.setRow(i, counter);
			SizeType ii = ipConnected[i];
			for (int k = reflection_.getRowPtr(ii); k < reflection_.getRowPtr(ii + 1); k++) {
				ComplexOrRealType val = reflection_.getValue(k);
				int col = ipConnectedInverse[reflection_.getCol(k)];
				if (col < 0)
					continue;
				A.pushCol(col);
				A.pushValue(val);
				counter++;
			}
		}
		A.setRow(A.rank(), counter);
		A.checkValidity();
	}

	void findIsolated()
	{
		for (SizeType i = 0; i < reflection_.rank(); i++) {
			SizeType nz = 0;
			bool hasDiagonal = false;
			for (int k = reflection_.getRowPtr(i); k < reflection_.getRowPtr(i + 1); k++) {
				ComplexOrRealType val = reflection_.getValue(k);
				SizeType col = reflection_.getCol(k);
				if (i == col) {
					val = val + 1.0;
					hasDiagonal = true;
				}
				if (isAlmostZero(val, 1e-4))
					continue;
				nz++;
			}
			if (!hasDiagonal)
				nz++;
			if (nz == 1)
				ipIsolated_.push_back(i);
			else
				ipConnected_.push_back(i);
		}
	}

	SizeType gencolorLabel(const SparseMatrixType& A, SizeType firstColor)
	{
		ilabel_.assign(A.rank(), AVAILABLE);

		SizeType ncolor = firstColor - 1;
		typename PsimagLite::Vector<SizeType>::Type ilist;
		RealType eps = 1e-3;

		// while (any( ilabel == unlabeled))
		for (SizeType ii = 0; ii < ilabel_.size(); ii++) {
			if (ilabel_[ii] != AVAILABLE)
				continue;
			ncolor++;
			SizeType icolor = ncolor;
			ilist.clear();
			findWithLabel(ilist, ilabel_, AVAILABLE);
			SizeType ifound = 0;
			for (SizeType i = 0; i < ilist.size(); i++) {
				//				if (ifound >= maxsize) break;

				SizeType ni = ilist[i];
				if (ilabel_[ni] == AVAILABLE) {
					ilabel_[ni] = icolor;
					ifound++;
					// --------------
					// mark neighbors
					// --------------
					typename PsimagLite::Vector<SizeType>::Type jlist;
					findConnected(jlist, ni, A, eps);
					//					jlist = find( A(ni,:) );
					for (SizeType j = 0; j < jlist.size(); j++) {
						SizeType nj = jlist[j];
						if (ilabel_[nj] == AVAILABLE) {
							ilabel_[nj] = NOT_AVAILABLE;
						}
					}
				}
			}
			ilist.clear();
			findWithLabel(ilist, ilabel_, NOT_AVAILABLE);
			//			ilist = find( ilabel == isseen);
			for (SizeType jj = 0; jj < ilist.size(); jj++) {
				ilabel_[ilist[jj]] = AVAILABLE;
			}
			//			if (ilist.size() >= 1) {
			//			   ilabel(ilist) = unlabeled;
			//			}
		}
		return ncolor;
	}

	void findConnected(typename PsimagLite::Vector<SizeType>::Type& jlist,
	    SizeType ni,
	    const SparseMatrixType& A,
	    const RealType& eps) const
	{
		bool hasDiagonal = false;
		for (int k = A.getRowPtr(ni); k < A.getRowPtr(ni + 1); k++) {
			ComplexOrRealType val = A.getValue(k);
			SizeType col = A.getCol(k);
			if (ni == col) {
				hasDiagonal = true;
				val += 1.0;
			}
			if (isAlmostZero(val, eps))
				continue;
			jlist.push_back(col);
		}
		if (!hasDiagonal)
			jlist.push_back(ni);
	}

	void gencolorCheck(const typename PsimagLite::Vector<SizeType>::Type& ilabel, SizeType ncolor) const
	{
		// ------------
		// RealType check
		// ------------
		typename PsimagLite::Vector<SizeType>::Type ilist;
		findWithLabel(ilist, ilabel, ncolor);
		//		isok = max( ilabel == ncolor);
		SizeType isok = *std::max(ilist.begin(), ilist.end());
		if (isok == 0) {
			PsimagLite::String s = "ncolor " + ttos(ncolor) + " max(ilabel) " + ttos(isok) + "\n";
			// throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void findWithLabel(typename PsimagLite::Vector<SizeType>::Type& ilist,
	    const typename PsimagLite::Vector<SizeType>::Type& ilabel,
	    SizeType icolor) const
	{
		for (SizeType i = 0; i < ilabel.size(); i++)
			if (ilabel[i] == icolor)
				ilist.push_back(i);
	}

	//	void gencolorPerm(const typename PsimagLite::Vector<SizeType>::Type& ilabel,
	//			  SizeType firstColor,
	//			  SizeType ncolor)
	//	{
	//		// find size:
	//		typename PsimagLite::Vector<SizeType>::Type ilist;
	//		SizeType ip = 0;
	//		iperm_.resize(reflection_.rank());

	//		for (SizeType icolor=firstColor;icolor<=ncolor;icolor++) {
	//			//			ilist = find( ilabel == icolor);
	//			ilist.clear();
	//			findWithLabel(ilist,ilabel,icolor);
	//			SizeType nc = ilist.size();
	//			//			nc = length(ilist);
	//			for (SizeType j=ip;j<ip+nc;j++) iperm_[j]=ilist[j-ip];

	//			// omitting check
	//			ip += nc;
	//		}
	//		if (!idebug_) return;

	//		//		% ----------------
	//		//		% print statistics
	//		//		% ----------------
	//		std::cout<<"ncolor="<<ncolor<<"\n";
	//		for (SizeType icolor=firstColor;icolor<ncolor;icolor++) {
	//			ilist.clear();
	//			findWithLabel(ilist,ilabel,icolor);
	//			std::cout<<"color="<<icolor<<" size="<<ilist.size()<<"\n";
	//		}

	//	}

	const SparseMatrixType& reflection_;
	bool idebug_;
	typename PsimagLite::Vector<SizeType>::Type ipIsolated_, ipConnected_;
	typename PsimagLite::Vector<SizeType>::Type ilabel_;
	//	typename PsimagLite::Vector<SizeType>::Type iperm_;
	typename PsimagLite::Vector<SizeType>::Type ipcolor_;
}; // class ReflectionColor

} // namespace Dmrg

/*@}*/
#endif // REFLECTION_COLOR_H
