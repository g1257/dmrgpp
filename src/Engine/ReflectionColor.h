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

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "LinearlyIndependentSet.h"
#include "LAPACK.h"
#include "Sort.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class ReflectionColor {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename LeftRightSuperType::SparseMatrixType
			SparseMatrixType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef std::vector<ComplexOrRealType> VectorType;
	typedef SparseVector<typename VectorType::value_type> SparseVectorType;
	typedef LinearlyIndependentSet<RealType,SparseMatrixType>  LinearlyIndependentSetType;

	enum {AVAILABLE,NOT_AVAILABLE,COLOR};

public:


	ReflectionColor(const SparseMatrixType& sSector)
	{
		std::vector<size_t> ipIsolated,ipConnected;
		findIsolated(ipIsolated,ipConnected,sSector);

		std::cerr<<"Isolated nodes="<<ipIsolated.size()<<"\n";
		std::cerr<<"Connected nodes="<<ipConnected.size()<<"\n";

		std::vector<size_t> ilabel(ipConnected.size(),AVAILABLE);
		size_t firstColor = 2;
		SparseMatrixType A;
		generateAmatrix(A,ipConnected,sSector);
		size_t ncolor = gencolorLabel(ilabel,A,firstColor);

		//		% print statistics
		//		% ----------------
		std::cout<<"ncolor="<<ncolor<<"\n";
		for (size_t icolor=firstColor;icolor<=ncolor;icolor++) {
			std::vector<size_t> ilist;
			findWithLabel(ilist,ilabel,icolor);
			std::cout<<"color="<<icolor<<" size="<<ilist.size()<<"\n";
		}

//		assert(ncolor-firstColor==1);

		gencolorCheck(ilabel,ncolor);
		std::vector<size_t> added,added2;
		for (size_t i=0;i<ipIsolated.size();i++) {
			size_t ind = ipIsolated[i];
			lis.fill(ind,sSector,1.0);
			added.push_back(ind);
		}
		std::vector<size_t> iperm(ipConnected.size());

		RealType sector = 1.0;

		size_t ip = 0;
		size_t icolor = firstColor;
		std::vector<size_t> ilist;
		findWithLabel(ilist,ilabel,icolor);
		size_t nc = ilist.size();
		for (size_t j=ip;j<ip+nc;j++) {
			iperm[j]=ilist[j-ip];
			size_t ind = ipConnected[iperm[j]];
			lis.fill(ind,sSector,sector);
		//	added2.push_back(ind);
			added.push_back(ind);
		}


//		for (size_t i=0;i<notadded.size();i++)
//			list.fill(notadded[i],sSector,sector,doTest);
		plusSector_ = lis.size();

		sector = -sector;
		for (size_t j=ip;j<ip+nc;j++) {
			lis.fill(ipConnected[iperm[j]],sSector,sector);
		}
//		for (size_t i=last;i<sSector.rank();i++) {
//			lis.fill(i,sSector,sector,doTest);
//		}

		size_t minuses = lis.size()-plusSector_;

		size_t extraSize = sSector.rank() - (minuses+plusSector_);

		addExtra(added,sSector,extraSize);
		std::ostringstream msg;
		msg<<plusSector_<<" +, "<<minuses<<" -.";
		progress_.printline(msg,std::cout);
		assert(minuses+plusSector_==sSector.rank());
		transform_ = lis.transform();
//		printFullMatrix(transform_,"transform");
	}

private:

	void generateAmatrix(SparseMatrixType& A,
			     const std::vector<size_t>& ipConnected,
			     const SparseMatrixType& sSector) const
	{
		A.resize(ipConnected.size());
		size_t counter=0;
		std::vector<int> ipConnectedInverse(sSector.rank(),-1);
		for (size_t i=0;i<ipConnected.size();i++)
			ipConnectedInverse[ipConnected[i]]=i;

		for (size_t i=0;i<ipConnected.size();i++) {
			A.setRow(i,counter);
			size_t ii = ipConnected[i];
			for (int k=sSector.getRowPtr(ii);k<sSector.getRowPtr(ii+1);k++) {
				ComplexOrRealType val = sSector.getValue(k);
				int col = ipConnectedInverse[sSector.getCol(k)];
				if (col<0) continue;
				A.pushCol(col);
				A.pushValue(val);
				counter++;
			}
		}
		A.setRow(A.rank(),counter);
		A.checkValidity();
	}

	void findIsolated(std::vector<size_t>& ipIsolated,
			  std::vector<size_t>& ipConnected,
			  const SparseMatrixType& sSector) const
	{
		for (size_t i=0;i<sSector.rank();i++) {
			size_t nz = 0;
			bool hasDiagonal = false;
			for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
				ComplexOrRealType val = sSector.getValue(k);
				size_t col = sSector.getCol(k);
				if (i==col) {
					val = val + 1.0;
					hasDiagonal=true;
				}
				if (isAlmostZero(val,1e-4)) continue;
				nz++;
			}
			if (!hasDiagonal) nz++;
			if (nz==1) ipIsolated.push_back(i);
			else ipConnected.push_back(i);
		}
	}

	size_t gencolorLabel(std::vector<size_t>& ilabel,
			     const SparseMatrixType& A,
			     size_t firstColor) const
	{
		size_t ncolor = firstColor-1;
		std::vector<size_t> ilist;
		RealType eps = 1e-3;

		//while (any( ilabel == unlabeled))
		for (size_t ii=0;ii<ilabel.size();ii++) {
			if (ilabel[ii]!=AVAILABLE) continue;
			ncolor++;
			size_t icolor = ncolor;
			ilist.clear();
			findWithLabel(ilist,ilabel,AVAILABLE);
			size_t ifound = 0;
			for (size_t i=0;i<ilist.size();i++) {
//				if (ifound >= maxsize) break;

				size_t ni = ilist[i];
				if (ilabel[ni] == AVAILABLE) {
					ilabel[ni] = icolor;
					ifound++;
					// --------------
					// mark neighbors
					// --------------
					std::vector<size_t> jlist;
					findConnected(jlist,ni,A,eps);
					//					jlist = find( A(ni,:) );
					for (size_t j=0;j<jlist.size();j++) {
						size_t nj = jlist[j];
						if (ilabel[nj] == AVAILABLE) {
							ilabel[nj] = NOT_AVAILABLE;
						}
					}
				}
			}
			ilist.clear();
			findWithLabel(ilist,ilabel,NOT_AVAILABLE);
			//			ilist = find( ilabel == isseen);
			for (size_t jj=0;jj<ilist.size();jj++) {
				ilabel[ilist[jj]] = AVAILABLE;
			}
			//			if (ilist.size() >= 1) {
			//			   ilabel(ilist) = unlabeled;
			//			}
		}
		return ncolor;
	}

	void findConnected(std::vector<size_t>& jlist,
			   size_t ni,
			   const SparseMatrixType& A,
			   const RealType& eps) const
	{
		bool hasDiagonal = false;
		for (int k=A.getRowPtr(ni);k<A.getRowPtr(ni+1);k++) {
			ComplexOrRealType val = A.getValue(k);
			size_t col = A.getCol(k);
			if (ni==col) {
				hasDiagonal=true;
				val += 1.0;
			}
			if (isAlmostZero(val,eps)) continue;
			jlist.push_back(col);
		}
		if (!hasDiagonal) jlist.push_back(ni);
	}

	void gencolorCheck(const std::vector<size_t>& ilabel,size_t ncolor) const
	{
		// ------------
		// double check
		// ------------
		std::vector<size_t> ilist;
		findWithLabel(ilist,ilabel,ncolor);
//		isok = max( ilabel == ncolor);
		size_t isok = *std::max(ilist.begin(),ilist.end());
		if (isok==0) {
			std::string s = "ncolor " + ttos(ncolor) + " max(ilabel) " + ttos(isok) + "\n";
			//throw std::runtime_error(s.c_str());
		}
	}

	void findWithLabel(std::vector<size_t>& ilist,
			   const std::vector<size_t>& ilabel,
			   size_t icolor) const
	{
		for (size_t i=0;i<ilabel.size();i++)
			if (ilabel[i]==icolor) ilist.push_back(i);
	}

	void gencolorPerm(std::vector<size_t>& iperm,
			  const std::vector<size_t>& ilabel,
			  size_t firstColor,
			  size_t ncolor) const
	{
		// find size:
		std::vector<size_t> ilist;
		size_t ip = 0;

		for (size_t icolor=firstColor;icolor<=ncolor;icolor++) {
			//			ilist = find( ilabel == icolor);
			ilist.clear();
			findWithLabel(ilist,ilabel,icolor);
			size_t nc = ilist.size();
			//			nc = length(ilist);
			for (size_t j=ip;j<ip+nc;j++) iperm[j]=ilist[j-ip];

			// omitting check
			ip += nc;
		}

		//		% ----------------
		//		% print statistics
		//		% ----------------
		std::cout<<"ncolor="<<ncolor<<"\n";
		for (size_t icolor=firstColor;icolor<ncolor;icolor++) {
			ilist.clear();
			findWithLabel(ilist,ilabel,icolor);
			std::cout<<"color="<<icolor<<" size="<<ilist.size()<<"\n";
		}

	}


}; // class ReflectionColor

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_COLOR_H
