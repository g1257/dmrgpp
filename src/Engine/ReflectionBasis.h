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

/*! \file ReflectionBasis
 *
 *
 */
#ifndef REFLECTION_BASIS_H
#define REFLECTION_BASIS_H

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "LinearlyIndependentSet.h"
#include "LAPACK.h"
#include "Sort.h"
#include "ReflectionColor.h"

namespace Dmrg {

template<typename RealType,typename SparseMatrixType>
class ReflectionBasis {

	typedef ReflectionColor<RealType,SparseMatrixType> ReflectionColorOrDomType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef std::vector<ComplexOrRealType> VectorType;
	typedef SparseVector<typename VectorType::value_type> SparseVectorType;
	typedef LinearlyIndependentSet<RealType,SparseMatrixType>  LinearlyIndependentSetType;

	enum {AVAILABLE,NOT_AVAILABLE,COLOR};
	enum {GREATER_THAN_ZERO,LESS_THAN_ZERO};

public:

	ReflectionBasis(const SparseMatrixType& reflection)
	: reflection_(reflection),idebug_(true)
	{
		ReflectionColorOrDomType colorOrDom(reflection_);

		setIsolated(colorOrDom.isolated());
		addColor(colorOrDom.ipcolor());

		std::vector<size_t> iavail;
		prepareAvailable(iavail,colorOrDom.ipcolor(),colorOrDom.isolated());

		std::vector<size_t> ip(iavail.size());
		permute(iavail,ip);

		choleskyFactor(iavail);
	}

	const SparseMatrixType& transform() const
	{
		std::string s(__FILE__);
		s += " tranform() unimplemented\n";
		throw std::runtime_error(s.c_str());
	}

private:

	void prepareAvailable(std::vector<size_t>& iavail,
			      const std::vector<size_t>& ipcolor,
			      const std::vector<size_t>& ipIsolated) const
	{
		std::vector<size_t> tmp(reflection_.rank(),1);
		for (size_t i=0;i<ipcolor.size();i++) tmp[ipcolor[i]]=0;
		for (size_t i=0;i<ipIsolated.size();i++) tmp[ipIsolated[i]]=0;
		for (size_t i=0;i<tmp.size();i++)
			if (tmp[i]>0) iavail.push_back(i);
	}

	void addColor(const std::vector<size_t>& ipcolor)
	{
		for (size_t i=0;i<ipcolor.size();i++) {
			ipPos_.push_back(ipcolor[i]);
			ipNeg_.push_back(ipcolor[i]);
		}
	}

	void setIsolated(const std::vector<size_t>& ipIsolated)
	{
		if (ipIsolated.size()==0) return;
		std::vector<ComplexOrRealType> dd;
		setDiagonal(dd,reflection_);
		findPermuted(ipPos_,dd,ipIsolated,GREATER_THAN_ZERO);
		findPermuted(ipNeg_,dd,ipIsolated,LESS_THAN_ZERO);
	}

	void setDiagonal(std::vector<ComplexOrRealType>& dd,const SparseMatrixType& m) const
	{
		for (size_t i=0;i<m.rank();i++) {
			ComplexOrRealType val = 0;
			for (int k=m.getRowPtr(i);k<m.getRowPtr(i+1);k++) {
				if (size_t(m.getCol(k))==i) {
					val = m.getValue(k);
					break;
				}
			}
			dd.push_back(val);
		}
	}

	void setDiagonal(SparseMatrixType& R1,
			 const std::vector<ComplexOrRealType>& dr,
			 const RealType& sector) const
	{
		const std::vector<size_t>& ipPosOrNeg = (sector>0) ? ipPos_ : ipNeg_;

		R1.resize(reflection_.rank());
		size_t counter = 0;
		for (size_t i=0;i<R1.rank();i++) {
			R1.setRow(i,counter);
			R1.pushCol(i);
			ComplexOrRealType val = 1 + sector*dr[ipPosOrNeg[i]];
			RealType val2 = std::norm(val);
			R1.pushValue(sqrt(2*val2));
			counter++;
		}
		R1.setRow(R1.rank(),counter);
		R1.checkValidity();
	}

	void findPermuted(std::vector<size_t>& x,
			  const std::vector<ComplexOrRealType>& dd,
			  const std::vector<size_t>& perm,
			  size_t lessOrGreater)
	{
		for (size_t i=0;i<perm.size();i++) {
			if (lessOrGreaterCondition(dd[perm[i]],lessOrGreater))
				x.push_back(perm[i]);
		}
	}

	bool lessOrGreaterCondition(const ComplexOrRealType& a,size_t lessOrGreater) const
	{
		if (lessOrGreater==GREATER_THAN_ZERO) {
			return (a>0);
		} else {
			return (a<0);
		}
	}

	void permute(std::vector<size_t>& iavail,const std::vector<size_t>& ip) const
	{
		for (size_t i=0;i<iavail.size();i++) iavail[i] = iavail[ip[i]];
	}

	void choleskyFactor(const std::vector<size_t>& iavail)
	{
		std::vector<ComplexOrRealType> dr(reflection_.rank());
		SparseMatrixType R1,Rm;

		setDiagonal(dr,reflection_);
		setDiagonal(R1,dr,1.0);
		setDiagonal(Rm,dr,-1.0);

		for (size_t i=0;i<iavail.size();i++) {
//			size_t ilast = i;
			size_t j = iavail[i];
			if (doneOneSector(i,j,R1,1.0)) break;
			if (doneOneSector(i,j,Rm,-1.0)) break;
		}
//		std::ostringstream msg;
//		msg<<plureflection__<<" +, "<<minuses<<" -.";
//		progress_.printline(msg,std::cout);
	}

	bool doneOneSector(size_t i,size_t j,SparseMatrixType& R1,const RealType& sector)
	{
		bool done = (ipPos_.size()+ipNeg_.size() >= reflection_.rank());
		if (done) return true;

		RealType tol = 1;

		// try to add vector to (sector) space, where sector= + or -
		std::vector<ComplexOrRealType> w;
		setColumn(w,sector,j);
		std::vector<ComplexOrRealType> T1w(w.size(),0);
		std::vector<size_t>& ipPosOrNeg = (sector>0) ? ipPos_ : ipNeg_;
		setT1w(T1w,ipPosOrNeg,w,sector);
		std::vector<ComplexOrRealType> r;
		linearSolverTriangular(r,R1,T1w);
		ComplexOrRealType rkk2 = scalarProduct(w,w) - scalarProduct(r,r);
		if (rkk2>tol) {
			// accept this column
			if (idebug_) {
				std::cerr<<__FILE__<<" "<<__LINE__<<" sector="<<sector;
				std::cerr<<" i="<<i<<" j="<<j<<" #pos="<<(ipPosOrNeg.size()-1)<<"\n";
			}
			grow(R1,r,sqrt(rkk2));
			ipPosOrNeg.push_back(j);
		}
		return false;
	}

	void setT1w(std::vector<ComplexOrRealType>& T1w,
		    const std::vector<size_t>& ipPosOrNeg,
		    const std::vector<ComplexOrRealType>& w,
		    const RealType& sector) const
	{
		for (size_t ii=0;ii<ipPosOrNeg.size();ii++) {
			size_t i = ipPosOrNeg[ii];
			for (int k = reflection_.getRowPtr(i);k<reflection_.getRowPtr(i+1);k++) {
				size_t col = reflection_.getCol(k);
				ComplexOrRealType val = reflection_.getVal(k);
				if (col==i) val += sector;
				val *= sector;
				T1w[col] += val*w[i];
			}
		}
	}

	void setColumn(std::vector<ComplexOrRealType>& w,const RealType& sector,size_t j) const
	{
		for (size_t i=0;i<reflection_.rank();i++) {
			ComplexOrRealType val = 0;
			for (int k = reflection_.getRowPtr(i);k<reflection_.getRowPtr(i+1);k++) {
				size_t col = reflection_.getCol(k);
				if (col==j) {
					val = reflection_.getVal(k);
					break;
				}
			}
			if (j==i) val += sector;
			val *= sector;
			w.push_back(val);
		}
	}

	void findIsolated(std::vector<size_t>& ipIsolated,
			  std::vector<size_t>& ipConnected,
			  const SparseMatrixType& reflection_) const
	{
		for (size_t i=0;i<reflection_.rank();i++) {
			size_t nz = 0;
			bool hasDiagonal = false;
			for (int k=reflection_.getRowPtr(i);k<reflection_.getRowPtr(i+1);k++) {
				ComplexOrRealType val = reflection_.getValue(k);
				size_t col = reflection_.getCol(k);
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

	const SparseMatrixType& reflection_;
	bool idebug_;
	std::vector<size_t> ipPos_,ipNeg_;

}; // class ReflectionBasis

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_BASIS_H
