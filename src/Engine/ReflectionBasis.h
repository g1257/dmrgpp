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

	enum {AVAILABLE,NOT_AVAILABLE,COLOR};
	enum {GREATER_THAN_ZERO,LESS_THAN_ZERO};

public:

	ReflectionBasis(const SparseMatrixType& reflection,bool idebug)
	: progress_("ReflectionBasis",0),reflection_(reflection),idebug_(idebug)
	{
		ReflectionColorOrDomType colorOrDom(reflection_);

		setIsolated(colorOrDom.isolated());
		addColor(colorOrDom.ipcolor());

		std::vector<size_t> iavail;
		prepareAvailable(iavail,colorOrDom.ipcolor(),colorOrDom.isolated());

		std::vector<size_t> ip(iavail.size());
		permute(iavail,ip);

		choleskyFactor(iavail);
		if (!idebug_) return;

		printFullMatrix(reflection,"reflection");
		std::cout<<"ipPos ";
		for (size_t i=0;i<ipPos_.size();i++)
			std::cout<<ipPos_[i]<<" ";
		std::cout<<"\n";

		std::cout<<"ipNeg ";
		for (size_t i=0;i<ipNeg_.size();i++)
			std::cout<<ipNeg_[i]<<" ";
		std::cout<<"\n";

		printFullMatrix(R1_,"R1_");
		printFullMatrix(Rm_,"Rm_");
	}

	const SparseMatrixType& R(const RealType& sector) const
	{
		return (sector>0) ? R1_ : Rm_;
	}

	const std::vector<size_t>& ipPosOrNeg(const RealType& sector) const
	{
		return (sector>0) ? ipPos_ : ipNeg_;
	}

	//! Invert triangular matrix R into Rinverse
	//! Hack due to not having rectangular CRS implemented
	void inverseTriangular(SparseMatrixType& R1Inverse,
			       const SparseMatrixType& R1,
			       const RealType& sector) const
	{
		SparseMatrixType R1t;
		transposeConjugate(R1t,R1);
		std::vector<ComplexOrRealType> r(R1t.rank());
		SparseMatrixType tmpMatrix(r.size(),r.size());
		size_t counter=0;
		for (size_t i=0;i<R1t.rank();i++) {
			tmpMatrix.setRow(i,counter);
			std::vector<ComplexOrRealType> rhs(R1t.rank(),0);
			rhs[i]=1;
			linearSolverTriangular(r,R1t,rhs);
			for (size_t i=0;i<r.size();i++) {
				if (isAlmostZero(r[i],1e-10)) continue;
				assert(i<((sector>0) ? ipPos_.size() : ipNeg_.size()));
				tmpMatrix.pushCol(i);
				tmpMatrix.pushValue(r[i]);
				counter++;
			}
		}
		tmpMatrix.setRow(R1.rank(),counter);
		tmpMatrix.checkValidity();
		transposeConjugate(R1Inverse,tmpMatrix);
#ifndef NDEBUG
		// check
		SparseMatrixType tmpMatrix2;
		multiply(tmpMatrix2,R1,R1Inverse);
		bool b = isTheIdentity(tmpMatrix2);
		if (b) return;
		printFullMatrix(R1,"R1");
		printFullMatrix(R1Inverse,"R1Inverse");
		printFullMatrix(tmpMatrix2,"tmpMatrix2");
		assert(b);
#endif
	}

	const SparseMatrixType reflection() const { return reflection_; }

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

		R1.resize(ipPosOrNeg.size());
		size_t counter = 0;
		for (size_t i=0;i<R1.rank();i++) {
			R1.setRow(i,counter);

			ComplexOrRealType val = 1.0 + sector*dr[ipPosOrNeg[i]];
			RealType val2 = sqrt(2.0*std::norm(val));
			if (isAlmostZero(val2,1e-10)) continue;

			R1.pushValue(val2);
			R1.pushCol(i);
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
			if (lessOrGreaterCondition(std::real(dd[perm[i]]),lessOrGreater))
				x.push_back(perm[i]);
		}
	}

	bool lessOrGreaterCondition(const RealType& a,size_t lessOrGreater) const
	{
		if (lessOrGreater==GREATER_THAN_ZERO) {
			return (a>0.0);
		} else {
			return (a<0.0);
		}
	}

	void permute(std::vector<size_t>& iavail,const std::vector<size_t>& ip) const
	{
		for (size_t i=0;i<iavail.size();i++) iavail[i] = iavail[ip[i]];
	}

	void choleskyFactor( const std::vector<size_t>& iavail)
	{
		std::vector<ComplexOrRealType> dr(reflection_.rank());

		setDiagonal(dr,reflection_);
		setDiagonal(R1_,dr,1.0);
		setDiagonal(Rm_,dr,-1.0);

		std::ostringstream msg2;
		msg2<<"needs extra churn, iavail="<<iavail.size();
		progress_.printline(msg2,std::cout);

		SparseMatrixType reflectionT;
		transposeConjugate(reflectionT,reflection_);
		for (size_t i=0;i<iavail.size();i++) {
//			size_t ilast = i;
			size_t j = iavail[i];
			if (doneOneSector(i,j,R1_,1.0,reflectionT)) break;
			if (doneOneSector(i,j,Rm_,-1.0,reflectionT)) break;
		}
		std::ostringstream msg;
		msg<<"R1.rank="<<R1_.rank()<<" Rm.rank="<<Rm_.rank();
		progress_.printline(msg,std::cout);
	}

	bool doneOneSector(size_t i,size_t j,SparseMatrixType& R1,const RealType& sector,const SparseMatrixType& reflectionT)
	{
		bool done = (ipPos_.size()+ipNeg_.size() >= reflection_.rank());
		if (done) return true;

		RealType tol = 1e-5;

		// try to add vector to (sector) space, where sector= + or -
		std::vector<ComplexOrRealType> w(reflection_.rank(),0.0);
		setColumn(w,sector,j,reflectionT);
		std::vector<ComplexOrRealType> T1w(R1.rank(),0);
		std::vector<size_t>& ipPosOrNeg = (sector>0) ? ipPos_ : ipNeg_;
		setT1w(T1w,ipPosOrNeg,w,sector,reflectionT);
		std::vector<ComplexOrRealType> r(R1.rank());
		SparseMatrixType R1t;
		transposeConjugate(R1t,R1);
		linearSolverTriangular(r,R1t,T1w);
		ComplexOrRealType rkk2 =w*w - r*r; // note: operator* will conjugate if needed
		if (std::norm(rkk2)>tol) {
			// accept this column
			if (idebug_) {
				std::cerr<<__FILE__<<" "<<__LINE__<<" sector="<<sector;
				std::cerr<<" i="<<i<<" j="<<j<<" #pos="<<(ipPosOrNeg.size()-1)<<"\n";
			}
			growOneRowAndOneColumn(R1,r,sqrt(rkk2),sector);
			ipPosOrNeg.push_back(j);
		}
		return false;
	}



	void growOneRowAndOneColumn(SparseMatrixType& R1,
				    const std::vector<ComplexOrRealType>& r,
				    const ComplexOrRealType& addedValue,
				    const RealType& sector) const
	{
		size_t n2 = (sector>0) ? ipPos_.size() : ipNeg_.size();
		size_t n = R1.rank();
		SparseMatrixType R1new(n+1,n+1);
		size_t counter = 0;
		for (size_t i=0;i<n;i++) {
			R1new.setRow(i,counter);
			for (int k = R1.getRowPtr(i);k<R1.getRowPtr(i+1);k++) {
				R1new.pushCol(R1.getCol(k));
				R1new.pushValue(R1.getValue(k));
				counter++;
			}
			// add extra column
			if (isAlmostZero(r[i],1e-10)) continue;
			R1new.pushCol(n2);
			R1new.pushValue(r[i]);
			counter++;
		}
		// add extra row and value
		R1new.setRow(n,counter);
		R1new.pushCol(n2);
		R1new.pushValue(addedValue);
		counter++;
		R1new.setRow(n+1,counter);
		R1new.checkValidity();
		R1 = R1new;
	}

	/**
	   Let R1t = transpose(R1), solve   L * r = rhs
	*/
	void linearSolverTriangular(std::vector<ComplexOrRealType>& r,
				    const SparseMatrixType& R1t,
				    const std::vector<ComplexOrRealType>& rhs) const
	{

		for(size_t irow=0; irow < R1t.rank(); irow++) {
			ComplexOrRealType dsum = 0.0;
			ComplexOrRealType diag = 0.0;
			for(int k=R1t.getRowPtr(irow); k < R1t.getRowPtr(irow+1); k++) {
				size_t j = R1t.getCol(k);
				ComplexOrRealType lij = R1t.getValue(k);
				if (j==irow) { // treat diagonal different
					diag = lij; // save diagonal in diag
					continue; // and don't sum it
				}
				dsum += lij * r[j];
			};
			r[irow] = (rhs[irow] - dsum) / diag; //<<<< you might store the inverse if you wish to avoid costly divisions
		};
	}

	void setT1wOld(std::vector<ComplexOrRealType>& T1w,
		    const std::vector<size_t>& ipPosOrNeg,
		    const std::vector<ComplexOrRealType>& w,
		    const RealType& sector) const
	{
		for (size_t ii=0;ii<ipPosOrNeg.size();ii++) {
			size_t i = ipPosOrNeg[ii];
			for (int k = reflection_.getRowPtr(i);k<reflection_.getRowPtr(i+1);k++) {
				size_t col = reflection_.getCol(k);
				ComplexOrRealType val = reflection_.getValue(k);
				//if (col==i) val += sector;
				val *= sector;
				T1w[col] += val*w[ii];
			}
			T1w[i] += w[i];
		}
	}

	void setT1w(std::vector<ComplexOrRealType>& T1w,
		     const std::vector<size_t>& ipPosOrNeg,
		     const std::vector<ComplexOrRealType>& w,
		     const RealType& sector,
		    const SparseMatrixType& reflectionT) const
	{
		size_t n = reflectionT.rank();
		std::vector<int> inverseP(n,-1);
		for (size_t ii=0;ii<ipPosOrNeg.size();ii++)
			inverseP[ipPosOrNeg[ii]]=ii;

		for (size_t col=0;col<T1w.size();col++) {
			ComplexOrRealType sum = 0.0;
			for (int k = reflectionT.getRowPtr(col);k<reflectionT.getRowPtr(col+1);k++) {
				int x = inverseP[reflectionT.getCol(k)];
				if (x<0) continue;
				sum += reflectionT.getValue(k)*w[x];
			}
			T1w[col] = w[col] + sector * sum;
		}

	}

	void setColumn(std::vector<ComplexOrRealType>& w,const RealType& sector,size_t j,const SparseMatrixType& reflectionT) const
	{
		for (int k = reflectionT.getRowPtr(j);k<reflectionT.getRowPtr(j+1);k++) {
			size_t col = reflectionT.getCol(k);
			w[col] = reflection_.getValue(k);
			if (col==j) w[col]+=sector;
			w[col] *= sector;
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

	PsimagLite::ProgressIndicator progress_;
	const SparseMatrixType& reflection_;
	bool idebug_;
	std::vector<size_t> ipPos_,ipNeg_;
	SparseMatrixType R1_,Rm_;

}; // class ReflectionBasis

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_BASIS_H
