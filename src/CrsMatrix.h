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

/*! \file CrsMatrix.h
 *
 *  A class to represent a sparse matrix in Compressed Row Storage
 *
 */

#ifndef CRSMATRIX_HEADER_H
#define CRSMATRIX_HEADER_H
#include <algorithm>
#include "BLAS.h" 
#include "Matrix.h" 
#include "Complex.h"
#include <cassert>

namespace PsimagLite {

	//! A Sparse Matrix in Compressed Row Storage (CRS) format.
	/** 
	The CRS format puts the subsequent nonzero elements of the matrix rows
	in contiguous memory locations. We create 3 vectors: one for complex numbers containing the values of the
	matrix entries 
	and the other two for integers ($colind$ and $rowptr$).
	The vector $values$ stores the values of the non-zero elements of the matrix,
	as they are traversed in a row-wise fashion.
	The $colind$ vector stores the column indices of the elements of the $values$
	vector. That is, if $values[k] = a[i][j]$ then $colind[k] = j$.
	The $rowptr$ vector stores the locations in the $values$ vector that start
	a row, that is $values[k] = a[i][j]$ if $rowptr[i] \le i < rowptr[i + 1]$.
	By convention, we define $rowptr[N_{dim}]$ to be equal to the number of non-zero elements,
	$n_z$, in the matrix. The storage savings of this approach are significant since instead of
	storing $N_{dim}^2$ elements, we need only $2n_z + N_{dim} + 1$ storage locations.\\
	To illustrate how the CRS format works, consider the non-symmetric matrix defined by
	\begin{equation}
		A=\left[\begin{tabular}{llllll}
  
		10 &  0 & 0 & 0  & -2 & 0 \\
		3 &  9 &  0 &  0 &  0 &  3 \\
		0 &  7 &  8 &  7 &  0 &  0 \\
		3 &  0 &  8 &  7  & 5 &  0 \\
		0 &   8 &  0 &  9 &  9 & 13 \\
		0 &  4 &  0 &  0 &  2&  -1 \\
	\end{tabular}\right]\end{equation}
	The CRS format for this matrix is then specified by the arrays:\\
	\begin{tt}
		values = [10 -2  3  9  3  7  8  7  3 ... 9 13  4  2 -1 ]\\
		colind = [ 0  4  0  1  5  1  2  3  0 ... 4  5  1  4  5 ]\\
		rowptr = [ 0  2  5  8 12 16 19 ]\\
	\end{tt}
	*/
	template<class T>
	class CrsMatrix {
    public:
    	typedef T MatrixElementType;
		typedef T value_type;
		
		CrsMatrix() { } 

		~CrsMatrix() {  }

		CrsMatrix(int nrow,int ncol)
		{
			if (nrow!=ncol) throw std::runtime_error("CrsMatrix: only works for square matrices.\n");
			resize(nrow);
		}

		template<typename S>
		CrsMatrix(const CrsMatrix<S>& a)
		{
			colind_=a.colind_;
			rowptr_=a.rowptr_;
			values_=a.values_;
			size_ = a.size_;
			/*
			resize(a.getSize(),a.nonZero());
			for (int i=0;i<colind_.size();i++) colind_[i]=a.getCol(i);
			for (int i=0;i<rowptr_.size();i++) rowptr_[i]=a.getRowPtr(i);
			for (int i=0;i<values_.size();i++) values_[i]=a.getValue(i);
			*/
		}

		template<typename S>
		CrsMatrix(const CrsMatrix<std::complex<S> >& a)
		{
			colind_=a.colind_;
			rowptr_=a.rowptr_;
			values_=a.values_;
			size_ = a.size_;
			/*resize(a.getSize(),a.nonZero());
			for (int i=0;i<colind_.size();i++) colind_[i]=a.getCol(i);
			for (int i=0;i<rowptr_.size();i++) rowptr_[i]=a.getRowPtr(i);
			for (int i=0;i<values_.size();i++) values_[i]=real(a.getValue(i));
			*/
		}

		CrsMatrix(const Matrix<T>& a)
		{
			//int i,j;
			//int n=a.n_row();
			int counter=0;
			double eps = 0;
			//if (n!=int(a.n_col())) throw std::runtime_error("CrsMatrix: only works for square matrices.\n");

			resize(a.n_row());	

			for (size_t i = 0; i < a.n_row(); i++) {
				setRow(i,counter);
				for (size_t j=0;j<a.n_col();j++) {
					if (std::norm(a(i,j))<=eps) continue;
					pushValue(a(i,j));
					pushCol(j);
					counter++;
				}

			}
			setRow(a.n_row(),counter);
		}

		void resize(int nsize)
		{
			colind_.clear();
			values_.clear();
			rowptr_.resize(nsize+1);
			size_=nsize;
		}

		void clear()
		{
			colind_.clear();
			values_.clear();
			rowptr_.clear();
			size_=0;
		}

		void resize(int nsize,int nonzero) 
		{
			resize(nsize);
			colind_.resize(nonzero);
			values_.resize(nonzero);
		}

		void setRow(int n,int v) 
		{
			rowptr_[n]=v;
		}

		void setCol(int n,int v) {
			colind_[n]=v;
		}

		void setValues(int n,const T &v) {
			values_[n]=v;
		}

		void operator*=(T x)
		{
			for (size_t i=0;i<values_.size();i++) values_[i] *= x;
		}

		template<typename VerySparseMatrixType>
		void operator=(const VerySparseMatrixType& m)
		{
			resize(m.rank());
			size_t counter=0;

			for (size_t i=0;i<m.rank();i++) {
				setRow(i,counter);
				size_t counter2=0;
				for (size_t j=counter;j<m.nonZero();j++) {
					if (m.getRow(j)!=i) break;
					pushCol(m.getColumn(j));
					pushValue(m.getValue(j));
					counter2++;
				}
				counter+=counter2;
				
			}
			setRow(m.rank(),counter);
		}

		void operator+=(CrsMatrix<T> const &m) 
		{
			CrsMatrix<T> c;
			if (size_t(this->size_)>=m.rank()) operatorPlus(c,*this,m);
			else operatorPlus(c,m,*this);
			*this =c;
		}

		T operator()(int i,int j) const 
		{
			for (int k=rowptr_[i];k<rowptr_[i+1];k++) if (colind_[k]==j) return values_[k];
			return static_cast<T>(0.0);
		}

		int nonZero() const { return colind_.size(); }		

		/** performs x = x + A * y
		 ** where x and y are vectors and A is a sparse matrix in
		 ** row-compressed format */
		template<typename S>
		void matrixVectorProduct (std::vector<S> &x, std::vector<S> const &y) const
		{ 
			for (size_t i = 0; i < y.size(); i++)
				for (int j = rowptr_[i]; j < rowptr_[i + 1]; j++)
					x[i] += values_[j] * y[colind_[j]];
		}
		

		size_t rank() const { return size_; }

		void pushCol(int i) { colind_.push_back(i); }

		void pushValue(T const &value) { values_.push_back(value); }

		//! Make a diagonal CRS matrix with value "value"
		void makeDiagonal(int rank,T const &value=0) 
		{
			size_=rank;
			rowptr_.resize(rank+1);
			values_.resize(rank);
			colind_.resize(rank);
			
			for (int i=0;i<rank;i++) {
				values_[i]=value;
				colind_[i]=i;
				rowptr_[i]=i;
			}
			rowptr_[rank]=rank;
		}

		const int& getRowPtr(size_t i) const { assert(i<rowptr_.size()); return rowptr_[i]; }

		const int& getCol(int i) const { return colind_[i]; }

		const T& getValue(int i) const { return values_[i]; }

		/*bool operator==(const CrsMatrix<T>& B) const
		{
			if (!utils::vectorEqual(values_,B.values_)) return false;
			if (!utils::vectorEqual(colind_,B.colind_)) return false;
			if (!utils::vectorEqual(rowptr_,B.rowptr_)) return false;
			return true;
		}*/

		void set(const std::vector<int> &rowptr,const std::vector<int>& colind,const std::vector<T>& values)
		{
			rowptr_=rowptr;
			colind_=colind;
			values_=values;
			size_=rowptr.size()-1;
		}
		
		template<typename S>
		friend std::ostream &operator<<(std::ostream &os,const CrsMatrix<S> &m);
		
		template<typename S,typename S2>
		friend void multiplyScalar(CrsMatrix<S> &ret,CrsMatrix<S> const &s,S2 const &v);
		
		template<class S>
		friend void difference(const CrsMatrix<S>& A,const CrsMatrix<S>& B);

		template<typename S>
		friend void MpiBroadcast(CrsMatrix<S> *v,int rank);
	
		template<typename S>
		friend void MpiSend(CrsMatrix<S>  *v,int iproc,int i);
	
		template<typename S>
		friend void MpiRecv(CrsMatrix<S> *v,int iproc,int i);
		
		template<typename CrsMatrixType>
		friend std::istream &operator>>(std::istream &is,CrsMatrix<CrsMatrixType>& m);

	private:
		std::vector<int> rowptr_;
		std::vector<int> colind_;
		std::vector<T> values_;
		int size_;
	}; // class CrsMatrix

	// Companion functions below:

	template<typename T>
	std::ostream &operator<<(std::ostream &os,const CrsMatrix<T> &m)
	{
		size_t n=m.rank();
		if (n==0) return os;
		os<<n<<"\n";
		for (size_t i=0;i<n+1;i++) os<<m.rowptr_[i]<<" ";
		os<<"\n";
		
		size_t nonzero=m.nonZero();
		os<<nonzero<<"\n";
		for (size_t i=0;i<nonzero;i++) os<<m.colind_[i]<<" ";
		os<<"\n";
		
		os<<nonzero<<"\n";
		for (size_t i=0;i<nonzero;i++) os<<m.values_[i]<<" ";
		os<<"\n";
		
		return os;
	}

	template<typename T>
	std::istream &operator>>(std::istream &is,CrsMatrix<T>& m)
	{
		int n;
		is>>n;
		if (n<0) throw std::runtime_error(
				"is>>CrsMatrix(...): Rank must be positive\n");
		m.resize(n);
		for (size_t i=0;i<size_t(n)+1;i++) is>>m.rowptr_[i];

		size_t nonzero;
		is>>nonzero;
		m.colind_.resize(nonzero);
		for (size_t i=0;i<nonzero;i++) is>>m.colind_[i];

		is>>nonzero;
		m.values_.resize(nonzero);
		for (size_t i=0;i<nonzero;i++) is>>m.values_[i];

		return is;
	}

	//! Transforms a Compressed-Row-Storage (CRS) into a full Matrix (Fast version)
	template<typename T,class FullMatrixTemplate>
	inline void crsMatrixToFullMatrix(FullMatrixTemplate& m,const CrsMatrix<T>& crsMatrix)
	{
		int i,k;
		int size=crsMatrix.rank();

		m.reset(size,size);
		for (i = 0; i < size ; i++) {
			for (k=0;k<size;k++) m(i,k)=0;
			for (k=crsMatrix.getRowPtr(i);k<crsMatrix.getRowPtr(i+1);k++) 
				m(i,crsMatrix.getCol(k))=crsMatrix.getValue(k);
		}

	}

	//! Transforms a full matrix into a Compressed-Row-Storage (CRS) Matrix
	// Use the constructor if possible
 	template<typename T,class FullMatrixTemplate>
	void fullMatrixToCrsMatrix(CrsMatrix<T>& crsMatrix, const FullMatrixTemplate& a)
	{
		int i,j;
		int n=a.n_row();
		int counter=0;
		
		if (n!=int(a.n_col())) throw std::runtime_error("fullMatrixToCrsMatrix: only works for square matrices.\n");
			
		crsMatrix.resize(n);	
		
		for (i = 0; i < n; i++) {
			crsMatrix.setRow(i,counter);
			for (j=0;j<n;j++) {
				if (a(i,j)==static_cast<T>(0)) continue;
				crsMatrix.pushValue(a(i,j));
				crsMatrix.pushCol(j);
				counter++;
			}
			
		}
		crsMatrix.setRow(n,counter);
	} 

	/** If order==false then 
	    creates B such that B_{i1+j1*nout,i2+j2*nout)=A(j1,j2)\delta_{i1,i2}
	    if order==true then
	    creates B such that B_{i1+j1*na,i2+j2*na)=A(i1,i2)\delta_{j1,j2}
	    where na=rank(A)
	  */
	template<class T>
	void externalProduct(CrsMatrix<T>  &B,CrsMatrix<T> const &A,int nout,
		std::vector<double> const &signs,bool order=true)
	{
		int na=A.rank();
		T tmp;

		B.resize(na*nout); //,nonzero*nout);

		int i,ii,jj,alpha,k,j,beta;
		int counter=0;
		for (ii=0;ii<na*nout;ii++) {
			if (order) {
				// ii = i+alpha*na;
				alpha = int(ii/na);
				i = ii-alpha*na;
			} else {
				//ii = alpha + i*nout;
				i = int(ii/nout);
				alpha = ii - i*nout;
			}
			B.setRow(ii,counter);
			for (k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
				j = A.getCol(k);
				beta=alpha;
				if (order) jj = j+beta*na;
				else       jj = beta+j*nout;
				//B.setCol(counter,jj);
				B.pushCol(jj);
				tmp = A.getValue(k);
				if (!order) tmp*=signs[alpha];
				//B.setValues(counter,tmp);
				B.pushValue(tmp);
				counter++;
			}
		}
		B.setRow(na*nout,counter);
	}

	//! Computes C = A external product B
	template<class T>
	void externalProduct(CrsMatrix<T>  &C,CrsMatrix<T> const &A,CrsMatrix<T> const &B)
	{
		int n=A.getSize()*B.getSize();
		C.resize(n);
		int na = A.getSize();
		T tmp;
		int i,k,kk,alpha,beta,j,counter=0;

		for (i=0;i<n;i++) {
			C.setRow(i,counter);
			// i = alpha + beta * na
			beta = int(i/na);
			alpha = i - beta * na;
			for (k=A.getRowPtr(alpha);k<A.getRowPtr(alpha+1);k++) {
				for (kk=B.getRowPtr(beta);kk<B.getRowPtr(beta+1);kk++) {
					j = A.getCol(k) + B.getCol(kk) *na;
					C.pushCol(j);
					tmp = A.getValue(k) * B.getValue(kk);
					C.pushValue(tmp);
					counter++;
				}
			}
		}
		C.setRow(n,counter);	
	}

	//! Computes C = A external product B (with signs)
	template<class T>
	void externalProduct(CrsMatrix<T>  &C,CrsMatrix<T> const &A,CrsMatrix<T> const &B,const std::vector<int>& signs,bool option=false)
	{
		int n=A.getSize()*B.getSize();
		C.resize(n);
		int na = A.getSize();
		T tmp;
		int i,k,kk,alpha,beta,j,counter=0;

		for (i=0;i<n;i++) {
			C.setRow(i,counter);
			// i = alpha + beta * na
			beta = int(i/na);
			alpha = i - beta * na;
			for (k=A.getRowPtr(alpha);k<A.getRowPtr(alpha+1);k++) {
				for (kk=B.getRowPtr(beta);kk<B.getRowPtr(beta+1);kk++) {
					j = A.getCol(k) + B.getCol(kk) *na;
					C.pushCol(j);
					int sign = signs[alpha];
					if (option) sign=signs[beta];
					tmp = A.getValue(k) * B.getValue(kk)*sign;
					C.pushValue(tmp);
					counter++;
				}
			}
		}
		C.setRow(n,counter);	
	}

	//! Sets ret = s * v where ret and s are CRS matrices and v is a scalar number
	template<typename S,typename T>
	void multiplyScalar(CrsMatrix<S> &ret,CrsMatrix<S> const &s,T const &v)
	{
		ret = s;
		for (size_t ii=0;ii<s.values_.size();ii++) {
			ret.values_[ii] *= v;
		}
	}

	//! C = A*B,  all matrices are CRS matrices
	//! idea is from http://web.maths.unsw.edu.au/~farid/Papers/Hons/node23.html
	template<typename S,typename S2>
	void multiply(CrsMatrix<S> &C,CrsMatrix<S> const &A,CrsMatrix<S2> const &B)
	{
		int j,k,s,mlast,itemp,jbk;
		size_t n = A.rank();
		std::vector<int> ptr(n),index(n);
		std::vector<S> temp(n);
		S tmp;

		if (n!=B.rank()) throw std::runtime_error("multiply: matrices must have the same rank.\n");
		
		C.resize(n);
		
		// mlast pointer to the last place we updated in the C vector 
		mlast = 0;
		for (size_t l=0;l<n;l++) ptr[l] = -1;
		// over the rows of A
		for (size_t i=0;i<n;i++) {
			C.setRow(i,mlast);
			// start calculations for row 
			itemp = 0;
			for(j = A.getRowPtr(i);j< A.getRowPtr(i+1);j++) { 
				for (k = B.getRowPtr(A.getCol(j)); k< B.getRowPtr(A.getCol(j)+1);k++) {
					jbk=B.getCol(k);
					tmp = A.getValue(j)*B.getValue(k);
					if( ptr[jbk]== -1) {
						ptr[jbk] = itemp;
						temp[ptr[jbk]] = tmp ;
						index[ptr[jbk]] = jbk;
						itemp++;
					} else  {
						temp[ptr[jbk]]= temp[ptr[jbk]]+tmp;
			   		}
				}
			}
			// before you leave this row update array c , jc 	
			for (s=0;s<itemp;s++) {
				C.pushValue(temp[s]);
				C.pushCol(index[s]);
				ptr[index[s]]= -1;
			}
			mlast += itemp;
		}
		C.setRow(n,mlast);
	}

	// vector2 = sparseMatrix * vector1
	template<class S>
	void multiply(std::vector<S>& v2, const CrsMatrix<S>& m, const std::vector<S>& v1)
	{
		int i,j,n=m.getSize();
		for (i=0;i<n;i++) {
			v2[i] =0;
			for (j=m.getRowPtr(i);j<m.getRowPtr(i+1);j++) {
				v2[i] += m.getValue(j)*v1[m.getCol(j)];
			}
		}
	}

	//! Sets B=transpose(conjugate(A))	
	template<typename S,typename S2>
	inline void transposeConjugate(CrsMatrix<S>  &B,CrsMatrix<S2> const &A)
	{
		int n=A.rank();
		int i,k;
		std::vector<std::vector<int> > col(n);
		std::vector<std::vector<S2> > value(n);

		// B(j,i) = conj(A(i,j))
		for (i=0;i<n;i++) {
			for (k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
				col[A.getCol(k)].push_back(i);
				S2 w = A.getValue(k);
				value[A.getCol(k)].push_back(w);
			}
		}

		B.resize(n); 

		int counter=0;
		for (i=0;i<n;i++) {
			
			B.setRow(i,counter);
			for (size_t j=0;j<col[i].size();j++) {
				if (value[i][j]==static_cast<S>(0.0)) continue;
				B.pushCol(col[i][j]);
				B.pushValue(std::conj(value[i][j]));
				counter++;
			}
		}
		B.setRow(n,counter);
	}

// 	template<class S>
// 	inline CrsMatrix<S> transposeConjugate(const CrsMatrix<S>& m)
// 	{
// 		CrsMatrix<S> temp;
// 		transposeConjugate(temp,m);
// 		return temp;
// 	}	

	//! Sets A = B(i,perm(j)), A and B CRS matrices	
	template<class S>
	void permute(CrsMatrix<S> &A,CrsMatrix<S> const &B,std::vector<size_t> const &perm)
	{
		int i,n = B.rank();

		A.resize(n);
		S tmp;

		std::vector<int> permInverse(n);
		for (i=0;i<n;i++) permInverse[perm[i]]=i;

		int k,counter=0;
		for (i=0;i<n;i++) {
			A.setRow(i,counter);
			for (k=B.getRowPtr(i);k<B.getRowPtr(i+1);k++) {
				A.pushCol(permInverse[B.getCol(k)]);
				tmp = B.getValue(k);
				A.pushValue(tmp);
				counter++;
			}
		}
		A.setRow(n,counter);
	}

	//! Sets A = B(perm(i),j), A and B CRS matrices		
	template<class S>
	void permuteInverse(CrsMatrix<S> &A,CrsMatrix<S> const &B,std::vector<size_t> const &perm)
	{
		int ii,n = B.rank();
		A.resize(n); // B.nonZero());
		S tmp;

		int k,counter=0;
		for (int i=0;i<n;i++) {
			ii = perm[i];
			A.setRow(i,counter);
			for (k=B.getRowPtr(ii);k<B.getRowPtr(ii+1);k++) {
				A.pushCol(B.getCol(k));
				tmp = B.getValue(k);
				A.pushValue(tmp);
				counter++;
			}
		}
		A.setRow(n,counter);
	}

	//! Sets A = B^\dagger * S * B
	template<class T>
	inline Matrix<T> transformFullFast(CrsMatrix<T> const &S,Matrix<T> const &fmB)
	{
		int nBig = S.rank();
		int nSmall = fmB.n_col();
		double alpha=1.0;
		double beta=0.0;
		Matrix<T> fmS,fmTmp(nBig,nSmall);
		
		crsMatrixToFullMatrix(fmS,S);

		psimag::BLAS::GEMM('N','N',nBig,nSmall,nBig,alpha,&(fmS(0,0)),nBig,&(fmB(0,0)),nBig,beta,&(fmTmp(0,0)),nBig);
		fmS.reset(nSmall,nSmall);
		psimag::BLAS::GEMM('C','N',nSmall,nSmall,nBig,alpha,&(fmB(0,0)),nBig,&(fmTmp(0,0)),nBig,beta,&(fmS(0,0)),nSmall);
		return fmS;
	}

	//! Sets A=B+C, restriction: B.size has to be larger or equal than C.size
	template<class T>
	void operatorPlus(CrsMatrix<T> &A,CrsMatrix<T> const &B,CrsMatrix<T> const &C)
	{
		size_t n = B.rank();
		T tmp;

		if (n<C.rank()) {
			std::cerr<<"B.size="<<B.rank()<<" C.size="<<C.rank()<<"\n";
			throw std::runtime_error("CrsMatrix<T>: operatorPlus: B.size must be larger or equal than C.size.\n");
		}

		std::vector<T>  valueTmp(n);
		std::vector<int> index;
		A.resize(n);

		size_t counter=0;
		for (size_t k2=0;k2<n;k2++) valueTmp[k2]= static_cast<T>(0.0);
		
		for (size_t i = 0; i < n; i++) {
			int k;
			A.setRow(i,counter);

			if (i<C.rank()) {
				// inspect this
				index.clear();
				for (k=B.getRowPtr(i);k<B.getRowPtr(i+1);k++) {
					if (B.getCol(k)<0 || size_t(B.getCol(k))>=n) throw std::runtime_error("operatorPlus (1)\n");
					valueTmp[B.getCol(k)]=B.getValue(k);
					index.push_back(B.getCol(k));
				}

				// inspect C 
				for (k=C.getRowPtr(i);k<C.getRowPtr(i+1);k++) {
					tmp = C.getValue(k);
					if (C.getCol(k)>=int(valueTmp.size()) || C.getCol(k)<0) throw std::runtime_error("operatorPlus (2)\n");

					valueTmp[C.getCol(k)] += tmp;
					index.push_back(C.getCol(k));
				}
				std::sort(index.begin(),index.end());
				k= -1;
				for (size_t kk=0;kk<index.size();kk++) {
					if (k==index[kk]) continue;
					k=index[kk];
					if (k<0 || size_t(k)>=n) throw std::runtime_error("operatorPlus (3)\n");
					tmp = valueTmp[k];
					if (tmp!=static_cast<T>(0.0)) {
						A.pushCol(k);
						A.pushValue(tmp);
						counter++;
						valueTmp[k]=static_cast<T>(0.0);
					}
				}
			} else {
				for (k=B.getRowPtr(i);k<B.getRowPtr(i+1);k++) {
					tmp = B.getValue(k);
					A.pushCol(B.getCol(k));
					A.pushValue(tmp);
					counter++;
				}
			}
		}
		A.setRow(n,counter);
	}

	template<typename T>
	bool isHermitian(const CrsMatrix<T>& A,bool doThrow=false)
	{
		for (size_t i=0;i<A.rank();i++) {
			for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
				if (std::norm(A.getValue(k)-std::conj(A(A.getCol(k),i)))<1e-6) continue;
				assert(false);
				return false;
			}
		}
		return true;
	}

	template<class T>
	void sumBlock(CrsMatrix<T> &A,CrsMatrix<T> const &B,size_t offset)
	{
		int counter=0;
		CrsMatrix<T> Bfull;
		
		Bfull.resize(A.rank());
		for (size_t i=0;i<offset;i++) Bfull.setRow(i,counter);
		for (size_t ii=0;ii<B.rank();ii++) {
			size_t i=ii+offset;
			Bfull.setRow(i,counter);
			for (int jj=B.getRowPtr(ii);jj<B.getRowPtr(ii+1);jj++) {
				size_t j = B.getCol(jj)+offset;
				T tmp  = B.getValue(jj);
				Bfull.pushCol(j);
				Bfull.pushValue(tmp);
				counter++;
			}
		}
		//if (i>A.rank()) throw std::runtime_error("sumBlock\n");
		for (size_t i=B.rank()+offset;i<A.rank();i++) Bfull.setRow(i,counter);
		Bfull.setRow(A.rank(),counter);
		A += Bfull;
	}
	
	template<class T>
	bool isDiagonal(const CrsMatrix<T>& A) 
	{
		size_t n = A.getSize();
		for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) if (i!=j && fabs(A(i,j))>1e-6) return false;	
		return true;
	}
	
// 	template<class T>
// 	void difference(const CrsMatrix<T>& A,const CrsMatrix<T>& B)
// 	{
// 		std::cerr<<"rowptr differences:\n";
// 		utils::difference(A.rowptr_,B.rowptr_);
// 		
// 		std::cerr<<"colind differences:\n";
// 		utils::difference(A.colind_,B.colind_);
// 		
// 		std::cerr<<"values differences:\n";
// 		utils::difference(A.values_,B.values_);
// 	}
	
	template<class T>
	bool isTheIdentity(const CrsMatrix<T>& A) 
	{
		T eps = 1e-6;
		size_t n = A.getSize();
		for (size_t i=0;i<n;i++) {
			if (fabs(A(i,i)-1.0)>eps) {
				std::cerr<<"Diagonal is "<<A(i,i)<<" at i="<<i<<"\n";
				return false; 
			}
			for (size_t j=0;j<n;j++) {
				if (i!=j && fabs(A(i,j))>eps) {
					std::cerr<<"A("<<i<<","<<j<<")="<<A(i,j)<<"\n";
					return false;
				}
			}
		}
		return true;
	}	

	template<typename T>
	Matrix<T> multiplyTc(const CrsMatrix<T>& a,const CrsMatrix<T>& b)
	{
		
		CrsMatrix<T> bb,c;
		transposeConjugate(bb,b);
		multiply(c,a,bb);
		Matrix<T> cc;
		crsMatrixToFullMatrix(cc,c);
		return cc;
	}

	template<typename T1,typename T2>
	CrsMatrix<T2> operator*(const T1& t1,const CrsMatrix<T2>& a)
	{
		CrsMatrix<T2> b;
		multiplyScalar(b,a,t1);
		return b;
	}

	template<typename T1,typename T2>
	CrsMatrix<T1> operator*(const CrsMatrix<T1>& a,const CrsMatrix<T2>& b)
	{
		CrsMatrix<T2> c;
		multiply(c,a,b);
		return c;
	}
		
} // namespace PsimagLite 
/*@}*/	
#endif
