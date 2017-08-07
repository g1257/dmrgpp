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
#include "loki/TypeTraits.h"
#include "Mpi.h"

namespace PsimagLite {

//! A Sparse Matrix in Compressed Row Storage (CRS) format.
/**
	The CRS format puts the subsequent nonzero elements of the matrix rows
	in contiguous memory locations. We create 3 vectors: one for complex numbers
	containing the values of the
	matrix entries
	and the other two for integers ($colind$ and $rowptr$).
	The vector $values$ stores the values of the non-zero elements of the matrix,
	as they are traversed in a row-wise fashion.
	The $colind$ vector stores the column indices of the elements of the $values$
	vector. That is, if $values[k] = a[i][j]$ then $colind[k] = j$.
	The $rowptr$ vector stores the locations in the $values$ vector that start
	a row, that is $values[k] = a[i][j]$ if $rowptr[i] \le i < rowptr[i + 1]$.
	By convention, we define $rowptr[N_{dim}]$ to be equal to the number of non-zero elements,
	$n_z$, in the matrix. The storage savings of this approach are significant
	because instead of
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

	CrsMatrix() : nrow_(0),ncol_(0) { }

	~CrsMatrix() {  }

	CrsMatrix(SizeType nrow,SizeType ncol)
	    : nrow_(nrow),ncol_(ncol)
	{
		resize(nrow,ncol);
	}

	template<typename S>
	CrsMatrix(const CrsMatrix<S>& a)
	{
		colind_=a.colind_;
		rowptr_=a.rowptr_;
		values_=a.values_;
		nrow_ = a.nrow_;
		ncol_ = a.ncol_;
	}

	template<typename S>
	CrsMatrix(const CrsMatrix<std::complex<S> >& a)
	{
		colind_=a.colind_;
		rowptr_=a.rowptr_;
		values_=a.values_;
		nrow_ = a.nrow_;
		ncol_ = a.ncol_;
	}

	explicit CrsMatrix(const Matrix<T>& a)
	{
		int counter=0;
		double eps = 0;

		resize(a.n_row(),a.n_col());

		for (SizeType i = 0; i < a.n_row(); i++) {
			setRow(i,counter);
			for (SizeType j=0;j<a.n_col();j++) {
				if (PsimagLite::norm(a(i,j))<=eps) continue;
				pushValue(a(i,j));
				pushCol(j);
				counter++;
			}

		}
		setRow(a.n_row(),counter);
	}

	// start closure ctors

	CrsMatrix(const std::ClosureOperator<CrsMatrix,
	          CrsMatrix,
	          std::ClosureOperations::OP_MULT>& c)
	{
		CrsMatrix& x = *this;
		const CrsMatrix& y = c.r1;
		const CrsMatrix& z = c.r2;
		multiply(x,y,z);
	}

	CrsMatrix(const std::ClosureOperator<T,CrsMatrix,std::ClosureOperations::OP_MULT>& c)
	{
		*this = c.r2;
		this->values_ *= c.r1;
	}

	// end all ctors

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType,
	                   String msg = "") const
	{
		String str = msg;
		str += "CrsMatrix";

		const char* start = reinterpret_cast<const char *>(this);
		const char* end = reinterpret_cast<const char *>(&colind_);
		SizeType total = mres.memResolv(&rowptr_, end-start, str + " rowptr");

		start = end;
		end = reinterpret_cast<const char *>(&values_);
		total += mres.memResolv(&colind_, end-start, str + " colind");

		start = end;
		end = reinterpret_cast<const char *>(&nrow_);
		total += mres.memResolv(&values_, end-start, str + " values");

		start = end;
		end = reinterpret_cast<const char *>(&ncol_);
		total += mres.memResolv(&nrow_, end-start, str + " nrow");

		total += mres.memResolv(&ncol_,
		                        sizeof(*this) - total,
		                        str + " ncol");

		return total;
	}

	void resize(SizeType nrow,SizeType ncol)
	{
		colind_.clear();
		values_.clear();
		rowptr_.clear();
		rowptr_.resize(nrow+1);
		nrow_ = nrow;
		ncol_ = ncol;
	}

	void clear()
	{
		colind_.clear();
		values_.clear();
		rowptr_.clear();
		nrow_=ncol_=0;
	}

	void resize(SizeType nrow,SizeType ncol,SizeType nonzero)
	{
		resize(nrow,ncol);
		colind_.resize(nonzero);
		values_.resize(nonzero);
	}

	void setRow(SizeType n,SizeType v)
	{
		assert(n<rowptr_.size());
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
		values_ *= x;
	}


	bool operator==(const CrsMatrix<T>& op) const
	{
		return (nrow_ == op.nrow_ &&
		        ncol_ == op.ncol_ &&
		        rowptr_ == op.rowptr_ &&
		        colind_ == op.colind_ &&
		        values_ == op.values_);
	}

	template<typename VerySparseMatrixType>
	typename EnableIf<!std::IsClosureLike<VerySparseMatrixType>::True,void>::Type
	operator=(const VerySparseMatrixType& m)
	{
		if (!m.sorted())
			throw RuntimeError("CrsMatrix: VerySparseMatrix must be sorted\n");

		clear();
		SizeType nonZeros = m.nonZero();
		resize(m.rows(),m.cols(),nonZeros);

		SizeType counter=0;
		for (SizeType i=0;i<m.rows();++i) {
			setRow(i,counter);

			while (counter < nonZeros && m.getRow(counter) == i) {
				colind_[counter] = m.getColumn(counter);
				values_[counter] = m.getValue(counter);
				counter++;
			}
		}

		setRow(m.rows(),counter);
		checkValidity();
	}

	void operator+=(const CrsMatrix& m)
	{
		CrsMatrix c;
		const T f1 = 1.0;
		add(c,m,f1);
		*this = c;
	}

	T element(int i,int j) const
	{
		for (int k=rowptr_[i];k<rowptr_[i+1];k++) if (colind_[k]==j) return values_[k];
		return static_cast<T>(0.0);
	}

	SizeType nonZeros() const { return colind_.size(); }

	/** performs x = x + A * y
		 ** where x and y are vectors and A is a sparse matrix in
		 ** row-compressed format */
	template<typename VectorLikeType>
	void matrixVectorProduct(VectorLikeType& x, const VectorLikeType& y) const
	{
		assert(x.size()==y.size());
		for (SizeType i = 0; i < y.size(); i++) {
			assert(i+1<rowptr_.size());
			for (int j = rowptr_[i]; j < rowptr_[i + 1]; j++) {
				assert(SizeType(j)<values_.size());
				assert(SizeType(j)<colind_.size());
				assert(SizeType(colind_[j])<y.size());
				x[i] += values_[j] * y[colind_[j]];
			}
		}
	}

#ifndef NO_DEPRECATED_ALLOWED
	int nonZero() const { return colind_.size(); } // DEPRECATED, use nonZeros()
#endif

	SizeType rows() const { return nrow_; }

	SizeType cols() const { return ncol_; }

	void pushCol(SizeType i) { colind_.push_back(i); }

	void pushValue(T const &value) { values_.push_back(value); }

	//! Make a diagonal CRS matrix with value "value"
	void makeDiagonal(SizeType row,T const &value=0)
	{
		nrow_=row;
		ncol_=row;
		rowptr_.resize(row+1);
		values_.resize(row);
		colind_.resize(row);

		for (SizeType i=0;i<row;i++) {
			values_[i]=value;
			colind_[i]=i;
			rowptr_[i]=i;
		}
		rowptr_[row]=row;
	}

	const int& getRowPtr(SizeType i) const { assert(i<rowptr_.size()); return rowptr_[i]; }

	const int& getCol(SizeType i) const { assert(i<colind_.size()); return colind_[i]; }

	const T& getValue(SizeType i) const { assert(i<values_.size()); return values_[i]; }

	Matrix<T> toDense() const
	{
		Matrix<T> m;
		crsMatrixToFullMatrix(m,*this);
		return m;
	}

	void checkValidity() const
	{
#ifndef NDEBUG
		SizeType n = nrow_;
		assert(n+1==rowptr_.size());
		assert(nrow_>0 && ncol_>0);
		for (SizeType i=0;i<n;i++) {
			typename Vector<SizeType>::Type p(ncol_,0);
			for (int k=rowptr_[i];k<rowptr_[i+1];k++) {
				SizeType col = colind_[k];
				assert(p[col]==0);
				p[col] = 1;
			}
		}
#endif
	}

	// closures operators start

	template<typename T1>
	typename EnableIf<Loki::TypeTraits<T1>::isArith || IsComplexNumber<T1>::True,
	CrsMatrix>::Type
	operator=(const std::ClosureOperator<T1,
	          CrsMatrix,
	          std::ClosureOperations::OP_MULT>& c)
	{
		*this = c.r2;
		this->values_ *= c.r1;
		return *this;
	}

	template<typename T1>
	typename EnableIf<Loki::TypeTraits<T1>::isArith || IsComplexNumber<T1>::True,
	CrsMatrix>::Type
	operator+=(const std::ClosureOperator<T1,
	           CrsMatrix,
	           std::ClosureOperations::OP_MULT>& c)
	{
		CrsMatrix s;
		add(s,c.r2,c.r1);
		*this = s;
		return *this;
	}

	template<typename T1>
	typename EnableIf<Loki::TypeTraits<T1>::isArith || IsComplexNumber<T1>::True,
	CrsMatrix>::Type
	operator+=(const std::ClosureOperator<std::ClosureOperator<T1,
	           CrsMatrix,
	           std::ClosureOperations::OP_MULT>,
	           CrsMatrix,
	           std::ClosureOperations::OP_MULT>& c)
	{
		CrsMatrix s;
		multiply(s,c.r1.r2,c.r2);
		CrsMatrix s2;
		add(s2,s,c.r1.r1);
		*this = s2;
		return *this;
	}

	// closures operators end

	void send(int root,int tag,MPI::CommType mpiComm)
	{
		MPI::send(nrow_,root,tag,mpiComm);
		MPI::send(ncol_,root,tag+1,mpiComm);
		MPI::send(rowptr_,root,tag+2,mpiComm);
		MPI::send(colind_,root,tag+3,mpiComm);
		MPI::send(values_,root,tag+4,mpiComm);
	}

	void recv(int root,int tag,MPI::CommType mpiComm)
	{
		MPI::recv(nrow_,root,tag,mpiComm);
		MPI::recv(ncol_,root,tag+1,mpiComm);
		MPI::recv(rowptr_,root,tag+2,mpiComm);
		MPI::recv(colind_,root,tag+3,mpiComm);
		MPI::recv(values_,root,tag+4,mpiComm);
	}

	friend bool isZero(const CrsMatrix& A, double eps = 0.0)
	{
		SizeType n = A.values_.size();
		for (SizeType i = 0; i < n; ++i) {
			if (std::abs(A.values_[i]) > eps)
				return false;
		}

		return true;
	}

	template<typename S>
	friend typename Real<S>::Type norm2(const CrsMatrix<S>& m);

	template<typename S>
	friend std::ostream &operator<<(std::ostream &os,const CrsMatrix<S> &m);

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

	template<typename S>
	friend void bcast(CrsMatrix<S>& m);

private:

	template<typename T1>
	void add(CrsMatrix<T>& c, const CrsMatrix<T>& m, const T1& t1) const
	{
		assert(m.rows() == m.cols());
		const T1 one = 1.0;
		if (nrow_>=m.rows()) operatorPlus(c,*this,one,m,t1);
		else operatorPlus(c,m,t1,*this,one);
	}

	//serializr start class CrsMatrix
	//serializr normal rowptr_
	typename Vector<int>::Type rowptr_;
	//serializr normal colind_
	typename Vector<int>::Type colind_;
	//serializr normal values_
	typename Vector<T>::Type values_;
	//serializr normal nrow_
	SizeType nrow_;
	//serializr normal ncol_
	SizeType ncol_;
}; // class CrsMatrix

// Companion functions below:

template<typename T>
std::ostream &operator<<(std::ostream &os,const CrsMatrix<T> &m)
{
	SizeType n=m.rows();
	if (n==0) {
		os<<"0 0\n";
		return os;
	}

	os<<n<<" "<<m.cols()<<"\n";
	for (SizeType i=0;i<n+1;i++) os<<m.rowptr_[i]<<" ";
	os<<"\n";

	SizeType nonzero=m.nonZero();
	os<<nonzero<<"\n";
	for (SizeType i=0;i<nonzero;i++) os<<m.colind_[i]<<" ";
	os<<"\n";

	os<<nonzero<<"\n";
	for (SizeType i=0;i<nonzero;i++) os<<m.values_[i]<<" ";
	os<<"\n";

	return os;
}

template<typename T>
std::istream &operator>>(std::istream &is,CrsMatrix<T>& m)
{
	int n;
	is>>n;
	if (n<0) throw RuntimeError(
	            "is>>CrsMatrix(...): Rows must be positive\n");

	int ncol=0;
	if (ncol<0) throw RuntimeError(
	            "is>>CrsMatrix(...): Cols must be positive\n");
	is>>ncol;

	if (n == 0 || ncol == 0) return is;

	m.resize(n,ncol);
	for (SizeType i=0;i<m.rowptr_.size();i++) is>>m.rowptr_[i];

	SizeType nonzero;
	is>>nonzero;
	m.colind_.resize(nonzero);
	for (SizeType i=0;i<m.colind_.size();i++) is>>m.colind_[i];

	is>>nonzero;
	m.values_.resize(nonzero);
	for (SizeType i=0;i<m.values_.size();i++) is>>m.values_[i];

	return is;
}

template<typename T>
class IsMatrixLike<CrsMatrix<T> > {
public:
	enum { True = true};
};

template<typename S>
void bcast(CrsMatrix<S>& m)
{
	MPI::bcast(m.rowptr_);
	MPI::bcast(m.colind_);
	MPI::bcast(m.values_);
	MPI::bcast(m.nrow_);
	MPI::bcast(m.ncol_);
}

//! Transforms a Compressed-Row-Storage (CRS) into a full Matrix (Fast version)
template<typename T>
void crsMatrixToFullMatrix(Matrix<T>& m,const CrsMatrix<T>& crsMatrix)
{
	m.reset(crsMatrix.rows(),crsMatrix.cols());
	for (SizeType i = 0; i < crsMatrix.rows() ; i++) {
		for (SizeType k=0;k<crsMatrix.cols();k++) m(i,k)=0;
		for (int k=crsMatrix.getRowPtr(i);k<crsMatrix.getRowPtr(i+1);k++)
			m(i,crsMatrix.getCol(k))=crsMatrix.getValue(k);
	}

}

//! Transforms a full matrix into a Compressed-Row-Storage (CRS) Matrix
// Use the constructor if possible
template<typename T>
void fullMatrixToCrsMatrix(CrsMatrix<T>& crsMatrix, const Matrix<T>& a)
{
	crsMatrix.resize(a.n_row(),a.n_col());

	SizeType counter = 0;
	for (SizeType i = 0; i < a.n_row(); i++) {
		crsMatrix.setRow(i,counter);
		for (SizeType j=0;j<a.n_col();j++) {
			if (a(i,j)==static_cast<T>(0)) continue;
			crsMatrix.pushValue(a(i,j));
			crsMatrix.pushCol(j);
			counter++;
		}

	}
	crsMatrix.setRow(crsMatrix.rows(),counter);
	crsMatrix.checkValidity();
}

/** If order==false then
		creates B such that B_{i1+j1*nout,i2+j2*nout)=A(j1,j2)\delta_{i1,i2}
		if order==true then
		creates B such that B_{i1+j1*na,i2+j2*na)=A(i1,i2)\delta_{j1,j2}
		where na=rank(A)
	  */

template<typename T,typename VectorLikeType>
typename EnableIf<IsVectorLike<VectorLikeType>::True &&
Loki::TypeTraits<typename VectorLikeType::value_type>::isFloat,
void>::Type
externalProduct(CrsMatrix<T>& B,
                const CrsMatrix<T>& A,
                int nout,
                const VectorLikeType& signs,
                bool order=true)
{
	int na=A.rows();
	T tmp;
	assert(A.rows()==A.cols());
	B.resize(na*nout,na*nout);

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
	assert(A.row()==A.col());
	assert(B.row()==B.col());
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
void externalProduct(CrsMatrix<T>  &C,
                     CrsMatrix<T> const &A,
                     CrsMatrix<T> const &B,
                     const typename Vector<int>::Type& signs,
                     bool option=false)
{
	assert(A.row()==A.col());
	assert(B.row()==B.col());
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

template<typename T>
void printFullMatrix(const CrsMatrix<T>& s,
                     const String& name,
                     SizeType how=0,
                     double eps = 1e-20)
{
	Matrix<T> fullm(s.rows(),s.cols());
	crsMatrixToFullMatrix(fullm,s);
	std::cout<<"--------->   "<<name;
	std::cout<<" rank="<<s.rows()<<"x"<<s.cols()<<" <----------\n";
	try {
		if (how==1) mathematicaPrint(std::cout,fullm);
		if (how==2) symbolicPrint(std::cout,fullm);
	} catch (std::exception& e) {

	}

	if (how==0) fullm.print(std::cout,eps);

}

//! C = A*B,  all matrices are CRS matrices
//! idea is from http://web.maths.unsw.edu.au/~farid/Papers/Hons/node23.html
template<typename S,typename S3,typename S2>
void multiply(CrsMatrix<S> &C,
              CrsMatrix<S3> const &A,
              CrsMatrix<S2> const &B)
{
	int j,s,mlast,itemp,jbk;
	SizeType n = A.rows();
	typename Vector<int>::Type ptr(B.cols(),-1),index(B.cols(),0);
	typename Vector<S>::Type temp(B.cols(),0);
	S tmp;

	C.resize(n,B.cols());

	// mlast pointer to the last place we updated in the C vector
	mlast = 0;
	// for (SizeType l=0;l<n;l++) ptr[l] = -1;
	// over the rows of A
	for (SizeType i=0;i<n;i++) {
		C.setRow(i,mlast);
		// start calculations for row
		itemp = 0;
		for (j = A.getRowPtr(i);j< A.getRowPtr(i+1);j++) {
			SizeType istart = B.getRowPtr(A.getCol(j));
			SizeType iend = B.getRowPtr(A.getCol(j)+1);
			for (SizeType k = istart; k< iend;k++) {
				jbk=B.getCol(k);
				tmp = A.getValue(j)*B.getValue(k);
				if (ptr[jbk]<0) {
					ptr[jbk] = itemp;
					temp[ptr[jbk]] = tmp ;
					index[ptr[jbk]] = jbk;
					itemp++;
				} else  {
					temp[ptr[jbk]]+= tmp;
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
	C.checkValidity();
}

// vector2 = sparseMatrix * vector1
template<class S>
void multiply(typename Vector<S>::Type& v2,
              const CrsMatrix<S>& m,
              const typename Vector<S>::Type& v1)
{
	SizeType n = m.rows();
	v2.resize(n);
	for (SizeType i=0;i<n;i++) {
		v2[i]=0;
		for (int j=m.getRowPtr(i);j<m.getRowPtr(i+1);j++) {
			v2[i] += m.getValue(j)*v1[m.getCol(j)];
		}
	}
}

//! Sets B=transpose(conjugate(A))
template<typename S,typename S2>
void transposeConjugate(CrsMatrix<S>& B, const CrsMatrix<S2>& A)
{
	SizeType n=A.rows();
	Vector<Vector<int>::Type >::Type col(n);
	typename Vector<typename Vector<S2>::Type>::Type value(n);

	// B(j,i) = conj(A(i,j))
	for (SizeType i=0;i<n;i++) {
		for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
			col[A.getCol(k)].push_back(i);
			S2 w = A.getValue(k);
			value[A.getCol(k)].push_back(w);
		}
	}

	B.resize(A.cols(),A.rows());

	SizeType counter=0;
	for (SizeType i=0; i<B.rows(); ++i) {

		B.setRow(i,counter);
		for (SizeType j=0;j<col[i].size();j++) {
			if (value[i][j]==static_cast<S>(0.0)) continue;
			B.pushCol(col[i][j]);
			B.pushValue(PsimagLite::conj(value[i][j]));
			counter++;
		}
	}

	B.setRow(B.rows(),counter);
}

//! Sets B=transpose(conjugate(A))
template<typename S,typename S2>
void transposeConjugate(CrsMatrix<S>& B,
                        const CrsMatrix<S2>& A,
                        typename Vector<typename Vector<int>::Type >::Type& col,
                        typename Vector<typename Vector<S2>::Type >::Type& value)
{
	SizeType n=A.rows();
	assert(col.size()==n);
	assert(value.size()==n);
	for (SizeType i=0;i<n;i++) {
		col[i].clear();
		value[i].clear();
	}

	// B(j,i) = conj(A(i,j))
	SizeType counter=0;
	for (SizeType i=0;i<n;i++) {
		for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
			col[A.getCol(k)].push_back(i);
			S2 w = A.getValue(k);
			value[A.getCol(k)].push_back(w);
			counter++;
		}
	}

	B.clear();
	B.resize(A.cols(),A.rows(),counter);

	counter=0;
	for (SizeType i=0;i<B.rows();i++) {

		B.setRow(i,counter);
		for (SizeType j=0;j<col[i].size();j++) {
			B.setCol(counter,col[i][j]);
			B.setValues(counter,PsimagLite::conj(value[i][j]));
			counter++;
		}
	}

	B.setRow(B.rows(),counter);
}

//! Sets A = B(i,perm(j)), A and B CRS matrices
template<class S>
void permute(CrsMatrix<S>& A,
             const CrsMatrix<S>& B,
             const Vector<SizeType>::Type& perm)
{
	SizeType  n = B.rows();

	assert(B.rows()==B.cols());
	A.resize(B.rows(),B.cols());

	typename Vector<int>::Type permInverse(n);
	assert(perm.size() == permInverse.size());
	for (SizeType i=0;i<n;i++) permInverse[perm[i]]=i;

	SizeType counter=0;
	for (SizeType i=0;i<n;i++) {
		A.setRow(i,counter);
		for (int k=B.getRowPtr(i);k<B.getRowPtr(i+1);k++) {
			A.pushCol(permInverse[B.getCol(k)]);
			S tmp = B.getValue(k);
			A.pushValue(tmp);
			counter++;
		}
	}
	A.setRow(n,counter);
	A.checkValidity();
}

//! Sets A = B(perm(i),j), A and B CRS matrices
template<class S>
void permuteInverse(CrsMatrix<S>& A,
                    const CrsMatrix<S>& B,
                    const Vector<SizeType>::Type& perm)
{
	SizeType n = B.rows();
	A.resize(n,B.cols());
	assert(B.rows()==B.cols());

	SizeType counter=0;
	for (SizeType i=0;i<n;i++) {
		SizeType ii = perm[i];
		A.setRow(i,counter);
		for (int k=B.getRowPtr(ii);k<B.getRowPtr(ii+1);k++) {
			A.pushCol(B.getCol(k));
			S tmp = B.getValue(k);
			A.pushValue(tmp);
			counter++;
		}
	}
	A.setRow(n,counter);
}

//! Sets A=B*b1+C*c1, restriction: B.size has to be larger or equal than C.size
template<typename T, typename T1>
void operatorPlus(CrsMatrix<T>& A,
                  const CrsMatrix<T>& B,
                  T1& b1,
                  const CrsMatrix<T>& C,
                  T1& c1)
{
	SizeType n = B.rows();
	assert(B.rows()==B.cols());
	assert(C.rows()==C.cols());

	T tmp;

	assert(n>=C.rows());

	typename Vector<T>::Type  valueTmp(n);
	typename Vector<int>::Type index;
	A.resize(n,B.cols());

	SizeType counter=0;
	for (SizeType k2=0;k2<n;k2++) valueTmp[k2]= static_cast<T>(0.0);

	for (SizeType i = 0; i < n; i++) {
		int k;
		A.setRow(i,counter);

		if (i<C.rows()) {
			// inspect this
			index.clear();
			for (k=B.getRowPtr(i);k<B.getRowPtr(i+1);k++) {
				if (B.getCol(k)<0 || SizeType(B.getCol(k))>=n)
					throw RuntimeError("operatorPlus (1)\n");
				valueTmp[B.getCol(k)]=B.getValue(k)*b1;
				index.push_back(B.getCol(k));
			}

			// inspect C
			for (k=C.getRowPtr(i);k<C.getRowPtr(i+1);k++) {
				tmp = C.getValue(k)*c1;
				if (C.getCol(k)>=int(valueTmp.size()) || C.getCol(k)<0)
					throw RuntimeError("operatorPlus (2)\n");

				valueTmp[C.getCol(k)] += tmp;
				index.push_back(C.getCol(k));
			}
			std::sort(index.begin(),index.end());
			k= -1;
			for (SizeType kk=0;kk<index.size();kk++) {
				if (k==index[kk]) continue;
				k=index[kk];
				if (k<0 || SizeType(k)>=n)
					throw RuntimeError("operatorPlus (3)\n");
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
				tmp = B.getValue(k)*b1;
				A.pushCol(B.getCol(k));
				A.pushValue(tmp);
				counter++;
			}
		}
	}

	A.setRow(n,counter);
}

template<typename T>
bool isHermitian(const CrsMatrix<T>& A,bool=false)
{
	if (A.rows()!=A.cols()) return false;
	for (SizeType i=0;i<A.rows();i++) {
		for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
			if (PsimagLite::norm(A.getValue(k)-PsimagLite::conj(A.element(A.getCol(k),i)))<1e-6)
				continue;
			assert(false);
			return false;
		}
	}
	return true;
}

template<class T>
void sumBlock(CrsMatrix<T> &A,CrsMatrix<T> const &B,SizeType offset)
{
	int counter=0;
	CrsMatrix<T> Bfull(A.rows(),A.cols());

	for (SizeType i=0;i<offset;i++) Bfull.setRow(i,counter);

	for (SizeType ii=0;ii<B.rows();ii++) {
		SizeType i=ii+offset;
		Bfull.setRow(i,counter);
		for (int jj=B.getRowPtr(ii);jj<B.getRowPtr(ii+1);jj++) {
			SizeType j = B.getCol(jj)+offset;
			T tmp  = B.getValue(jj);
			Bfull.pushCol(j);
			Bfull.pushValue(tmp);
			counter++;
		}
	}

	for (SizeType i=B.rows()+offset;i<A.rows();i++) Bfull.setRow(i,counter);
	Bfull.setRow(A.rows(),counter);
	Bfull.checkValidity();
	A += Bfull;
	A.checkValidity();
}

template<class T>
bool isDiagonal(const CrsMatrix<T>& A,double eps=1e-6,bool checkForIdentity=false)
{
	if (A.rows()!=A.cols()) return false;
	SizeType n = A.rows();
	const T f1 = (-1.0);
	for (SizeType i=0;i<n;i++) {
		for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
			SizeType col = A.getCol(k);
			const T& val = A.getValue(k);
			if (checkForIdentity && col==i && PsimagLite::norm(val + f1)>eps) {
				return false;
			}
			if (col!=i && PsimagLite::norm(val)>eps) {
				return false;
			}
		}
	}
	return true;
}

template<class T>
bool isTheIdentity(const CrsMatrix<T>& A,double eps=1e-6)
{
	return isDiagonal(A,eps,true);
}

template<typename T>
typename Real<T>::Type norm2(const CrsMatrix<T>& m)
{
	T val = 0;
	for (SizeType i=0;i<m.values_.size();i++)
		val += PsimagLite::conj(m.values_[i])*m.values_[i];

	return PsimagLite::real(val);
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

} // namespace PsimagLite
/*@}*/
#endif

