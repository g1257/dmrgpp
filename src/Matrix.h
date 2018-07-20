/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
#ifndef MATRIX_H_
#define MATRIX_H_

#include <unistd.h>
#include "Vector.h"
#include <stdexcept>
#include <iostream>
#include "LapackExtra.h"
#include "LAPACK.h"
#include "Complex.h"
#include <cassert>
#include "TypeToString.h"
#include "Mpi.h"
#include "Io/IoSerializerStub.h"

namespace PsimagLite {

template<typename T>
struct SvdCmplxFuncPtr {

	typedef std::complex<T> ComplexType;

	typedef int (*F)(char*,
	                 int*,
	                 int*,
	                 ComplexType*,
	                 int*,
	                 T*,
	                 ComplexType*,
	                 int*,
	                 ComplexType*,
	                 int*,
	                 ComplexType*,
	                 int*,
	                 T*,
	                 int*,
	                 int*);
};

template<typename T>
struct SvdRealFuncPtr {
	typedef int (*F)(char*,
	                 int*,
	                 int*,
	                 T*,
	                 int*,
	                 T*,
	                 T*,
	                 int*,
	                 T*,
	                 int*,
	                 T*,
	                 int*,
	                 int*,
	                 int*);
};

template<typename RealType>
void expComplexOrReal(RealType& x,const RealType& y)
{
	x = exp(y);
}

template<typename RealType>
void expComplexOrReal(std::complex<RealType>& x,const RealType& y)
{
	x = std::complex<RealType>(cos(y),sin(y));
}

template<typename T2>
class MatrixNonOwned;

template<typename T>
class  Matrix  {
public:
	typedef T value_type; // legacy name

	Matrix()
	    : nrow_(0), ncol_(0)
	{}

	Matrix(SizeType nrow,SizeType ncol)
	    : nrow_(nrow),
	      ncol_(ncol),
	      data_(nrow*ncol, 0) // the 0 should not be here, FIXME TODO
	{}

	Matrix(const typename Vector<T>::Type& data,
	       SizeType nrow,
	       SizeType ncol)
	    : nrow_(nrow),ncol_(ncol),data_(data)
	{
		if (data.size() < nrow*ncol)
			throw RuntimeError("Matrix::ctor failed\n");

		data_.resize(nrow*ncol);
	}

	template<typename RealType>
	Matrix(const Matrix<RealType>& m,
	       typename EnableIf<!IsComplexNumber<RealType>::True,int>::Type = 0)
	{
		nrow_=m.rows();
		ncol_=m.cols();
		data_.resize(nrow_*ncol_);
		for (SizeType i=0;i<nrow_;i++)
			for (SizeType j=0;j<ncol_;j++)
				data_[i+j*nrow_] = m(i,j);
	}

	// ctor closures
	Matrix(const std::ClosureOperator<Matrix,Matrix,std::ClosureOperations::OP_MULT>& c)
	{
		const T f1 = 1.0;
		matrixMatrix(c.r1,c.r2,f1);
	}

	template<typename T1>
	Matrix(const std::ClosureOperator<
	       std::ClosureOperator<T1,Matrix,std::ClosureOperations::OP_MULT>,
	       Matrix,
	       std::ClosureOperations::OP_MULT>& c,
	       typename EnableIf<Loki::TypeTraits<T1>::isArith,int>::Type = 0)
	{
		*this = c.r1.r1*c.r1.r2*c.r2;
	}

	template<typename T1>
	Matrix(const std::ClosureOperator<T1,
	       std::ClosureOperator<Matrix,Matrix,std::ClosureOperations::OP_PLUS>,
	       std::ClosureOperations::OP_MULT>& c,
	       typename EnableIf<Loki::TypeTraits<T1>::isArith,int>::Type = 0)
	{
		*this = c.r1*(c.r2.r1+c.r2.r2);
	}
	// end all ctors

	void read(int fd)
	{
		::read(fd,&ncol_,sizeof(ncol_));
		::read(fd,&nrow_,sizeof(nrow_));
		data_.resize(nrow_*ncol_);
		::read(fd,&(data_[0]),sizeof(T)*nrow_*ncol_);
	}

	void read(String label, IoSerializer& ioSerializer)
	{
		ioSerializer.read(nrow_, label + "/nrow_");
		ioSerializer.read(ncol_, label + "/ncol_");
		if (nrow_ == 0 || ncol_ == 0) return;
		ioSerializer.read(data_, label + "/data_");
	}

	void write(String label, IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);
		ioSerializer.write(label + "/nrow_", nrow_);
		ioSerializer.write(label + "/ncol_", ncol_);
		if (nrow_ == 0 || ncol_ == 0) return;
		ioSerializer.write(label + "/data_", data_);
	}

	void print(int fd) const
	{
		::write(fd,(const void*)&ncol_,sizeof(ncol_));
		::write(fd,(const void*)&nrow_,sizeof(nrow_));
		::write(fd,(const void*)&(data_[0]),sizeof(T)*nrow_*ncol_);

	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType,
	                   String msg) const
	{
		String str = msg;
		str += "Matrix";
		const char* start = (const char *)&nrow_;
		const char* end = (const char*)&ncol_;
		SizeType total = mres.memResolv(&nrow_,end-start,str + " nrow_");

		start = end;
		end = (const char*)&data_;
		total += mres.memResolv(&ncol_,end-start,str + " ncol_");

		total += mres.memResolv(&data_,sizeof(*this)-total, str+" data_");

		return total;
	}

	void clear()
	{
		nrow_ = ncol_ = 0;
		data_.clear();
	}

	// default assigment operator is fine

	SizeType nonZeros() const
	{
		const T zval = 0.0;
		SizeType n = nrow_*ncol_;
		assert(data_.size() >= n);
		SizeType count = 0;
		for (SizeType i = 0; i < n; ++i)
			if (data_[i] != zval) ++count;

		return count;
	}

#ifndef NO_DEPRECATED_ALLOWED
	SizeType n_row() const { return nrow_; } // legacy name

	SizeType n_col() const { return ncol_; } // legacy name
#endif

	SizeType rows() const { return nrow_; }

	SizeType cols() const { return ncol_; }

	const typename Vector<T>::Type data() const
	{
		return data_;
	}

	const T& operator()(SizeType i,SizeType j) const
	{
		assert(i<nrow_ && j<ncol_);
		assert(i+j*nrow_<data_.size());
		return data_[i+j*nrow_];
	}

	T& operator()(SizeType i,SizeType j)
	{
		assert(i<nrow_ && j<ncol_);
		assert(i+j*nrow_<data_.size());
		return data_[i+j*nrow_];
	}

	bool operator==(const Matrix<T>& other) const
	{
		return (nrow_ == other.nrow_ &&
		        ncol_ == other.ncol_ &&
		        data_ == other.data_);
	}

	void resize(SizeType nrow, SizeType ncol, const T& val)
	{
		nrow_ = nrow;
		ncol_ = ncol;
		data_.resize(nrow*ncol, val);
	}

	void resize(SizeType nrow, SizeType ncol)
	{
		// this throw should not be here, FIXME TODO
		if (nrow_ != 0 || ncol_ != 0) throw
			RuntimeError("Matrix::resize(...): only applies when Matrix is empty\n");

		nrow_ = nrow;
		ncol_ = ncol;
		data_.resize(nrow*ncol, 0); // the 0 should not be here, FIXME TODO
	}

	Matrix<T>& operator+=(const Matrix<T>& other)
	{
		assert(data_.size() == ncol_*nrow_);
		SizeType total = std::min(data_.size(), other.data_.size());
		for (SizeType i = 0; i < total; ++i)
			data_[i] += other.data_[i];
		return *this;
	}

	Matrix<T>& operator-=(const Matrix<T>& other)
	{
		// domain checking ???
		for (SizeType i=0;i<ncol_*nrow_;i++)
			data_[i] -= other.data_[i];
		return *this;
	}

	Matrix<T>& operator*=(const T& value)
	{
		for (SizeType i=0;i<ncol_*nrow_;i++)
			data_[i] *= value;
		return *this;
	}

	void print(std::ostream& os,const double& eps) const
	{
		os<<nrow_<<" "<<ncol_<<"\n";
		for (SizeType i=0;i<nrow_;i++) {
			for (SizeType j=0;j<ncol_;j++) {
				T val = data_[i+j*nrow_];
				if (PsimagLite::norm(val)<eps) val=0.0;
				os<<val<<" ";
			}
			os<<"\n";
		}
	}

	void setTo(const T& val)
	{
		for (SizeType i=0;i<data_.size();i++) data_[i]=val;
	}

	void send(int root,int tag,MPI::CommType mpiComm)
	{
		MPI::send(nrow_,root,tag,mpiComm);
		MPI::send(ncol_,root,tag+1,mpiComm);
		MPI::send(data_,root,tag+2,mpiComm);
	}

	void recv(int root,int tag,MPI::CommType mpiComm)
	{
		MPI::recv(nrow_,root,tag,mpiComm);
		MPI::recv(ncol_,root,tag+1,mpiComm);
		MPI::recv(data_,root,tag+2,mpiComm);
	}

	void bcast(int root,MPI::CommType mpiComm)
	{
		MPI::bcast(nrow_,root,mpiComm);
		MPI::bcast(ncol_,root,mpiComm);
		MPI::bcast(data_,root,mpiComm);
	}

	// start closure members

	template<typename T1>
	Matrix& operator+=(const std::ClosureOperator<T1,Matrix,std::ClosureOperations::OP_MULT>& c)
	{
		nrow_ = c.r2.nrow_;
		ncol_ = c.r2.ncol_;
		this->data_ += c.r1*c.r2.data_;
		return *this;
	}

	template<typename T1>
	Matrix& operator+=(const std::ClosureOperator<Matrix,T1,std::ClosureOperations::OP_MULT>& c)
	{
		nrow_ = c.r1.nrow_;
		ncol_ = c.r1.ncol_;
		this->data_ += c.r2*c.r1.data_;
		return *this;
	}

	template<typename T1>
	Matrix& operator=(const std::ClosureOperator<T1,Matrix,std::ClosureOperations::OP_MULT>& c)
	{
		nrow_ = c.r2.nrow_;
		ncol_ = c.r2.ncol_;
		this->data_ <= c.r1*c.r2.data_;
		return *this;
	}

	template<typename T1>
	Matrix& operator=(const std::ClosureOperator<Matrix,T1,std::ClosureOperations::OP_MULT>& c)
	{
		nrow_ = c.r1.nrow_;
		ncol_ = c.r1.ncol_;
		this->data_ <= c.r2*c.r1.data_;
		return *this;
	}

	template<typename T1>
	Matrix& operator=(const std::ClosureOperator<
	                  std::ClosureOperator<T1,Matrix,std::ClosureOperations::OP_MULT>,
	                  Matrix,
	                  std::ClosureOperations::OP_MULT>& c)
	{
		matrixMatrix(c.r1.r2,c.r2,c.r1.r1);
		return *this;
	}

	template<typename T1>
	Matrix& operator=(const std::ClosureOperator<T1,
	                  std::ClosureOperator<Matrix,Matrix,std::ClosureOperations::OP_PLUS>,
	                  std::ClosureOperations::OP_MULT>& c)
	{
		this->nrow_ = c.r2.r1.nrow_;
		this->ncol_ = c.r2.r1.ncol_;
		this->data_ <= c.r1*(c.r2.r1.data_ + c.r2.r2.data_);
		return *this;
	}

	template<typename T1>
	Matrix& operator=(const std::ClosureOperator<T1,
	                  std::ClosureOperator<Matrix,Matrix,std::ClosureOperations::OP_MINUS>,
	                  std::ClosureOperations::OP_MULT>& c)
	{
		this->nrow_ = c.r2.r1.nrow_;
		this->ncol_ = c.r2.r1.ncol_;
		this->data_ <= c.r1*(c.r2.r1.data_ - c.r2.r2.data_);
		return *this;
	}

	Matrix& operator=
	(const std::ClosureOperator<Matrix,Matrix,std::ClosureOperations::OP_PLUS>& c)
	{
		nrow_ = c.r1.nrow_;
		ncol_ = c.r1.ncol_;
		assert(nrow_ == c.r2.nrow_);
		assert(ncol_ == c.r2.ncol_);
		this->data_ <= c.r1.data_ + c.r2.data_;
		return *this;
	}

	Matrix& operator=
	(const std::ClosureOperator<Matrix,Matrix,std::ClosureOperations::OP_MINUS>& c)
	{
		nrow_ = c.r1.nrow_;
		ncol_ = c.r1.ncol_;
		assert(nrow_ == c.r2.nrow_);
		assert(ncol_ == c.r2.ncol_);
		this->data_ <= c.r1.data_ - c.r2.data_;
		return *this;
	}

	template<typename T1>
	Matrix& operator=
	(const std::ClosureOperator<Matrix,
	 std::ClosureOperator<T1,Matrix,std::ClosureOperations::OP_MULT>,
	 std::ClosureOperations::OP_PLUS>& c)
	{
		nrow_ = c.r1.nrow_;
		ncol_ = c.r1.ncol_;
		assert(nrow_ == c.r2.r2.nrow_);
		assert(ncol_ == c.r2.r2.ncol_);
		this->data_ <= c.r1.data_ + c.r2.r1*c.r2.r2.data_;
		return *this;
	}

	template<typename T1>
	Matrix& operator=
	(const std::ClosureOperator<
	 std::ClosureOperator<T1,Matrix,std::ClosureOperations::OP_MULT>,
	 Matrix,
	 std::ClosureOperations::OP_PLUS>& c)
	{
		nrow_ = c.r2.nrow_;
		ncol_ = c.r2.ncol_;
		assert(nrow_ == c.r1.r2.nrow_);
		assert(ncol_ == c.r1.r2.ncol_);
		this->data_ <= c.r2.data_ + c.r1.r1*c.r1.r2.data_;
		return *this;
	}

	Matrix& operator=
	(const std::ClosureOperator<Matrix<T>,Matrix,std::ClosureOperations::OP_MULT>& c)
	{
		const Matrix<T>& a = c.r1;
		const Matrix<T>& b = c.r2;
		const T f1 = 1.0;
		matrixMatrix(a,b,f1);

		return *this;
	}

	// end closure members

	template<typename T2>
	friend class MatrixNonOwned;

private:

	template<typename T1>
	void matrixMatrix(const Matrix<T>& a, const Matrix<T>& b, const T1& t1)
	{
		Matrix<T> m(a.rows(), b.cols());
		assert(a.cols()==b.rows());
		for (SizeType i=0;i<a.rows();i++) {
			for (SizeType j=0;j<b.cols();j++) {
				T sum = 0.0;
				for (SizeType k=0;k<a.cols();k++) {
					sum += a(i,k) * b(k,j);
				}

				m(i,j) = sum*t1;
			}
		}

		*this = m; // copy is needed in case a or b are *this matrix
	}

	SizeType nrow_;
	SizeType ncol_;
	typename Vector<T>::Type data_;
}; // class Matrix

// start in Matrix.cpp
void geev(char jobvl,
          char jobvr,
          Matrix<std::complex<double> >& a,
          Vector<std::complex<double> >::Type & w,
          Matrix<std::complex<double> >& vl,
          Matrix<std::complex<double> >& vr);

void diag(Matrix<double> &m,Vector<double> ::Type& eigs,char option);

void diag(Matrix<std::complex<double> > &m,Vector<double> ::Type&eigs,char option);

void diag(Matrix<float> &m,Vector<float> ::Type& eigs,char option);

void diag(Matrix<std::complex<float> > &m,Vector<float> ::Type& eigs,char option);

void inverse(Matrix<std::complex<double> >& m);

void inverse(Matrix<double>& m);
// end in Matrix.cpp

template<typename T>
class IsMatrixLike {
public:
	enum { True = false};
};

template<typename T>
class IsMatrixLike<Matrix<T> > {
public:
	enum { True = true};
};

template<typename T>
std::ostream &operator<<(std::ostream &os,Matrix<T> const &A)
{
	SizeType i,j;
	os<<A.rows()<<" "<<A.cols()<<"\n";
	for (i=0;i<A.rows();i++) {
		for (j=0;j<A.cols();j++) os<<A(i,j)<<" ";
		os<<"\n";
	}
	return os;
}

template<typename T>
void mathematicaPrint(std::ostream& os,const Matrix<T>& A)
{
	os<<"m:={";
	for (SizeType i=0;i<A.rows();i++) {
		os<<"{";
		for (SizeType j=0;j<A.cols();j++) {
			os<<A(i,j);
			if (j+1<A.cols()) os<<",";
		}
		os<<"}";
		if (i+1<A.rows()) os<<",\n";
	}
	os<<"}\n";
}

template<typename T>
void symbolicPrint(std::ostream& os,const Matrix<T>& A)
{
	SizeType i,j;
	os<<A.rows()<<" "<<A.cols()<<"\n";
	typename Vector<T>::Type values;
	String s = "symbolicPrint: Not enough characters\n";
	SizeType maxCharacters = 25;
	for (i=0;i<A.rows();i++) {
		for (j=0;j<A.cols();j++) {

			const T& val = A(i,j);
			if (PsimagLite::norm(val)<1e-6) {
				os<<" 0 ";
				continue;
			}

			SizeType k=0;
			for (;k<values.size();k++)
				if (PsimagLite::norm(values[k]-val)<1e-6) break;
			bool b1 = (k==values.size());

			SizeType k2 = 0;
			for (;k2<values.size();k2++)
				if (PsimagLite::norm(values[k2]+val)<1e-6) break;

			bool b2 = (k2==values.size());

			if (b1) {
				if (b2) {
					values.push_back(val);
					if (values.size()>maxCharacters)
						throw RuntimeError(s.c_str());
					char chark = k + 65;
					os<<" "<<chark<<" ";
				} else {
					char chark = k2 + 65;
					os<<"-"<<chark<<" ";
				}
			} else {
				char chark = k + 65;
				os<<" "<<chark<<" ";
			}
		}
		os<<"\n";
	}
	os<<"---------------------------\n";
	for (SizeType i=0;i<values.size();i++) {
		char chari = i + 65;
		os<<chari<<"="<<values[i]<<" ";
	}
	os<<"\n";
}

template<typename T>
void printNonZero(const Matrix<T>& m,std::ostream& os)
{
	os<<"#size="<<m.rows()<<"x"<<m.cols()<<"\n";
	for (SizeType i=0;i<m.rows();i++) {
		SizeType nonzero = 0;
		for (SizeType j=0;j<m.cols();j++) {
			const T& val = m(i,j);
			if (val!=static_cast<T>(0)) {
				os<<val<<" ";
				nonzero++;
			}
		}
		if (nonzero>0) os<<"\n";
	}
	os<<"#diagonal non-zero\n";
	for (SizeType i=0;i<m.rows();i++) {
		const T& val = m(i,i);
		if (val!=static_cast<T>(0)) {
			os<<val<<" ";
		}
	}
	os<<"\n";
}

template<typename T>
std::istream& operator>>(std::istream& is, Matrix<T>& A)
{
	SizeType nrow=0,ncol=0;
	is >> nrow;
	is >> ncol;

	if (is) {
		A.reset(nrow,ncol);
		for (SizeType j=0; j<A.rows(); j++) for (SizeType i=0; i<A.cols(); i++) {
			is >> A(j,i);
		}
	}
	if (!is) {
		String str("ERROR istream& operator >> ");
		str += "(std::istream&, Matrix<T>&): read past end stream";
		throw RangeError(str);
	}
	return is;
}

template<typename T>
bool isHermitian(Matrix<T> const &A,bool verbose=false)
{
	SizeType n=A.rows();
	double eps=1e-6;
	if (n!=A.cols())
		throw RuntimeError("isHermitian called on a non-square matrix.\n");
	for (SizeType i=0;i<n;i++) for (SizeType j=0;j<n;j++) {
		if (PsimagLite::norm(A(i,j)-PsimagLite::conj(A(j,i)))>eps) {
			if (verbose) {
				std::cerr<<"A("<<i<<","<<j<<")="<<A(i,j);
				std::cerr<<" A("<<j<<","<<i<<")="<<A(j,i)<<"\n";
			}
			return false;
		}
	}
	return true;
}

template<typename T>
bool isTheIdentity(Matrix<T> const &a)
{

	for (SizeType i=0;i<a.rows();i++) {
		for (SizeType j=0;j<a.cols();j++) {
			if (i!=j && PsimagLite::norm(a(i,j))>1e-6)  {
				std::cerr<<"a("<<i<<","<<j<<")="<<a(i,j)<<"\n";
				return false;
			}
		}
	}

	for (SizeType i=0;i<a.rows();i++)
		if (PsimagLite::norm(a(i,i)-static_cast<T>(1.0))>1e-6) return false;

	return true;
}

template<typename T>
bool isZero(Matrix<std::complex<T> > const &a)
{

	for (SizeType i=0;i<a.rows();i++)
		for (SizeType j=0;j<a.cols();j++)
			if (PsimagLite::norm(a(i,j))>0) return false;
	return true;
}

template<typename T>
bool isZero(const Matrix<T>& m)
{
	for (SizeType i=0;i<m.rows();i++)
		for (SizeType j=0;j<m.cols();j++)
			if (fabs(m(i,j))>0) return false;
	return true;
}

// closures start

template<typename T1,typename T2>
typename EnableIf<(IsMatrixLike<T1>::True || IsMatrixLike<T2>::True)
&& !std::IsClosureLike<T1>::True && !std::IsClosureLike<T2>::True,
std::ClosureOperator<T1,T2,std::ClosureOperations::OP_MULT> >::Type operator*(const T1& a,
                                                                              const T2& b)
{
	return std::ClosureOperator<T1,T2,std::ClosureOperations::OP_MULT>(a,b);
}

template<typename T>
std::ClosureOperator<Matrix<T>,Matrix<T>,std::ClosureOperations::OP_PLUS>
operator+(const Matrix<T>& a,const Matrix<T>& b)
{
	return std::ClosureOperator<Matrix<T>,Matrix<T>,std::ClosureOperations::OP_PLUS>(a,b);
}

template<typename T>
std::ClosureOperator<Matrix<T>,Matrix<T>,std::ClosureOperations::OP_MINUS>
operator-(const Matrix<T>& a,const Matrix<T>& b)
{
	return std::ClosureOperator<Matrix<T>,Matrix<T>,std::ClosureOperations::OP_MINUS>(a,b);
}

template<typename T,typename A>
void operator<=(std::vector<T,A>& v, const std::ClosureOperator<Matrix<T>,
                std::vector<T,A>,
                std::ClosureOperations::OP_MULT>& c)
{
	const std::vector<T,A>& b = c.r2;
	const Matrix<T>& a = c.r1;
	assert(a.cols()==b.size());
	v.resize(a.rows());
	for (SizeType i=0;i<a.rows();i++) {
		T sum = 0;
		for (SizeType j=0;j<b.size();j++) sum += a(i,j)*b[j];
		v[i] = sum;
	}
}

template<typename T,typename A>
void operator<=(std::vector<T,A>& v, const std::ClosureOperator<std::vector<T,A>,
                Matrix<T>,
                std::ClosureOperations::OP_MULT>& c)
{
	const std::vector<T,A>& b = c.r1;
	const Matrix<T>& a = c.r2;
	assert(a.rows()==b.size());
	v.resize(a.cols());
	for (SizeType i=0;i<a.cols();i++) {
		T sum = 0;
		for (SizeType j=0;j<b.size();j++) sum += b[j] * a(j,i);
		v[i] = sum;
	}
}

template<typename T>
Matrix<T> multiplyTransposeConjugate(const Matrix<T>& O1,
                                     const Matrix<T>& O2,
                                     char modifier='C')
{
	SizeType n=O1.rows();
	Matrix<T> ret(n,n);
	if (modifier=='C') {
		for (SizeType s=0;s<n;s++)
			for (SizeType t=0;t<n;t++)
				for (SizeType w=0;w<n;w++)
					ret(s,t) += PsimagLite::conj(O1(w,s))*O2(w,t);
	} else {
		for (SizeType s=0;s<n;s++)
			for (SizeType t=0;t<n;t++)
				for (SizeType w=0;w<n;w++)
					ret(s,t) += O1(w,s)*O2(w,t);
	}
	return ret;
}

template<class T>
void transposeConjugate(Matrix<T>& m2,const Matrix<T>& m)
{
	m2.resize(m.cols(),m.rows());
	for (SizeType i=0;i<m2.rows();i++)
		for (SizeType j=0;j<m2.cols();j++)
			m2(i,j)=PsimagLite::conj(m(j,i));

}

template<class T>
void transpose(Matrix<T>& m2,const Matrix<T>& m)
{
	m2.resize(m.cols(),m.rows());
	for (SizeType i=0;i<m2.rows();++i)
		for (SizeType j=0;j<m2.cols();++j)
			m2(i,j) = m(j,i);
}

template<typename T>
Matrix<T> transposeConjugate(const Matrix<T>& A)
{
	Matrix<T> ret(A.cols(),A.rows());
	for (SizeType i=0;i<A.cols();i++)
		for (SizeType j=0;j<A.rows();j++)
			ret(i,j)=PsimagLite::conj(A(j,i));
	return ret;
}

template<typename T>
void exp(Matrix<T>& m)
{
	SizeType n = m.rows();
	typename Vector<typename Real<T>::Type>::Type eigs(n);
	diag(m,eigs,'V');
	Matrix<T> expm(n,n);
	for (SizeType i=0;i<n;i++) {
		for (SizeType j=0;j<n;j++) {
			T sum = 0;
			for (SizeType k=0;k<n;k++) {
				typename Real<T>::Type alpha = eigs[k];
				T tmp = 0.0;
				expComplexOrReal(tmp,alpha);
				sum+= PsimagLite::conj(m(i,k))*m(j,k)*tmp;
			}
			expm(i,j) = sum;
		}
	}
	m = expm;

}

template<typename T>
T norm2(const Matrix<T>& m)
{
	T sum = 0;
	for (SizeType i=0;i<m.rows();i++)
		for (SizeType j=0;j<m.cols();j++)
			sum += m(i,j) * m(i,j);
	return sum;
}

template<typename T>
T norm2(const Matrix<std::complex<T> >& m)
{
	T sum = 0;
	for (SizeType i=0;i<m.rows();i++)
		for (SizeType j=0;j<m.cols();j++)
			sum += PsimagLite::real(PsimagLite::conj(m(i,j)) * m(i,j));
	return sum;
}

template<typename T>
void outerProduct(Matrix<T>& A,const Matrix<T>& B,const Matrix<T>& C)
{
	SizeType ni = B.rows();
	SizeType nj = B.cols();

	A.resize(B.rows()*C.rows(),B.cols()*C.cols());

	for (SizeType i1 = 0; i1 < B.rows(); ++i1)
		for (SizeType j1 = 0; j1 < B.cols(); ++j1)
			for (SizeType i2 = 0; i2 < C.rows(); ++i2)
				for (SizeType j2 = 0; j2 < C.cols(); ++j2)
					A(i1+i2*ni,j1+j2*nj) = B(i1,j1) * C(i2,j2);
}

// closures end

void svd(char jobz,
         Matrix<double>& a,
         Vector<double>::Type& s,
         Matrix<double>& vt,
         SvdRealFuncPtr<double>::F = psimag::LAPACK::dgesdd_,
         int = 0);
void svd(char jobz,
         Matrix<std::complex<double> >& a,
         Vector<double>::Type& s,
         Matrix<std::complex<double> >& vt,
         SvdCmplxFuncPtr<double>::F = psimag::LAPACK::zgesdd_,
         int = 0);
void svd(char jobz, Matrix<float>& a,
         Vector<float>::Type& s,
         Matrix<float>& vt,
         SvdRealFuncPtr<float>::F backend = psimag::LAPACK::sgesdd_,
         int = 0);
void svd(char jobz,
         Matrix<std::complex<float> >& a,
         Vector<float>::Type& s,
         Matrix<std::complex<float> >& vt,
         SvdCmplxFuncPtr<float>::F backend = psimag::LAPACK::cgesdd_,
         int = 0);

#ifdef USE_MPI
namespace MPI {

template<typename SomeMatrixType>
typename EnableIf<IsMatrixLike<SomeMatrixType>::True &
Loki::TypeTraits<typename SomeMatrixType::value_type>::isArith,
void>::Type allReduce(SomeMatrixType& v,
                      MPI_Op op = MPI_SUM,
                      CommType mpiComm = MPI_COMM_WORLD)
{
	SomeMatrixType recvbuf = v;
	MPI_Datatype datatype = MpiData<typename SomeMatrixType::value_type>::Type;
	SizeType total = v.rows() * v.cols();
	int errorCode = MPI_Allreduce(&(v(0,0)),&(recvbuf(0,0)),total,datatype,op,mpiComm);
	checkError(errorCode,"MPI_Allreduce",mpiComm);
	v = recvbuf;
}

template<typename SomeMatrixType>
typename EnableIf<IsMatrixLike<SomeMatrixType>::True &
IsComplexNumber<typename SomeMatrixType::value_type>::True,
void>::Type allReduce(SomeMatrixType& v,
                      MPI_Op op = MPI_SUM,
                      CommType mpiComm = MPI_COMM_WORLD)
{
	SomeMatrixType recvbuf = v;
	MPI_Datatype datatype = MpiData<typename SomeMatrixType::value_type::value_type>::Type;
	SizeType total = v.rows() * v.cols();
	int errorCode = MPI_Allreduce(&(v(0,0)),&(recvbuf(0,0)),2*total,datatype,op,mpiComm);
	checkError(errorCode,"MPI_Allreduce",mpiComm);
	v = recvbuf;
}

} // namespace MPI
#endif // USE_MPI

} // namespace PsimagLite
#endif
