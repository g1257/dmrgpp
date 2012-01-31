// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
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
// END LICENSE BLOCK
#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <stdexcept>
#include <iostream>
#include "Lapack.h"
#include "Complex.h"
#include <cassert>

namespace PsimagLite {
	template<typename T>
	class  Matrix  {
	public:
		typedef T value_type; // legacy name
	
		Matrix()
		: nrow_(0), ncol_(0)
		{}

		Matrix(size_t nrow,size_t ncol) 
		: nrow_(nrow),ncol_(ncol),data_(nrow*ncol)
		{}

		// copy constructor
		Matrix(const Matrix<T>& m)
		{
			nrow_=m.nrow_;
			ncol_=m.ncol_;
			data_=m.data_;
		}

		template<typename RealType>
		Matrix(const Matrix<RealType>& m)
		{
			nrow_=m.n_row();
			ncol_=m.n_col();
			data_.resize(nrow_*ncol_);
			for (size_t i=0;i<nrow_;i++)
				for (size_t j=0;j<ncol_;j++) 
					data_[i+j*nrow_] = m(i,j);
		}
		
		template<typename SomeMatrixType>
		Matrix(const SomeMatrixType& m)
		{
			nrow_=m.rank();
			ncol_=m.rank();
			data_.resize(nrow_*ncol_);
			for (size_t i=0;i<nrow_;i++)
				for (size_t j=0;j<ncol_;j++) 
					data_[i+j*nrow_] = m(i,j);
		}

		// default assigment operator is fine

		size_t n_row() const { return nrow_; } // legacy name
		
		size_t n_col() const { return ncol_; } // legacy name

		const T& operator()(size_t i,size_t j) const
		{
			return data_[i+j*nrow_];
		}

		T& operator()(size_t i,size_t j)
		{
			assert(i+j*nrow_<data_.size());
			return data_[i+j*nrow_];
		}

		void resize(size_t nrow,size_t ncol)
		{
			if (nrow_!=0 || ncol_!=0) throw
				std::runtime_error("Matrix::resize(...): only applies when Matrix is empty\n");
			reset(nrow,ncol);
		}

		void reset(size_t nrow,size_t ncol)
		{
			nrow_=nrow; ncol_=ncol;
			data_.resize(nrow*ncol);
			//for (size_t i=0;i<data_.size();i++) data_[i] = 0;
		}
		
		Matrix<T>& operator += (const Matrix<T>& other)
		{
			// domain checking ??? 
			for(size_t i=0;i<ncol_*nrow_;i++) 
				data_[i] += other.data_[i];
			return *this;
		}
		
		Matrix<T>& operator -= (const Matrix<T>& other)
		{
		  // domain checking ??? 
		  for(size_t i=0;i<ncol_*nrow_;i++) 
		    data_[i] -= other.data_[i];
		  return *this;
		}

	private:
		size_t nrow_,ncol_;
		std::vector<T> data_;
	}; // class Matrix

	template<typename T>
	std::ostream &operator<<(std::ostream &os,Matrix<T> const &A)
	{
		size_t i,j;
		os<<A.n_row()<<" "<<A.n_col()<<"\n";
		for (i=0;i<A.n_row();i++) {
			for (j=0;j<A.n_col();j++) os<<A(i,j)<<" ";
			os<<"\n";
		}
		return os;
	}

	template<typename T>
	void mathematicaPrint(std::ostream& os,const Matrix<T>& A)
	{
		os<<"m:={";
		for (size_t i=0;i<A.n_row();i++) {
			os<<"{";
			for (size_t j=0;j<A.n_col();j++) {
				os<<A(i,j);
				if (j+1<A.n_col()) os<<",";
			}
			os<<"}";
			if (i+1<A.n_row()) os<<",\n";
		}
		os<<"}\n";
	}

	template<typename T>
	void symbolicPrint(std::ostream& os,const Matrix<T>& A)
	{
		size_t i,j;
		os<<A.n_row()<<" "<<A.n_col()<<"\n";
		std::vector<T> values;
		std::string s = "symbolicPrint: Not enough characters\n";
		size_t maxCharacters = 25;	
		for (i=0;i<A.n_row();i++) {
			for (j=0;j<A.n_col();j++) {

				const T& val = A(i,j);
				if (std::norm(val)<1e-6) {
					os<<" 0 ";
					continue;
				}

				size_t k=0;
				for (;k<values.size();k++)
					if (std::norm(values[k]-val)<1e-6) break;
				bool b1 = (k==values.size());

				size_t k2 = 0;
				for (;k2<values.size();k2++)
					if (std::norm(values[k2]+val)<1e-6) break;

				bool b2 = (k2==values.size());

				if (b1) {
					if (b2) {
						values.push_back(val);
						if (values.size()>maxCharacters)
							throw std::runtime_error(s.c_str());
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
		for (size_t i=0;i<values.size();i++) {
			char chari = i + 65;
			os<<chari<<"="<<values[i]<<" ";
		}
		os<<"\n";
	}

	template <class T>
	std::istream& operator >> (std::istream& is, Matrix<T>& A)
	{
		size_t nrow=0,ncol=0;
		is >> nrow;
		is >> ncol;
		if(is) {
			A.reset(nrow,ncol);
			for (size_t j=0; j<A.n_row(); j++) for (size_t i=0; i<A.n_col(); i++) {
				is >> A(j,i);
			}
		}
		if(!is) {
			throw std::range_error("ERROR istream& operator >> (std::istream&, Matrix<T>&): read past end stream");
		}
		return is;
	}
	
	template<typename T>
	Matrix<T> operator+(const Matrix<T>& a,const Matrix<T>& b)
	{
		Matrix<T> c(a.n_row(),a.n_col());
		for (size_t i=0;i<a.n_row();i++) for (size_t j=0;j<a.n_col();j++) c(i,j) = a(i,j) + b(i,j);
		return c;
	}
	
	template<typename T>
	Matrix<T> operator-(const Matrix<T>& a,const Matrix<T>& b)
	{
		Matrix<T> c(a.n_row(),a.n_col());
		for (size_t i=0;i<a.n_row();i++) for (size_t j=0;j<a.n_col();j++) c(i,j) = a(i,j) - b(i,j);
		return c;
	}
	

	void diag(Matrix<double> &m,std::vector<double> &eigs,char option)
	{
		char jobz=option;
		char uplo='U';
		int n=m.n_row();
		int lda=m.n_col();
		std::vector<double> work(3);
		int info,lwork= -1;

		if (lda<=0) throw std::runtime_error("lda<=0\n");

		eigs.resize(n);

		// query:
		dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
		if (info!=0) {
			std::cerr<<"info="<<info<<"\n";
			throw std::runtime_error("diag: dsyev_: failed with info!=0.\n");
		}
		lwork = int(work[0])+1;
		work.resize(lwork+1);
		// real work:
		dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
		if (info!=0) {
			std::cerr<<"info="<<info<<"\n";
			throw std::runtime_error("diag: dsyev_: failed with info!=0.\n");
		}

	}

	void diag(Matrix<std::complex<double> > &m,std::vector<double> &eigs,char option)
	{
		char jobz=option;
		char uplo='U';
		int n=m.n_row();
		int lda=m.n_col();
		std::vector<std::complex<double> > work(3);
		std::vector<double> rwork(3*n);
		int info,lwork= -1;

		eigs.resize(n);

		// query:
		zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
		lwork = int(real(work[0]))+1;
		work.resize(lwork+1);
		// real work:
		zheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
		if (info!=0) {
			std::cerr<<"info="<<info<<"\n";
			throw std::runtime_error("diag: zheev: failed with info!=0.\n");
		}

	}

	void diag(Matrix<std::complex<float> > &m,std::vector<float> &eigs,char option)
	{
		char jobz=option;
		char uplo='U';
		int n=m.n_row();
		int lda=m.n_col();
		std::vector<std::complex<float> > work(3);
		std::vector<float> rwork(3*n);
		int info,lwork= -1;

		eigs.resize(n);

		// query:
		cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
		lwork = int(real(work[0]))+1;
		work.resize(lwork+1);
		// real work:
		cheev_(&jobz,&uplo,&n,&(m(0,0)),&lda,&(eigs[0]),&(work[0]),&lwork,&(rwork[0]),&info);
		if (info!=0) {
			std::cerr<<"info="<<info<<"\n";
			throw std::runtime_error("diag: cheev: failed with info!=0.\n");
		}

	}

	template<typename T>
	bool isHermitian(Matrix<T> const &A,bool verbose=false)
	{
		size_t n=A.n_row();
		double eps=1e-6;
		if (n!=A.n_col()) throw std::runtime_error
			("isHermitian called on a non-square matrix.\n");
		for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
			if (std::norm(A(i,j)-std::conj(A(j,i)))>eps) {
				if (verbose) std::cerr<<"A("<<i<<","<<j<<")="<<A(i,j)<<" A("<<j<<","<<i<<")="<<A(j,i)<<"\n";
			return false;
		}
		return true;
	}

	template<typename T>
	void printNonZero(const Matrix<T>& m,std::ostream& os)
	{
		os<<"#size="<<m.n_row()<<"x"<<m.n_col()<<"\n";
		for (size_t i=0;i<m.n_row();i++) {
			size_t nonzero = 0;
			for (size_t j=0;j<m.n_col();j++) {
				const T& val = m(i,j);
				if (val!=static_cast<T>(0)) {
					os<<val<<" ";
					nonzero++;
				}
			}
			if (nonzero>0) os<<"\n";
		}
		os<<"#diagonal non-zero\n";
		for (size_t i=0;i<m.n_row();i++) {
			const T& val = m(i,i);
			if (val!=static_cast<T>(0)) {
				os<<val<<" ";
			}
		}
		os<<"\n";
	}
	
	template<typename T>
	bool isTheIdentity(Matrix<T> const &a)
	{
		
		for (size_t i=0;i<a.n_row();i++) { 
			for (size_t j=0;j<a.n_col();j++) { 
				if (i!=j && std::norm(a(i,j))>0)  {
					std::cerr<<"a("<<i<<","<<j<<")="<<a(i,j)<<"\n";
					return false;
				}
			}
		}
		
		for (size_t i=0;i<a.n_row();i++) 
			if (std::norm(a(i,i)-1.0)>0) return false;
			 
		return true;
	}
	
	template<typename T>
	bool isZero(Matrix<std::complex<T> > const &a)
	{
		
		for (size_t i=0;i<a.n_row();i++) 
			for (size_t j=0;j<a.n_col();j++) 
				if (norm(a(i,j))>0) return false;
		return true;
	}
	
	template<typename T>
	bool isZero(const PsimagLite::Matrix<T>& m)
	{
		bool eps=1e-5;
		for (size_t i=0;i<m.n_row();i++)
			for (size_t j=0;j<m.n_row();j++)
				if (fabs(m(i,j))>eps) return false;
		return true;
	}

	
	template<typename T>
	Matrix<T> multiplyTransposeConjugate(
		const Matrix<T>& O1,
		const Matrix<T>& O2,char modifier='C')
	{
		size_t n=O1.n_row();
		Matrix<T> ret(n,n);
		if (modifier=='C') {
			for (size_t s=0;s<n;s++) 
				for (size_t t=0;t<n;t++) 
					for (size_t w=0;w<n;w++) 
						ret(s,t) += std::conj(O1(w,s))*O2(w,t);
		} else {
			for (size_t s=0;s<n;s++) 
				for (size_t t=0;t<n;t++) 
					for (size_t w=0;w<n;w++) 
						ret(s,t) += O1(w,s)*O2(w,t);
		}
		return ret;
	}
	
	template<class T>
	void transposeConjugate(Matrix<T>& m2,const Matrix<T>& m)
	{
		m2.resize(m.n_row(),m.n_col());
		for (size_t i=0;i<m2.n_row();i++)
			for (size_t j=0;j<m2.n_col();j++) 
				m2(i,j)=std::conj(m(j,i));
		
	}
	
	template<typename T>
	PsimagLite::Matrix<T> transposeConjugate(const Matrix<T>& A)
	{
		PsimagLite::Matrix<T> ret(A.n_col(),A.n_row());
		for (size_t i=0;i<A.n_col();i++)
			for (size_t j=0;j<A.n_row();j++)
				ret(i,j)=std::conj(A(j,i));
		return ret;
	}

} // namespace PsimagLite
#endif

