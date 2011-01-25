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

		// default assigment operator is fine

		size_t n_row() const { return nrow_; } // legacy name
		
		size_t n_col() const { return ncol_; } // legacy name

		const T& operator()(size_t i,size_t j) const
		{
			return data_[i+j*nrow_];
		}

		T& operator()(size_t i,size_t j)
		{
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


	void diag(Matrix<double> &m,std::vector<double> &eigs,char option)
	{
		char jobz=option;
		char uplo='U';
		int n=m.n_row();
		int lda=m.n_col();
		std::vector<double> work(3);
		int info,lwork= -1;

		eigs.resize(n);

		// query:
		dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
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

	template<class T>
	bool isHermitian(Matrix<T> const &A,bool verbose=false)
	{
		size_t n=A.n_row();
		double eps=1e-6;
		if (n!=A.n_col()) throw std::runtime_error
			("isHermitian called on a non-square matrix.\n");
		for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
			if (norm(A(i,j)-std::conj(A(j,i)))>eps) {
				if (verbose) std::cerr<<"A("<<i<<","<<j<<")="<<A(i,j)<<" A("<<j<<","<<i<<")="<<A(j,i)<<"\n";
			return false;
		}
		return true;
	}


} // namespace PsimagLite
#endif

