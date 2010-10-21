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
#ifndef UTILS_HEADER_H
#define UTILS_HEADER_H

//! Utility functions
#include <vector>
#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include "Matrix.h" // psimag
#include "Sort.h"

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *, 
	std::complex<double> *,int *, double *, int *);
extern "C" void dsyev_(char *,char *,int *,double *,int *, double *,double *,int *,int *);

namespace std {
	inline double conj(double const &v) { return v; }
	inline double real(double const &v) { return v; }
	inline double imag(double const &v) { return 0.0; }

	template<class T1,class T2>
	ostream &operator<<(std::ostream &os,const pair<T1,T2>& p)
	{
		os<<p.first<<" "<<p.second;
		return os;
	}
	template<class X>
	std::ostream &operator<<(std::ostream &s,std::vector<X> const &v)
	{
		for (size_t i=0;i<v.size();i++) s<<i<<" "<<v[i]<<"\n";
		return s;
	}
	
	template<class X>
	inline X operator*(std::vector<X> const &v,std::vector<X> const &w)
	{
		X result=0;
		for (size_t i=0;i<v.size();i++) result += v[i]*conj(w[i]);
		return result;
	}
	
	
	template<class X>
	X norm(std::vector<X> const &v)
	{
		return sqrt(v*v);
	}
	
	template<class X>
	X norm(std::vector<std::complex<X> > const &v)
	{
		std::complex<X> x = v*v;
		if (fabs(imag(x))>1e-5) throw std::runtime_error("Norm isn't real\n");
		return sqrt(real(x));
	}
	
	template<class X>
	void vectorNormalize(std::vector<X>& v,X x=0)
	{
		if (x==0) x=norm(v);
		for (size_t i=0;i<v.size();++i) v[i] /= x;
	}

	template<typename T>
	void vectorInvert(std::vector<T>& dest,const std::vector<T>& src)
	{
		int n = src.size();
		dest.clear();
		for (int i=n-1;i>=0;i--) dest.push_back(src[i]);
	}
	
	template<typename T>
	bool operator!=(const std::vector<T>& v1,const std::vector<T>& v2)
	{
		bool eps = 1e-6;
		if (v1.size()!=v2.size()) return false;
		for (size_t i=0;i<v1.size();i++) if (fabs(v1[i]-v2[i])>eps) return false;
		return true;
	}
	
	template<typename FieldType>
	void vectorConvert(std::vector<FieldType>& dest,const std::vector<FieldType>&src)
	{
		dest = src;
	}
	
	template<typename FieldType>
	void vectorConvert(std::vector<std::complex<FieldType> >& dest,const std::vector<FieldType>&src)
	{
		dest.resize(src.size());
		for (size_t i=0;i<dest.size();i++) dest[i] = src[i];
	}
	
	template<typename FieldType,typename FieldType2>
	inline std::vector<FieldType2> operator*=(const std::vector<FieldType2>& v,FieldType value)
	{
		std::vector<FieldType2> w = v;
		for (size_t i=0;i<w.size();i++) w[i] *= value;
		return w;
	}
	
	template<typename FieldType>
	inline std::vector<FieldType> operator+=(std::vector<FieldType>& v,const std::vector<FieldType>& w)
	{
		for (size_t i=0;i<w.size();i++) v[i] += w[i];
		return v;
	}
}
	
namespace utils {


	template<typename ContainerType>
	void sort(ContainerType& x,std::vector<size_t>& iperm)
	{
		Sort<ContainerType> s;
		s.sort(x,iperm);
	}
	
	/*template<template<typename,typename> class ContainerTemplate,typename Field,typename A>
	void sort(ContainerTemplate<Field,A>& x,std::vector<size_t>& iperm)
	{
		sort<double,std::vector<Field,A>,Field>(x,iperm);				
	}
	
	template<template<typename,typename,typename,typename> class ContainerTemplate,typename Field,typename Cmp,typename A>
	void sort(ContainerTemplate<size_t,Field,Cmp,A>& x,std::vector<size_t>& iperm)
	{
		sort<double,ContainerTemplate<size_t,Field,Cmp,A>,Field>(x,iperm);				
	}*/
	
	std::string getTimeDate()
	{
		struct timeval tv;
		time_t tt;
	
		gettimeofday(&tv,0);
		tt=tv.tv_sec; /* seconds since 1970 */
		return asctime(localtime(&tt));
	}
	
	template<typename SomeType>
	void truncateVector(std::vector<SomeType> &v,std::vector<SomeType> const &removedIndices)
	{
		std::vector<SomeType> tmpVector;
		for (size_t i=0;i<v.size();i++) {
			if (isInVector(removedIndices,i)>=0) continue;
			tmpVector.push_back(v[i]);
		}
		v=tmpVector;
	}
	
	template<typename SomeType>
	void reorder(std::vector<SomeType> &v,std::vector<size_t> const &permutation)
	{
		std::vector<SomeType> tmpVector(v.size());
		for (size_t i=0;i<v.size();i++) tmpVector[i]=v[permutation[i]]; 
		v = tmpVector;
	}
	
	template<typename SomeType>
	void reorder(psimag::Matrix<SomeType>& v,std::vector<size_t> const &permutation)
	{
		psimag::Matrix<SomeType> tmpVector(v.n_row(),v.n_col());
		for (size_t i=0;i<v.n_row();i++) 
			for (size_t j=0;j<v.n_col();j++)
				tmpVector(i,j)=v(permutation[i],permutation[j]); 
		v = tmpVector;
	}
	
	template<typename T,class  Op >
	size_t vectorMax(const std::vector<T>& v,T max)
	{
		size_t isaved=0;
		Op f;
		for (size_t i=0;i<v.size();i++) {
			if (f(v[i],max)) {
				max=v[i];
				isaved=i;
			}
		}
		return isaved;
	}
		
	template<typename FieldType>
	bool isAnInteger(const FieldType& t)
	{
		if (t==int(t)) return true;
		return false;
	}		
	
	template<class T>
	std::string ttos(T t)
	{
		std::stringstream ss;
		std::string str;
		ss.precision(10);
		ss<<t;
		ss>>str;
		return str;
	}

	template<class T>
	inline T square(T x)
	{
		return (x*x);
	}
	
	
	//! given ind and n, get x and y such that ind = x + y*n
	inline void getCoordinates(size_t &x,size_t &y,size_t ind,size_t n)
	{
		//y  = ind/n;
		//x = ind - y*n;
		//x= ind % n;
		div_t q = div(ind, n);
		y = q.quot;
		x = q.rem;
	}


	template<typename X,typename Y>
	inline int isInVector(std::vector<X> const &natBasis,Y const &v)
	{
// 		if (natBasis.size()==0) return -1;
// 		for (size_t ii=0;ii<natBasis.size();ii++) if (natBasis[ii]==v) return ii;
// 		return -1;
		typename std::vector<X>::const_iterator x = find(natBasis.begin(),natBasis.end(),v);
		if (x==natBasis.end()) return -1;
		return x-natBasis.begin();
		
	}
	
	//! Find smax and emin such that:
	//! smax = maximum(B intersection sBlock)
	//! emin = minimum(B intersection envBlock)
	void findExtremes(size_t& smax,size_t& emin,std::vector<size_t> const &B,std::vector<size_t> const &sBlock)
	{
		smax = 0;
		emin = 1000;
		int r;
		for (size_t i=0;i<B.size();i++) {
			r = utils::isInVector(sBlock,B[i]);
			if (r<0) { //it's in the environment
				if (B[i]<emin) emin=B[i];
			} else { //it's in the system
				if (B[i]>smax) smax=B[i];
			}
		}
		
	}

	template<class T>
	void transposeConjugate(psimag::Matrix<T>& m2,const psimag::Matrix<T>& m)
	{
		size_t i,j;
		m2.resize(m.n_row(),m.n_col());
		for (i=0;i<m2.n_row();i++) for (j=0;j<m2.n_col();j++) m2(i,j)=std::conj(m(j,i));
		
	}
	
	template<class T>
	void vectorPrint(std::vector<T> const &v,char const *name,std::ostream &s)
	{
		unsigned int i;
		for (i=0;i<v.size();i++) s<<name<<"["<<i<<"]="<<v[i]<<std::endl;
	}
	
	template<typename T1,typename T2>
	void vectorPrint(std::map<T1,T2>  &v,char const *name,std::ostream &s)
	{
		// Stroustrup page 483
		typedef typename std::map<T1,T2>::const_iterator CI;
		
		for (CI p=v.begin();p!=v.end();++p) {
			s<<name<<"["<<p->first<<"]="<<p->second<<"\n";
		}
	}
	
	
	
	template<class T>
	bool vectorEqual(std::vector<T> const &a,std::vector<T> const &b)
	{
		size_t n=a.size();
		if (n!=b.size()) return false;
		T eps=1e-6;
		for (size_t i=0;i<n;i++) if (fabs(a[i]-b[i])>eps) return false;
		return true;
	}
	
	template<class T>
	bool vectorEqual(std::vector<std::pair<T,T> > const &a,std::vector<std::pair<T,T> > const &b)
	{
		size_t n=a.size();
		if (n!=b.size()) return false;
		
		for (size_t i=0;i<n;i++) if (a[i]!=b[i]) return false;
		return true;
	}
	
	template<class T>
	void difference	(std::vector<T> const &a,std::vector<T> const &b)	
	{
		size_t n=a.size();
		if (n!=b.size()) std::cerr<<"Vector length different\n";
		T eps=1e-6;
		for (size_t i=0;i<n;i++) if (fabs(a[i]-b[i])>eps) std::cerr<<"a["<<i<<"]="<<a[i]<<" but b="<<b[i]<<"\n";
	}
	
	template<class T>
	void matrixPrint(psimag::Matrix<T> const &a,std::ostream &s,char separator='\t')
	{
		T eps=1e-8;
		for (size_t i=0;i<a.n_row();i++) {
			for (size_t j=0;j<a.n_col();j++) {
				T tmp= a(i,j);
				if (fabs(tmp)<eps) tmp=0;
				s<<tmp<<separator;
			}
			s<<"\n";
		}
	}
	
	
	template<class T>
	void matrixPrint(T const &a,std::ostream &s)
	{
		for (int i=0;i<a.getSize();i++) {
			for (int j=0;j<a.getSize();j++) s<<a(i,j)<<"\t";
			s<<"\n";
		}
	}
	
	void diag(psimag::Matrix<double> &m,std::vector<double> &eigs,char option)
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
		lwork = int(std::real(work[0]))+1;
		work.resize(lwork+1);	
		// real work:
		dsyev_(&jobz,&uplo,&n,&(m(0,0)),&lda, &(eigs[0]),&(work[0]),&lwork, &info);
		if (info!=0) {
			std::cerr<<"info="<<info<<"\n";
			throw std::runtime_error("diag: dsyev_: failed with info!=0.\n");
		}	
		
	}
	
	void diag(psimag::Matrix<std::complex<double> > &m,std::vector<double> &eigs,char option)
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
	
	template<typename T,typename CrsMatrixType>
	void diagTest(const CrsMatrixType& m,const std::string& label,bool option=false)
	{
		psimag::Matrix<T> fullm;
		crsMatrixToFullMatrix(fullm,m);
		
		std::cerr<<"MAAATRIX"<<label<<" ";
		for (size_t i=0;i<fullm.n_col();i++) 
			std::cerr<<fullm(0,i)<<" ";
		std::cerr<<"\n";			
		std::vector<T> eigs(fullm.n_row());
		utils::diag(fullm,eigs,'V');
		std::cerr<<label<<" eigs[0]="<<eigs[0]<<"\n";
					
		
		
		if (!option) return;
		std::cout<<"MAAAAAAAAAAAATRIX\n";
		mathematicaPrint(std::cout,fullm);
		std::cout<<"---------------------\n";
					
	}
	
	template<typename T>
	bool isZero(const psimag::Matrix<T>& m)
	{
		bool eps=1e-5;
		for (size_t i=0;i<m.n_row();i++)
			for (size_t j=0;j<m.n_row();j++)
				if (fabs(m(i,j))>eps) return false;
		return true;
	}
	
	template<typename T>
	void transform(psimag::Matrix<T>& m,const psimag::Matrix<T>& transform)
	{		
		int nBig = m.n_row();
		int nSmall = transform.n_col();
		double alpha=1.0;
		double beta=0.0;
		psimag::Matrix<T> fmS,fmTmp;
		
		
		
		fmTmp.resize(nBig,nSmall);
		
		
		psimag::BLAS::GEMM('N','N',nBig,nSmall,nBig,alpha,&(m(0,0)),nBig,&(transform(0,0)),nBig,beta,&(fmTmp(0,0)),nBig);
		fmS.resize(nSmall,nSmall);
		psimag::BLAS::GEMM('C','N',nSmall,nSmall,nBig,alpha,&(transform(0,0)),nBig,&(fmTmp(0,0)),nBig,beta,&(fmS(0,0)),nSmall);
		
		m=fmS;
	}
	
	//! Sets A = B(i,perm(j)), A and B CRS matrices	
	template<typename S>
	void permute(psimag::Matrix<S>& A,const psimag::Matrix<S>& B,const std::vector<size_t>& perm)
	{
		size_t n = B.n_row();
		A.resize(n,n);
	
		for (size_t i=0;i<n;i++) 
			for (size_t k=0;k<B.n_col();k++) 
				A(i,k)=B(i,perm[k]);
				
				
	}
	
	//! Sets A = B(perm(i),j), A and B CRS matrices		
	template<typename S>
	void permuteInverse(psimag::Matrix<S>& A,const psimag::Matrix<S>& B,const std::vector<size_t>& perm)
	{
		size_t n = B.n_row();
		A.resize(n,n);
		
		
		for (size_t i=0;i<n;i++) 
			for (size_t k=0;k<B.n_col();k++) 
				A(i,k)=B(perm[i],k);
		
				
	}
	
	
	
	//! swap column i and column j of matrix m
	template<class Field>
	void swapMatrix(psimag::Matrix<Field> &m,int i,int j)
	{
		if (i==j) return;
		size_t k;
		std::vector<int> temp(m.n_row());
		
		for (k=0;k<m.n_row();k++) {
			temp[k] = m(k,i);
			m(k,i)=m(k,j);
		}
		for (k=0;k<m.n_row();k++) {
			m(k,j)=temp[i];
		}
		
	}
	
	template<class Field>
	void sortTransform(psimag::Matrix<int> &m,std::vector<Field> const &v)
	{
		size_t i,j;
		for (i=0;i<m.n_row();i++) { 
			for (j=0;j<m.n_col();j++) {
				m(i,j)=0;
				if (i==j) m(i,j)=1;
			}
		}
		 
		for (i=0;i<v.size();i++) {
			for (j=v.size()-1;j>=i+1;j--) {
				if (v[j-1]>v[j]) {
					swapMatrix(m,j-1,j);
				}
			}
		}
	}

	
	
       // A = B union C
        template<typename Block>
        void blockUnion(Block &A,Block const &B,Block const &C)
        {
		A=B;
		for (size_t i=0;i<C.size();i++) A.push_back(C[i]);
        }
	
	template<typename Field>
	void myRandomT(std::complex<Field> &value)
	{
		value = std::complex<Field>(drand48 () - 0.5, drand48 () - 0.5);
	}

	template<typename Field>
	void myRandomT(Field &value)
	{
		value = drand48 () - 0.5;
	}
	
	template<typename Field>
	Field myProductT(Field const &value1,Field const &value2)
	{
		return std::real(value1*std::conj(value2));
	}
	
	template<typename Field>
	Field myProductT(std::complex<Field> const &value1,std::complex<Field> const &value2)
	{
		return real(value1*conj(value2));
	}
	
	void getrusageForLinux(std::ostream& os)
	{
		
		std::ifstream fin("/proc/self/stat");
		std::string s;
		std::vector<std::string> elems;
		while(!fin.eof()) {
			fin>>s;
			elems.push_back(s);
		}
		fin.close();
		if (elems.size()<24) return;
		os<<"MU: proc/self/stat[23]="<<atof(elems[23].c_str())<<"\n";
	}

	
	void memoryUsage(std::ostream& os)
	{
		struct rusage usage;
		
		int ret = getrusage(RUSAGE_SELF,&usage); 
		if (ret!=0) {
			os<<"MU: error code="<<ret<<" returning...\n";
			return;
		}
		if (usage.ru_maxrss==0) {
			getrusageForLinux(os);
			return;
		}
		os<<"MU: maximum resident set size ="<<usage.ru_maxrss<<"\n";
		os<<"MU: integral shared memory size="<<usage.ru_ixrss<<"\n";
		os<<"MU: integral unshared data size="<<usage.ru_idrss<<"\n";
		os<<"MU: integral unshared stack size="<<usage.ru_isrss<<"\n";
		os<<"MU: page reclaims="<<usage.ru_minflt<<"\n";
		os<<"MU: page faults="<<usage.ru_majflt<<"\n";
		os<<"MU: page swaps="<<usage.ru_nswap<<"\n";
		// etc, etc
	}
	
	
	
}




namespace psimag {
	
	
	template<typename T>
	void transposeConjugate(psimag::Matrix<T>& dest,const psimag::Matrix<T>& src)
	{
		size_t n = src.n_row();
		if (n!=src.n_col()) throw std::runtime_error("transposeConjugate: only for square matrices\n");
		dest.resize(n,n);
		for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) dest(i,j)=conj(src(j,i));
	}
	
	template<typename T>
	psimag::Matrix<T> multiplyTransposeConjugate(const psimag::Matrix<T>& O1,const psimag::Matrix<T>& O2,char modifier='C')
	{
		size_t n=O1.n_row();
		psimag::Matrix<T> ret(n,n);
		if (modifier=='C') {
			for (size_t s=0;s<n;s++) for (size_t t=0;t<n;t++) for (size_t w=0;w<n;w++) ret(s,t) += conj(O1(w,s))*O2(w,t);
		} else {
			for (size_t s=0;s<n;s++) for (size_t t=0;t<n;t++) for (size_t w=0;w<n;w++) ret(s,t) += O1(w,s)*O2(w,t);
		}
		return ret;
	}
		
	template<typename T>
	psimag::Matrix<T> operator+(const psimag::Matrix<T>& a,const psimag::Matrix<T>& b)
	{
		psimag::Matrix<T> c(a.n_row(),a.n_col());
		for (size_t i=0;i<a.n_row();i++) for (size_t j=0;j<a.n_col();j++) c(i,j) = a(i,j) + b(i,j);
		return c;
	}
	
	template<typename T>
	psimag::Matrix<T> operator-(const psimag::Matrix<T>& a,const psimag::Matrix<T>& b)
	{
		psimag::Matrix<T> c(a.n_row(),a.n_col());
		for (size_t i=0;i<a.n_row();i++) for (size_t j=0;j<a.n_col();j++) c(i,j) = a(i,j) - b(i,j);
		return c;
	}
	
	inline double norm(double const &v) { return fabs(v); }
	
	inline double norm(std::complex<double> const &v) 
	{
		double x = real(v)*real(v)+imag(v)*imag(v);
		return sqrt(x); 
	}
	
	template<typename T>
	void enforcePhase(T* v,size_t n)
	{
		T sign1=0;
		for (size_t j=0;j<n;j++) {
			if (psimag::norm(v[j])>1e-6) {
				if (real(v[j])>0) sign1=1;
				else sign1= -1;
				break;
			}
		}
		// get a consistent phase
		for (size_t j=0;j<n;j++) v[j] *= sign1;
	}
	
	template<typename T>
	void enforcePhase(std::vector<T>& v)
	{
		enforcePhase(&(v[0]),v.size());	
	}
	
	template<typename T>
	void enforcePhase(psimag::Matrix<T>& a)
	{
		T* vpointer = &(a(0,0));
		for (size_t i=0;i<a.n_col();i++) 
			enforcePhase(&(vpointer[i*a.n_row()]),a.n_row());
	}
	
	template<typename T,typename T2>
	bool almostEqual(const psimag::Matrix<T>& a,const psimag::Matrix<T>& b,const T2& eps)
	{
		for (size_t i=0;i<a.n_row();i++) {
			for (size_t j=0;j<a.n_col();j++) { 
				if (psimag::norm(a(i,j)-b(i,j))>eps) {
					std::cerr<<"a("<<i<<","<<j<<")="<<a(i,j);
					std::cerr<<" b("<<i<<","<<j<<")="<<b(i,j)<<"\n";
					throw std::runtime_error("almostEqual\n");
				}
			}
		}
		return true;
	}
	
	template<typename T,typename T2>
	bool isTheIdentity(psimag::Matrix<T> const &a,T2 eps = 1e-7)
	{
		
		for (size_t i=0;i<a.n_row();i++) { 
			for (size_t j=0;j<a.n_col();j++) { 
				if (i!=j && psimag::norm(a(i,j))>eps)  {
					std::cerr<<"a("<<i<<","<<j<<")="<<a(i,j)<<"\n";
					return false;
				}
			}
		}
		
		for (size_t i=0;i<a.n_row();i++) if (norm(a(i,i)-1.0)>eps) return false;
			 
		return true;
	}	
	
	template<typename T>
	bool isZero(psimag::Matrix<T> const &a,T eps = 1e-7)
	{
		
		for (size_t i=0;i<a.n_row();i++) 
			for (size_t j=0;j<a.n_col();j++) 
				if (psimag::norm(a(i,j))>eps) return false;
		return true;
	}	
	
	template<typename T>
	bool isZero(psimag::Matrix<std::complex<T> > const &a,T eps = 1e-7)
	{
		
		for (size_t i=0;i<a.n_row();i++) 
			for (size_t j=0;j<a.n_col();j++) 
				if (psimag::norm(a(i,j))>eps) return false;
		return true;
	}	
	
	template<typename T>
	void matrixIdentity(psimag::Matrix<T>& identity,size_t n)
	{
		identity.resize(n,n);
		for (size_t i=0;i<n;i++) {
			for (size_t j=0;j<n;j++) {
				identity(i,j)=0;
				if (i==j) identity(i,i)=1.0;
			}
		}
	}
	template<typename T>
	T maxElement(const psimag::Matrix<T>& S)
	{
		size_t i0=0;
		size_t j0=0;
		T maxe = 0.0;
		for (size_t i=0;i<S.n_row();i++) {
			for (size_t j=0;j<S.n_col();j++) {
				if (psimag::norm(maxe)<psimag::norm(S(i,j))) {
					maxe = S(i,j);
					i0 = i;
					j0 = j;
				}
			}
		}
		std::cerr<<"MAX ELEMENT S("<<i0<<","<<j0<<")="<<maxe<<"\n";
		return maxe;		
	}
	
	template<typename T,typename T2>
	bool isUnitary(const psimag::Matrix<T>& S,const T2& eps)
	{
		size_t n = S.n_row();
		size_t n2 = S.n_col();
		T zzero = 0.0;
		T zone = 1.0;
		psimag::Matrix<T> C(n,n);
		psimag::BLAS::GEMM('N','C',n,n,n2,zone,&(S(0,0)),n,&(S(0,0)),n,zzero,&(C(0,0)),n);
		for (size_t i=0;i<n;i++) C(i,i) -= 1.0;
		T maxe = maxElement(C);
		if (psimag::norm(maxe)<eps) return true;
		throw std::runtime_error("isUnitary: false!!\n");
		return false;
		
	}
	
	template<typename T>
	psimag::Matrix<T> transposeConjugate(const psimag::Matrix<T>& A)
	{
		psimag::Matrix<T> ret(A.n_col(),A.n_row());
		for (size_t i=0;i<A.n_col();i++) for (size_t j=0;j<A.n_row();j++) ret(i,j)=conj(A(j,i));
		return ret;
	}
	
	template<typename T>
	psimag::Matrix<T> multiply(const psimag::Matrix<T>& A,const psimag::Matrix<T>& B)
	{
		psimag::Matrix<T> ret(A.n_row(),B.n_col());
		for (size_t i=0;i<A.n_row();i++) {
			for (size_t j=0;j<B.n_col();j++) {
				ret(i,j)=0;
				for (size_t k=0;k<A.n_col();k++) 
					ret(i,j) += A(i,k)*B(k,j);
			}
		}
		
		return ret;
	}
			
	
	template<class T>
	void operator*=(Matrix<T> &A,T const &v)
	{
		int i,j;
		for (i=0;i<A.n_row();i++) for (j=0;j<A.n_col();j++) A(i,j)*=v;
	}
	
	template<class T>
	void accumulate(Matrix<T> &A,Matrix<T> const &B)
	{
		size_t i,j;
		int nrow = B.n_row();
		int ncol = B.n_col();
		
		for (i=0;i<nrow;i++) for (j=0;j<ncol;j++) A(i,j) += B(i,j);
	}
	
	template<class T>
	bool isUnitary(Matrix<T> const &A)
	{
		size_t i,j,k;
		bool flag;
		Matrix<T> m(A.n_col(),A.n_col());
		double eps=1e-6;
		
		for (i=0;i<A.n_col();i++) {
			for (j=0;j<A.n_col();j++) {
				m(i,j)=static_cast<T>(0.0);
				for (k=0;k<A.n_row();k++) m(i,j) += A(k,i)*conj(A(k,j));
			}
		}
		flag=true;
		for (i=0;i<m.n_row();i++) {
			for (j=0;j<m.n_col();j++) {
				if (i==j) {
					if (norm(m(i,j)-1.0)>eps) {
						flag=false;
						break;
					}
				} else {
					if (norm(m(i,j)-0.0)>eps) {
						flag=false;
						break;
					}
				}
			}
		}
		return flag;
		
	}
	
	
	
	template<class T>
	void setValue(Matrix<T>  &A,int n,T const &value)
	{
		A.resize(n,n);
		for (size_t i=0;i<n;i++) A(i,i)=value;
		
	}
	
	template<class T>
	int matrixRank(Matrix<T> const &A)
	{
		int n = A.n_row();
		if (n!=A.n_col()) throw std::runtime_error("matrixRank: matrix must be square.\n");
		return n;
	}
	
	template<class T>
	void invert(psimag::Matrix<T> &A)
	{
		int n = A.n_row();
		int info;
		std::vector<int> ipiv(n);
		Matrix<T> identity;
		int i,j;

		// the identity
		identity.resize(n,n);
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				identity(i,j)=0;
				if (i==j) identity(i,i)=1.0;
			}
		}
		// AA^(-1) = I
		LAPACK::GESV(n,n,&(A(0,0)),n,&(ipiv[0]),&(identity(0,0)),n,info);
		if (info<0) {
			std::cerr<<"solver Linear: the "<<(-info)<<"-th argument had an illegal value.\n";
			throw std::runtime_error("psimag::invert\n");
		}  else if (info>0) {
			std::cerr<<"solve Linear: U(i,i) is exactly zero for i="<<info<<".  The factorization "
			"has been completed, but the factor U is exactly  singular,  so "
			"the solution could not be computed.";
			throw std::runtime_error("psimag::invert\n");
		}
        
        	for (i=0;i<n;i++) for (j=0;j<n;j++) A(i,j)=identity(i,j);


	}

	template<class T>
	void truncate(Matrix<T> &A,std::vector<size_t> const &removed,bool rowOption)
	{
		size_t j;
		int x=removed.size();
		if (x<=0) return;
		size_t nrow = A.n_row();
		size_t ncol = A.n_col();
		
		size_t n = ncol;
		if (rowOption)  n = nrow;
		
		if (int(n)<=x) {
			std::cerr<<"psimag::truncate: n="<<n<<" must be larger than x="<<x<<" rowoption="<<rowOption<<"\n";
			throw std::runtime_error("psimag::truncated\n");
		}
		
		std::vector<int> remap(n);
		
		
		//! find remapping
		j=0;
		for (size_t i=0;i<n;i++) {
			remap[i] = -1;
			if (utils::isInVector(removed,i)>=0) continue;
			remap[i]=j;
			j++;
		}
		if (j!=n-x) throw std::runtime_error("truncate: psimag::Matrix is throwing...\n");
		
		//! truncate
		if (rowOption) {
			Matrix<T> B(nrow-x,ncol);
			for (size_t i=0;i<ncol;i++) {
				for (j=0;j<nrow;j++) {
					if (remap[j]<0) continue;
					B(remap[j],i)=A(j,i);
				}
			}
			A=B; 
		} else {
			Matrix<T> B(nrow,ncol-x);
			for (size_t i=0;i<nrow;i++) {
				for (j=0;j<ncol;j++) {
					if (remap[j]<0) continue;
					B(i,remap[j])=A(i,j);
				}
			}
			A=B; 
		}
			
		
	
	}
	
	
	template<class T>
	bool isHermitian(Matrix<T> const &A,bool verbose=false)
	{
		size_t n=A.n_row();
		double eps=1e-6;
		
		if (n!=A.n_col()) throw std::runtime_error("isHermitian called on a non-square matrix.\n");
		for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) 
			if (psimag::norm(A(i,j)-conj(A(j,i)))>eps) {
				if (verbose) std::cerr<<"A("<<i<<","<<j<<")="<<A(i,j)<<" A("<<j<<","<<i<<")="<<A(j,i)<<"\n";
				return false;
			}
		return true;	
	}
	
	template<class T>
	void mathematicaPrint(std::ostream& os,Matrix<T> const &a)
	{
		
		os<<"{";
		for (size_t i=0;i<a.n_row();i++) {
			os<<"{";
			for (size_t j=0;j<a.n_col();j++) {
				if (fabs(a(i,j))>1e-6) os<<a(i,j);
				else os<<"0";
				if (j<a.n_col()-1) os<<",";
			}
			os<<"}";
			if (i<a.n_row()-1) os<<",";
			os<<"\n";
		}
		os<<"}\n";
	}
	
	
} //namespace Dmrg
/*@}*/
#endif


