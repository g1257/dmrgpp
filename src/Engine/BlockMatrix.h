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

/*! \file BlockMatrix.h
 *
 *  A class to represent a block diagonal matrix
 *
 */
#ifndef BLOCKMATRIX_HEADER_H
#define BLOCKMATRIX_HEADER_H
#include <vector>
#include <iostream>
#include "Matrix.h" // in psimag

namespace Dmrg {

	//! A block matrix class
	//! Blocks can be of any type and are templated with the type MatrixInBlockTemplate
	//
	template<class T,class MatrixInBlockTemplate>
	class BlockMatrix {
	
	public:
		typedef MatrixInBlockTemplate BuildingBlockType;
		
		BlockMatrix(int rank,int blocks) : rank_(rank),offsets_(blocks+1),data_(blocks) 
		{
			offsets_[blocks]=rank;
		}
		//! operator+= 
		void operator+=(BlockMatrix<T,MatrixInBlockTemplate> const &m)
		{
			BlockMatrix<T,MatrixInBlockTemplate> c;
			if (offsets_.size()<m.blocks()) operatorPlus(c,*this,m);
			else operatorPlus(c,m,*this);
			*this = c;
		}
		
		//! Sets all blocks of size 1 and value value
		void makeDiagonal(int n,T const &value)
		{
			rank_ = n;
			MatrixInBlockTemplate m;
			
			setValue(m,1,value);
			
			offsets_.clear();
			data_.clear();
			for (int i=0;i<n;i++) {
				offsets_.push_back(i);
				data_.push_back(m);
			}
			offsets_.push_back(n);
		}
		
		//! Set block to the zero matrix, 
		void setToZero(int n)
		{
			T value=static_cast<T>(0.0);
			makeDiagonal(n,value);
		}
		
		//! set block to the identity matrix
		void setToIdentity(int n)
		{
			T value=static_cast<T>(1.0);
			makeDiagonal(n,value);
		}
		
		void setBlock(int i,int offset,MatrixInBlockTemplate const &m) { data_[i]=m; offsets_[i]=offset; }
		
		void sumBlock(int i,MatrixInBlockTemplate const &m) { data_[i] += m; }
		
		int rank() const { return rank_; }
		
		int offsets(int i) const { return offsets_[i]; }
		
		size_t blocks() const { return data_.size(); }
		
		T operator()(int i,int j) const 
		{	
			int k;
			for (k=data_.size()-1;k>=0;k--) if (i>=offsets_[k]) break;
			
			if (j<offsets_[k] || j>=offsets_[k+1]) return static_cast<T>(0.0);
			return data_[k](i-offsets_[k],j-offsets_[k]);
			
		}
		
		
		MatrixInBlockTemplate operator()(int i) const { return data_[i]; }
		
		template<class S,class MatrixInBlockTemplate2>
		friend void operatorPlus(BlockMatrix<S,MatrixInBlockTemplate2>  &C,BlockMatrix<S,MatrixInBlockTemplate2>  const &A,BlockMatrix<S,MatrixInBlockTemplate2> const &B);
		
		template<class S,class MatrixInBlockTemplate2>
		friend void operatorPlus(BlockMatrix<S,MatrixInBlockTemplate2>   &A,BlockMatrix<S,MatrixInBlockTemplate2> const &B);
	
		template<class S,class MatrixInBlockTemplate2>
		friend std::ostream &operator<<(std::ostream &s,BlockMatrix<S,MatrixInBlockTemplate2> const &A);
		
		template<typename S,typename Field,typename ConcurrencyTemplate>
		friend void diagonalise(BlockMatrix<S,PsimagLite::Matrix<S> >  &C,std::vector<Field> &eigs,char option,ConcurrencyTemplate &concurrency);
		
	private:
		int rank_; //the rank of this matrix
		std::vector<int> offsets_; //starting of diagonal offsets for each block
		std::vector<MatrixInBlockTemplate> data_; // data on each block	
	
	}; // class BlockMatrix

	// Companion Functions
	template<class S,class MatrixInBlockTemplate>
	std::ostream &operator<<(std::ostream &s,BlockMatrix<S,MatrixInBlockTemplate> const &A)
	{
		for (size_t m=0;m<A.blocks();m++) {
			int nrank = A.offsets(m+1)-A.offsets(m);
			s<<"block number "<<m<<" has rank "<<nrank<<"\n";
			s<<A.data_[m];
		}
		return s;
	}

	//C=A+ B where A.offsets is contained in B.offsets
	template<class S,class MatrixInBlockTemplate>
	void operatorPlus(BlockMatrix<S,MatrixInBlockTemplate>  &C,
			  BlockMatrix<S,MatrixInBlockTemplate>  const &A,BlockMatrix<S,MatrixInBlockTemplate> const &B)
	{
		size_t i;
		int counter=0;
		C.rank_ = A.rank_;
		C.offsets_=A.offsets_;
		C.data_.resize(A.data_.size());
		for (i=0;i<A.data_.size();i++) {
			C.data_[i] = A.data_[i];
		
			while(B.offsets_[counter] <A.offsets_[i+1]) {
				accumulate(C.data_[i],B.data_[counter++]);
				if (counter>=B.offsets_.size()) break;
			}
			if (counter>=B.offsets_.size() && i<A.offsets_.size()-1) throw std::runtime_error("operatorPlus: restriction not met.\n");	
		}
	}

	//A+= B where A.offsets is contained in B.offsets
	template<class S,class MatrixInBlockTemplate>
	void operatorPlus(BlockMatrix<S,MatrixInBlockTemplate>   &A,BlockMatrix<S,MatrixInBlockTemplate> const &B)
	{
		size_t i;
		int counter=0;

		for (i=0;i<A.data_.size();i++) {
			while(B.offsets_[counter] < A.offsets_[i+1]) {
				accumulate(A.data_[i],B.data_[counter++]);
				if (counter>=B.offsets_.size()) break;
			}
			if (counter>=B.offsets_.size() && i<A.offsets_.size()-1) throw std::runtime_error("operatorPlus: restriction not met.\n");	
		}
	 	
	}

	//! Parallel version of the diagonalization of a block diagonal matrix
	template<typename S,typename Field,typename ConcurrencyTemplate>
	void diagonalise(BlockMatrix<S,PsimagLite::Matrix<S> >  &C,std::vector<Field> &eigs,char option,ConcurrencyTemplate &concurrency)
	{
		std::vector<Field> eigsTmp;
		std::vector<std::vector<Field> > eigsForGather;
		std::vector<size_t> weights(C.blocks());
		
		eigsForGather.resize(C.blocks());
		size_t m;
		for (m=0;m<C.blocks();m++) {
			eigsForGather[m].resize(C.offsets(m+1)-C.offsets(m));
			weights[m] =  C.offsets(m+1)-C.offsets(m);
		}
		
		eigs.resize(C.rank());
		
		concurrency.loopCreate(C.blocks(),weights);
		
		while(concurrency.loop(m)) {
			utils::diag(C.data_[m],eigsTmp,option);
			utils::enforcePhase(C.data_[m]);
			for (int j=C.offsets(m);j< C.offsets(m+1);j++) eigsForGather[m][j-C.offsets(m)] = eigsTmp[j-C.offsets(m)];
		}
		
		concurrency.gather(C.data_);
		concurrency.gather(eigsForGather);
		
		for (m=0;m<C.blocks();m++) {
			for (int j=C.offsets(m);j< C.offsets(m+1);j++) eigs[j]=eigsForGather[m][j-C.offsets(m)];
		}

		concurrency.broadcast(eigs);
		concurrency.broadcast(C.data_);
	}

	template<class S,class MatrixInBlockTemplate>
	bool isUnitary(BlockMatrix<S,MatrixInBlockTemplate> const &B)
	{
		bool flag=true;
		MatrixInBlockTemplate matrixTmp;
		
		for (size_t m=0;m<B.blocks();m++) {
			matrixTmp = B(m);
			if (!isUnitary(matrixTmp)) flag=false;
		}
		return flag;
	}

	template<class S>
	void blockMatrixToFullMatrix(PsimagLite::Matrix<S> &fm,BlockMatrix<S,PsimagLite::Matrix<S> > const &B)
	{
		int n = B.rank();
		int i,j;
		fm.reset(n,n);
		for (j=0;j<n;j++) for (i=0;i<n;i++) fm(i,j)=B(i,j);
	}

} // namespace Dmrg
/*@}*/	

#endif
