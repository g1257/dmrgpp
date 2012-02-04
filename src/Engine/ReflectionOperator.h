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

/*! \file ReflectionOperator
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
#ifndef REFLECTION_OPERATOR_H
#define REFLECTION_OPERATOR_H

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "LinearlyIndependentSet.h"
#include "LAPACK.h"
#include "Sort.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class ReflectionOperator {

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

	ReflectionOperator(const LeftRightSuperType& lrs,size_t n0,bool isEnabled,size_t expandSys)
	: lrs_(lrs),
	  n0_(n0),
	  progress_("ReflectionOperator",0),
	  plusSector_(0),
	  isEnabled_(isEnabled),
	  expandSys_(expandSys),
	  reflectedLeft_(n0_,n0_),
	  reflectedRight_(n0_,n0_)
	{
		size_t counter=0;
		for (size_t i=0;i<reflectedLeft_.rank();i++) {
			reflectedLeft_.setRow(i,counter);
			reflectedLeft_.pushCol(i);
			reflectedLeft_.pushValue(1);
			counter++;
		}

		reflectedLeft_.setRow(reflectedLeft_.rank(),counter);
		reflectedRight_=reflectedLeft_;
	}

	void check(const std::vector<size_t>& sectors)
	{
		if (!isEnabled_) return;
//		SparseMatrixType sSuper;
//		setS(sSuper);
		SparseMatrixType sSector;
		setSsector(sSector,sectors);
		updateReflected();
//		extractCurrentSector(sSector,sSuper,sectors);
		computeItems(sSector);
		checkTransform(sSector);
	}

	void changeBasis(const PsimagLite::Matrix<ComplexOrRealType>& transform,size_t direction)
	{
		if (!isEnabled_) return;
		SparseMatrixType newreflected;
		static PsimagLite::Matrix<ComplexOrRealType> cachedTransform;
		if (direction==expandSys_) {
//			changeBasis(newreflected,reflectedLeft_,transform);
//			reflectedLeft_ = newreflected;
			cachedTransform = transform;
		} else {
			changeBasis(newreflected,reflectedLeft_,cachedTransform,transform);
			reflectedLeft_ = newreflected;

			changeBasis(newreflected,reflectedRight_,transform,cachedTransform);
			reflectedRight_ = newreflected;
		}
//		SparseMatrixType one1;
//		multiply(one1,reflectedLeft_,reflectedLeft_);
//		printFullMatrix(one1,"Should be I for ReflectedLeft:");
//		assert(isTheIdentity(one1));

//		multiply(one1,reflectedRight_,reflectedRight_);
//		printFullMatrix(one1,"Should be I for ReflectedRight:");
//		assert(isTheIdentity(one1));
	}

	void transform(SparseMatrixType& matrixA,
		       SparseMatrixType& matrixB,
		       const SparseMatrixType& matrix) const
	 {
		assert(isEnabled_);

		SparseMatrixType rT;
		transposeConjugate(rT,transform_);

		SparseMatrixType tmp;
		multiply(tmp,matrix,rT);

		SparseMatrixType matrix2;
//		printFullMatrix(matrix,"OriginalHam");
		multiply(matrix2,transform_,tmp);
//		printFullMatrix(transform_,"transform");
//		printFullMatrix(matrix2,"FinalHam");
		split(matrixA,matrixB,matrix2);
	 }

	template<typename SomeVectorType>
	void setInitState(const SomeVectorType& initVector,
			  SomeVectorType& initVector1,
			  SomeVectorType& initVector2) const
	{
		assert(isEnabled_);
		size_t minusSector = initVector.size()-plusSector_;
		initVector1.resize(plusSector_);
		initVector2.resize(minusSector);
		for (size_t i=0;i<initVector.size();i++) {
			if (i<plusSector_) initVector1[i]=initVector[i];
			else  initVector2[i-plusSector_]=initVector[i];
		}
	}

	RealType setGroundState(VectorType& gs,
				const RealType& gsEnergy1,
				const VectorType& gsVector1,
				const RealType& gsEnergy2,
				const VectorType& gsVector2) const
	{
		assert(isEnabled_);
		size_t rank = gsVector1.size() + gsVector2.size();
		if (gsEnergy1<=gsEnergy2) {
			setGs(gs,gsVector1,rank,0);
			return gsEnergy1;
		}
		setGs(gs,gsVector2,rank,gsVector1.size());
		return gsEnergy2;
	}

	const LeftRightSuperType& leftRightSuper() const { return lrs_; }

	bool isEnabled() const { return isEnabled_; }

private:

	void checkTransform(const SparseMatrixType& sSector)
	{
		SparseMatrixType rT;
		transposeConjugate(rT,transform_);
		SparseMatrixType tmp3;
		multiply(tmp3,transform_,rT);

//		printFullMatrix(rT,"Transform");

		SparseMatrixType tmp4;
		multiply(tmp3,sSector,rT);
		multiply(tmp4,transform_,tmp3);
//		printFullMatrix(tmp4,"R S R^\\dagger");
	}

	void changeBasis(SparseMatrixType& newreflected,
			 const SparseMatrixType& reflected,
			 const PsimagLite::Matrix<ComplexOrRealType>& transform1,
			 const PsimagLite::Matrix<ComplexOrRealType>& transform2)
	{
		newreflected.resize(reflected.rank());
		size_t counter = 0;
		assert(reflected.rank()==transform1.n_row());
		assert(reflected.rank()==transform2.n_row());

		size_t total = reflected.rank();
		std::vector<int> ptr(total,-1);
		std::vector<size_t> index(total,0);
		std::vector<ComplexOrRealType> temp(total,0);
		std::cerr<<"transform1="<<transform1.n_row()<<"x"<<transform1.n_col()<<"\n";
		std::cerr<<"transform2="<<transform2.n_row()<<"x"<<transform2.n_col()<<"\n";

		for (size_t x=0;x<total;x++) {
			newreflected.setRow(x,counter);

			size_t itemp = 0;
			for (size_t xprime=0;xprime<transform1.n_row();xprime++) {
				ComplexOrRealType wl1 =  transform1(xprime,x);
//				ComplexOrRealType wl1 =  transform(x,xprime);
				for (int k=reflected.getRowPtr(xprime);k<reflected.getRowPtr(xprime+1);k++) {
					size_t xsecond = reflected.getCol(k);
					ComplexOrRealType r = reflected.getValue(k);
					for (size_t xthird=0;xthird<transform2.n_row();xthird++) {
//						ComplexOrRealType val = wl1 * r * transform(xthird,xsecond);
						ComplexOrRealType val = wl1 * r * transform2(xsecond,xthird);
						//if (isAlmostZero(val)) continue;
						if (ptr[xthird]<0) {
							ptr[xthird] = itemp;
							temp[ptr[xthird]] = val;
							index[ptr[xthird]] = xthird;
							itemp++;
						} else {
							temp[ptr[xthird]] += val;
						}
//						newreflected.pushCol(xthird);
//						newreflected.pushValue(val);
//						counter++;
					}
				}
			}
			for (size_t s=0;s<itemp;s++) {
				newreflected.pushCol(index[s]);
				newreflected.pushValue(temp[s]);
				ptr[index[s]] = -1;
			}
			counter += itemp;
		}
		newreflected.setRow(reflected.rank(),counter);
		newreflected.checkValidity();
	}

	void setSsector(SparseMatrixType& sSector,const std::vector<size_t>& sectors) const
	{
		assert(sectors.size()==1);
		size_t m = sectors[0];
		size_t offset = lrs_.super().partition(m);
		size_t total = lrs_.super().partition(m+1)-offset;
		sSector.resize(total);
		size_t counter = 0;
		size_t ns = lrs_.left().size();
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		PackIndicesType pack1(ns);
		assert(reflectedLeft_.rank()==ns/n0_);
		assert(reflectedRight_.rank()==ns/n0_);
		std::vector<int> ptr(total,-1);
		std::vector<size_t> index(total,0);
		std::vector<ComplexOrRealType> temp(total,0);

		for (size_t i=0;i<total;i++) {
			sSector.setRow(i,counter);
			size_t x = 0, y = 0;
			pack1.unpack(x,y,lrs_.super().permutation(i+offset));

			size_t x0=0,x1=0;
			pack2.unpack(x0,x1,lrs_.left().permutation(x));

			size_t y0=0,y1=0;
			pack3.unpack(y0,y1,lrs_.right().permutation(y));

			size_t itemp = 0;

			for (int k=reflectedLeft_.getRowPtr(x0);k<reflectedLeft_.getRowPtr(x0+1);k++) {
				ComplexOrRealType val1 = reflectedLeft_.getValue(k);
				for (int k2=reflectedRight_.getRowPtr(y1);k2<reflectedRight_.getRowPtr(y1+1);k2++) {
					ComplexOrRealType val2 = val1 * reflectedRight_.getValue(k2);
					if (isAlmostZero(val2)) continue;
					size_t x0prime = reflectedLeft_.getCol(k);
					size_t xprime = pack3.pack(x1,x0prime,lrs_.right().permutationInverse());

					size_t y1prime = reflectedRight_.getCol(k2);
					size_t yprime = pack2.pack(y1prime,y0,lrs_.left().permutationInverse());

					size_t iprime = pack1.pack(yprime,xprime,lrs_.super().permutationInverse());
					assert(iprime>=offset && iprime<offset+total);
					size_t col = iprime - offset;
					if (ptr[col]<0) {
						ptr[col] = itemp;
						temp[ptr[col]] = val2;
						index[ptr[col]] = col;
						itemp++;
					} else {
						temp[ptr[col]] += val2;
					}
//					sSector.pushCol(iprime-offset);
//					sSector.pushValue(val2);
//					counter++;
				}
			}
			for (size_t s=0;s<itemp;s++) {
				sSector.pushCol(index[s]);
				sSector.pushValue(temp[s]);
				ptr[index[s]] = -1;
				counter++;
			}
		}
		sSector.setRow(total,counter);
		sSector.checkValidity();

//		printSparseMatrix(sSector,"sSector",1e-3);

		SparseMatrixType tmp;
		multiply(tmp,sSector,sSector);
		//printFullMatrix(tmp,"sSector*sSector");
		assert(isTheIdentity(tmp,1e-12));
	}


	void printSparseMatrix(SparseMatrixType& s1,const std::string& label,const RealType& eps) const
	{
		std::cout<<label<<"\n";
		for (size_t i=0;i<s1.rank();i++) {
			for (int k=s1.getRowPtr(i);k<s1.getRowPtr(i+1);k++) {
				ComplexOrRealType val = s1.getValue(k);
				if (isAlmostZero(val,eps)) continue;
				std::cout<<(i+1)<<" "<<(1+s1.getCol(k))<<" "<<val<<"\n";
			}
		}
	}

	void updateReflected()
	{
		size_t ns = lrs_.left().size();
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		SparseMatrixType reflectedLeft(ns,ns);

		size_t counter = 0;
		for (size_t x=0;x<ns;x++) {
			reflectedLeft.setRow(x,counter);
			size_t x0 = 0, x1=0;
			pack2.unpack(x0,x1,lrs_.left().permutation(x));
			for (int k=reflectedLeft_.getRowPtr(x0);k<reflectedLeft_.getRowPtr(x0+1);k++) {
				size_t x0r = reflectedLeft_.getCol(k);
				size_t col = pack3.pack(x1,x0r,lrs_.right().permutationInverse());
				ComplexOrRealType val = reflectedLeft_.getValue(k);
				//if (isAlmostZero(val)) continue;
				reflectedLeft.pushCol(col);
				reflectedLeft.pushValue(val);
				counter++;
			}
		}
		reflectedLeft.setRow(ns,counter);
		reflectedLeft.checkValidity();
		reflectedLeft_ = reflectedLeft;

		SparseMatrixType reflectedRight(ns,ns);
		counter=0;
		for (size_t x=0;x<ns;x++) {
			size_t x0 = 0, x1=0;
			reflectedRight.setRow(x,counter);
			pack3.unpack(x0,x1,lrs_.right().permutation(x));
			for (int k=reflectedRight_.getRowPtr(x1);k<reflectedRight_.getRowPtr(x1+1);k++) {
				size_t x1r = reflectedRight_.getCol(k);
				ComplexOrRealType val = reflectedRight_.getValue(k);
				//if (isAlmostZero(val)) continue;
				size_t col = pack2.pack(x1r,x0,lrs_.left().permutationInverse());
				reflectedRight.pushCol(col);
				reflectedRight.pushValue(val);
				counter++;
			}
		}
		reflectedRight.setRow(ns,counter);
		reflectedRight.checkValidity();
		reflectedRight_ = reflectedRight;
	}

	void setGs(VectorType& gs,const VectorType& v,size_t rank,size_t offset) const
	{
		VectorType gstmp(rank,0);

		for (size_t i=0;i<v.size();i++) {
			gstmp[i+offset]=v[i];
		}
		SparseMatrixType rT;
		transposeConjugate(rT,transform_);
		multiply(gs,rT,gstmp);
	}

	void computeItems(const SparseMatrixType& sSector)
	{
		std::vector<size_t> ipIsolated,ipConnected;
		findIsolated(ipIsolated,ipConnected,sSector);

		std::cerr<<"Isolated nodes="<<ipIsolated.size()<<"\n";
		std::cerr<<"Connected nodes="<<ipConnected.size()<<"\n";

		LinearlyIndependentSetType lis(sSector.rank());

		std::vector<size_t> ilabel(ipConnected.size(),AVAILABLE);
		size_t firstColor = 2;
		SparseMatrixType A;
		generateAmatrix(A,ipConnected,sSector);
		size_t ncolor = gencolorLabel(ilabel,A,firstColor);

		gencolorCheck(ilabel,ncolor);

		std::vector<size_t> iperm;
		gencolorPerm(iperm,ilabel,firstColor,ncolor);

		for (size_t i=0;i<iperm.size();i++) {
			size_t ind = ipConnected[iperm[i]];
			ipIsolated.push_back(ind);
		}
		// FIXME THERE IS NO NEED TO SORT HERE
		std::vector<size_t> iperm2(ipIsolated.size());
		Sort<std::vector<size_t> > sort;
		sort.sort(ipIsolated,iperm2);
		for (size_t i=0;i<ipIsolated.size();i++)
			lis.fill(ipIsolated[i],sSector,1.0);

		plusSector_ = lis.size();

		std::vector<size_t> xminus;
		for (size_t i=0;i<iperm.size();i++) {
			size_t ind = ipConnected[iperm[i]];
			xminus.push_back(ind);
		}

		iperm2.resize(xminus.size());
		sort.sort(xminus,iperm2);
		for (size_t i=0;i<xminus.size();i++)
			lis.fill(xminus[i],sSector,-1.0);

		size_t minuses = lis.size()-plusSector_;

		std::ostringstream msg;
		msg<<plusSector_<<" +, "<<minuses<<" -.";
		progress_.printline(msg,std::cout);
		assert(minuses+plusSector_==sSector.rank());
		transform_ = lis.transform();
//		printFullMatrix(transform_,"transform");
	}

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
		for (size_t icolor=firstColor;icolor<ncolor;icolor++) {
			ilist.clear();
			findWithLabel(ilist,ilabel,icolor);
			ip += ilist.size();
		}
		//		 --------------------
		//		 generate permutation
		//		 --------------------

		iperm.resize(ip);
		ip = 0;
		for (size_t icolor=firstColor;icolor<ncolor;icolor++) {
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

	void classify(std::vector<size_t>& p,const SparseMatrixType& sSector,const RealType& eps) const
	{
		for (size_t i=0;i<sSector.rank();i++) {
			if (p[i]==AVAILABLE) p[i]=COLOR;
			if (p[i]==NOT_AVAILABLE) continue;
			for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
				ComplexOrRealType val = sSector.getValue(k);
				size_t col = sSector.getCol(k);
				if (i==col) continue; //val -= sector;
				if (isAlmostZero(val,eps)) continue;
				p[col]=NOT_AVAILABLE;
			}
		}
	}

//	size_t computeOneSector(LinearlyIndependentSetType& lis,
//				const SparseMatrixType& sSector,
//				RealType sector)
//	{
//		std::vector<size_t> p(sSector.rank(),AVAILABLE);
//		classify(p,sSector,1e-1);

//		std::vector<SparseVectorType> v(sSector.rank(),sSector.rank());

//		size_t r = 0;
//		for (size_t i=0;i<sSector.rank();i++) {
//			assert(p[i]==NOT_AVAILABLE || p[i]==COLOR);
//			if (sector<0) {
//				for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
//					size_t col = sSector.getCol(k);
//					if (isAlmostZero(sSector.getValue(k),1e-8)) continue;
//					if (i==col) {
//						p[i]=NOT_AVAILABLE;
//						break;
//					}
//				}
//			}
//			if (p[i]!=COLOR) continue;
//			r++;
//		}
//		std::cerr<<"sector="<<sector<<" SIZE="<<r<<"\n";
//		r = 0;

//		for (size_t i=0;i<sSector.rank();i++) {
//			if (p[i]!=COLOR) continue;
//			bool hasDiagonal = false;
//			for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
//				size_t col = sSector.getCol(k);
//				//if (isAlmostZero(sSector.getValue(k))) continue;
//				if (i==col) {
//					v[r].add(col,sSector.getValue(k)+sector);
//					hasDiagonal = true;
//				} else {
//					v[r].add(col,sSector.getValue(k));
//				}
//			}
//			if (!hasDiagonal) {
//				v[r].add(i,sector);
//			}
//			r++;
//		}
//		size_t rnew = 0;
//		for (size_t i=0;i<r;i++) {
//			v[i].sort();
//			RealType norma = PsimagLite::norm(v[i]);
//			if (isAlmostZero(norma)) continue;

//			v[i] *= (1.0/norma);

//			lis.fill(v[i]);
//			v[i].clear();
//			rnew++;
//		}

//		assert(r==rnew);
//		lis.clear();


//		return lis.size();
//	}

	void printFullMatrix(const SparseMatrixType& s,const std::string& name) const
	{
		PsimagLite::Matrix<ComplexOrRealType> fullm(s.rank(),s.rank());
		crsMatrixToFullMatrix(fullm,s);
		std::cout<<"--------->   "<<name<<" <----------\n";
		try {
//			mathematicaPrint(std::cout,fullm);
//			symbolicPrint(std::cout,fullm);
		} catch (std::exception& e) {
			std::cout<<fullm;
		}

		std::cout<<fullm;

	}

	void split(SparseMatrixType& matrixA,SparseMatrixType& matrixB,const SparseMatrixType& matrix) const
	{
		size_t counter = 0;
		matrixA.resize(plusSector_);
		for (size_t i=0;i<plusSector_;i++) {
			matrixA.setRow(i,counter);
			for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
				size_t col = matrix.getCol(k);
				ComplexOrRealType val = matrix.getValue(k);
				if (col<plusSector_) {
					matrixA.pushCol(col);
					matrixA.pushValue(val);
					counter++;
					continue;
				}
				if (!isAlmostZero(val,1e-5)) {
					std::string s(__FILE__);
					s += " Hamiltonian has no reflection symmetry.";
					throw std::runtime_error(s.c_str());
				}
			}
		}
		matrixA.setRow(plusSector_,counter);

		size_t rank = matrix.rank();
		size_t minusSector=rank-plusSector_;
		matrixB.resize(minusSector);
		counter=0;
		for (size_t i=plusSector_;i<rank;i++) {
			matrixB.setRow(i-plusSector_,counter);
			for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
				size_t col = matrix.getCol(k);
				ComplexOrRealType val = matrix.getValue(k);
				if (col>=plusSector_) {
					matrixB.pushCol(col-plusSector_);
					matrixB.pushValue(val);
					counter++;
					continue;
				}
				if (!isAlmostZero(val,1e-5)) {
					std::string s(__FILE__);
					s += " Hamiltonian has no reflection symmetry.";
					throw std::runtime_error(s.c_str());
				}
			}
		}
		matrixB.setRow(minusSector,counter);
	}

	//SparseMatrixType& operator()() const { return s_; }

	const LeftRightSuperType& lrs_;
	size_t n0_; // hilbert size of one site
	PsimagLite::ProgressIndicator progress_;
	size_t plusSector_;
	bool isEnabled_;
	size_t expandSys_;
	SparseMatrixType reflectedLeft_,reflectedRight_;
	SparseMatrixType transform_;
}; // class ReflectionOperator

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_OPERATOR_H
