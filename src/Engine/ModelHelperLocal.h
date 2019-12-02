/*
Copyright (c) 2009-2016-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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
#ifndef DMRG_MODELHELPER_H
#define DMRG_MODELHELPER_H

#include "PackIndices.h" // in PsimagLite
#include "Link.h"
#include "Concurrency.h"
#include "Vector.h"

/** \ingroup DMRG */
/*@{*/

/*! \file ModelHelperLocal.h
 *
 *  A class to contain state information about the Hamiltonian
 *  to help with the calculation of x+=Hy
 *
 */

namespace Dmrg {
template<typename LeftRightSuperType_>
class ModelHelperLocal {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef std::pair<SizeType,SizeType> PairType;

public:

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef typename LeftRightSuperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef typename OperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BasisType::RealType RealType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef Link<SparseElementType> LinkType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorSparseElementType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename BasisType::QnType QnType;
	typedef typename PsimagLite::Vector<OperatorStorageType*>::Type VectorOperatorStorageType;
	typedef typename PsimagLite::Vector<VectorOperatorStorageType>::Type
	VectorVectorOperatorStorageType;
	typedef PsimagLite::Concurrency ConcurrencyType;

	class Aux {

	public:

		Aux(SizeType m, const LeftRightSuperType& lrs) :
		    m_(m), buffer_(lrs.left().size())

		{
			createBuffer(lrs);
			createAlphaAndBeta(lrs);
		}

		SizeType m() const { return m_; }

		int buffer(SizeType i, SizeType j) const
		{
			assert(i < buffer_.size());
			assert(j < buffer_[i].size());
			return buffer_[i][j];
		}

		const PsimagLite::Vector<int>::Type& buffer(SizeType i) const
		{
			assert(i < buffer_.size());
			return buffer_[i];
		}

		SizeType alpha(SizeType i) const
		{
			assert(i < alpha_.size());
			return alpha_[i];
		}

		SizeType beta(SizeType i) const
		{
			assert(i < beta_.size());
			return beta_[i];
		}

		bool fermionSigns(SizeType i) const
		{
			assert(i < fermionSigns_.size());
			return fermionSigns_[i];
		}

	private:

		void createBuffer(const LeftRightSuperType& lrs)
		{
			SizeType ns = lrs.left().size();
			SizeType ne = lrs.right().size();
			int offset = lrs.super().partition(m_);
			int total = lrs.super().partition(m_+1) - offset;

			typename PsimagLite::Vector<int>::Type  tmpBuffer(ne);
			for (SizeType alphaPrime=0;alphaPrime<ns;alphaPrime++) {
				for (SizeType betaPrime=0;betaPrime<ne;betaPrime++) {
					tmpBuffer[betaPrime] = lrs.super().
					        permutationInverse(alphaPrime + betaPrime*ns) -offset;
					if (tmpBuffer[betaPrime]>=total) tmpBuffer[betaPrime]= -1;
				}
				buffer_[alphaPrime]=tmpBuffer;
			}
		}

		void createAlphaAndBeta(const LeftRightSuperType& lrs)
		{
			SizeType ns = lrs.left().size();
			int offset = lrs.super().partition(m_);
			int total = lrs.super().partition(m_+1) - offset;

			PackIndicesType pack(ns);
			alpha_.resize(total);
			beta_.resize(total);
			fermionSigns_.resize(total);
			for (int i=0;i<total;i++) {
				// row i of the ordered product basis
				pack.unpack(alpha_[i],beta_[i],lrs.super().permutation(i+offset));
				int fs = lrs.left().fermionicSign(alpha_[i],-1);
				fermionSigns_[i] = (fs < 0) ? true : false;
			}
		}

		SizeType m_;
		typename PsimagLite::Vector<PsimagLite::Vector<int>::Type>::Type buffer_;
		VectorSizeType alpha_;
		VectorSizeType beta_;
		typename PsimagLite::Vector<bool>::Type fermionSigns_;
	};

	ModelHelperLocal(const LeftRightSuperType& lrs)
	    : lrs_(lrs),
	      garbage_(ConcurrencyType::codeSectionParams.npthreads),
	      seen_(ConcurrencyType::codeSectionParams.npthreads)
	{
		ConcurrencyType::mutexInit(&mutex_);
	}

	~ModelHelperLocal()
	{
		const SizeType n = garbage_.size();
		for (SizeType i = 0; i < n; ++i) {
			const SizeType m = garbage_[i].size();
			for (SizeType j = 0; j < m; ++j) {
				delete garbage_[i][j];
				garbage_[i][j] = 0;
			}
		}

		ConcurrencyType::mutexDestroy(&mutex_);
	}

	void clearThreadSelves() const
	{
		threadSelves_.clear();
	}

	const OperatorStorageType& reducedOperator(char modifier,
	                                           SizeType i,
	                                           SizeType sigma,
	                                           const ProgramGlobals::SysOrEnvEnum type) const
	{

		assert(!BasisType::useSu2Symmetry());

		const OperatorStorageType* m = 0;
		PairType ii;
		if (type == ProgramGlobals::SysOrEnvEnum::SYSTEM) {
			ii = lrs_.left().getOperatorIndices(i,sigma);
			m = &(lrs_.left().getOperatorByIndex(ii.first).getStorage());
		} else {
			assert(type == ProgramGlobals::SysOrEnvEnum::ENVIRON);
			ii = lrs_.right().getOperatorIndices(i,sigma);
			m =&(lrs_.right().getOperatorByIndex(ii.first).getStorage());
		}

		m->checkValidity();
		if (modifier == 'N') return *m;

		assert(modifier == 'C');
		SizeType typeIndex = (type == ProgramGlobals::SysOrEnvEnum::SYSTEM) ? 0 : 1;
		SizeType packed = typeIndex + ii.first*2;
		const ConcurrencyType::PthreadtType threadSelf = ConcurrencyType::threadSelf();
		const SizeType threadNum = threadNumberFromSelf(threadSelf);

		if (garbage_.size() != seen_.size())
			err("reducedOperator: FATAL: internal error\n");

		if (garbage_.size() <= threadNum || seen_.size() <= threadNum)
			err("reducedOperator: FATAL: " + ttos(threadNum) + " >= " +
			    ttos(garbage_.size()) + "\n");

		int indexOfSeen = PsimagLite::indexOrMinusOne(seen_[threadNum], packed);
		if (indexOfSeen >= 0) {
			assert(static_cast<SizeType>(indexOfSeen) < garbage_[threadNum].size());
			return *(garbage_[threadNum][indexOfSeen]);
		}

		OperatorStorageType* mc = new OperatorStorageType;
		transposeConjugate(*mc, *m);
		garbage_[threadNum].push_back(mc);
		seen_[threadNum].push_back(packed);
		mc->checkValidity();

		return *mc;
	}

	static bool isSu2() { return false; }

	int size(SizeType mm) const
	{
		int tmp = lrs_.super().partition(mm + 1) - lrs_.super().partition(mm);
		return tmp; //reflection_.size(tmp);
	}

	const QnType& quantumNumber(SizeType mm) const
	{
		return lrs_.super().qnEx(mm);
	}

	//! Does matrixBlock= (AB), A belongs to pSprime and B
	// belongs to pEprime or viceversa (inter)
	void fastOpProdInter(SparseMatrixType const &A,
	                     SparseMatrixType const &B,
	                     SparseMatrixType &matrixBlock,
	                     const LinkType& link,
	                     const Aux& aux) const
	{
		RealType fermionSign = (link.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION)
		        ? -1 : 1;

		//! work only on partition m
		if (link.type==ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM)  {
			LinkType link2 = link;
			link2.value *= fermionSign;
			link2.type = ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON;
			fastOpProdInter(B, A, matrixBlock, link2, aux);
			return;
		}

		SizeType m = aux.m();
		int offset = lrs_.super().partition(m);
		int total = lrs_.super().partition(m + 1) - offset;
		int counter=0;
		matrixBlock.resize(total,total);

		int i;
		for (i=0;i<total;i++) {
			// row i of the ordered product basis
			matrixBlock.setRow(i,counter);
			int alpha = aux.alpha(i);
			int beta = aux.beta(i);

			SparseElementType fsValue = (fermionSign < 0 && aux.fermionSigns(i))
			        ? -link.value
			        : link.value;

			for (int k=A.getRowPtr(alpha);k<A.getRowPtr(alpha+1);k++) {
				int alphaPrime = A.getCol(k);
				for (int kk=B.getRowPtr(beta);kk<B.getRowPtr(beta+1);kk++) {
					int betaPrime = B.getCol(kk);
					int j = aux.buffer(alphaPrime, betaPrime);
					if (j<0) continue;
					/* fermion signs note:
					here the environ is applied first and has to "cross"
					the system, hence the sign factor pSprime.fermionicSign(alpha,tmp)
					*/
					SparseElementType tmp = A.getValue(k) * B.getValue(kk)*fsValue;
					//if (tmp==static_cast<MatrixElementType>(0.0)) continue;
					matrixBlock.pushCol(j);
					matrixBlock.pushValue(tmp);
					counter++;
				}
			}
		}

		matrixBlock.setRow(i,counter);
	}

	// Does x+= (AB)y, where A belongs to pSprime and B  belongs to pEprime or
	// viceversa (inter)
	// Has been changed to accomodate for reflection symmetry
	void fastOpProdInter(VectorSparseElementType& x,
	                     const VectorSparseElementType& y,
	                     const SparseMatrixType& A,
	                     const SparseMatrixType& B,
	                     const LinkType& link,
	                     const Aux& aux) const
	{
		RealType fermionSign =  (link.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION)
		        ? -1 : 1;

		if (link.type==ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM)  {
			LinkType link2 = link;
			link2.value *= fermionSign;
			link2.type = ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON;
			fastOpProdInter(x,y,B,A,link2);
			return;
		}

		//! work only on partition m
		int m = aux.m();
		int offset = lrs_.super().partition(m);
		int total = lrs_.super().partition(m + 1) - offset;

		for (int i=0;i<total;++i) {
			// row i of the ordered product basis
			int alpha = aux.alpha(i);
			int beta = aux.beta(i);
			SparseElementType sum = 0.0;
			int startkk = B.getRowPtr(beta);
			int endkk = B.getRowPtr(beta+1);
			int startk = A.getRowPtr(alpha);
			int endk = A.getRowPtr(alpha+1);
			/* fermion signs note:
			 * here the environ is applied first and has to "cross"
			 * the system, hence the sign factor pSprime.fermionicSign(alpha,tmp)
			 */

			SparseElementType fsValue = (fermionSign < 0 && aux.fermionSigns(i))
			        ? -link.value
			        : link.value;

			for (int k=startk;k<endk;++k) {
				int alphaPrime = A.getCol(k);
				SparseElementType tmp2 = A.getValue(k) *fsValue;
				const typename PsimagLite::Vector<int>::Type& bufferTmp =
				        aux.buffer(alphaPrime);

				for (int kk=startkk;kk<endkk;++kk) {
					int betaPrime= B.getCol(kk);
					int j = bufferTmp[betaPrime];
					if (j<0) continue;

					SparseElementType tmp = tmp2 * B.getValue(kk);
					sum += tmp * y[j];
				}
			}

			x[i] += sum;
		}
	}

	// Let H_{alpha,beta; alpha',beta'} =
	// basis2.hamiltonian_{alpha,alpha'} \delta_{beta,beta'}
	// Let H_m be  the m-th block (in the ordering of basis1) of H
	// Then, this function does x += H_m * y
	// This is a performance critical function
	// Has been changed to accomodate for reflection symmetry
	void hamiltonianLeftProduct(VectorSparseElementType& x,
	                            const VectorSparseElementType& y,
	                            const Aux& aux) const
	{
		int m = aux.m();
		int offset = lrs_.super().partition(m);
		int i,k,alphaPrime;
		int bs = lrs_.super().partition(m+1)-offset;
		const SparseMatrixType& hamiltonian = lrs_.left().hamiltonian().getCRS();
		SizeType ns = lrs_.left().size();
		SparseElementType sum = 0.0;
		PackIndicesType pack(ns);
		for (i=0;i<bs;i++) {
			SizeType r,beta;
			pack.unpack(r,beta,lrs_.super().permutation(i+offset));

			// row i of the ordered product basis
			for (k=hamiltonian.getRowPtr(r);k<hamiltonian.getRowPtr(r+1);k++) {
				alphaPrime = hamiltonian.getCol(k);
				int j = aux.buffer(alphaPrime, beta);
				if (j<0) continue;
				sum += hamiltonian.getValue(k)*y[j];
			}

			x[i] += sum;
			sum = 0.0;
		}
	}

	// Let  H_{alpha,beta; alpha',beta'} =
	// basis2.hamiltonian_{beta,beta'} \delta_{alpha,alpha'}
	// Let H_m be  the m-th block (in the ordering of basis1) of H
	// Then, this function does x += H_m * y
	// This is a performance critical function
	void hamiltonianRightProduct(VectorSparseElementType& x,
	                             const VectorSparseElementType& y,
	                             const Aux& aux) const
	{
		int m = aux.m();
		int offset = lrs_.super().partition(m);
		int i,k;
		int bs = lrs_.super().partition(m+1)-offset;
		const SparseMatrixType& hamiltonian = lrs_.right().hamiltonian().getCRS();
		SizeType ns = lrs_.left().size();
		SparseElementType sum = 0.0;
		PackIndicesType pack(ns);
		for (i=0;i<bs;i++) {
			SizeType alpha,r;
			pack.unpack(alpha,r,lrs_.super().permutation(i+offset));

			// row i of the ordered product basis
			for (k=hamiltonian.getRowPtr(r);k<hamiltonian.getRowPtr(r+1);k++) {
				int j = aux.buffer(alpha, hamiltonian.getCol(k));
				if (j < 0) continue;
				sum += hamiltonian.getValue(k)*y[j];
			}

			x[i] += sum;
			sum = 0.0;
		}
	}

	// if option==true let H_{alpha,beta; alpha',beta'} =
	// basis2.hamiltonian_{alpha,alpha'} \delta_{beta,beta'}
	// if option==false let  H_{alpha,beta; alpha',beta'} =
	// basis2.hamiltonian_{beta,beta'} \delta_{alpha,alpha'}
	// returns the m-th block (in the ordering of basis1) of H
	// Note: USed only for debugging
	void calcHamiltonianPart(SparseMatrixType &matrixBlock,
	                         bool option,
	                         const Aux& aux) const
	{
		int m  = aux.m();
		SizeType offset = lrs_.super().partition(m);
		int k,alphaPrime=0,betaPrime=0;
		int bs = lrs_.super().partition(m+1)-offset;
		SizeType ns=lrs_.left().size();
		SparseMatrixType hamiltonian;
		if (option) {
			hamiltonian = lrs_.left().hamiltonian().getCRS();
		} else {
			hamiltonian = lrs_.right().hamiltonian().getCRS();
		}

		matrixBlock.resize(bs,bs);

		int counter=0;
		PackIndicesType pack(ns);
		for (SizeType i=offset;i<lrs_.super().partition(m+1);i++) {
			matrixBlock.setRow(i-offset,counter);
			SizeType alpha,beta;
			pack.unpack(alpha,beta,lrs_.super().permutation(i));
			SizeType r=beta;
			if (option) {
				betaPrime=beta;
				r = alpha;
			} else {
				alphaPrime=alpha;
				r=beta;
			}

			assert(r<hamiltonian.rows());
			// row i of the ordered product basis
			for (k=hamiltonian.getRowPtr(r);k<hamiltonian.getRowPtr(r+1);k++) {

				if (option) alphaPrime = hamiltonian.getCol(k);
				else 	    betaPrime  = hamiltonian.getCol(k);
				SizeType j = lrs_.super().permutationInverse(alphaPrime + betaPrime*ns);
				if (j<offset || j>=lrs_.super().partition(m+1)) continue;
				SparseElementType tmp = hamiltonian.getValue(k);
				matrixBlock.pushCol(j-offset);
				matrixBlock.pushValue(tmp);
				counter++;
			}
		}

		matrixBlock.setRow(lrs_.super().partition(m+1)-offset,counter);

#ifndef NDEBUG
		if (!isHermitian(matrixBlock)) {
			std::cerr<<matrixBlock.toDense();
			PsimagLite::String lOrR = (option) ? "Left" : "Right";
			err(lOrR + " Hamiltonian Matrix not Hermitian\n");
		}
#endif
	}

	const LeftRightSuperType& leftRightSuper() const
	{
		return lrs_;
	}

private:

	SizeType threadNumberFromSelf(ConcurrencyType::PthreadtType threadSelf) const
	{
		ConcurrencyType::mutexLock(&mutex_);

		int threadPreNum = PsimagLite::indexOrMinusOne(threadSelves_, threadSelf);
		if (threadPreNum < 0) {
			threadPreNum = threadSelves_.size();
			threadSelves_.push_back(threadSelf);
		}

		ConcurrencyType::mutexUnlock(&mutex_);

		return threadPreNum;
	}

	const LeftRightSuperType& lrs_;
	mutable VectorVectorOperatorStorageType garbage_;
	mutable typename PsimagLite::Vector<BlockType>::Type seen_;
	mutable ConcurrencyType::MutexType mutex_;
	mutable PsimagLite::Vector<ConcurrencyType::PthreadtType>::Type threadSelves_;
}; // class ModelHelperLocal
} // namespace Dmrg
/*@}*/

#endif

