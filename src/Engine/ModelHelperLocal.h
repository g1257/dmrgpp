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
#ifndef MODELHELPER_LOC_HEADER_H
#define MODELHELPER_LOC_HEADER_H

#include "PackIndices.h" // in PsimagLite
#include "Link.h"

/** \ingroup DMRG */
/*@{*/

/*! \file ModelHelperLocal.h
 *
 *  A class to contain state information about the Hamiltonian to help with the calculation of x+=Hy
 *
 */

namespace Dmrg {
	template<typename LeftRightSuperType_,typename ConcurrencyType_>
	class ModelHelperLocal {

		typedef PsimagLite::PackIndices PackIndicesType;
		typedef std::pair<SizeType,SizeType> PairType;

	public:
		typedef LeftRightSuperType_ LeftRightSuperType;
		typedef typename LeftRightSuperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef typename OperatorsType::BasisType BasisType;
		typedef typename BasisType::BlockType BlockType;
		typedef typename BasisType::RealType RealType;
		typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
		typedef Link<SparseElementType,RealType> LinkType;

		enum { System=0,Environ=1 };

		ModelHelperLocal(SizeType m,
		                 const LeftRightSuperType& lrs,
		                 bool useReflection=false)
		: m_(m),
		  lrs_(lrs),
		  buffer_(lrs_.left().size()),
		  basis2tc_(lrs_.left().numberOfOperators()),
		  basis3tc_(lrs_.right().numberOfOperators())
		{
			createBuffer();
			createTcOperators(basis2tc_,lrs_.left());
			createTcOperators(basis3tc_,lrs_.right());
			createAlphaAndBeta();
		}

		SizeType m() const { return m_; }

		static bool isSu2() { return false; }

		const SparseMatrixType& getReducedOperator(char modifier,
		                                          SizeType i,
		                                          SizeType sigma,
		                                          SizeType type) const
		{
			if (modifier=='N') {
				if (type==System) {
					PairType ii =lrs_.left().getOperatorIndices(i,sigma); 
					return lrs_.left().getOperatorByIndex(ii.first).data;
				} else {
					PairType ii =lrs_.right().getOperatorIndices(i,sigma);
					return lrs_.right().getOperatorByIndex(ii.first).data;
				}
			}
			return getTcOperator(i,sigma,type);
		}

		int size() const
		{
			int tmp = lrs_.super().partition(m_+1)-lrs_.super().partition(m_);
			return tmp; //reflection_.size(tmp);
		}

		int quantumNumber() const
		{
			int state = lrs_.super().partition(m_);
			return lrs_.super().qn(state);
		}

		//! Does matrixBlock= (AB), A belongs to pSprime and B
		// belongs to pEprime or viceversa (inter)
		void fastOpProdInter(SparseMatrixType const &A,
		                     SparseMatrixType const &B,
		                     SparseMatrixType &matrixBlock,
		                     const LinkType& link,
		                     bool flipped = false) const
		{
			RealType fermionSign =(link.fermionOrBoson==ProgramGlobals::FERMION) ? -1 : 1;

			//! work only on partition m
			if (link.type==ProgramGlobals::ENVIRON_SYSTEM)  {
				LinkType link2 = link;
				link2.value *= fermionSign;
				link2.type = ProgramGlobals::SYSTEM_ENVIRON;
				fastOpProdInter(B,A,matrixBlock,link2,true);
				return;
			}

			int m = m_;
			int offset = lrs_.super().partition(m);
			int total = lrs_.super().partition(m+1) - offset;
			int counter=0;
			matrixBlock.resize(total,total);

			int i;
			for (i=0;i<total;i++) {
				// row i of the ordered product basis
				matrixBlock.setRow(i,counter);
				int alpha=alpha_[i];
				int beta=beta_[i];

				for (int k=A.getRowPtr(alpha);k<A.getRowPtr(alpha+1);k++) {
					int alphaPrime = A.getCol(k);
					for (int kk=B.getRowPtr(beta);kk<B.getRowPtr(beta+1);kk++) {
						int betaPrime= B.getCol(kk);
						int j = buffer_[alphaPrime][betaPrime];
						if (j<0) continue;
						/* fermion signs note:
						   here the environ is applied first and has to "cross"
						   the system, hence the sign factor pSprime.fermionicSign(alpha,tmp)
						  */
						SparseElementType tmp = A.getValue(k) * B.getValue(kk)*link.value;
						if (link.fermionOrBoson == ProgramGlobals::FERMION)
							tmp *= lrs_.left().fermionicSign(alpha,int(fermionSign));
						//if (tmp==static_cast<MatrixElementType>(0.0)) continue;
						matrixBlock.pushCol(j);
						matrixBlock.pushValue(tmp);
						counter++;
					}
				}
			}
			matrixBlock.setRow(i,counter);
		}

		//! Does x+= (AB)y, where A belongs to pSprime and B  belongs to pEprime or viceversa (inter)
		//! Has been changed to accomodate for reflection symmetry
		void fastOpProdInter(typename PsimagLite::Vector<SparseElementType>::Type&x,
		                     const typename PsimagLite::Vector<SparseElementType>::Type&y,
		                     SparseMatrixType const &A,
		                     SparseMatrixType const &B,
		                     const LinkType& link,
		                     bool flipped = false) const
		{
			RealType fermionSign =  (link.fermionOrBoson==ProgramGlobals::FERMION) ? -1 : 1;

			if (link.type==ProgramGlobals::ENVIRON_SYSTEM)  {
				LinkType link2 = link;
				link2.value *= fermionSign;
				link2.type = ProgramGlobals::SYSTEM_ENVIRON;
				fastOpProdInter(x,y,B,A,link2,true);
				return;
			}

			//! work only on partition m
			int m = m_;
			int offset = lrs_.super().partition(m);
			int total = lrs_.super().partition(m+1) - offset;

			for (int i=0;i<total;i++) {
				// row i of the ordered product basis
				int alpha=alpha_[i];
				int beta=beta_[i];
				SparseElementType& xSubI = x[i];
				int startkk = B.getRowPtr(beta);
				int endkk = B.getRowPtr(beta+1);
				int startk = A.getRowPtr(alpha);
				int endk = A.getRowPtr(alpha+1);
				/* fermion signs note:
				 *   here the environ is applied first and has to "cross"
				 *   the system, hence the sign factor pSprime.fermionicSign(alpha,tmp)
				 */
				RealType fs = lrs_.left().fermionicSign(alpha,int(fermionSign));
				SparseElementType fsValue = (link.fermionOrBoson == ProgramGlobals::FERMION)
				        ? fs*link.value
				        : link.value;

				for (int k=startk;k<endk;k++) {
					int alphaPrime = A.getCol(k);
					SparseElementType tmp2 = A.getValue(k) *fsValue;
					const typename PsimagLite::Vector<int>::Type& bufferTmp = buffer_[alphaPrime];

					for (int kk=startkk;kk<endkk;kk++) {
						int betaPrime= B.getCol(kk);
						int j = bufferTmp[betaPrime];
						if (j<0) continue;

						SparseElementType tmp = tmp2 * B.getValue(kk);
						xSubI += tmp * y[j];
					}
				}
			}
		}

		//! Let H_{alpha,beta; alpha',beta'} = basis2.hamiltonian_{alpha,alpha'} \delta_{beta,beta'}
		//! Let H_m be  the m-th block (in the ordering of basis1) of H
		//! Then, this function does x += H_m * y
		//! This is a performance critical function
		//! Has been changed to accomodate for reflection symmetry
		void hamiltonianLeftProduct(typename PsimagLite::Vector<SparseElementType> ::Type& x,
		                            const typename PsimagLite::Vector<SparseElementType>::Type& y) const
		{
			int m = m_;
			int offset = lrs_.super().partition(m);
			int i,k,alphaPrime;
			int bs = lrs_.super().partition(m+1)-offset;
			const SparseMatrixType& hamiltonian = lrs_.left().hamiltonian();
			SizeType ns = lrs_.left().size();

			PackIndicesType pack(ns);
			for (i=0;i<bs;i++) {
				SizeType r,beta;
				pack.unpack(r,beta,lrs_.super().permutation(i+offset));

				// row i of the ordered product basis
				for (k=hamiltonian.getRowPtr(r);k<hamiltonian.getRowPtr(r+1);k++) {
					alphaPrime = hamiltonian.getCol(k);
					int j = buffer_[alphaPrime][beta];
					if (j<0) continue;
					x[i]+= hamiltonian.getValue(k)*y[j];
				}
			}
		}

		//! Let  H_{alpha,beta; alpha',beta'} = basis2.hamiltonian_{beta,beta'} \delta_{alpha,alpha'}
		//! Let H_m be  the m-th block (in the ordering of basis1) of H
		//! Then, this function does x += H_m * y
		//! This is a performance critical function
		void hamiltonianRightProduct(typename PsimagLite::Vector<SparseElementType>::Type& x,
		                             const typename PsimagLite::Vector<SparseElementType>::Type& y) const
		{
			int m = m_;
			int offset = lrs_.super().partition(m);
			int i,k;
			int bs = lrs_.super().partition(m+1)-offset;
			const SparseMatrixType& hamiltonian = lrs_.right().hamiltonian();
			SizeType ns = lrs_.left().size();

			PackIndicesType pack(ns);
			for (i=0;i<bs;i++) {
				SizeType alpha,r;
				pack.unpack(alpha,r,lrs_.super().permutation(i+offset));

				// row i of the ordered product basis
				for (k=hamiltonian.getRowPtr(r);k<hamiltonian.getRowPtr(r+1);k++) {
					int j = buffer_[alpha][hamiltonian.getCol(k)];
					if (j<0) continue;
					x[i]+= hamiltonian.getValue(k)*y[j];
				}
			}
		}

		//! if option==true let H_{alpha,beta; alpha',beta'} = basis2.hamiltonian_{alpha,alpha'} \delta_{beta,beta'}
		//! if option==false let  H_{alpha,beta; alpha',beta'} = basis2.hamiltonian_{beta,beta'} \delta_{alpha,alpha'}
		//! returns the m-th block (in the ordering of basis1) of H
		//! Note: USed only for debugging
		void calcHamiltonianPart(SparseMatrixType &matrixBlock,
		                         bool option) const
		{
			int m  = m_;
			SizeType offset = lrs_.super().partition(m);
			int k,alphaPrime=0,betaPrime=0;
			int bs = lrs_.super().partition(m+1)-offset;
			SizeType ns=lrs_.left().size();
			SparseMatrixType hamiltonian;
			if (option) {
				hamiltonian = lrs_.left().hamiltonian();
			} else {
				hamiltonian = lrs_.right().hamiltonian();
			}
			PsimagLite::Matrix<SparseElementType> fullm;
			crsMatrixToFullMatrix(fullm,hamiltonian);
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
				assert(r<hamiltonian.row());
				// row i of the ordered product basis
				for (k=hamiltonian.getRowPtr(r);k<hamiltonian.getRowPtr(r+1);k++) {

					if (option) alphaPrime = hamiltonian.getCol(k);
					else 	    betaPrime  = hamiltonian.getCol(k);
					SizeType j = lrs_.super().permutationInverse(alphaPrime + betaPrime * ns);
					if (j<offset || j>=lrs_.super().partition(m+1)) continue;
					SparseElementType tmp = hamiltonian.getValue(k);
					matrixBlock.pushCol(j-offset);
					matrixBlock.pushValue(tmp);
					counter++;
				}
			}
			matrixBlock.setRow(lrs_.super().partition(m+1)-offset,counter);
		}

		const LeftRightSuperType& leftRightSuper() const
		{
			return lrs_;
		}

	private:
		int m_;
		const LeftRightSuperType&  lrs_;
		typename PsimagLite::Vector<PsimagLite::Vector<int>::Type>::Type buffer_;
		typename PsimagLite::Vector<SparseMatrixType>::Type basis2tc_,basis3tc_;
		typename PsimagLite::Vector<SizeType>::Type alpha_,beta_;

		const SparseMatrixType& getTcOperator(int i,SizeType sigma,SizeType type) const
		{
			if (type==System) {
				PairType ii =lrs_.left().getOperatorIndices(i,sigma);
				assert(ii.first<basis2tc_.size());
				return basis2tc_[ii.first];
			}
			PairType ii =lrs_.right().getOperatorIndices(i,sigma);
			assert(ii.first<basis3tc_.size());
			return basis3tc_[ii.first];
		}

		void createBuffer()
		{
			SizeType ns=lrs_.left().size();
			SizeType ne=lrs_.right().size();
			int offset = lrs_.super().partition(m_);
			int total = lrs_.super().partition(m_+1) - offset;

			typename PsimagLite::Vector<int>::Type  tmpBuffer(ne);
			for (SizeType alphaPrime=0;alphaPrime<ns;alphaPrime++) {
				for (SizeType betaPrime=0;betaPrime<ne;betaPrime++) {
					tmpBuffer[betaPrime] =lrs_.super().permutationInverse(alphaPrime + betaPrime*ns) -offset;
					if (tmpBuffer[betaPrime]>=total) tmpBuffer[betaPrime]= -1;
				}
				buffer_[alphaPrime]=tmpBuffer;
			}
		}

		void createTcOperators(typename PsimagLite::Vector<SparseMatrixType>::Type& basistc,
							   const BasisWithOperatorsType& basis)
		{
			if (basistc.size()==0) return;
			SizeType n=basis.getOperatorByIndex(0).data.row();
			bool b = true;
			for (SizeType i=0;i<basistc.size();i++) {
				if (basis.getOperatorByIndex(i).data.row()!=n) {
					b=false;
					break;
				}
			}
			if (b) createTcOperatorsCached(basistc,basis);
			else createTcOperatorsSimple(basistc,basis);
		}

		void createTcOperatorsSimple(typename PsimagLite::Vector<SparseMatrixType>::Type& basistc,
		                       const BasisWithOperatorsType& basis)
		{
			for (SizeType i=0;i<basistc.size();i++)
				transposeConjugate(basistc[i],basis.getOperatorByIndex(i).data);
		}

		void createTcOperatorsCached(typename PsimagLite::Vector<SparseMatrixType>::Type& basistc,
							   const BasisWithOperatorsType& basis)
		{
			if (basistc.size()==0) return;
			SizeType n=basis.getOperatorByIndex(0).data.row();
			typename PsimagLite::Vector<PsimagLite::Vector<int>::Type>::Type col(n);
			typename PsimagLite::Vector<typename PsimagLite::Vector<typename SparseMatrixType::value_type>::Type>::Type value(n);
			for (SizeType i=0;i<basistc.size();i++) {
				const SparseMatrixType& tmp = basis.getOperatorByIndex(i).data;
				assert(tmp.row()==n);
				transposeConjugate(basistc[i],tmp,col,value);

			}
		}

		void createAlphaAndBeta()
		{
			SizeType ns=lrs_.left().size();
			int offset = lrs_.super().partition(m_);
			int total = lrs_.super().partition(m_+1) - offset;

			PackIndicesType pack(ns);
			alpha_.resize(total);
			beta_.resize(total);
			for (int i=0;i<total;i++) {
				// row i of the ordered product basis
				pack.unpack(alpha_[i],beta_[i],lrs_.super().permutation(i+offset));
			}
		}
	}; // class ModelHelperLocal
} // namespace Dmrg
/*@}*/

#endif

