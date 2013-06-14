/*
Copyright (c) 2012, UT-Battelle, LLC
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

/*! \file KronMatrix.h
 *
 *
 */

#ifndef KRON_MATRIX_HEADER_H
#define KRON_MATRIX_HEADER_H

#include "Matrix.h"
#include "KronConnections.h"

namespace Dmrg {

template<typename InitKronType>
class KronMatrix {

	typedef typename InitKronType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename InitKronType::GenGroupType GenGroupType;

public:

	KronMatrix(const InitKronType& initKron)
	: initKron_(initKron)
	{
		std::cout<<"KronMatrix: EXPERIMENTAL: preparation done for size="<<initKron.size()<<"\n";
	}

	void matrixVectorProduct(typename PsimagLite::Vector<ComplexOrRealType>::Type& vout,
				 const typename PsimagLite::Vector<ComplexOrRealType>::Type& vin) const
	{
		const typename PsimagLite::Vector<SizeType>::Type& permInverse = initKron_.lrs().super().permutationInverse();
		const typename PsimagLite::Vector<SizeType>::Type& perm = initKron_.lrs().super().permutationVector();
		const SparseMatrixType& left = initKron_.lrs().left().hamiltonian();
		const SparseMatrixType& right = initKron_.lrs().right().hamiltonian();
		SizeType nl = left.row();
		SizeType nr = right.row();
		SizeType nq = initKron_.size();
		SizeType offset = initKron_.offset();

		MatrixType V(nl,nr);
		for (SizeType i=0;i<nl;i++) {
			for (SizeType j=0;j<nr;j++) {
				SizeType r = permInverse[i+j*nl];
				if (r<offset || r>=offset+nq) continue;
				V(i,j) = vin[r-offset];
			}
		}

		MatrixType W(nl,nr);

		computeRight(W,V);
		computeLeft(W,V);
		computeConnections(W,V);

		for (SizeType r=0;r<vout.size();r++) {
			div_t divresult = div(perm[r+offset],nl);
			SizeType i = divresult.rem;
			SizeType j = divresult.quot;
			vout[r] += W(i,j);
		}
	}

private:

	void computeRight(MatrixType& W,const MatrixType& V) const
	{
		SizeType npatches = initKron_.patch();
		const GenGroupType& istartLeft = initKron_.istartLeft();
		const GenGroupType& istartRight = initKron_.istartRight();
		const ArrayOfMatStructType& artStruct = initKron_.aRt();

		for (SizeType ipatch = 0;ipatch<npatches;ipatch++) {
			SizeType i = initKron_.patch(GenIjPatchType::LEFT,ipatch);
			SizeType j = initKron_.patch(GenIjPatchType::RIGHT,ipatch);

			SizeType i1 = istartLeft(i);
			SizeType i2 = istartLeft(i+1);

			SizeType j1 = istartRight(j);
			//SizeType j2 = istartRight(j+1);

			const SparseMatrixType& tmp = artStruct(j,j);
			for (SizeType ii=i1;ii<i2;ii++) {
				for (SizeType mr=0;mr<tmp.row();mr++) {
					for (int kk=tmp.getRowPtr(mr);kk<tmp.getRowPtr(mr+1);kk++) {
						SizeType col = tmp.getCol(kk) + j1;
//						if (col>=i2) continue;
						W(ii,col) += V(ii, mr+j1) * tmp.getValue(kk);
					}
				}
			}
		}
	}

	void computeLeft(MatrixType& W,const MatrixType& V) const
	{
		SizeType npatches = initKron_.patch();
		const GenGroupType& istartLeft = initKron_.istartLeft();
		const GenGroupType& istartRight = initKron_.istartRight();
		const ArrayOfMatStructType& alStruct = initKron_.aL();

		for (SizeType ipatch = 0;ipatch<npatches;ipatch++) {
			SizeType i = initKron_.patch(GenIjPatchType::LEFT,ipatch);
			SizeType j = initKron_.patch(GenIjPatchType::RIGHT,ipatch);

			SizeType i1 = istartLeft(i);
//			SizeType i2 = istartLeft(i+1);

			SizeType j1 = istartRight(j);
			SizeType j2 = istartRight(j+1);

			const SparseMatrixType& tmp = alStruct(i,i);

			for (SizeType jj=j1;jj<j2;jj++) {
				for (SizeType mr=0;mr<tmp.row();mr++) {
					for (int kk=tmp.getRowPtr(mr);kk<tmp.getRowPtr(mr+1);kk++) {
						SizeType col = tmp.getCol(kk) + i1;
//						if (col>=i2) continue;
						W(mr+i1,jj) += tmp.getValue(kk) *  V(col, jj);
					}
				}
			}
		}
	}

	// ATTENTION: MPI is not supported, only pthreads
	void computeConnections(MatrixType& W,const MatrixType& V) const
	{
		typedef KronConnections<InitKronType> KronConnectionsType;
		KronConnectionsType kc(initKron_,W,V);
		typedef typename InitKronType::template ParallelConnectionsInner<KronConnectionsType> ParallelConnectionsInnerType;
		typedef typename ParallelConnectionsInnerType::Type ParallelConnectionsInnerTypeType;
		ParallelConnectionsInnerTypeType parallelConnections(PsimagLite::Concurrency::npthreads,
		                                                     PsimagLite::MPI::COMM_WORLD);

		SizeType npatches = initKron_.patch();
		parallelConnections.loopCreate(npatches,kc);
		//hc.sync(parallelConnections,concurrency_);
	}

	const InitKronType& initKron_;

}; //class KronMatrix

} // namespace PsimagLite

/*@}*/

#endif // KRON_MATRIX_HEADER_H
