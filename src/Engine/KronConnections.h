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

/*! \file KronConnections.h
 *
 *
 */

#ifndef KRON_CONNECTIONS_H
#define KRON_CONNECTIONS_H

#include "Matrix.h"

namespace Dmrg {

template<typename InitKronType>
class KronConnections {

	typedef typename InitKronType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename InitKronType::GenGroupType GenGroupType;

public:

	typedef typename InitKronType::RealType RealType;

	KronConnections(const InitKronType& initKron,MatrixType& W,const MatrixType& V)
	: initKron_(initKron),W_(W),V_(V),hasMpi_(PsimagLite::Concurrency::hasMpi())
	{
	}

	//! ATTENTION: ONLY VALID ON X86 AND X86_64 WHERE += IS ATOMIC
	void thread_function_(SizeType threadNum,SizeType blockSize,SizeType total,pthread_mutex_t* myMutex)
	{
		SizeType nC = initKron_.connections();
		const GenGroupType& istartLeft = initKron_.istartLeft();
		const GenGroupType& istartRight = initKron_.istartRight();
		MatrixType intermediate(W_.n_row(),W_.n_col());

		SizeType mpiRank = (hasMpi_) ? PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD) : 0;
		SizeType npthreads = PsimagLite::Concurrency::npthreads;

		for (SizeType p=0;p<blockSize;p++) {
			SizeType outPatch = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (outPatch>=total) break;

			for (SizeType inPatch=0;inPatch<total;++inPatch) {
				for (SizeType ic=0;ic<nC;++ic) {
					const ComplexOrRealType& val = initKron_.value(ic);
					const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
					const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

					SizeType i = initKron_.patch(GenIjPatchType::LEFT,inPatch);
					SizeType j = initKron_.patch(GenIjPatchType::RIGHT,inPatch);

					SizeType ip = initKron_.patch(GenIjPatchType::LEFT,outPatch);
					SizeType jp = initKron_.patch(GenIjPatchType::RIGHT,outPatch);

					SizeType i1 = istartLeft(i);
//					SizeType i2 = istartLeft(i+1);

					SizeType j1 = istartRight(j);
					SizeType j2 = istartRight(j+1);

					SizeType ip1 = istartLeft(ip);
//					SizeType ip2 = istartLeft(ip+1);

					SizeType jp1 = istartRight(jp);
//					SizeType jp2 = istartLeft(jp+1);

					const SparseMatrixType& tmp1 =  xiStruct(ip,i);
					const SparseMatrixType& tmp2 =  yiStruct(j,jp);

					SizeType colsize = j2 -j1;
					for (SizeType mr2=0;mr2<colsize;++mr2)
						for (SizeType mr=0;mr<tmp1.row();++mr)
							intermediate(mr,mr2)=0.0;

					for (SizeType mr=0;mr<tmp1.row();++mr) {
						for (int k3=tmp1.getRowPtr(mr);k3<tmp1.getRowPtr(mr+1);++k3) {
							SizeType col3 = tmp1.getCol(k3)+i1;
							ComplexOrRealType valtmp = val * tmp1.getValue(k3);
							for (SizeType mr2=j1;mr2<j2;++mr2) {
								intermediate(mr,mr2-j1) += valtmp * V_(col3,mr2);
							}
						}
					}

					for (SizeType mr2=0;mr2<colsize;++mr2) {
						SizeType start = tmp2.getRowPtr(mr2);
						SizeType end = tmp2.getRowPtr(mr2+1);
						for (SizeType k4=start;k4<end;++k4) {
							ComplexOrRealType value1 = tmp2.getValue(k4);
							SizeType col4 = tmp2.getCol(k4)+jp1;
							for (SizeType mr=0;mr<tmp1.row();++mr) {
								W_(mr+ip1,col4) += intermediate(mr,mr2) * value1;
							}
						}
					}
				}
			}
		}
	}

	void sync()
	{
		if (!hasMpi_) return;
		PsimagLite::MPI::allReduce(W_);
	}

	const InitKronType& initKron_;
	MatrixType& W_;
	const MatrixType& V_;
	bool hasMpi_;
}; //class KronConnections

} // namespace PsimagLite

/*@}*/

#endif // KRON_CONNECTIONS_H
