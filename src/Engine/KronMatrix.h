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

	KronMatrix(const InitKronType& initKron,const std::string& options)
	: initKron_(initKron),options_(options)
	{
		std::cout<<"KronMatrix: preparation done for size="<<initKron.size()<<"\n";
	}

	void matrixVectorProduct(std::vector<ComplexOrRealType>& vout,
				 const std::vector<ComplexOrRealType>& vin) const
	{
		const std::vector<size_t>& permInverse = initKron_.lrs().super().permutationInverse();
		const std::vector<size_t>& perm = initKron_.lrs().super().permutationVector();
		const SparseMatrixType& left = initKron_.lrs().left().hamiltonian();
		const SparseMatrixType& right = initKron_.lrs().right().hamiltonian();
		size_t nl = left.row();
		size_t nr = right.row();
		size_t nq = initKron_.size();
		size_t offset = initKron_.offset();

		MatrixType V(nl,nr);
		for (size_t i=0;i<nl;i++) {
			for (size_t j=0;j<nr;j++) {
				size_t r = permInverse[i+j*nl];
				if (r<offset || r>=offset+nq) continue;
				V(i,j) = vin[r-offset];
			}
		}

		MatrixType W(nl,nr);

		if (options_.find("noright")==std::string::npos) computeRight(W,V);
		if (options_.find("noleft")==std::string::npos) computeLeft(W,V);
		if (options_.find("noconnections")==std::string::npos) computeConnections(W,V);

		for (size_t r=0;r<vout.size();r++) {
			div_t divresult = div(perm[r+offset],nl);
			size_t i = divresult.rem;
			size_t j = divresult.quot;
			vout[r] += W(i,j);
		}
	}

private:

	void computeRight(MatrixType& W,const MatrixType& V) const
	{
		size_t npatches = initKron_.patch();
		const GenGroupType& istartLeft = initKron_.istartLeft();
		const GenGroupType& istartRight = initKron_.istartRight();
		const ArrayOfMatStructType& artStruct = initKron_.aRt();

		for (size_t ipatch = 0;ipatch<npatches;ipatch++) {
			size_t i = initKron_.patch(GenIjPatchType::LEFT,ipatch);
			size_t j = initKron_.patch(GenIjPatchType::RIGHT,ipatch);

			size_t i1 = istartLeft(i);
			size_t i2 = istartLeft(i+1);

			size_t j1 = istartRight(j);
			//size_t j2 = istartRight(j+1);

			const SparseMatrixType& tmp = artStruct(j,j);
			for (size_t ii=i1;ii<i2;ii++) {
				for (size_t mr=0;mr<tmp.row();mr++) {
					for (int kk=tmp.getRowPtr(mr);kk<tmp.getRowPtr(mr+1);kk++) {
						size_t col = tmp.getCol(kk) + j1;
//						if (col>=i2) continue;
						W(ii,col) += V(ii, mr+j1) * tmp.getValue(kk);
					}
				}
			}
		}
	}

	void computeLeft(MatrixType& W,const MatrixType& V) const
	{
		size_t npatches = initKron_.patch();
		const GenGroupType& istartLeft = initKron_.istartLeft();
		const GenGroupType& istartRight = initKron_.istartRight();
		const ArrayOfMatStructType& alStruct = initKron_.aL();

		for (size_t ipatch = 0;ipatch<npatches;ipatch++) {
			size_t i = initKron_.patch(GenIjPatchType::LEFT,ipatch);
			size_t j = initKron_.patch(GenIjPatchType::RIGHT,ipatch);

			size_t i1 = istartLeft(i);
			size_t i2 = istartLeft(i+1);

			size_t j1 = istartRight(j);
			size_t j2 = istartRight(j+1);

			const SparseMatrixType& tmp = alStruct(i,i);

			for (size_t jj=j1;jj<j2;jj++) {
				for (size_t mr=0;mr<tmp.row();mr++) {
					for (int kk=tmp.getRowPtr(mr);kk<tmp.getRowPtr(mr+1);kk++) {
						size_t col = tmp.getCol(kk) + i1;
						if (col>=i2) continue;
						W(mr+i1,jj) += tmp.getValue(kk) *  V(col, jj);
					}
				}
			}
		}
	}

	void computeConnections(MatrixType& W,const MatrixType& V) const
	{
		size_t npatches = initKron_.patch();
		size_t nC = initKron_.connections();
		const GenGroupType& istartLeft = initKron_.istartLeft();
		const GenGroupType& istartRight = initKron_.istartRight();
		MatrixType intermediate(W.n_row(),W.n_col());

		for (size_t outPatch=0;outPatch<npatches;outPatch++) {
			for (size_t inPatch=0;inPatch<npatches;inPatch++) {
				for (size_t ic=0;ic<nC;ic++) {
					const ComplexOrRealType& val = initKron_.value(ic);
					const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
					const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

					size_t i = initKron_.patch(GenIjPatchType::LEFT,inPatch);
					size_t j = initKron_.patch(GenIjPatchType::RIGHT,inPatch);

					size_t ip = initKron_.patch(GenIjPatchType::LEFT,outPatch);
					size_t jp = initKron_.patch(GenIjPatchType::RIGHT,outPatch);

					size_t i1 = istartLeft(i);
//					size_t i2 = istartLeft(i+1);

					size_t j1 = istartRight(j);
					size_t j2 = istartRight(j+1);

					size_t ip1 = istartLeft(ip);
//					size_t ip2 = istartLeft(ip+1);

					size_t jp1 = istartRight(jp);
//					size_t jp2 = istartLeft(jp+1);

					const SparseMatrixType& tmp1 =  xiStruct(ip,i);
					const SparseMatrixType& tmp2 =  yiStruct(j,jp);

					size_t colsize = j2 -j1;
					for (size_t mr2=0;mr2<colsize;mr2++)
						for (size_t mr=0;mr<tmp1.row();mr++)
							intermediate(mr,mr2)=0.0;

					for (size_t mr=0;mr<tmp1.row();mr++) {
						for (int k3=tmp1.getRowPtr(mr);k3<tmp1.getRowPtr(mr+1);k3++) {
							size_t col3 = tmp1.getCol(k3)+i1;
							ComplexOrRealType valtmp = val * tmp1.getValue(k3);
							for (size_t mr2=j1;mr2<j2;mr2++) {
								intermediate(mr,mr2-j1) += valtmp * V(col3,mr2);
							}
						}
					}

					for (size_t mr=0;mr<tmp1.row();mr++) {
						for (size_t mr2=0;mr2<colsize;mr2++) {
							size_t start = tmp2.getRowPtr(mr2);
							size_t end = tmp2.getRowPtr(mr2+1);
							for (size_t k4=start;k4<end;k4++) {
								size_t col4 = tmp2.getCol(k4)+jp1;
								W(mr+ip1,col4) += intermediate(mr,mr2) * tmp2.getValue(k4) ;
							}
						}
					}
				}
			}
		}
	}

	const InitKronType& initKron_;
	std::string options_;
}; //class KronMatrix

} // namespace PsimagLite

/*@}*/

#endif // KRON_MATRIX_HEADER_H
