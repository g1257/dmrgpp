/*
Copyright (c) 2009,-2015 UT-Battelle, LLC
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
/** \ingroup DMRG */
/*@{*/
/** \file Parallel4PointDs.h
*/

#ifndef PARALLEL_4POINT_DS_H
#define PARALLEL_4POINT_DS_H

#include "Matrix.h"
#include "Mpi.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename ModelType,typename FourPointCorrelationsType>
class Parallel4PointDs {

	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename FourPointCorrelationsType::MatrixType MatrixType;
	typedef typename FourPointCorrelationsType::BraketType BraketType;
	typedef typename MatrixType::value_type FieldType;
	typedef typename FourPointCorrelationsType::SparseMatrixType SparseMatrixType;
	typedef PsimagLite::Concurrency ConcurrencyType;

public:

	enum FourPointModeEnum {MODE_NORMAL, MODE_THIN, MODE_THINupdn};

	typedef typename ModelType::RealType RealType;

	Parallel4PointDs(MatrixType& fpd,
	                 const FourPointCorrelationsType& fourpoint,
	                 const ModelType& model,
	                 const typename PsimagLite::Vector<SizeType>::Type& gammas,
	                 const typename PsimagLite::Vector<PairType>::Type& pairs,
	                 FourPointModeEnum mode)
	    : fpd_(fpd),
	      fourpoint_(fourpoint),
	      model_(model),
	      gammas_(gammas),
	      pairs_(pairs),
	      mode_(mode)
	{}

	void doTask(SizeType taskNumber, SizeType)
	{
		SizeType i = pairs_[taskNumber].first;
		SizeType j = pairs_[taskNumber].second;

		fpd_(i,j) = (mode_ == MODE_NORMAL) ? fourPointDelta(2*i, 2*j, gammas_)
		                                   : fourPointThin(i, j);
		if (mode_ == MODE_NORMAL) {
			fpd_(i,j) = fourPointDelta(2*i,2*j,gammas_);
		} else if (mode_ == MODE_THIN) {
			fpd_(i,j) = fourPointThin(i, j);
		} else if (mode_ == MODE_THINupdn) {
			fpd_(i,j) = fourPointThinupdn(i, j);
		} else {
			throw PsimagLite::RuntimeError("Parallel4PointDs: No matching mode_ found \n");
		}
	}

	SizeType tasks() const { return pairs_.size(); }

private:

	FieldType fourPointDelta(SizeType i,
	                         SizeType j,
	                         const typename PsimagLite::Vector<SizeType>::Type& gammas) const
	{
		SizeType hs = model_.hilbertSize(0);
		SizeType nx = 0;
		while (hs) {
			hs>>=1;
			nx++;
		}

		nx /= 2;
		SizeType site = 0;
		// C_{gamma0,up}
		PsimagLite::String str("<gs|c[" + ttos(site) + "]?" + ttos(gammas[0] + 0*nx) + "';");
		// const SparseMatrixType& opC0 = model.naturalOperator("c",site,gammas[0] + 0*nx).data;

		// C_{gamma1,down}
		str += "<gs|c[" + ttos(site) + "]?" + ttos(gammas[1] + 1*nx) + "';";
		// const SparseMatrixType& opC1 = model.naturalOperator("c",site,gammas[1] + 1*nx).data;

		// C_{gamma2,down}
		str += "<gs|c[" + ttos(site) + "]?" + ttos(gammas[2] + 1*nx) + ";";
		//const SparseMatrixType& opC2 = model.naturalOperator("c",site,gammas[2] + 1*nx).data;

		// C_{gamma3,up}
		str += "<gs|c[" + ttos(site) + "]?" + ttos(gammas[3] + 0*nx) + "|gs>";
		//const SparseMatrixType& opC3 = model.naturalOperator("c",site,gammas[3] + 0*nx).data;

		BraketType braket(model_, str);
		return fourpoint_(i,i+1,j,j+1,braket);
	}

	FieldType fourPointThin(SizeType i, SizeType j) const
	{
		SizeType number1 = fpd_.n_row()/2;
		SizeType spin0 = i/number1;
		SizeType tmp = i % number1;
		SizeType number2 = sqrt(number1);
		SizeType thini2 = tmp/number2;
		SizeType thini1 = tmp % number2;

		SizeType spin1 = j/number1;
		tmp = j % number1;
		SizeType thinj2 = tmp/number2;
		SizeType thinj1 = tmp % number2;
		//int sign = gammas[0] - 1;

		SizeType site = 0;

		// c(i1,orb1,spin0)
		PsimagLite::String str = "<gs|c?"+ ttos(spin0) + "[" + ttos(site) + "];";
		// SparseMatrixType O1 = model_.naturalOperator("c",site,spin0).data;
		// c(i2,orb2,1-spin0)

		str += "c?"+ ttos(1- spin0) + "[" + ttos(site) + "];";
		//SparseMatrixType O2 = model_.naturalOperator("c",site,1-spin0).data;

		// c(i2,orb2,spin1)
		str += "c?"+ ttos(spin1) + "[" + ttos(site) + "];";
		// SparseMatrixType O3 = model_.naturalOperator("c",site,spin1).data;

		// c(i3,orb1,1-spin1)
		str += "c?"+ ttos(1 - spin1) + "[" + ttos(site) + "]|gs>";
		//SparseMatrixType O4 = model_.naturalOperator("c",site,1-spin1).data;


//		SizeType val = spin0 + spin1 + 1;
//		int signTerm = (val & 1) ? sign : 1;
//		FieldType fourval = fourpoint_(thini1,thini2,thinj1,thinj2,braket);
//		return signTerm*fourval;

		BraketType braket(model_, str);
		FieldType fourval = fourpoint_(thini1,thini2,thinj1,thinj2,braket);
		return fourval;
	}

	FieldType fourPointThinupdn(SizeType i, SizeType j) const
	{
		SizeType number1 = fpd_.n_row()/2;
		SizeType spin0 = i/number1;
		SizeType tmp = i % number1;
		SizeType number2 = sqrt(number1);
		SizeType thini2 = tmp/number2;
		SizeType thini1 = tmp % number2;

		SizeType spin1 = j/number1;
		tmp = j % number1;
		SizeType thinj2 = tmp/number2;
		SizeType thinj1 = tmp % number2;
		//int sign = gammas[0] - 1;

		SizeType site = 0;

		// c(i1,orb1,spin0)
		PsimagLite::String str = "<gs|c?"+ ttos(spin0) + "[" + ttos(site) + "];";

		// c(i2,orb2,spin0)
		str += "c?"+ ttos(spin0) + "[" + ttos(site) + "];";

		// c(i2,orb2,spin1)
		str += "c?"+ ttos(spin1) + "[" + ttos(site) + "];";

		// c(i3,orb1,spin1)
		str += "c?"+ ttos(spin1) + "[" + ttos(site) + "]|gs>";

		BraketType braket(model_, str);
		FieldType fourval = fourpoint_(thini1,thini2,thinj1,thinj2,braket);
		return fourval;
	}

	MatrixType& fpd_;
	const FourPointCorrelationsType& fourpoint_;
	const ModelType& model_;
	const typename PsimagLite::Vector<SizeType>::Type& gammas_;
	const typename PsimagLite::Vector<PairType>::Type& pairs_;
	const FourPointModeEnum mode_;
}; // class Parallel4PointDs
} // namespace Dmrg

/*@}*/
#endif // PARALLEL_4POINT_DS_H

