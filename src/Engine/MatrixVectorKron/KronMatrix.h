/*
Copyright (c) 2012-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file KronMatrix.h
 *
 *
 */

#ifndef KRON_MATRIX_HEADER_H
#define KRON_MATRIX_HEADER_H

#include "Matrix.h"
#include "KronConnections.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "PsimagLite.h"
#include "ProgressIndicator.h"
#ifdef PLUGIN_SC
#include "DmrgppPluginSC.h"
#endif

namespace Dmrg {

template<typename InitKronType>
class KronMatrix {

	typedef typename InitKronType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef KronConnections<InitKronType> KronConnectionsType;
	typedef typename KronConnectionsType::MatrixType MatrixType;
	typedef typename KronConnectionsType::VectorType VectorType;
	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename GenIjPatchType::BasisType BasisType;

public:

	KronMatrix(InitKronType& initKron, PsimagLite::String name)
	    : initKron_(initKron),
	      progress_("KronMatrix")
	{
		PsimagLite::String str((initKron.loadBalance()) ? "true" : "false");
		PsimagLite::OstringStream msg;
		msg<<"KronMatrix: "<<name<<" sizes="<<initKron.size(InitKronType::NEW);
		msg<<" "<<initKron.size(InitKronType::OLD);
		msg<<" loadBalance "<<str;
		progress_.printline(msg, std::cout);
	}

	void matrixVectorProduct(VectorType& vout, const VectorType& vin) const
	{
//		if (usePluginSc_)
//			return matrixVectorPluginSc(vout, vin);

		initKron_.copyIn(vout, vin);

		KronConnectionsType kc(initKron_);

		typedef PsimagLite::Parallelizer<KronConnectionsType> ParallelizerType;
		ParallelizerType parallelConnections(PsimagLite::Concurrency::npthreads,
		                                     PsimagLite::MPI::COMM_WORLD);

		if (initKron_.loadBalance())
			parallelConnections.loopCreate(kc, initKron_.weightsOfPatchesNew());
		else
			parallelConnections.loopCreate(kc);

		kc.sync();

		initKron_.copyOut(vout);
	}

private:

	KronMatrix(const KronMatrix&);

	const KronMatrix& operator=(const KronMatrix&);

	InitKronType& initKron_;
	PsimagLite::ProgressIndicator progress_;
}; //class KronMatrix

} // namespace PsimagLite

/*@}*/

#endif // KRON_MATRIX_HEADER_H
