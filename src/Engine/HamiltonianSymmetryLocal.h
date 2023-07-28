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

/*! \file HamiltonianSymmetryLocal.h
 *
 *  This class contains the implementation of local symmetries
 *  An object of this class is contained by DmrgBasisImplementation
 *
 */
#ifndef HAM_SYMM_LOCAL_H
#define HAM_SYMM_LOCAL_H
#include "Io/IoSelector.h"
#include "ProgramGlobals.h"
#include "Sort.h"

namespace Dmrg
{
template <typename SparseMatrixType>
class HamiltonianSymmetryLocal
{

	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename PsimagLite::Real<SparseElementType>::Type RealType;
	typedef PsimagLite::CrsMatrix<RealType> FactorsType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	SizeType getFlavor(SizeType) const
	{
		return 0; // meaningless
	}

	template <typename SolverParametersType>
	void calcRemovedIndices(VectorSizeType& removedIndices,
	    const VectorSizeType& perm,
	    SizeType kept,
	    const SolverParametersType&) const
	{
		const SizeType permSize = perm.size();
		if (permSize <= kept)
			return;

		SizeType target = permSize - kept;

		removedIndices.clear();
		for (SizeType i = 0; i < target; ++i) {
			if (removedIndices.size() >= target)
				break;
			if (PsimagLite::indexOrMinusOne(removedIndices, perm[i]) >= 0)
				continue;
			removedIndices.push_back(perm[i]);
		}
	}

	template <typename IoInputter>
	void read(IoInputter&,
	    PsimagLite::String,
	    bool,
	    typename PsimagLite::EnableIf<
		PsimagLite::IsInputLike<IoInputter>::True,
		int>::Type
	    = 0)
	{
	}

	void write(PsimagLite::IoSelector::Out& io,
	    PsimagLite::String label,
	    PsimagLite::IoNgSerializer::WriteMode mode) const
	{
		io.write(0, label + "FACTORSSIZE", mode);
	}
}; // class HamiltonianSymmetryLocal
} // namespace Dmrg

/*@}*/
#endif
