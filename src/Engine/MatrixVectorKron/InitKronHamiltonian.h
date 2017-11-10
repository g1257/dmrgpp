/*
Copyright (c) 2009-2017, UT-Battelle, LLC
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

/*! \file InitKronHamiltonian.h
 *
 *
 */
#ifndef INITKRON_HAMILTONIAN_H
#define INITKRON_HAMILTONIAN_H
#include "ProgramGlobals.h"
#include "InitKronBase.h"
#include "Vector.h"

namespace Dmrg {

template<typename ModelType_>
class InitKronHamiltonian : public InitKronBase<typename ModelType_::LeftRightSuperType> {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef ModelType_ ModelType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef InitKronBase<LeftRightSuperType> BaseType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BaseType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename PsimagLite::Vector<ArrayOfMatStructType*>::Type VectorArrayOfMatStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;

	InitKronHamiltonian(const ModelType& model,
	                       const ModelHelperType& modelHelper)
	    : BaseType(modelHelper_.leftRightSuper()),
	      model_(model),
	      modelHelper_(modelHelper)
	{
		addHlAndHr();
		convertXcYcArrays();
	}

	bool loadBalance() const
	{
			return (model_.params().options.find("KronLoadBalance") != PsimagLite::String::npos);
	}

private:

	void addHlAndHr()
	{
		const RealType value = 1.0;
		const SparseMatrixType& aL = modelHelper_.leftRightSuper().left().hamiltonian();
		const SparseMatrixType& aR = modelHelper_.leftRightSuper().right().hamiltonian();
		identityL_.makeDiagonal(aL.rows(), value);
		identityR_.makeDiagonal(aR.rows(), value);
		std::pair<SizeType, SizeType> ops(0,0);
		std::pair<char, char> mods('n', 'n');
		LinkType link(0,
		              0,
		              ProgramGlobals::SYSTEM_SYSTEM,
		              value,
		              0,
		              ProgramGlobals::BOSON,
		              ops,
		              mods,
		              1,
		              value,
		              0);
		addOneConnection(aL,identityR_,link);
		addOneConnection(identityL_,aR,link);
	}

	void convertXcYcArrays()
	{
		SizeType total = model_.getLinkProductStruct(modelHelper_);

		for (SizeType ix=0;ix<total;ix++) {
			SparseMatrixType const* A = 0;
			SparseMatrixType const* B = 0;

			LinkType link2 = model_.getConnection(&A,&B,ix,modelHelper_);
			if (link2.type==ProgramGlobals::ENVIRON_SYSTEM)  {
				LinkType link3 = link2;
				link3.type = ProgramGlobals::SYSTEM_ENVIRON;
				if (link3.fermionOrBoson == ProgramGlobals::FERMION)
					link3.value *= -1.0;
				addOneConnection(*B,*A,link3);
				continue;
			}

			addOneConnection(*A,*B,link2);
		}
	}

	InitKronHamiltonian(const InitKronHamiltonian&);

	InitKronHamiltonian& operator=(const InitKronHamiltonian&);

	const ModelType& model_;
	const ModelHelperType& modelHelper_;
	SparseMatrixType identityL_;
	SparseMatrixType identityR_;
};
} // namespace Dmrg

#endif // INITKRON_HAMILTONIAN_H
