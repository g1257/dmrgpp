/*
Copyright (c) 2009-2012-2018, UT-Battelle, LLC
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

/*! \file ModelCommon.h
 *
 *  An abstract class to represent the strongly-correlated-electron models that
 *  can be used with the DmrgSolver
 *
 */

#ifndef MODEL_COMMON_H
#define MODEL_COMMON_H
#include <iostream>

#include "HamiltonianConnection.h"
#include "InputCheck.h"
#include "InputNg.h"
#include "LabeledOperators.h"
#include "ModelLinks.h"
#include "NoPthreads.h"
#include "Profiling.h"
#include "ProgressIndicator.h"
#include "Sort.h"
#include "Su2SymmetryGlobals.h"

namespace Dmrg {

template <typename ParametersType, typename SuperGeometryType, typename ModelHelperType>
class ModelCommon {

	using SparseMatrixType  = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType = typename SparseMatrixType::value_type;
	using LinkType          = typename ModelHelperType::LinkType;

public:

	using VectorType             = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using InputValidatorType     = PsimagLite::InputNg<InputCheck>::Readable;
	using OperatorsType          = typename ModelHelperType::OperatorsType;
	using OperatorType           = typename OperatorsType::OperatorType;
	using Block                  = typename ModelHelperType::BlockType;
	using RealType               = typename ModelHelperType::RealType;
	using MyBasis                = typename ModelHelperType::BasisType;
	using BasisWithOperatorsType = typename ModelHelperType::BasisWithOperatorsType;
	using LabeledOperatorsType   = LabeledOperators<OperatorType>;
	using ModelLinksType         = ModelLinks<LabeledOperatorsType, SuperGeometryType>;
	using HamiltonianConnectionType
	    = HamiltonianConnection<ModelLinksType, ModelHelperType, ParametersType>;
	using VectorLinkType       = typename HamiltonianConnectionType::VectorLinkType;
	using LeftRightSuperType   = typename ModelHelperType::LeftRightSuperType;
	using VectorOperatorType   = typename PsimagLite::Vector<OperatorType>::Type;
	using VectorVectorLinkType = typename PsimagLite::Vector<VectorLinkType>::Type;
	using VectorSizeType       = typename HamiltonianConnectionType::VectorSizeType;
	using VerySparseMatrixType = typename HamiltonianConnectionType::VerySparseMatrixType;

	ModelCommon(const ParametersType& params, const SuperGeometryType& superGeometry)
	    : params_(params)
	    , superGeometry_(superGeometry)
	    , progress_("ModelCommon")
	{
		Su2SymmetryGlobals<RealType>::init(ModelHelperType::isSu2());
		MyBasis::useSu2Symmetry(ModelHelperType::isSu2());
		if (params.options.isSet("OperatorsChangeAll"))
			OperatorsType::setChangeAll(true);
	}

	const ParametersType& params() const { return params_; }

	const SuperGeometryType& superGeometry() const { return superGeometry_; }

	void addConnectionsInNaturalBasis(SparseMatrixType&         hmatrix,
	                                  const VectorOperatorType& cm,
	                                  const Block&              block,
	                                  RealType                  time) const
	{
		if (block.size() != 1)
			err("addConnectionsInNaturalBasis(): unimplemented\n");
	}

	const ParametersType&         params_;
	const SuperGeometryType&      superGeometry_;
	PsimagLite::ProgressIndicator progress_;
}; // class ModelCommon
} // namespace Dmrg
/*@}*/
#endif
