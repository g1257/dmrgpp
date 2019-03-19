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

/*! \file DmrgSerializer.h
 *
 *  Serialize dmrg data
 */
#ifndef DMRG_SERIAL_H
#define DMRG_SERIAL_H

#include "BLAS.h"
#include "Io/IoSelector.h"
#include "FermionSign.h"
#include "ProgramGlobals.h"
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"

namespace Dmrg {
// Move also checkpointing from DmrgSolver to here (FIXME)
template<typename LeftRightSuperType,typename VectorType>
class DmrgSerializer {

	typedef DmrgSerializer<LeftRightSuperType,VectorType> ThisType;
	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

public:

	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::VectorQnType VectorQnType;
	typedef typename BasisType::VectorSizeType VectorSizeType;
	typedef FermionSign FermionSignType;
	typedef typename BasisType::RealType RealType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType> BlockOffDiagMatrixType;

	DmrgSerializer(const FermionSignType& fS,
	               const FermionSignType& fE,
	               const LeftRightSuperType& lrs,
	               const VectorType& wf,
	               const BlockDiagonalMatrixType& transform,
	               ProgramGlobals::DirectionEnum direction)
	    : fS_(fS),
	      fE_(fE),
	      lrs_(lrs),
	      wavefunction_(wf),
	      transform_(transform),
	      direction_(direction)
	{}

	// used only by IoNg:
	template<typename IoInputType>
	DmrgSerializer(IoInputType& io,
	               PsimagLite::String prefix,
	               bool bogus,
	               bool isObserveCode,
	               typename PsimagLite::EnableIf<
	               PsimagLite::IsInputLike<IoInputType>::True, int>::Type = 0)
	    : fS_(io, prefix + "/fS", bogus),
	      fE_(io, prefix + "/fE", bogus),
	      lrs_(io, prefix, isObserveCode),
	      transform_(io, prefix + "/transform")
	{
		if (bogus) return;

		wavefunction_.read(io, prefix + "/WaveFunction");
		io.read(direction_, prefix + "/direction");
	}

	template<typename SomeIoOutType>
	void write(SomeIoOutType& io,
	           PsimagLite::String prefix,
	           typename BasisWithOperatorsType::SaveEnum option,
	           SizeType numberOfSites,
	           SizeType counter,
	           typename PsimagLite::EnableIf<
	           PsimagLite::IsOutputLike<SomeIoOutType>::True, int>::Type = 0) const
	{
		if (counter == 0) io.createGroup(prefix);

		io.write(counter + 1,
		         prefix + "/Size",
		         (counter == 0) ? SomeIoOutType::Serializer::NO_OVERWRITE :
		                          SomeIoOutType::Serializer::ALLOW_OVERWRITE);

		prefix += ("/" + ttos(counter));

		io.createGroup(prefix);

		fS_.write(io, prefix + "/fS");
		fE_.write(io, prefix + "/fE");
		bool minimizeWrite = (lrs_.super().block().size() == numberOfSites);
		lrs_.write(io, prefix, option, minimizeWrite);

		wavefunction_.write(io, prefix + "/WaveFunction");

		transform_.write(prefix + "/transform", io);
		io.write(direction_, prefix + "/direction");
	}

	const FermionSignType& fermionicSignLeft() const
	{
		return fS_;
	}

	const FermionSignType& fermionicSignRight() const
	{
		return fE_;
	}

	const LeftRightSuperType& leftRightSuper() const
	{
		return lrs_;
	}

	const VectorType& wavefunction() const { return wavefunction_; }

	SizeType cols() const
	{
		return transform_.cols();
	}

	SizeType rows() const
	{
		return transform_.rows();
	}

	ProgramGlobals::DirectionEnum direction() const { return direction_; }

	SizeType site() const
	{
		return (direction_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ?
		            lrs_.right().block()[0] - 1 : lrs_.right().block()[0];
	}

	void transform(SparseMatrixType& ret, const SparseMatrixType& O) const
	{
		BlockOffDiagMatrixType m(O, transform_.offsetsRows());
		m.transform(transform_);
		m.toSparse(ret);
	}

private:

	void fillOffsets(VectorSizeType& v, const BasisType& basis) const
	{
		SizeType n = basis.partition();
		if (n == 0) return;
		v.resize(n);
		for (SizeType i = 0; i < n; ++i)
			v[i] = basis.partition(i);
	}

	// Disallowing copy and assignment here:
	DmrgSerializer(const ThisType& ds);

	ThisType& operator=(const ThisType& ds);

	FermionSignType fS_;
	FermionSignType fE_;
	LeftRightSuperType lrs_;
	VectorType wavefunction_;
	BlockDiagonalMatrixType transform_;
	ProgramGlobals::DirectionEnum direction_;
}; // class DmrgSerializer
} // namespace Dmrg 

/*@}*/
#endif
