/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file DmrgSerializer.h
 *
 *  Serialize dmrg data
 */
#ifndef DMRG_SERIAL_H
#define DMRG_SERIAL_H

#include "BLAS.h"
#include "IoSelector.h"
#include "FermionSign.h"
#include "ProgramGlobals.h"

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
	typedef FermionSign FermionSignType;
	typedef typename BasisType::RealType RealType;

	DmrgSerializer(const FermionSignType& fS,
		       const FermionSignType& fE,
		       const LeftRightSuperType& lrs,
		       const VectorType& wf,
		       const SparseMatrixType& transform,
		       ProgramGlobals::DirectionEnum direction)
		: fS_(fS),
		  fE_(fE),
		  lrs_(lrs),
		  wavefunction_(wf),
		  transform_(transform),
		  direction_(direction)
	{
		transposeConjugate(transformC_,transform_);
	}


	DmrgSerializer(typename PsimagLite::IoSelector::In& io,
	               bool bogus,
	               bool isObserveCode)
		: fS_(io,bogus),
		  fE_(io,bogus),
		  lrs_(io, isObserveCode)
	{
		if (bogus) return;
		PsimagLite::String s = "#WAVEFUNCTION_sites=";
		wavefunction_.load(io,s);
		s = "#TRANSFORM_sites=";
		io.readMatrix(transform_,s);
		transposeConjugate(transformC_,transform_);
		s = "#DIRECTION=";
		io.readline(direction_,s);
	}

	// Save to disk everything needed to compute any observable (OBSOLETE!!)
	template<typename IoOutputter>
	void save(IoOutputter& io,
	          SizeType option,
	          SizeType numberOfSites,
	          typename PsimagLite::EnableIf<
	          PsimagLite::IsOutputLike<IoOutputter>::True, int>::Type = 0) const
	{
		fS_.save(io);
		fE_.save(io);
		lrs_.save(io,option,numberOfSites);

		// save wavefunction
		PsimagLite::String label = "#WAVEFUNCTION_sites=";
		for (SizeType i=0;i<lrs_.super().block().size();i++) {
			label += ttos(lrs_.super().block()[i])+",";
		}
		//SparseVector<typename TargettingType::TargetVectorType::value_type> psiSparse(target.gs());
		wavefunction_.save(io,label);

		label = "#TRANSFORM_sites=";
		for (SizeType i=0;i<lrs_.left().block().size();i++) {
			label += ttos(lrs_.left().block()[i])+",";
		}
		io.printMatrix(transform_,label);
		PsimagLite::String s = "#DIRECTION="+ttos(direction_);
		io.printline(s);
//		io.print("#DIRECTION=",direction_);
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

	SizeType columns() const { return transform_.cols(); }

	SizeType rows() const { return transform_.rows(); }

	ProgramGlobals::DirectionEnum direction() const { return direction_; }

	SizeType site() const
	{
		if (direction_==ProgramGlobals::EXPAND_SYSTEM) return lrs_.right().block()[0]-1;
		else return lrs_.right().block()[0];
	}

	void transform(SparseMatrixType& ret,const SparseMatrixType& O) const
	{
		SparseMatrixType ret2;
		multiply(ret2,transformC_,O);
		multiply(ret,ret2,transform_);
	}

private:

	// Disallowing copy and assignment here:
	DmrgSerializer(const ThisType& ds);
	ThisType& operator=(const ThisType& ds);


	FermionSignType fS_,fE_;
	LeftRightSuperType lrs_;
	VectorType wavefunction_;
	SparseMatrixType transform_;
	SparseMatrixType transformC_;
	ProgramGlobals::DirectionEnum direction_;
}; // class DmrgSerializer
} // namespace Dmrg 

/*@}*/
#endif
