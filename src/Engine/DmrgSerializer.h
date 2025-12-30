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

#include "BasisTraits.hh"
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"
#include "FermionSign.h"
#include "Io/IoSelector.h"
#include "ProgramGlobals.h"

namespace Dmrg {
// Move also checkpointing from DmrgSolver to here (FIXME)
template <typename LeftRightSuperType, typename VectorWithOffsetType> class DmrgSerializer {

public:

	typedef DmrgSerializer<LeftRightSuperType, VectorWithOffsetType> ThisType;
	typedef typename LeftRightSuperType::SparseMatrixType            SparseMatrixType;
	typedef typename SparseMatrixType::value_type                    ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                    MatrixType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType      BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType               BasisType;
	typedef typename BasisType::VectorQnType                         VectorQnType;
	typedef typename BasisType::VectorSizeType                       VectorSizeType;
	typedef FermionSign                                              FermionSignType;
	typedef typename BasisType::RealType                             RealType;
	typedef BlockDiagonalMatrix<MatrixType>                          BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType>                           BlockOffDiagMatrixType;
	typedef
	    typename PsimagLite::Vector<typename PsimagLite::Vector<VectorWithOffsetType*>::Type>::
	        Type VectorVectorVectorWithOffsetType;

	DmrgSerializer(const FermionSignType&                  fS,
	               const FermionSignType&                  fE,
	               const LeftRightSuperType&               lrs,
	               const VectorVectorVectorWithOffsetType& wf,
	               const BlockDiagonalMatrixType&          transform,
	               ProgramGlobals::DirectionEnum           direction)
	    : fS_(fS)
	    , fE_(fE)
	    , wavefunction_(wf)
	    , ownWf_(false)
	    , transform_(transform)
	    , direction_(direction)
	{
		lrs_ = new LeftRightSuperType(lrs);
	}

	// used only by IoNg:
	template <typename IoInputType>
	DmrgSerializer(
	    IoInputType&       io,
	    PsimagLite::String prefix,
	    bool               bogus,
	    const BasisTraits& basisTraits,
	    bool               readOnDemand,
	    typename PsimagLite::EnableIf<PsimagLite::IsInputLike<IoInputType>::True, int>::Type
	    = 0)
	    : fS_(io, prefix + "/fS", bogus)
	    , fE_(io, prefix + "/fE", bogus)
	    , lrs_(new LeftRightSuperType(io, prefix, basisTraits))
	    , ownWf_(true)
	    , transform_(io, prefix + "/transform", readOnDemand)
	{
		if (bogus)
			return;

		try {
			SizeType nsectors = 0;
			io.read(nsectors, prefix + "/WaveFunction/Size");
			wavefunction_.resize(nsectors);
			for (SizeType i = 0; i < nsectors; ++i) {
				SizeType nexcited = 0;
				io.read(nexcited, prefix + "/WaveFunction/" + ttos(i) + "/Size");
				wavefunction_[i].resize(nexcited);
				for (SizeType j = 0; j < nexcited; ++j) {
					wavefunction_[i][j] = new VectorWithOffsetType;
					wavefunction_[i][j]->read(io,
					                          prefix + "/WaveFunction/"
					                              + ttos(i) + "/" + ttos(j));
				}
			}
		} catch (...) {
			wavefunction_.resize(1);
			wavefunction_[0].resize(1);
			wavefunction_[0][0] = new VectorWithOffsetType;
			wavefunction_[0][0]->read(io, prefix + "/WaveFunction");
			std::cerr << "WARNING: Outdated branch of execution?!\n";
		}

		io.read(direction_, prefix + "/direction");
	}

	~DmrgSerializer()
	{
		freeLrs();

		if (!ownWf_)
			return; // <<--- EARLY EXIT HERE

		const SizeType nsectors = wavefunction_.size();
		for (SizeType i = 0; i < nsectors; ++i) {
			const SizeType nexcited = wavefunction_[i].size();
			for (SizeType j = 0; j < nexcited; ++j) {
				delete wavefunction_[i][j];
				wavefunction_[i][j] = nullptr;
			}
		}
	}

	void freeLrs()
	{
		delete lrs_;
		lrs_ = nullptr;
	}

	template <typename SomeIoOutType>
	void write(
	    SomeIoOutType&                            io,
	    PsimagLite::String                        prefix,
	    typename BasisWithOperatorsType::SaveEnum option,
	    SizeType                                  numberOfSites,
	    SizeType                                  counter,
	    typename PsimagLite::EnableIf<PsimagLite::IsOutputLike<SomeIoOutType>::True, int>::Type
	    = 0) const
	{
		if (counter == 0)
			io.createGroup(prefix);

		io.write(counter + 1,
		         prefix + "/Size",
		         (counter == 0) ? SomeIoOutType::Serializer::NO_OVERWRITE
		                        : SomeIoOutType::Serializer::ALLOW_OVERWRITE);

		prefix += ("/" + ttos(counter));

		io.createGroup(prefix);

		fS_.write(io, prefix + "/fS");
		fE_.write(io, prefix + "/fE");
		assert(lrs_);
		bool minimizeWrite = (lrs_->super().block().size() == numberOfSites);
		lrs_->write(io, prefix, option, minimizeWrite);

		io.createGroup(prefix + "/WaveFunction");
		const SizeType nsectors = wavefunction_.size();
		io.write(nsectors, prefix + "/WaveFunction/Size");
		for (SizeType i = 0; i < nsectors; ++i) {
			const SizeType nexcited = wavefunction_[i].size();
			io.createGroup(prefix + "/WaveFunction/" + ttos(i));
			io.write(nexcited, prefix + "/WaveFunction/" + ttos(i) + "/Size");
			for (SizeType j = 0; j < nexcited; ++j)
				wavefunction_[i][j]->write(
				    io, prefix + "/WaveFunction/" + ttos(i) + "/" + ttos(j));
		}

		transform_.write(prefix + "/transform", io);
		io.write(direction_, prefix + "/direction");
	}

	const FermionSignType& fermionicSignLeft() const { return fS_; }

	const FermionSignType& fermionicSignRight() const { return fE_; }

	const LeftRightSuperType& leftRightSuper() const
	{
		assert(lrs_);
		return *lrs_;
	}

	const VectorWithOffsetType& psiConst(SizeType sectorIndex, SizeType levelIndex) const
	{
		if (sectorIndex >= wavefunction_.size())
			err(PsimagLite::String(__FILE__) + "::wavefunction(): sectorIndex = "
			    + ttos(sectorIndex) + ">=" + ttos(wavefunction_.size()) + "\n");

		if (levelIndex >= wavefunction_[sectorIndex].size())
			err(PsimagLite::String(__FILE__)
			    + "::wavefunction(): levelIndex = " + ttos(levelIndex)
			    + ">=" + ttos(wavefunction_[sectorIndex].size()) + "\n");

		return *(wavefunction_[sectorIndex][levelIndex]);
	}

	SizeType cols() const { return transform_.cols(); }

	SizeType rows() const { return transform_.rows(); }

	ProgramGlobals::DirectionEnum direction() const { return direction_; }

	void transform(SparseMatrixType& ret, const SparseMatrixType& O) const
	{
		transform(ret, O, transform_);
	}

	static void transform(SparseMatrixType&              ret,
	                      const SparseMatrixType&        O,
	                      const BlockDiagonalMatrixType& transformExternal)
	{
		BlockOffDiagMatrixType m(O, transformExternal.offsetsRows());
		static const SizeType  gemmRnb         = 0; // disable GemmR
		static const SizeType  threadsForGemmR = 1; // disable GemmR parallel
		m.transform(transformExternal, gemmRnb, threadsForGemmR);
		m.toSparse(ret);
	}

	// This is only called from DMRG and so lrs_ has never been freed
	SizeType centerOfOrthogonality() const
	{
		assert(lrs_);
		SizeType max = lrs_->left().block().size();
		assert(max > 0);
		const SizeType site = lrs_->left().block()[max - 1];
		return (direction_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? site
		                                                                    : site + 1;
	}

private:

	void fillOffsets(VectorSizeType& v, const BasisType& basis) const
	{
		SizeType n = basis.partition();
		if (n == 0)
			return;
		v.resize(n);
		for (SizeType i = 0; i < n; ++i)
			v[i] = basis.partition(i);
	}

	// Disallowing copy and assignment here:
	DmrgSerializer(const ThisType& ds);

	ThisType& operator=(const ThisType& ds);

	FermionSignType                  fS_;
	FermionSignType                  fE_;
	LeftRightSuperType*              lrs_;
	VectorVectorVectorWithOffsetType wavefunction_;
	bool                             ownWf_;
	BlockDiagonalMatrixType          transform_;
	ProgramGlobals::DirectionEnum    direction_;
}; // class DmrgSerializer
} // namespace Dmrg

/*@}*/
#endif
