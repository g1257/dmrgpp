/*
Copyright (c) 2009-2017-2018 UT-Battelle, LLC
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

/*! \file Kitaev.h
 *
 *  An implementation of the Kitaev model (started March 2018)
 *
 */

#ifndef DMRG_KITAEV_H
#define DMRG_KITAEV_H

#include <algorithm>
#include "ModelBase.h"
#include "ParametersKitaev.h"
#include "LinkProductKitaev.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramGlobals.h"
#include "ModelCommon.h"
#include "Utils.h"
#include "Complex.h"

namespace Dmrg {

template<typename ModelBaseType>
class Kitaev : public ModelBaseType {

public:

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::TargetQuantumElectronsType TargetQuantumElectronsType;
	typedef	typename ModelBaseType::VectorType VectorType;

private:

	typedef typename ModelBaseType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef unsigned int long WordType;
	typedef LinkProductKitaev<ModelHelperType> LinkProductType;
	typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelBaseType::VectorRealType VectorRealType;

	static const int TWICE_THE_SPIN = 1;

	enum InternalDir {DIR_X, DIR_Y, DIR_Z};

public:

	typedef typename PsimagLite::Vector<unsigned int long>::Type HilbertBasisType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::PairType PairType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;

	Kitaev(const SolverParamsType& solverParams,
	       InputValidatorType& io,
	       const GeometryType& geometry)
	    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
	      modelParameters_(io),
	      geometry_(geometry)
	{

		SizeType n = geometry_.numberOfSites();
		SizeType m = modelParameters_.magneticField.size();

		if (m > 0 && m != n) {
			PsimagLite::String msg("Kitaev: If provided, ");
			msg += " MagneticField must be a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (BasisType::useSu2Symmetry()) {
			PsimagLite::String msg("Kitaev does not have SU(2) symmetry\n");
			throw PsimagLite::RuntimeError(msg);
		}
	}

	SizeType memResolv(PsimagLite::MemResolv& mres,
	                   SizeType,
	                   PsimagLite::String msg = "") const
	{
		return 0;
	}

	void print(std::ostream& os) const { os<<modelParameters_; }

	SizeType hilbertSize(SizeType) const { return TWICE_THE_SPIN + 1; }

	void setQuantumNumbers(SymmetryElectronsSzType& q, const BlockType& block) const
	{
		VectorSizeType qns;
		HilbertBasisType basis;
		setNaturalBasis(basis, qns, block);
		setSymmetryRelated(q, basis, block.size());
	}

	//! set operator matrices for sites in block
	void setOperatorMatrices(VectorOperatorType& operatorMatrices,
	                         const BlockType& block) const
	{
		HilbertBasisType natBasis;
		SparseMatrixType tmpMatrix;

		VectorSizeType qvector;
		setNaturalBasis(natBasis,qvector,block);

		typename MatrixType::value_type dummy = 0.0;

		operatorMatrices.clear();
		for (SizeType i=0;i<block.size();i++) {
			// Set the operators S^x_i in the natural basis
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_X, dummy);

			typename OperatorType::Su2RelatedType su2related;

			OperatorType myOp(tmpMatrix, 1, PairType(0, 0), 1.0, su2related);
			operatorMatrices.push_back(myOp);

			// Set the operators S^y_i in the natural basis
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_Y, dummy);
			OperatorType myOp2(tmpMatrix, 1, PairType(0, 0), 1.0, su2related);
			operatorMatrices.push_back(myOp2);

			// Set the operators S^z_i in the natural basis
			tmpMatrix = findSdirMatrices(i, natBasis, DIR_Z, dummy);
			OperatorType myOp3(tmpMatrix, 1, PairType(0, 0), 1.0, su2related);
			operatorMatrices.push_back(myOp3);
		}
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatrices(creationMatrix,block);
		assert(creationMatrix.size()>0);

		if (what == "sx") // S^x
			return creationMatrix[0];

		if (what == "sy") // S^y
			return creationMatrix[1];

		if (what == "sz") // S^z
			return creationMatrix[2];

		PsimagLite::String str("Kitaev: naturalOperator: no label ");
		str += what + "\n";
		throw PsimagLite::RuntimeError(str);
	}

	//! find all states in the natural basis for a block of n sites
	void setNaturalBasis(HilbertBasisType& basis,
	                     VectorSizeType& q,
	                     const VectorSizeType& block) const
	{
		SizeType total = utils::powUint(TWICE_THE_SPIN + 1, block.size());

		for (SizeType i = 0; i < total; ++i) basis.push_back(i);

		SymmetryElectronsSzType qq;
		setSymmetryRelated(qq, basis, block.size());
		qq.findQuantumNumbers(q, MyBasis::useSu2Symmetry());
	}

	//! Dummy since this model has no fermion sign
	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType) const
	{
		electrons.resize(basis.size());
		std::fill(electrons.begin(), electrons.end(), 0);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                const VectorOperatorType& cm,
	                                const BlockType& block,
	                                RealType,
	                                RealType factorForDiagonals=1.0)  const
	{
		SizeType linSize = geometry_.numberOfSites();
		if (modelParameters_.magneticField.size() != linSize)
			return; // <<---- PLEASE NOTE EARLY EXIT HERE

		SizeType n=block.size();

		for (SizeType i = 0; i < n; ++i) {
			// magnetic field
			RealType tmp = modelParameters_.magneticField[block[i*2]]*factorForDiagonals;
			hmatrix += tmp*cm[1+i*2].data;
		}
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return modelParameters_.targetQuantum;
	}

private:

	SizeType logBase2(SizeType x) const
	{
		SizeType counter = 0;
		while (x > 0) {
			x >>= 1;
			counter++;
		}

		return (counter == 0) ? counter : counter - 1;
	}

	SparseMatrixType findSdirMatrices(SizeType,
	                                  const HilbertBasisType&,
	                                  InternalDir,
	                                  RealType) const
	{
		err("Kitaev needs useComplex in SolverOptions in the input file\n");
		throw PsimagLite::RuntimeError("FATAL\n");
	}

	SparseMatrixType findSdirMatrices(SizeType,// site,
	                                  const HilbertBasisType& natBasis,
	                                  InternalDir dir,
	                                  std::complex<RealType>) const
	{
		SizeType total = natBasis.size();
		MatrixType cm(total,total);
		cm.setTo(0.0);
		assert(total == TWICE_THE_SPIN + 1);
		assert(TWICE_THE_SPIN == 1);

		if (dir == DIR_X) {
			cm(0,1) = cm(1,0) = 1.0;
		} else if (dir == DIR_Y) {
			cm(0, 1) = std::complex<RealType>(0.0, 1.0);  // FIXME TODO CHECK SIGN HERE
			cm(1, 0) = std::complex<RealType>(0.0, -1.0); // FIXME TODO CHECK SIGN HERE
		} else if (dir == DIR_Z) {
			cm(0, 0) = 1.0;
			cm(1, 1) = -1.0;
		} else {
			assert(false);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	void setSymmetryRelated(SymmetryElectronsSzType& q,
	                        const HilbertBasisType& basis,
	                        int n) const
	{
		assert(n == 1);
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		typedef std::pair<SizeType,SizeType> PairType;
		SizeType nbasis = basis.size();
		typename PsimagLite::Vector<PairType>::Type jmvalues(nbasis, PairType(0, 0));
		VectorSizeType flavors(nbasis, 1);

		VectorSizeType szPlusConst(nbasis, 0);
		VectorSizeType electrons(basis.size(), 0);

		q.set(jmvalues,flavors,electrons,szPlusConst);
	}

	ParametersKitaev<RealType>  modelParameters_;
	const GeometryType& geometry_;
}; // class Kitaev

} // namespace Dmrg
/*@}*/
#endif //DMRG_KITAEV_H

