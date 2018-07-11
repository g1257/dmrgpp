/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file SuperExtendedHubbard1Orb.h
 *
 *  Hubbard + V^n_{ij} n_i n_j + V^S_{ij} S_i S_j + V^P_{ij} P^dag_i P_j
 *
 */
#ifndef EXTENDED_SUPER_HUBBARD_1ORB_H
#define EXTENDED_SUPER_HUBBARD_1ORB_H
#include "../Models/ExtendedHubbard1Orb/ExtendedHubbard1Orb.h"
#include "LinkProdExtendedSuperHubbard1Orb.h"

namespace Dmrg {
//! Super Extended Hubbard for DMRG solver, uses ModelHubbard by containment
template<typename ModelBaseType>
class ExtendedSuperHubbard1Orb : public ModelBaseType {

public:

	typedef typename ModelBaseType::VectorSizeType VectorSizeType;
	typedef ExtendedHubbard1Orb<ModelBaseType> ExtendedHubbard1OrbType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
	typedef typename ModelBaseType::LinkType LinkType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelBaseType::TargetQuantumElectronsType TargetQuantumElectronsType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef LinkProdExtendedSuperHubbard1Orb<ModelHelperType, GeometryType> LinkProductType;
	typedef	typename ModelBaseType::MyBasis MyBasis;
	typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
	typedef typename MyBasis::SymmetryElectronsSzType SymmetryElectronsSzType;
	typedef typename ExtendedHubbard1OrbType::HilbertBasisType HilbertBasisType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename ExtendedHubbard1OrbType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
	typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelBaseType::VectorOperatorType VectorOperatorType;
	typedef typename PsimagLite::Vector<HilbertState>::Type VectorHilbertStateType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;

	ExtendedSuperHubbard1Orb(const SolverParamsType& solverParams,
	                         InputValidatorType& io,
	                         GeometryType const &geometry)
	    : ModelBaseType(solverParams,
	                    geometry,
	                    new LinkProductType((geometry_.orbitals(LinkProductType::TERM_HOPPING,
	                                                           0) == 2))),
	      modelParameters_(io),
	      geometry_(geometry),
	      extendedHubbard_(solverParams,io,geometry)
	{}

	SizeType memResolv(PsimagLite::MemResolv&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType hilbertSize(SizeType) const
	{
		return (1<<2);
	}

	void setOperatorMatrices(VectorOperatorType& creationMatrix,
	                         const BlockType& block) const
	{
		extendedHubbard_.setOperatorMatrices(creationMatrix,block);
		addAditionalOperatorMatrices(creationMatrix,block);
		assert(creationMatrix.size() == 6);
	}

	OperatorType naturalOperator(const PsimagLite::String& what,
	                             SizeType site,
	                             SizeType dof) const
	{
		BlockType block;
		block.resize(1);
		block[0]=site;
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatrices(creationMatrix,block);

		if (what=="n") {
			VectorSizeType allowed(1,0);
			ModelBaseType::checkNaturalOperatorDof(dof,what,allowed);
			return creationMatrix[2];
		}

		return extendedHubbard_.naturalOperator(what,site,dof);
	}

	//! find total number of electrons for each state in the basis
	void findElectrons(typename PsimagLite::Vector<SizeType> ::Type&electrons,
	                   const VectorHilbertStateType& basis,
	                   SizeType site) const
	{
		extendedHubbard_.findElectrons(electrons,basis,site);
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		if (!io.doesGroupExist(label1))
		        io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		extendedHubbard_.write(label, io);
	}

	virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
	                                        const VectorOperatorType& cm,
	                                        const BlockType& block,
	                                        RealType time,
	                                        RealType factorForDiagonals=1.0)  const
	{
		// remove the other matrices before sending to Hubbard
		VectorOperatorType cmCorrected;
		SizeType k = 0;
		for (SizeType i = 0; i < block.size(); ++i) {
			for (SizeType j = 0; j < 2; ++j) {
				cmCorrected.push_back(cm[k++]);
			}

			cmCorrected.push_back(cm[k++]);
			k += 3;
		}

		extendedHubbard_.addDiagonalsInNaturalBasis(hmatrix,
		                                            cmCorrected,
		                                            block,
		                                            time,
		                                            factorForDiagonals);
	}

	virtual const TargetQuantumElectronsType& targetQuantum() const
	{
		return extendedHubbard_.targetQuantum();
	}

	void setBasis(HilbertBasisType& basis,
	              SymmetryElectronsSzType& qq,
	              const VectorSizeType& block) const
	{
		extendedHubbard_.setBasis(basis, qq, block);
	}

private:

	void addAditionalOperatorMatrices(VectorOperatorType& creationMatrix,
	                                  const BlockType& block) const
	{
		HilbertBasisType natBasis;
		SymmetryElectronsSzType qq;
		extendedHubbard_.setBasis(natBasis, qq, block);

		for (SizeType i=0;i<block.size();i++) {
			setSplusi(creationMatrix,i);
			setSzi(creationMatrix,i);
			setPairi(creationMatrix,i);
		}

		assert(creationMatrix.size() == 6);
	}

	void setSplusi(VectorOperatorType& cm,
	               SizeType) const
	{
		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType m;
		transposeConjugate(m,cm[1].data);
		SparseMatrixType splus;
		multiply(splus,cm[0].data,m);

		OperatorType myOp(splus,
		                  1,
		                  typename OperatorType::PairType(0,0),
		                  -1,
		                  su2related);
		cm.push_back(myOp);
	}

	void setSzi(VectorOperatorType& cm,
	            SizeType) const
	{
		typename OperatorType::Su2RelatedType su2related;

		SparseMatrixType n1 = n(cm[0].data);
		SparseMatrixType n2 = n(cm[1].data);
		MatrixType dn1;
		MatrixType dn2;
		crsMatrixToFullMatrix(dn1,n1);
		crsMatrixToFullMatrix(dn2,n2);

		SizeType n = dn1.rows();
		MatrixType szmatrix(n,n);


		for (SizeType i = 0; i < n; ++i)
			szmatrix(i,i) = static_cast<RealType>(0.5)*(dn1(i,i) - dn2(i,i));

		OperatorType sz(SparseMatrixType(szmatrix),
		                1,
		                typename OperatorType::PairType(0,0),
		                1,
		                su2related);

		cm.push_back(sz);
	}

	void setPairi(VectorOperatorType& cm,
	              SizeType) const
	{
		typename OperatorType::Su2RelatedType su2related;
		SparseMatrixType pair;
		multiply(pair,cm[0].data,cm[1].data);

		OperatorType myOp(pair,
		                  1,
		                  typename OperatorType::PairType(0,0),
		                  1,
		                  su2related);
		cm.push_back(myOp);
	}

	SparseMatrixType n(const SparseMatrixType& c) const
	{
		SparseMatrixType tmpMatrix;
		SparseMatrixType cdagger;
		transposeConjugate(cdagger,c);
		multiply(tmpMatrix,c,cdagger);

		return tmpMatrix;
	}

	//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
	RealType sign(typename HilbertSpaceHubbardType::HilbertState const &ket,
	              int i,
	              int sigma) const
	{
		int value=0;
		value += HilbertSpaceHubbardType::calcNofElectrons(ket,0,i,0);
		value += HilbertSpaceHubbardType::calcNofElectrons(ket,0,i,1);
		int tmp1 = HilbertSpaceHubbardType::get(ket,0) &1;
		int tmp2 = HilbertSpaceHubbardType::get(ket,0) &2;
		if (i>0 && tmp1>0) value++;
		if (i>0 && tmp2>0) value++;

		if (sigma==1) { // spin down
			if ((HilbertSpaceHubbardType::get(ket,i) &1)) value++;

		}
		if (value%2==0) return 1.0;

		return -1;
	}


	//serializr start class ExtendedSuperHubbard1Orb
	//serializr vptr
	//serializr normal modelParameters_
	ParametersModelHubbard<RealType>  modelParameters_;
	//serializr ref geometry_ start
	const GeometryType &geometry_;
	//serializr normal modelExtHubbard_
	ExtendedHubbard1OrbType extendedHubbard_;
};	//class ExtendedSuperHubbard1Orb

} // namespace Dmrg
/*@}*/
#endif // EXTENDED_SUPER_HUBBARD_1ORB_H

