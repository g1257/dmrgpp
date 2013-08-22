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

/*! \file ModelHeisenberg.h
 *
 *  An implementation of the Quantum Heisenberg Model to use with  DmrgSolver
 *
 */
 
#ifndef DMRG_MODEL_HEISENBERG_HEADER_H
#define DMRG_MODEL_HEISENBERG_HEADER_H

#include <algorithm>
#include "ModelBase.h"
#include "ParametersModelHeisenberg.h"
#include "LinkProductHeisenbergSpinOneHalf.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramGlobals.h"
#include "ModelCommon.h"

namespace Dmrg {	
	
	template<typename ModelBaseType>
	class ModelHeisenberg : public ModelBaseType {

	public:

		typedef typename ModelBaseType::ModelHelperType ModelHelperType;
		typedef typename ModelBaseType::GeometryType GeometryType;
		typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
		typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
		typedef typename ModelBaseType::LinkType LinkType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename ModelHelperType::RealType RealType;
		typedef	typename ModelBaseType::VectorType VectorType;

	private:

		typedef typename ModelBaseType::BlockType BlockType;
		typedef typename ModelBaseType::SolverParamsType SolverParamsType;
		typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef unsigned int long long WordType;
		typedef LinkProductHeisenbergSpinOneHalf<ModelHelperType> LinkProductType;
		typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
		typedef typename ModelBaseType::InputValidatorType InputValidatorType;

		static const int NUMBER_OF_ORBITALS=1;
		static const int DEGREES_OF_FREEDOM=2; // spin up and down

	public:

		typedef typename PsimagLite::Vector<unsigned int long long>::Type HilbertBasisType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
		typedef	typename ModelBaseType::MyBasis MyBasis;
		typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;

		ModelHeisenberg(const SolverParamsType& solverParams,
		                InputValidatorType& io,
		                GeometryType const &geometry)
		    : ModelBaseType(solverParams,io,geometry,new ModelCommonType(geometry)),
		      modelParameters_(io),
		      geometry_(geometry),
		      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM)
		{}

		void print(std::ostream& os) const { os<<modelParameters_; }

		SizeType hilbertSize(SizeType site) const { return modelParameters_.twiceTheSpin+1; }

		//! find  operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(VectorOperatorType& operatorMatrices,
		                     SparseMatrixType &hamiltonian,
		                     BasisDataType &q,
		                     const BlockType& block,
		                     const RealType& time) const
		{
			HilbertBasisType natBasis;
			
			typename PsimagLite::Vector<SizeType>::Type qvector;
			setNaturalBasis(natBasis,qvector,block);
			
			setOperatorMatrices(operatorMatrices,block);
			
			setSymmetryRelated(q,natBasis,block.size());
			
			this->calcHamiltonian(hamiltonian,operatorMatrices,block,time);
		}

		//! set operator matrices for sites in block
		void setOperatorMatrices(VectorOperatorType& operatorMatrices,
		                         const BlockType& block) const
		{
			HilbertBasisType natBasis;
			SparseMatrixType tmpMatrix;
			
			typename PsimagLite::Vector<SizeType>::Type qvector;
			setNaturalBasis(natBasis,qvector,block);

			operatorMatrices.clear();
			for (SizeType i=0;i<block.size();i++) {
				// Set the operators S^+_i in the natural basis
				tmpMatrix=findSplusMatrices(i,natBasis);

				typename OperatorType::Su2RelatedType su2related;
				su2related.source.push_back(i*DEGREES_OF_FREEDOM);
				su2related.source.push_back(i*DEGREES_OF_FREEDOM+NUMBER_OF_ORBITALS);	
				su2related.source.push_back(i*DEGREES_OF_FREEDOM);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(-1);
				su2related.transpose.push_back(1);
				su2related.offset = NUMBER_OF_ORBITALS;

				OperatorType myOp(tmpMatrix,1,typename OperatorType::PairType(2,2),-1,su2related);
				operatorMatrices.push_back(myOp);

				// Set the operators S^z_i in the natural basis
				tmpMatrix = findSzMatrices(i,natBasis);
				typename OperatorType::Su2RelatedType su2related2;
				OperatorType myOp2(tmpMatrix,1,typename OperatorType::PairType(2,1),1.0/sqrt(2.0),su2related2);
				operatorMatrices.push_back(myOp2);
			}
		}

		PsimagLite::Matrix<SparseElementType> naturalOperator(const PsimagLite::String& what,SizeType site,SizeType dof) const
		{
			BlockType block;
			block.resize(1);
			block[0]=site;
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="+" or what=="i") { // S^+
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[0].data);
				return tmp;
			} else if (what=="-") { // S^-
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[0].data);
				transposeConjugate(tmp);
				return tmp;
			} else if (what=="z") { // S^z
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[1].data);
				return tmp;
			}
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}
		
		//! find all states in the natural basis for a block of n sites
		void setNaturalBasis(HilbertBasisType& basis,
		                     typename PsimagLite::Vector<SizeType>::Type& q,
		                     const typename PsimagLite::Vector<SizeType>::Type& block) const
		{
			assert(block.size()==1);
			SizeType total = modelParameters_.twiceTheSpin + 1;
			for (SizeType i=0;i<total;i++) basis.push_back(i);
			BasisDataType qq;
			setSymmetryRelated(qq,basis,block.size());
			MyBasis::findQuantumNumbers(q,qq);
		}
		
		//! Dummy since this model has no fermion sign
		void findElectrons(typename PsimagLite::Vector<SizeType>::Type& electrons,
				   const HilbertBasisType& basis,
		                   SizeType site) const
		{
			electrons.resize(basis.size());
			for (SizeType i=0;i<electrons.size();i++)
				electrons[i] = 0;
		}

		virtual void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
		                                        const VectorOperatorType& cm,
		                                        const BlockType& block,
		                                        RealType time,
		                                        RealType factorForDiagonals=1.0)  const
		{}

	private:

		ParametersModelHeisenberg<RealType>  modelParameters_;
		GeometryType const &geometry_;
		SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
		SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;

		//! Find S^+_i in the natural basis natBasis
		SparseMatrixType findSplusMatrices(int i,const HilbertBasisType& natBasis) const
		{
			SizeType total = natBasis.size();
			PsimagLite::Matrix<SparseElementType> cm(total,total);
			RealType j = 0.5*modelParameters_.twiceTheSpin;

			for (SizeType ii=0;ii<total;ii++) {
				SizeType ket = natBasis[ii];
				SizeType bra = ket + 1;
				if (bra>=total) continue;
				RealType m = ket - j;
				RealType x = j*(j+1)-m*(m+1);
				assert(x>=0);
				cm(ket,bra) = sqrt(x);
			}

			SparseMatrixType operatorMatrix(cm);
			return operatorMatrix;
		}

		//! Find S^z_i in the natural basis natBasis
		SparseMatrixType findSzMatrices(int i,const HilbertBasisType& natBasis) const
		{
			SizeType total = natBasis.size();
			PsimagLite::Matrix<SparseElementType> cm(total,total);
			RealType j = 0.5*modelParameters_.twiceTheSpin;

			for (SizeType ii=0;ii<total;ii++) {
				SizeType ket = natBasis[ii];
				RealType m = ket - j;
				cm(ket,ket) = m;
			}

			SparseMatrixType operatorMatrix(cm);
			return operatorMatrix;
		}

		void setSymmetryRelated(BasisDataType& q,const HilbertBasisType& basis,int n) const
		{
			if (n!=1) PsimagLite::RuntimeError("ModelFeAs::setSymmetryRelated() implemented for n=1 only\n");
			
			// find j,m and flavors (do it by hand since we assume n==1)
			// note: we use 2j instead of j
			// note: we use m+j instead of m
			// This assures us that both j and m are SizeType
			typedef std::pair<SizeType,SizeType> PairType;
			typename PsimagLite::Vector<PairType>::Type jmvalues;
			typename PsimagLite::Vector<SizeType>::Type flavors; 
			PairType jmSaved;
			jmSaved.first = modelParameters_.twiceTheSpin;
			jmSaved.second = basis[0];
			jmSaved.first++;
			jmSaved.second++;

			typename PsimagLite::Vector<SizeType>::Type electronsUp(basis.size());
			typename PsimagLite::Vector<SizeType>::Type electronsDown(basis.size());
			for (SizeType i=0;i<basis.size();i++) {
				PairType jmpair;
				jmpair.first = modelParameters_.twiceTheSpin;
				jmpair.second = basis[i];
				SizeType ket = basis[i];
				jmvalues.push_back(jmpair);

				// nup
				electronsUp[i] = (ket==0) ?  0 : 1;
				// ndown
				electronsDown[i] = (ket==1) ?  0 : 1;

				flavors.push_back(electronsUp[i]+electronsDown[i]);
				jmSaved = jmpair;
			}
			q.jmValues=jmvalues;
			q.flavors = flavors;
			q.electrons = electronsUp + electronsDown;
			q.szPlusConst.resize(electronsUp.size());
			for (SizeType i=0;i<q.szPlusConst.size();i++) {
				q.szPlusConst[i] = (modelParameters_.twiceTheSpin + electronsUp[i]) - electronsDown[i];
				assert(!(q.szPlusConst[i] & 1));
				q.szPlusConst[i] = static_cast<SizeType>(q.szPlusConst[i]*0.5);
			}
		}
	}; // class ModelHeisenberg

} // namespace Dmrg
/*@}*/
#endif //DMRG_MODEL_HEISENBERG_HEADER_H
