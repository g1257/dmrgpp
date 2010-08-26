// BEGIN LICENSE BLOCK
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
// END LICENSE BLOCK
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
#include "HilbertSpaceHubbard.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "IoSimple.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramGlobals.h"

namespace Dmrg {	
	
	template<typename ModelHelperType_,
	typename SparseMatrixType,
	typename GeometryType,
	template<typename> class SharedMemoryTemplate>
	class ModelHeisenberg
		: public ModelBase<ModelHelperType_,SparseMatrixType,GeometryType,
  	LinkProductHeisenbergSpinOneHalf<ModelHelperType_>,SharedMemoryTemplate> {

	public:
		typedef ModelHelperType_ ModelHelperType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename ModelHelperType::RealType RealType;

	private:
		typedef typename ModelHelperType::BlockType Block;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename ModelHelperType::ReflectionSymmetryType ReflectionSymmetryType;
		typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
		typedef unsigned int long long WordType;
		typedef HilbertSpaceHubbard<WordType> HilbertSpaceType;
		typedef typename HilbertSpaceType::HilbertState HilbertStateType;
		typedef LinkProductHeisenbergSpinOneHalf<ModelHelperType> LinkProductType;
		typedef ModelBase<ModelHelperType,SparseMatrixType,GeometryType,LinkProductType,SharedMemoryTemplate>
				ModelBaseType;
		
		

		static const int NUMBER_OF_ORBITALS=OperatorsType::NUMBER_OF_ORBITALS;
		static const int DEGREES_OF_FREEDOM=2; // spin up and down
		static int const maxNumberOfSites=ProgramGlobals::MaxNumberOfSites;

	public:
		typedef std::vector<HilbertStateType> HilbertBasisType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef	typename ModelBaseType::MyBasis MyBasis;
		typedef	typename ModelBaseType::MyBasisWithOperators MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;

		ModelHeisenberg(ParametersModelHeisenberg<RealType> const &mp,GeometryType const &geometry) 
			: ModelBaseType(geometry),modelParameters_(mp), geometry_(geometry), 
			spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM),
					reinterpretX_(maxNumberOfSites),reinterpretY_(maxNumberOfSites)
		{
		}

		void print(std::ostream& os) const { os<<modelParameters_; }

		size_t orbitals() const { return NUMBER_OF_ORBITALS; }

		//! find  operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(
				std::vector<OperatorType> &operatorMatrices,
    				SparseMatrixType &hamiltonian,
				BasisDataType &q,
    				Block const &block) const
		{
			std::vector<HilbertStateType> natBasis;
			
			setNaturalBasis(natBasis,block.size());
			
			setOperatorMatrices(operatorMatrices,block);
			
			setSymmetryRelated(q,natBasis,block.size());
			
			calcHamiltonian(hamiltonian,operatorMatrices,block);
		}

		//! set operator matrices for sites in block
		void setOperatorMatrices(std::vector<OperatorType> &operatorMatrices,Block const &block) const
		{
			std::vector<HilbertStateType> natBasis;
			SparseMatrixType tmpMatrix;

			setNaturalBasis(natBasis,block.size());

			operatorMatrices.clear();
			for (size_t i=0;i<block.size();i++) {
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
		
		psimag::Matrix<SparseElementType> getOperator(const std::string& what,size_t gamma=0,size_t spin=0) const
		{
			Block block;
			block.resize(1);
			block[0]=0;
			std::vector<OperatorType> creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="+" or what=="i") { // S^+
				psimag::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[0].data);
				return tmp;
			} else if (what=="-") { // S^-
				psimag::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[0].data);
				transposeConjugate(tmp);
				return tmp;
			} else if (what=="z") { // S^z
				psimag::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[1].data);
				return tmp;
			}
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}
		
		//! find all states in the natural basis for a block of n sites
		void setNaturalBasis(std::vector<HilbertStateType>  &basis,int n) const
		{
			if (n!=1) throw std::runtime_error("setNaturalBasis: implemented only for blocks of size=1\n");
			basis.push_back(0);
			basis.push_back(1);
			basis.push_back(2);
		}
		
		//! Dummy since this model has no fermion sign
		void findElectrons(std::vector<size_t>& electrons,
				   std::vector<HilbertStateType>& basis) const
		{
			electrons.resize(basis.size());
			for (size_t i=0;i<electrons.size();i++) electrons[i] = 0;
		}

	private:
		ParametersModelHeisenberg<RealType>  modelParameters_;
		GeometryType const &geometry_;
		SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
		SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
		size_t reinterpretX_,reinterpretY_;

		//! Find S^+_i in the natural basis natBasis
		SparseMatrixType findSplusMatrices(int i,std::vector<HilbertStateType> const &natBasis) const
		{
			HilbertStateType bra,ket;
			int n = natBasis.size();
			psimag::Matrix<SparseElementType> cm(n,n);

			for (size_t ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				if (HilbertSpaceType::get(ket,i)==2) {
					// it is a down electron, then flip it:
					HilbertSpaceType::destroy(bra,i,1);
					HilbertSpaceType::create(bra,i,0);
					int jj = utils::isInVector(natBasis,bra);
					if (jj<0) throw std::runtime_error("findOperatorMatrices: internal error while"
								"creating.\n");
					cm(ii,jj)=1.0;
				}
			}
			SparseMatrixType operatorMatrix(cm);
			return operatorMatrix;
		}

		//! Find S^z_i in the natural basis natBasis
		SparseMatrixType findSzMatrices(int i,std::vector<HilbertStateType> const &natBasis) const
		{
			HilbertStateType bra,ket;
			int n = natBasis.size();
			psimag::Matrix<SparseElementType> cm(n,n);

			for (size_t ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				switch (HilbertSpaceType::get(ket,i)) {
					case 0:
						cm(ii,ii)=0;
						break;
					case 1:
						cm(ii,ii)=0.5;
						break;
					case 2:
						cm(ii,ii)= -0.5;
						break;
					default:
						throw std::runtime_error("findOperatorMatrices: internal error while"
								"creating.\n");
				}
			}
			SparseMatrixType operatorMatrix(cm);
			return operatorMatrix;
		}

		//! Full hamiltonian from operator matrices cm
		void calcHamiltonian(SparseMatrixType &hmatrix,std::vector<OperatorType> const &cm,Block const &block) const
		{
			size_t n=block.size();
			SparseMatrixType tmpMatrix,tmpMatrix2,niup,nidown;

			hmatrix.makeDiagonal(cm[0].data.rank());
			
			//! exchange
			for (size_t i=0;i<n;i++) {
				SparseMatrixType sPlusOperatorI = cm[i].data; //S^+_i
				SparseMatrixType szOperatorI =cm[i+n].data; //S^z_i
				for (size_t j=0;j<n;j++) {
					for (size_t term=0;term<geometry_.terms();term++) {
						for (size_t dofs=0;dofs<LinkProductType::dofs(term);dofs++) {
							std::pair<size_t,size_t> edofs = LinkProductType::connectorDofs(term,dofs);
							RealType tmp = geometry_(block[i],edofs.first,block[j],edofs.second,term);
						
							if (i==j || tmp==0.0) continue;
			
							SparseMatrixType sPlusOperatorJ = cm[j].data;//S^+_j
							SparseMatrixType tJ, tI;
							transposeConjugate(tJ,sPlusOperatorJ);
							transposeConjugate(tI,sPlusOperatorI);
							hmatrix += 0.5*tmp*(sPlusOperatorI*tJ);
							hmatrix += 0.5*tmp*(tI*sPlusOperatorJ);
		
							// S^z_i S^z_j
							SparseMatrixType szOperatorJ=cm[j+n].data; //S^z_j
							hmatrix += tmp*szOperatorI*szOperatorJ;
						}
					}
				}
			}
		}

		void setSymmetryRelated(BasisDataType& q,std::vector<HilbertStateType>  const &basis,int n) const
		{
			if (n!=1) std::runtime_error("ModelFeAs::setSymmetryRelated() implemented for n=1 only\n");
			
			// find j,m and flavors (do it by hand since we assume n==1)
			// note: we use 2j instead of j
			// note: we use m+j instead of m
			// This assures us that both j and m are size_t
			typedef std::pair<size_t,size_t> PairType;
			std::vector<PairType> jmvalues;
			std::vector<size_t> flavors; 
			PairType jmSaved = calcJmvalue<PairType>(basis[0]);
			jmSaved.first++;
			jmSaved.second++;

			std::vector<size_t> electronsUp(basis.size());
			std::vector<size_t> electronsDown(basis.size());
			for (size_t i=0;i<basis.size();i++) {
				PairType jmpair = calcJmvalue<PairType>(basis[i]);
				
				jmvalues.push_back(jmpair);
				// nup
				electronsUp[i] = HilbertSpaceType::getNofDigits(basis[i],HilbertSpaceType::SPIN_UP);
				// ndown
				electronsDown[i] = HilbertSpaceType::getNofDigits(basis[i],HilbertSpaceType::SPIN_DOWN);

				flavors.push_back(electronsUp[i]+electronsDown[i]);
				jmSaved = jmpair;
			}
			q.jmValues=jmvalues;
			q.flavors = flavors;
			q.electronsUp = electronsUp;
			q.electronsDown = electronsDown;
		}
		
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are size_t
		// Reinterprets 6 and 9
		template<typename PairType>
		PairType calcJmvalue(const HilbertStateType& ket) const
		{
			PairType jm(0,0);

			size_t x=reinterpretX_,y=reinterpretY_; // these states need reinterpretation
			if (ket==x) {
				jm=std::pair<size_t,size_t>(2,1);
			} else if (ket==y) {
				jm=std::pair<size_t,size_t>(0,0);
			} else jm=calcJmValueAux<PairType>(ket);
			//std::cerr<<jm.first<<" "<<jm.second<<" |--------------\n";
			return jm; 
		}

		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are size_t
		// does not work for 6 or 9
		template<typename PairType>
		PairType calcJmValueAux(const HilbertStateType& ket) const
		{
			size_t site0=0;
			size_t site1=0;

			spinSquared_.doOnePairOfSitesA(ket,site0,site1);
			spinSquared_.doOnePairOfSitesB(ket,site0,site1);
			spinSquared_.doDiagonal(ket,site0,site1);

			RealType sz = spinSquared_.spinZ(ket,site0);
			PairType jm= spinSquaredHelper_.getJmPair(sz);
			return jm;
		}
	}; // class ModelHeisenberg

	template<typename ModelHelperType,
	typename SparseMatrixType,
	typename GeometryType,
	template<typename> class SharedMemoryTemplate>
	std::ostream &operator<<(std::ostream &os,const ModelHeisenberg<ModelHelperType,SparseMatrixType,GeometryType,SharedMemoryTemplate>& model)
	{
		model.print(os);
		return os;
	}

	template<typename SparseMatrixType>
	SparseMatrixType operator*(const SparseMatrixType& a,const SparseMatrixType& b)
	{
		SparseMatrixType temp;
		multiply(temp,a,b);
		return temp;
	}

	template<typename RealTypeType,typename SparseMatrixType>
	SparseMatrixType operator*(const RealTypeType& a,const SparseMatrixType& b)
	{
		SparseMatrixType temp;
		multiplyScalar(temp,b,static_cast<typename SparseMatrixType::value_type>(a));
		return temp;
	}
} // namespace Dmrg
/*@}*/
#endif //DMRG_MODEL_HEISENBERG_HEADER_H
