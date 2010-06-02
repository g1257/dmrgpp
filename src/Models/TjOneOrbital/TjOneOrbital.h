// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

/*! \file TjOneOrbital.h
 *
 *  An implementation of the Quantum Tj one-orbital Model to use with  DmrgSolver
 *
 */
 
#ifndef TJ_ONEORBITAL_HEADER_H
#define TJ_ONEORBITAL_HEADER_H

#include <algorithm>
#include "ModelBase.h"
#include "ParametersTjOneOrbital.h"
#include "LinkProductTjOneOrbital.h"
#include "HilbertSpaceHubbard.h"
#include "CrsMatrix.h"
#include "VerySparseMatrix.h"
#include "IoSimple.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "ProgramLimits.h"

namespace Dmrg {	
	
	template<typename ModelHelperType_,
	typename SparseMatrixType,
	typename GeometryType,
	template<typename> class SharedMemoryTemplate>
	class TjOneOrbital
		: public ModelBase<ModelHelperType_,SparseMatrixType,GeometryType,
  	LinkProductTjOneOrbital<typename SparseMatrixType::MatrixElementType,ModelHelperType_>,SharedMemoryTemplate> {
	public:
		typedef ModelHelperType_ ModelHelperType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;

	private:
		typedef typename ModelHelperType::Block Block;
		typedef typename ModelHelperType::Field Field;
		typedef typename ModelHelperType::MatrixElementType MatrixElementType;
		typedef typename ModelHelperType::ReflectionSymmetryType ReflectionSymmetryType;
		typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
		typedef unsigned int long long WordType;
		typedef HilbertSpaceHubbard<WordType> HilbertSpaceType;
		typedef typename HilbertSpaceType::HilbertState HilbertStateType;
		typedef LinkProductTjOneOrbital<MatrixElementType,ModelHelperType> LinkProductType;
		typedef ModelBase<ModelHelperType,SparseMatrixType,GeometryType,LinkProductType,SharedMemoryTemplate>
				ModelBaseType;
		typedef	typename ModelBaseType::MyBasis MyBasis;
		typedef	typename ModelBaseType::MyBasisWithOperators MyBasisWithOperators;
		typedef typename OperatorsType::OperatorType OperatorType;

		static const int NUMBER_OF_ORBITALS=OperatorsType::NUMBER_OF_ORBITALS;
		static const int DEGREES_OF_FREEDOM=2; // spin up and down
		static int const maxNumberOfSites=ProgramLimits::MaxNumberOfSites;
		static const int FERMION_SIGN = -1;
		
	public:
		typedef SharedMemoryTemplate<LinkProductType> SharedMemoryType;
		typedef typename MyBasis::BasisDataType BasisDataType;

		TjOneOrbital(ParametersTjOneOrbital<Field> const &mp,GeometryType const &geometry) 
			: ModelBaseType(DEGREES_OF_FREEDOM,geometry),modelParameters_(mp), geometry_(geometry), 
			spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM),
					reinterpretX_(maxNumberOfSites),reinterpretY_(maxNumberOfSites)
		{
		}

		void print(std::ostream& os) const { os<<modelParameters_; }

		//! returns internal degrees of freedom for this model
		int dof() const { return DEGREES_OF_FREEDOM; } 

		size_t orbitals() const { return NUMBER_OF_ORBITALS; }

		int absoluteSystemSize() const {return modelParameters_.linSize; }

		//! site of the hamiltonian matrix
		int getSize(ModelHelperType const &modelHelper) const 
		{
			return modelHelper.size();
		}

		double density() const { return modelParameters_.density; }

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
				su2related.source.push_back(i*dof());
				su2related.source.push_back(i*dof()+NUMBER_OF_ORBITALS);	
				su2related.source.push_back(i*dof());
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
				
				// SEt c^\dagger_{i\sigma}
				size_t tdof = 2;
				for (int sigma=0;sigma<dof();sigma++) {
					tmpMatrix = findCreationMatrices(i,sigma,natBasis);
					int asign= 1;
					if (sigma>0) asign= 1;
					typename OperatorType::Su2RelatedType su2related;
					if (sigma==0) {
						su2related.source.push_back(i*dof()+tdof);
						su2related.source.push_back(i*dof()+1+tdof);	
						su2related.transpose.push_back(-1);
						su2related.transpose.push_back(-1);
						su2related.offset = NUMBER_OF_ORBITALS;
					}
					OperatorType myOp(tmpMatrix,-1,typename OperatorType::PairType(1,1-sigma),asign,su2related);
					
					operatorMatrices.push_back(myOp);
				}
				
				// FIXME: set n_i matrices
			}
		}

		//! set block of sites S X Y and E (see paper for geometry description)
		void setBlocksOfSites(Block &S,std::vector<Block> &X,std::vector<Block> &Y,Block &E) const 
		{
			geometry_.setBlocksOfSites(S,X,Y,E);
		}
		
		psimag::Matrix<MatrixElementType> getOperator(const std::string& what,size_t gamma=0,size_t spin=0)
		{
			Block block;
			block.resize(1);
			block[0]=0;
			std::vector<OperatorType> creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="+") { //S^+
				psimag::Matrix<MatrixElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[0].data);
				return tmp;
			}
			if (what=="-") { //S^-
				psimag::Matrix<MatrixElementType> tmp;
				size_t rank = creationMatrix[0].data.rank();
				for (size_t i=0;i<rank;i++) for (size_t j=0;j<rank;j++)
						tmp(j,i) = (creationMatrix[0].data)(i,j);
				//transposeConjugate(tmp,);
				return tmp;
			}
			if (what=="z") { //S^z
				psimag::Matrix<MatrixElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[1].data);
				return tmp;
			}
			if (what=="n") { // n
				psimag::Matrix<MatrixElementType> tmp =  multiplyTc(creationMatrix[2].data,creationMatrix[2].data)
						+ multiplyTc(creationMatrix[3].data,creationMatrix[3].data);
				return tmp;
			}
			if (what=="c") { //c^\dagger
				psimag::Matrix<MatrixElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[spin].data);
				return tmp;
			}
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}
	private:
		ParametersTjOneOrbital<Field>  modelParameters_;
		GeometryType const &geometry_;
		SpinSquaredHelper<Field,WordType> spinSquaredHelper_;
		SpinSquared<SpinSquaredHelper<Field,WordType> > spinSquared_;
		size_t reinterpretX_,reinterpretY_;
		
		//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
		Field sign(typename HilbertSpaceType::HilbertState const &ket, int i,int sigma) const
		{
			int value=0;
			value += HilbertSpaceType::calcNofElectrons(ket,0,i,0);
			value += HilbertSpaceType::calcNofElectrons(ket,0,i,1);
			int tmp1 = HilbertSpaceType::get(ket,0) &1;
			int tmp2 = HilbertSpaceType::get(ket,0) &2;
			if (i>0 && tmp1>0) value++;
			if (i>0 && tmp2>0) value++;

			if (sigma==1) { // spin down
				if ((HilbertSpaceType::get(ket,i) &1)) value++;
				
			}
			if (value==0 || value%2==0) return 1.0;

			return FERMION_SIGN;
		}

		//! find all states in the natural basis for a block of n sites
		void setNaturalBasis(std::vector<HilbertStateType>  &basis,int n) const
		{
			if (n!=1) throw std::runtime_error("setNaturalBasis: implemented only for blocks of size=1\n");
			basis.push_back(0);
			basis.push_back(1);
			basis.push_back(2);
		}

		//! Find S^+_i in the natural basis natBasis
		SparseMatrixType findSplusMatrices(int i,std::vector<HilbertStateType> const &natBasis) const
		{
			HilbertStateType bra,ket;
			int n = natBasis.size();
			psimag::Matrix<MatrixElementType> cm(n,n);

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
			CrsMatrix<MatrixElementType> operatorMatrix(cm);
			return operatorMatrix;
		}

		//! Find S^z_i in the natural basis natBasis
		SparseMatrixType findSzMatrices(int i,std::vector<HilbertStateType> const &natBasis) const
		{
			HilbertStateType bra,ket;
			int n = natBasis.size();
			psimag::Matrix<MatrixElementType> cm(n,n);

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
			CrsMatrix<MatrixElementType> operatorMatrix(cm);
			return operatorMatrix;
		}
		
		//! Find c^\dagger_isigma in the natural basis natBasis
		SparseMatrixType findCreationMatrices(int i,int sigma,std::vector<typename HilbertSpaceType::HilbertState> const &natBasis) const
		{
			typename HilbertSpaceType::HilbertState bra,ket;
			int n = natBasis.size();
			psimag::Matrix<MatrixElementType> cm(n,n);
			
			for (size_t ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				if (HilbertSpaceType::isNonZero(ket,i,sigma)) {
					
				} else {
					HilbertSpaceType::create(bra,i,sigma);
					int jj = utils::isInVector(natBasis,bra);
					if (jj<0) continue; // double occupied states are not allowed in this model
					cm(ii,jj) =sign(ket,i,sigma);
				}
			}

			CrsMatrix<MatrixElementType> creationMatrix(cm);
			return 	creationMatrix;
		}
		
		//! Full hamiltonian from operator matrices cm
		void calcHamiltonian(SparseMatrixType &hmatrix,std::vector<OperatorType> const &cm,Block const &block) const
		{
			int i,j,n=block.size();
			int type;
			SparseMatrixType tmpMatrix,tmpMatrix2,niup,nidown;
			int smax,emin;
			MatrixElementType tmp;
			//size_t of = 2;
			geometry_.findExtremes(smax,emin,block);
			hmatrix.makeDiagonal(cm[0].data.rank());
			size_t nOfOperatorsPerSite = OperatorsType::numberOfOperatorsPerSite();
			
			for (i=0;i<n;i++) {
				SparseMatrixType sPlusOperatorI = cm[i*nOfOperatorsPerSite].data; //S^+_i
				SparseMatrixType szOperatorI =cm[i*nOfOperatorsPerSite+1].data; //S^z_i
				for (j=0;j<n;j++) {
					
					// 0.5*(S^+_i S^-_j + S^-_i S^+_j)
					type = geometry_.calcConnectorType(block[i],block[j]);
					tmp = geometry_.calcConnectorValue(type,block[i],0,block[j],0,smax,emin,LinkProductType::JVALUES);
					
					if (tmp==0.0) continue;
					
					SparseMatrixType sPlusOperatorJ = cm[j*nOfOperatorsPerSite].data;//S^+_j
					SparseMatrixType tJ, tI;
					transposeConjugate(tJ,sPlusOperatorJ);
					transposeConjugate(tI,sPlusOperatorI);
					hmatrix += 0.5*(sPlusOperatorI*tJ);
					hmatrix += 0.5*(tI*sPlusOperatorJ);

					// S^z_i S^z_j
					SparseMatrixType szOperatorJ=cm[j*nOfOperatorsPerSite+1].data; //S^z_j
					hmatrix += szOperatorI*szOperatorJ;
				}

				//! hopping part
// 				for (j=0;j<n;j++) {
// 					type = geometry_.calcConnectorType(block[i],block[j]);
// 					tmp = geometry_.calcConnectorValue(type,block[i],0,block[j],0,smax,emin,HOPPINGS);
// 					
// 					if (i==j || tmp==0.0) continue;
// 			
// 					for (int sigma=0;sigma<dof();sigma++) {
// 						transposeConjugate(tmpMatrix2,cm[of+sigma+j*dof()].data);
// 						multiply(tmpMatrix,cm[of+sigma+i*dof()].data,tmpMatrix2);
// 						multiplyScalar(tmpMatrix2,tmpMatrix,tmp);
// 						hmatrix += tmpMatrix2;
// 					}
// 				}
				
				// FIXME ADD n_i n_j term
			} // loop over i
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

			Field sz = spinSquared_.spinZ(ket,site0);
			PairType jm= spinSquaredHelper_.getJmPair(sz);
			return jm;
		}
	}; // class TjOneOrbital

	template<typename ModelHelperType,
	typename SparseMatrixType,
	typename GeometryType,
	template<typename> class SharedMemoryTemplate>
	std::ostream &operator<<(std::ostream &os,const TjOneOrbital<ModelHelperType,SparseMatrixType,GeometryType,SharedMemoryTemplate>& model)
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

	template<typename FieldType,typename SparseMatrixType>
	SparseMatrixType operator*(const FieldType& a,const SparseMatrixType& b)
	{
		SparseMatrixType temp;
		multiplyScalar(temp,b,a);
		return temp;
	}
} // namespace Dmrg
/*@}*/
#endif //TJ_ONEORBITAL_HEADER_H
