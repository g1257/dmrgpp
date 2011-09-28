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

/*! \file ModelFeBasedSc.h
 *
 *  An implementation of a Hubbard model for Fe-based superconductors to use with the DmrgSolver
 *
 */
#ifndef MODEL_FEAS_DMRG
#define MODEL_FEAS_DMRG
#include "ModelBase.h"
#include "ParametersModelFeAs.h"
#include "HilbertSpaceFeAs.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "LinkProductFeAs.h"
#include "ProgramGlobals.h"

namespace Dmrg {
	template<
		typename ModelHelperType_,
		typename SparseMatrixType,
		typename GeometryType,
  		template<typename> class SharedMemoryTemplate>
	class ModelFeBasedSc : public ModelBase<ModelHelperType_,SparseMatrixType,GeometryType,
 	LinkProductFeAs<ModelHelperType_>,SharedMemoryTemplate> {
		
		typedef unsigned int long long WordType;
		typedef  HilbertSpaceFeAs<WordType> HilbertSpaceFeAsType;

	public:	
		typedef ModelHelperType_ ModelHelperType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename HilbertSpaceFeAsType::HilbertState HilbertState;
		
		typedef typename ModelHelperType::BlockType Block;
		typedef typename ModelHelperType::ReflectionSymmetryType ReflectionSymmetryType;
		
		static int const maxNumberOfSites=ProgramGlobals::MaxNumberOfSites;;
		static const int FERMION_SIGN = -1;
		static const int NUMBER_OF_ORBITALS=OperatorsType::NUMBER_OF_ORBITALS;
		static const size_t DEGREES_OF_FREEDOM=OperatorsType::DEGREES_OF_FREEDOM;
		static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
		static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

	public:
		typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
		typedef std::vector<HilbertState> HilbertBasisType;
		typedef LinkProductFeAs<ModelHelperType> LinkProductType;
		typedef   ModelBase<ModelHelperType,SparseMatrixType,GeometryType,LinkProductType,SharedMemoryTemplate> ModelBaseType;
		typedef	 typename ModelBaseType::MyBasis MyBasis;
		typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;

		static size_t const REINTERPRET  = 1;
		ModelFeBasedSc(ParametersModelFeAs<RealType> const &mp,GeometryType const &geometry) 
			: ModelBaseType(geometry),reinterpretX_(6),reinterpretY_(9),modelParameters_(mp), geometry_(geometry),
					   spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM)
		{
			setPauliMatrix();
		}

		size_t orbitals() const { return NUMBER_OF_ORBITALS; }

		size_t hilbertSize() const { return (size_t)pow(2,NUMBER_OF_ORBITALS*2); } 

		void print(std::ostream& os) const { operator<<(os,modelParameters_); }

		//! find creation operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(std::vector<OperatorType> &creationMatrix,SparseMatrixType &hamiltonian,
				BasisDataType &q,Block const &block)  const
		{
			std::vector<HilbertState> natBasis;
			std::vector<size_t> qvector;
			setNaturalBasis(natBasis,qvector,block.size());			

			setOperatorMatrices(creationMatrix,block);

			//! Set symmetry related
			setSymmetryRelated(q,natBasis,block.size());

			//! set hamiltonian
			calcHamiltonian(hamiltonian,creationMatrix,block);

			SparseMatrixType tmpMatrix2;
			tmpMatrix2.makeDiagonal(natBasis.size(),0.0);
		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(std::vector<OperatorType> &creationMatrix,Block const &block) const
		{
			std::vector<HilbertState> natBasis;
			SparseMatrixType tmpMatrix;
			std::vector<size_t> qvector;
			setNaturalBasis(natBasis,qvector,block.size());

			//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
			creationMatrix.clear();
			for (size_t i=0;i<block.size();i++) {
				for (size_t sigma=0;sigma<DEGREES_OF_FREEDOM;sigma++) {
					findOperatorMatrices(tmpMatrix,i,sigma,natBasis);
					size_t m=0;
					int asign=1;
					if (sigma>1) {
						m=1;
						asign= -1;
					}
					typename OperatorType::Su2RelatedType su2related;
					if (sigma <2) {
						su2related.source.push_back(i*DEGREES_OF_FREEDOM+sigma);
						su2related.source.push_back(i*DEGREES_OF_FREEDOM+sigma + NUMBER_OF_ORBITALS);	
						su2related.transpose.push_back(-1);
						su2related.transpose.push_back(-1);
						su2related.offset = NUMBER_OF_ORBITALS;
					}	
					OperatorType myOp(tmpMatrix,-1,typename OperatorType::PairType(1,m),asign,su2related);
					creationMatrix.push_back(myOp);
				}
			}
		}

		PsimagLite::Matrix<SparseElementType> getOperator(const std::string& what,size_t orbital=0,size_t spin=0) const
		{
			Block block;
			block.resize(1);
			block[0]=0;
			std::vector<OperatorType> creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="+" or what=="i") {
				PsimagLite::Matrix<SparseElementType> tmp = multiplyTc(creationMatrix[0].data,creationMatrix[2].data);
				PsimagLite::Matrix<SparseElementType> tmp2 = multiplyTc(creationMatrix[1].data,creationMatrix[3].data);
				return tmp+tmp2;
			}
			if (what=="-") {
				PsimagLite::Matrix<SparseElementType> tmp = multiplyTc(creationMatrix[2].data,creationMatrix[0].data);
				PsimagLite::Matrix<SparseElementType> tmp2 = multiplyTc(creationMatrix[3].data,creationMatrix[1].data);
				return tmp+tmp2;
			}
			if (what=="z") {
				PsimagLite::Matrix<SparseElementType> tmp =multiplyTc(creationMatrix[0].data,creationMatrix[0].data)
						+ multiplyTc(creationMatrix[1].data,creationMatrix[1].data);
				PsimagLite::Matrix<SparseElementType> tmp2 =multiplyTc(creationMatrix[2].data,creationMatrix[2].data)
						+ multiplyTc(creationMatrix[3].data,creationMatrix[3].data);
				return tmp-tmp2;
			}
			if (what=="n") {
				PsimagLite::Matrix<SparseElementType> tmp =  multiplyTc(creationMatrix[0].data,creationMatrix[0].data)
						+ multiplyTc(creationMatrix[1].data,creationMatrix[1].data)
				 		+ multiplyTc(creationMatrix[2].data,creationMatrix[2].data)
						+ multiplyTc(creationMatrix[3].data,creationMatrix[3].data);
				return tmp;
			} 
			if (what=="c") {
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[orbital + spin*NUMBER_OF_ORBITALS].data);
				return tmp;
			}
			if (what=="d") { // delta = c^\dagger * c^dagger
				SparseMatrixType atmp;
				multiply(atmp,creationMatrix[orbital+orbital+NUMBER_OF_ORBITALS].data,creationMatrix[orbital].data);
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,atmp);
				return tmp;
			}
			std::cerr<<"what="<<what<<"\n";
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}
		
		//! find all states in the natural basis for a block of n sites
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(std::vector<HilbertState>  &basis,
		                     std::vector<size_t>& q,
		                     int n) const
		{
			HilbertState a=0;
			int sitesTimesDof=n*DEGREES_OF_FREEDOM;
			HilbertState total = (1<<sitesTimesDof);

			std::vector<HilbertState>  basisTmp;
			for (a=0;a<total;a++) basisTmp.push_back(a);

			// reorder the natural basis (needed for MULTIPLE BANDS)
			findQuantumNumbers(q,basisTmp,n);
			std::vector<size_t> iperm(q.size());
			Sort<std::vector<size_t> > sort;
			sort.sort(q,iperm);
			basis.clear();
			for (a=0;a<total;a++) basis.push_back(basisTmp[iperm[a]]);
		}
		
		void findElectrons(std::vector<size_t>& electrons,const std::vector<HilbertState>  &basis) const
		{
			electrons.resize(basis.size());
			for (size_t i=0;i<basis.size();i++) {
				// nup
				size_t nup = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_UP);
				// ndown
				size_t ndown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_DOWN);
				electrons[i] = nup + ndown;
			}
		}
		

	private:
		HilbertState reinterpretX_,reinterpretY_;
		const ParametersModelFeAs<RealType>&  modelParameters_;
		GeometryType const &geometry_;
		SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
		SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
		std::vector<PsimagLite::Matrix<RealType> > pauliMatrix_;

		//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS 
		RealType sign(HilbertState const &ket, int i,size_t sigma) const
		{
			int value=0;
			for (size_t alpha=0;alpha<DEGREES_OF_FREEDOM;alpha++) 
				value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);
			// add electron on site 0 if needed
			if (i>0) value += HilbertSpaceFeAsType::electronsAtGivenSite(ket,0);

			//order for sign is: up a (sigma==0), down a (sigma==2), up b (sigma==1), down b(sigma==3)
			unsigned int x = HilbertSpaceFeAsType::get(ket,i);	
			switch (sigma) {
				case 0:
					break;
				case 1:
					if (x & 1) value++;
					if (x & 4) value++;
					break; 
				case 2:
					if (x&1) value++;
					break;
				case 3:
					if (x&1) value++;
					if (x&4) value++;
					if (x&2) value++;
					break;
				
			}
			if (value==0 || value%2==0) return 1.0;

			return FERMION_SIGN;
		}

		//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void findOperatorMatrices(SparseMatrixType& creationMatrix,int i,int sigma,std::vector<HilbertState> const &natBasis) const
		{
			HilbertState bra,ket;
			int n = natBasis.size();
			PsimagLite::Matrix<SparseElementType> cm(n,n);

			for (size_t ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				
				if (HilbertSpaceFeAsType::isNonZero(ket,i,sigma)) {
					
				} else {
					HilbertSpaceFeAsType::create(bra,i,sigma);
					int jj = PsimagLite::isInVector(natBasis,bra);
					if (jj<0) throw std::runtime_error("findOperatorMatrices: internal error while"
							"creating.\n");
					if (ii==size_t(jj)) {
						std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra<<" sigma="<<sigma<<"\n";
						throw std::runtime_error("Creation operator cannot be diagonal\n");
					}
					cm(ii,jj) =sign(ket,i,sigma);
				}
			}
			if (REINTERPRET) reinterpret(cm,natBasis);

			SparseMatrixType temp;
			fullMatrixToCrsMatrix(temp,cm);
			transposeConjugate(creationMatrix,temp);
		}

		void findQuantumNumbers(std::vector<size_t>& q,const std::vector<HilbertState>  &basis,int n) const
		{
			BasisDataType qq;
			setSymmetryRelated(qq,basis,n);
			MyBasis::findQuantumNumbers(q,qq);
		}

		void setSymmetryRelated(BasisDataType& q,std::vector<HilbertState>  const &basis,int n) const
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

				size_t na = HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,0) +
						HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,0+2);
				size_t nb = HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,1) +
						HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,1+2); //

				size_t flavor = na  + 3*nb;

				flavors.push_back(flavor);
				jmSaved = jmpair;

				// nup
				electronsUp[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_UP);
				// ndown
				electronsDown[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_DOWN);
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
		PairType calcJmvalue(const HilbertState& ket) const
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
		PairType calcJmValueAux(const HilbertState& ket) const
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

		//! SU(2) symmetry related block
		//! Let |9> = |up a down b> and
		//! Let |6> = |up b down a>  then
		void reinterpret(PsimagLite::Matrix<SparseElementType>& cm,std::vector<HilbertState>  const &basis) const
		{
			int n  = cm.n_row();
			if (n!=16) throw std::runtime_error("reinterpret (unimplemented case): "
						"blocks must be of size 1, and basis of size 16\n");
			PsimagLite::Matrix<SparseElementType> cmCopy(n,n);
			int i,j;
			int x=PsimagLite::isInVector(basis,reinterpretX_);
			int y=PsimagLite::isInVector(basis,reinterpretY_);

			RealType factor = 0.7071067811865475244;
			for (i=0;i<n;i++) {
				if (i==x || i==y) continue;
				for (j=0;j<n;j++) {
					if (j==x || j==y) continue;
					cmCopy(i,j)=cm(i,j);
				}
			}
			for (j=0;j<n;j++) {
				if (j==x || j==y) continue;
				cmCopy(x,j)=factor*(cm(x,j)+cm(y,j));
				cmCopy(y,j)=factor*(cm(x,j)-cm(y,j));
			}
			for (i=0;i<n;i++) {
				if (i==x || i==y) continue;
				cmCopy(i,x)=factor*(cm(i,x)+cm(i,y));
				cmCopy(i,y)=factor*(cm(i,x)-cm(i,y));
			}
			RealType sf= -1;
			cmCopy(x,x)=0.5*(cm(x,x)+cm(x,y)+cm(y,x)+cm(y,y));
			cmCopy(x,y)=0.5*(cm(x,x)-sf*cm(x,y)+sf*cm(y,x)-cm(y,y));
			cmCopy(y,x)=0.5*(cm(x,x)+sf*cm(x,y)-sf*cm(y,x)-cm(y,y));
			cmCopy(y,y)=0.5*(cm(x,x)-cm(x,y)-cm(y,x)+cm(y,y));

			cm = cmCopy;	
		}

		int nini(const HilbertState& state,size_t i,size_t gamma1,size_t spin1,size_t gamma2,size_t spin2) const
		{
			int tmp1 = HilbertSpaceFeAsType::calcNofElectrons(state,i,gamma1+spin1*NUMBER_OF_ORBITALS);
		 	int tmp2 = HilbertSpaceFeAsType::calcNofElectrons(state,i,gamma2+spin2*NUMBER_OF_ORBITALS);
			return tmp1*tmp2;
		}

		//! Full hamiltonian from creation matrices cm
		void calcHamiltonian(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,Block const &block) const
		{
			size_t n=block.size();
			SparseMatrixType tmpMatrix,tmpMatrix2;

			hmatrix.makeDiagonal(cm[0].data.rank());
			
			for (size_t i=0;i<n;i++) {
				//! hopping part
				for (size_t j=0;j<n;j++) {
					for (size_t term=0;term<geometry_.terms();term++) {
						for (size_t dofs=0;dofs<LinkProductType::dofs(term);dofs++) {
							std::pair<size_t,size_t> edofs = LinkProductType::connectorDofs(term,dofs);
							RealType tmp = geometry_(block[i],edofs.first,block[j],edofs.second,term);
						
							if (i==j || tmp==0.0) continue;

							size_t spin = dofs/4;
							size_t xtmp = (spin==0) ? 0 : 4;
							xtmp = dofs - xtmp;
							size_t orb1 = xtmp/2;
							size_t orb2 = (xtmp & 1);
							
							size_t dof1 = orb1 + spin*2;
							size_t dof2 = orb2 + spin*2;
							transposeConjugate(tmpMatrix2,cm[dof2+j*DEGREES_OF_FREEDOM].data);
							multiply(tmpMatrix,cm[dof1+i*DEGREES_OF_FREEDOM].data,tmpMatrix2);
							multiplyScalar(tmpMatrix2,tmpMatrix,tmp);
							hmatrix += tmpMatrix2;
						}
					}
				}
				addInteraction(hmatrix,cm,i);
			}
		}

		void addInteraction(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
		{
			addInteractionU1(hmatrix,cm,i);
			addInteractionU2(hmatrix,cm,i);
			addInteractionJ1(hmatrix,cm,i);
			addInteractionJ2(hmatrix,cm,i);
		}

		//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
		void addInteractionU1(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
		{
			int dof=DEGREES_OF_FREEDOM;
			SparseMatrixType tmpMatrix,tmpMatrix2;

			for (size_t alpha=0;alpha<size_t(NUMBER_OF_ORBITALS);alpha++) {
				SparseMatrixType m1=cm[alpha+SPIN_UP*NUMBER_OF_ORBITALS+i*dof].data;
				SparseMatrixType m2=cm[alpha+SPIN_DOWN*NUMBER_OF_ORBITALS+i*dof].data;

				multiply(tmpMatrix,n(m1),n(m2));
				multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[0]); // this is U
				hmatrix += tmpMatrix2;
			}
		}

		//! Term is U[1] n_{i BAND0 } n_{i BAND1}
		void addInteractionU2(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
		{
			size_t orbital0=0,orbital1=1;
			SparseMatrixType tmpMatrix,tmpMatrix2;

			multiply(tmpMatrix, nSummedOverSpin(cm,i,orbital0),nSummedOverSpin(cm,i,orbital1));
			multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[1]);// this is U'-J/2
			hmatrix += tmpMatrix2;
		}

		//! Term is U[2] S_{i BAND0 } S_{i BAND1}
		void addInteractionJ1(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
		{
			size_t orbital0=0,orbital1=1;
			SparseMatrixType tmpMatrix2,tmpMatrix;
			RealType val=0;
			RealType val2=2.0;
			RealType val3=4.0;

			multiply(tmpMatrix, spinOperator(cm,i,orbital0,0),spinOperator(cm,i,orbital1,1));
			val = modelParameters_.hubbardU[2]/val2;
			multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
			hmatrix += tmpMatrix2;

			multiply(tmpMatrix, spinOperator(cm,i,orbital0,1),spinOperator(cm,i,orbital1,0));
			val = modelParameters_.hubbardU[2]/val2;
			multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
			hmatrix += tmpMatrix2;

			multiply(tmpMatrix, spinOperator(cm,i,orbital0,2),spinOperator(cm,i,orbital1,2));
			val = modelParameters_.hubbardU[2]/val3;
			multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
			hmatrix += tmpMatrix2;
			return;
		}

		//! Term is U[3] \sum_{\alpha}\bar{n}_{i\alpha UP} \bar{n}_{i\alpha DOWN} 
		//! where \bar{n}_{i\alpha \spin} = c^\dagger_{i\alpha\spin} c_{i\bar{\alpha}\bar{spin}}
		void addInteractionJ2(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
		{
			SparseMatrixType tmpMatrix,tmpMatrix2;

			for (size_t alpha=0;alpha<size_t(NUMBER_OF_ORBITALS);alpha++) {
				multiply(tmpMatrix,nBar(cm,i,alpha,SPIN_UP),nBar(cm,i,alpha,SPIN_DOWN));
				multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[3]); // this is -J
				hmatrix += tmpMatrix2;
			}
		}

		SparseMatrixType n(const SparseMatrixType& c) const
		{
			SparseMatrixType tmpMatrix;
			SparseMatrixType cdagger;
			transposeConjugate(cdagger,c);
			multiply(tmpMatrix,c,cdagger);

			return tmpMatrix;
		}

		SparseMatrixType nBar(const std::vector<OperatorType>& cm,size_t i,size_t orbital,size_t spin) const
		{
			SparseMatrixType tmpMatrix,cdagger=cm[orbital+spin*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data;
			SparseMatrixType cbar;
			transposeConjugate(cbar,cm[1-orbital+(1-spin)*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
			multiply(tmpMatrix,cdagger,cbar);
			return tmpMatrix;
		}

		SparseMatrixType nSummedOverSpin(const std::vector<OperatorType>& cm,size_t i,size_t orbital) const
		{
			SparseMatrixType tmpMatrix = n(cm[orbital+SPIN_UP*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
			tmpMatrix += n(cm[orbital+SPIN_DOWN*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
			return tmpMatrix;
		}

		SparseMatrixType spinOperator(const std::vector<OperatorType>& cm,size_t i,size_t orbital,size_t component) const
		{
			switch (component) {
				case 0: // S^+
					return spinOperatorAux(cm,i,orbital,SPIN_UP,SPIN_DOWN);
					break;
				case 1: // S^-
					return spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_UP);
					break;
			}
			SparseMatrixType tmpMatrix=spinOperatorAux(cm,i,orbital,SPIN_UP,SPIN_UP);
			SparseMatrixType tmpMatrix2=spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_DOWN);
			SparseMatrixType tmpMatrix3;
			multiplyScalar(tmpMatrix3,tmpMatrix2,-1.0);
			tmpMatrix += tmpMatrix3;
			return tmpMatrix;
		}

		SparseMatrixType spinOperatorAux(const std::vector<OperatorType>& cm,size_t i,size_t orbital,size_t spin1,size_t spin2) const
		{
			SparseMatrixType result,temp;
			transposeConjugate(temp,cm[orbital+spin2*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
			multiply(
				result, // =
				cm[orbital+spin1*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data, // times
				temp
			);
					
			return result;
		}

		void setPauliMatrix()
		{
			PsimagLite::Matrix<RealType> matrixTmp(2,2);
			// x component
			matrixTmp(0,0)=matrixTmp(1,0)=matrixTmp(0,1)=matrixTmp(1,1)=0;
			matrixTmp(1,0)=matrixTmp(0,1)=1;
			pauliMatrix_.push_back(matrixTmp);
			// y component
			matrixTmp(0,0)=matrixTmp(1,0)=matrixTmp(0,1)=matrixTmp(1,1)=0;
			matrixTmp(1,0)=matrixTmp(0,1)=1;
			pauliMatrix_.push_back(matrixTmp);
			// z component
			matrixTmp(0,0)=matrixTmp(1,0)=matrixTmp(0,1)=matrixTmp(1,1)=0;
			matrixTmp(0,0)= 1; matrixTmp(1,1)= -1;
			pauliMatrix_.push_back(matrixTmp);
		}

		void diagTest(const SparseMatrixType& fullm,const std::string& str) const
		{
			if (fullm.rank()!=256) return;
			PsimagLite::Matrix<SparseElementType> fullm2;
			crsMatrixToFullMatrix(fullm2,fullm);
			std::vector<SparseElementType> eigs(fullm2.n_row());
			PsimagLite::diag(fullm2,eigs,'V');
			std::cout<<str<<" diagTest size="<<fullm.rank()<<" eigs[0]="<<eigs[0]<<"\n";
			std::cout<<fullm;
		}
	};     //class ModelFeBasedSc

	template<
		typename ModelHelperType,
		typename SparseMatrixType,
		typename GeometryType,
  		template<typename> class SharedMemoryTemplate
		>
	std::ostream &operator<<(std::ostream &os,const ModelFeBasedSc<
		ModelHelperType,
		SparseMatrixType,
		GeometryType,
  		SharedMemoryTemplate
		>& model)
	{
		model.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif
