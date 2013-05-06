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
//		typedef typename ModelHelperType::ReflectionSymmetryType ReflectionSymmetryType;
		
		static const int FERMION_SIGN = -1;
		static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
		static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

	public:
		typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
		typedef typename PsimagLite::Vector<HilbertState>::Type HilbertBasisType;
		typedef LinkProductFeAs<ModelHelperType> LinkProductType;
		typedef   ModelBase<ModelHelperType,SparseMatrixType,GeometryType,LinkProductType,SharedMemoryTemplate> ModelBaseType;
		typedef	 typename ModelBaseType::MyBasis MyBasis;
		typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef typename ModelBaseType::InputValidatorType InputValidatorType;

		static size_t const REINTERPRET  = 1;
		ModelFeBasedSc(InputValidatorType& io,GeometryType const &geometry,ConcurrencyType& concurrency)
			: ModelBaseType(geometry,concurrency),
			  reinterpretX_(6),
			  reinterpretY_(9),
			  modelParameters_(io),
			  geometry_(geometry),
			  spinSquared_(spinSquaredHelper_,modelParameters_.orbitals,2*modelParameters_.orbitals)
		{
			LinkProductType::setOrbitals(modelParameters_.orbitals);
			HilbertSpaceFeAsType::setOrbitals(modelParameters_.orbitals);
//			setPauliMatrix();
			if (modelParameters_.potentialV.size()!=2*modelParameters_.orbitals*geometry.numberOfSites()) {
				std::string str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "potentialV length must be 2*orbitals times the number of sites\n";
				throw std::runtime_error(str.c_str());
			}
		}

// 		size_t orbitals() const { return NUMBER_OF_ORBITALS; }

		size_t hilbertSize(size_t site) const { return (size_t)pow(2,modelParameters_.orbitals*2); }

		void print(std::ostream& os) const { operator<<(os,modelParameters_); }

		//! find creation operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
				     SparseMatrixType &hamiltonian,
				     BasisDataType &q,
				     Block const &block,
				     const RealType& time)  const
		{
			typename PsimagLite::Vector<HilbertState>::Type natBasis;
			typename PsimagLite::Vector<size_t>::Type qvector;
			setNaturalBasis(natBasis,qvector,block);			

			setOperatorMatrices(creationMatrix,block);

			//! Set symmetry related
			setSymmetryRelated(q,natBasis,block.size());

			//! set hamiltonian
			calcHamiltonian(hamiltonian,creationMatrix,block);

			SparseMatrixType tmpMatrix2;
			tmpMatrix2.makeDiagonal(natBasis.size(),0.0);
		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,Block const &block) const
		{
			typename PsimagLite::Vector<HilbertState>::Type natBasis;
			SparseMatrixType tmpMatrix;
			typename PsimagLite::Vector<size_t>::Type qvector;
			setNaturalBasis(natBasis,qvector,block);

			//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
			creationMatrix.clear();
			size_t dofs = 2*modelParameters_.orbitals;
			for (size_t i=0;i<block.size();i++) {
				for (size_t sigma=0;sigma<dofs;sigma++) {
					findOperatorMatrices(tmpMatrix,i,sigma,natBasis);
					size_t m=0;
					int asign=1;
					if (sigma>modelParameters_.orbitals-1) {
						m=1;
						asign= -1;
					}
					typename OperatorType::Su2RelatedType su2related;
					if (sigma <modelParameters_.orbitals) {
						su2related.source.push_back(i*dofs+sigma);
						su2related.source.push_back(i*dofs+sigma + modelParameters_.orbitals);
						su2related.transpose.push_back(-1);
						su2related.transpose.push_back(-1);
						su2related.offset = modelParameters_.orbitals;
					}	
					OperatorType myOp(tmpMatrix,-1,typename OperatorType::PairType(1,m),asign,su2related);
					creationMatrix.push_back(myOp);
				}
			}
		}

		PsimagLite::Matrix<SparseElementType> naturalOperator(const std::string& what,size_t site,size_t dof) const
		{
			Block block;
			block.resize(1);
			block[0]=site;
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			setOperatorMatrices(creationMatrix,block);
			size_t orbital = dof % modelParameters_.orbitals;
			size_t spin = dof/modelParameters_.orbitals;
			assert(creationMatrix.size()>0);
			size_t nrow = creationMatrix[0].data.row();

			if (what=="+" or what=="i") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (size_t x=0;x<modelParameters_.orbitals;x++)
					tmp += multiplyTc(creationMatrix[x].data,creationMatrix[x+modelParameters_.orbitals].data);
				return tmp;
			}
			if (what=="-") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (size_t x=0;x<modelParameters_.orbitals;x++)
					tmp += multiplyTc(creationMatrix[x+modelParameters_.orbitals].data,creationMatrix[x].data);
				return tmp;
			}
			if (what=="z") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				PsimagLite::Matrix<SparseElementType> tmp2(nrow,nrow);
				for (size_t x=0;x<modelParameters_.orbitals;x++) {
					tmp += multiplyTc(creationMatrix[x].data,creationMatrix[x].data);
					tmp2 += multiplyTc(creationMatrix[x+modelParameters_.orbitals].data,creationMatrix[x+modelParameters_.orbitals].data);
				}
				return tmp-tmp2;
			}
			if (what=="n") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (size_t x=0;x<2*modelParameters_.orbitals;x++)
					tmp += multiplyTc(creationMatrix[x].data,creationMatrix[x].data);
				return tmp;
			} 
			if (what=="c") {
				PsimagLite::Matrix<SparseElementType> tmp;
				SparseMatrixType cdagger;
				transposeConjugate(cdagger,creationMatrix[orbital + spin*modelParameters_.orbitals].data);
				crsMatrixToFullMatrix(tmp,cdagger);
				return tmp;
			}
			if (what=="d") { // delta = c^\dagger * c^dagger
				SparseMatrixType atmp;
				multiply(atmp,creationMatrix[orbital+orbital+modelParameters_.orbitals].data,creationMatrix[orbital].data);
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,atmp);
				return tmp;
			}
			if (what=="identity") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (size_t i=0;i<tmp.n_row();i++) tmp(i,i) = 1.0;
				return tmp;
			}
			std::cerr<<"what="<<what<<"\n";
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}
		
		//! find all states in the natural basis for a block of n sites
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(typename PsimagLite::Vector<HilbertState>  ::Type&basis,
		                     typename PsimagLite::Vector<size_t>::Type& q,
		                     const typename PsimagLite::Vector<size_t>::Type& block) const
		{
			assert(block.size()==1);
			HilbertState a=0;
			int sitesTimesDof=2*modelParameters_.orbitals;
			HilbertState total = (1<<sitesTimesDof);

			typename PsimagLite::Vector<HilbertState>::Type  basisTmp;
			for (a=0;a<total;a++) basisTmp.push_back(a);

			// reorder the natural basis (needed for MULTIPLE BANDS)
			findQuantumNumbers(q,basisTmp,block.size());
			typename PsimagLite::Vector<size_t>::Type iperm(q.size());
			PsimagLite::Sort<typename PsimagLite::Vector<size_t>::Type > sort;
			sort.sort(q,iperm);
			basis.clear();
			for (a=0;a<total;a++) basis.push_back(basisTmp[iperm[a]]);
		}
		
		void findElectrons(typename PsimagLite::Vector<size_t>::Type& electrons,
		                   const typename PsimagLite::Vector<HilbertState>::Type& basis,
		                   size_t site) const
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
		ParametersModelFeAs<RealType>  modelParameters_;
		GeometryType const &geometry_;
		SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
		SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
		//typename PsimagLite::Vector<PsimagLite::Matrix<RealType>::Type > pauliMatrix_;

		//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS 
		RealType sign(HilbertState const &ket, int i,size_t sigma) const
		{
			int value=0;
			size_t dofs=2*modelParameters_.orbitals;
			for (size_t alpha=0;alpha<dofs;alpha++)
				value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);
			// add electron on site 0 if needed
			if (i>0) value += HilbertSpaceFeAsType::electronsAtGivenSite(ket,0);

			//order for sign is: a up, a down, b up, b down, etc
			unsigned int x = HilbertSpaceFeAsType::get(ket,i);
			int spin = sigma/modelParameters_.orbitals;
			size_t orb = sigma % modelParameters_.orbitals;

			for (size_t j=0;j<orb;j++) {
				for (size_t k=0;k<2;k++) {
					size_t ind = j + k * modelParameters_.orbitals;
					int mask = (1<<ind);
					if (x & mask) value++;
				}
			}
			if (spin==SPIN_DOWN) {
				int mask = (1<<orb);
				if (x & mask) value++;
			}
			if (value==0 || value%2==0) return 1.0;

			return FERMION_SIGN;
		}

		//! Find c^\dagger_i\gamma\sigma in the natural basis natBasis
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void findOperatorMatrices(SparseMatrixType& creationMatrix,
		                          int i,
		                          int sigma,
		                          const typename PsimagLite::Vector<HilbertState>::Type& natBasis) const
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
			if (REINTERPRET && modelParameters_.orbitals==2) reinterpret(cm,natBasis);

			SparseMatrixType temp;
			fullMatrixToCrsMatrix(temp,cm);
			transposeConjugate(creationMatrix,temp);
		}

		void findQuantumNumbers(typename PsimagLite::Vector<size_t>::Type& q,const typename PsimagLite::Vector<HilbertState>  ::Type&basis,int n) const
		{
			BasisDataType qq;
			setSymmetryRelated(qq,basis,n);
			MyBasis::findQuantumNumbers(q,qq);
		}

		void setSymmetryRelated(BasisDataType& q,
		                        const typename PsimagLite::Vector<HilbertState>::Type& basis,
		                        int n) const
		{
			if (n!=1) std::runtime_error("ModelFeAs::setSymmetryRelated() implemented for n=1 only\n");

			// find j,m and flavors (do it by hand since we assume n==1)
			// note: we use 2j instead of j
			// note: we use m+j instead of m
			// This assures us that both j and m are size_t
			typedef std::pair<size_t,size_t> PairType;
			typename PsimagLite::Vector<PairType>::Type jmvalues;
			typename PsimagLite::Vector<size_t>::Type flavors; 
			PairType jmSaved = calcJmvalue<PairType>(basis[0]);
			jmSaved.first++;
			jmSaved.second++;

			typename PsimagLite::Vector<size_t>::Type electronsUp(basis.size());
			typename PsimagLite::Vector<size_t>::Type electronsDown(basis.size());
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
			if (modelParameters_.orbitals!=2) return jm;
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
		void reinterpret(PsimagLite::Matrix<SparseElementType>& cm,
		                 const typename PsimagLite::Vector<HilbertState>::Type& basis) const
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
			int tmp1 = HilbertSpaceFeAsType::calcNofElectrons(state,i,gamma1+spin1*modelParameters_.orbitals);
			int tmp2 = HilbertSpaceFeAsType::calcNofElectrons(state,i,gamma2+spin2*modelParameters_.orbitals);
			return tmp1*tmp2;
		}

		//! Full hamiltonian from creation matrices cm
		void calcHamiltonian(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,Block const &block) const
		{
			size_t n=block.size();
			SparseMatrixType tmpMatrix,tmpMatrix2;

			hmatrix.makeDiagonal(cm[0].data.row());
			
			size_t orbitalsSquared = (modelParameters_.orbitals*modelParameters_.orbitals);
			size_t degreesOfFreedom = 2*modelParameters_.orbitals;

			for (size_t i=0;i<n;i++) {
				//! hopping part
				for (size_t j=0;j<n;j++) {
					for (size_t term=0;term<geometry_.terms();term++) {
						typename GeometryType::AdditionalDataType additional;
						geometry_.fillAdditionalData(additional,term,block[i],block[j]);
						size_t dofsTotal = LinkProductType::dofs(term,additional);
						for (size_t dofs=0;dofs<dofsTotal;dofs++) {
							std::pair<size_t,size_t> edofs = LinkProductType::connectorDofs(term,dofs,additional);
							RealType tmp = geometry_(block[i],edofs.first,block[j],edofs.second,term);
						
							if (i==j || tmp==0.0) continue;

							size_t spin = dofs/orbitalsSquared;
							size_t xtmp = (spin==0) ? 0 : orbitalsSquared;
							xtmp = dofs - xtmp;
							size_t orb1 = xtmp/modelParameters_.orbitals;
							size_t orb2 = xtmp % modelParameters_.orbitals;
							
							size_t dof1 = orb1 + spin*modelParameters_.orbitals;
							size_t dof2 = orb2 + spin*modelParameters_.orbitals;

							transposeConjugate(tmpMatrix2,cm[dof2+j*degreesOfFreedom].data);
							multiply(tmpMatrix,cm[dof1+i*degreesOfFreedom].data,tmpMatrix2);
							multiplyScalar(tmpMatrix2,tmpMatrix,tmp);
							hmatrix += tmpMatrix2;
						}
					}
				}
				addInteraction(hmatrix,cm,i);
				addMagneticField(hmatrix,cm,i,block[i]);
				addPotentialV(hmatrix,cm,i,block[i]);
			}
		}

		void addInteraction(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i) const
		{
			addInteractionU1(hmatrix,cm,i);
			addInteractionU2(hmatrix,cm,i);
			addInteractionJ1(hmatrix,cm,i);
			addInteractionJ2(hmatrix,cm,i);
		}

		//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
		void addInteractionU1(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i) const
		{
			int dof=2*modelParameters_.orbitals;
			SparseMatrixType tmpMatrix,tmpMatrix2;

			for (size_t alpha=0;alpha<size_t(modelParameters_.orbitals);alpha++) {
				SparseMatrixType m1=cm[alpha+SPIN_UP*modelParameters_.orbitals+i*dof].data;
				SparseMatrixType m2=cm[alpha+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;

				multiply(tmpMatrix,n(m1),n(m2));
				multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[0]); // this is U
				hmatrix += tmpMatrix2;
			}
		}

		//! Term is U[1] n_{i BAND0 } n_{i BAND1}
		void addInteractionU2(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i) const
		{

			for (size_t orb1=0;orb1<modelParameters_.orbitals;orb1++) {
				for (size_t orb2=orb1+1;orb2<modelParameters_.orbitals;orb2++) {
				SparseMatrixType tmpMatrix,tmpMatrix2;
				multiply(tmpMatrix, nSummedOverSpin(cm,i,orb1),nSummedOverSpin(cm,i,orb2));
				multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[1]);// this is U'-J/2
				hmatrix += tmpMatrix2;
				}
			}
		}

		//! Term is U[2] S_{i BAND0 } S_{i BAND1}
		void addInteractionJ1(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i) const
		{
			RealType val=0;
			RealType val2=2.0;
			RealType val3=4.0;

			for (size_t orb1=0;orb1<modelParameters_.orbitals;orb1++) {
				for (size_t orb2=orb1+1;orb2<modelParameters_.orbitals;orb2++) {
					SparseMatrixType tmpMatrix2,tmpMatrix;


					multiply(tmpMatrix, spinOperator(cm,i,orb1,0),spinOperator(cm,i,orb2,1));
					val = modelParameters_.hubbardU[2]/val2;
					multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
					hmatrix += tmpMatrix2;

					multiply(tmpMatrix, spinOperator(cm,i,orb1,1),spinOperator(cm,i,orb2,0));
					val = modelParameters_.hubbardU[2]/val2;
					multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
					hmatrix += tmpMatrix2;

					multiply(tmpMatrix, spinOperator(cm,i,orb1,2),spinOperator(cm,i,orb2,2));
					val = modelParameters_.hubbardU[2]/val3;
					multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
					hmatrix += tmpMatrix2;
				}
			}
		}

		//! Term is U[3] \sum_{\alpha}\bar{n}_{i\alpha UP} \bar{n}_{i\alpha DOWN} 
		//! where \bar{n}_{i\alpha \spin} = c^\dagger_{i\alpha\spin} c_{i\bar{\alpha}\bar{spin}}
		void addInteractionJ2(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i) const
		{
			SparseMatrixType tmpMatrix,tmpMatrix2;

			for (size_t orb1=0;orb1<modelParameters_.orbitals;orb1++) {
				for (size_t orb2=0;orb2<modelParameters_.orbitals;orb2++) {
					if (orb1==orb2) continue;
					multiply(tmpMatrix,nBar(cm,i,orb1,orb2,SPIN_UP),nBar(cm,i,orb1,orb2,SPIN_DOWN));
					multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[3]); // this is -J
					hmatrix += tmpMatrix2;
				}
			}
		}

		void addMagneticField(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t actualIndexOfSite) const
		{
			if (modelParameters_.magneticField.n_row()<3) return;
			assert(actualIndexOfSite<modelParameters_.magneticField.n_col());
			for (size_t orb=0;orb<modelParameters_.orbitals;orb++)
				addMagneticField(hmatrix,cm,i,actualIndexOfSite,orb);
		}

		void addMagneticField(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t actualIndexOfSite,size_t orbital) const
		{
			int dof=2*modelParameters_.orbitals;
			SparseMatrixType cup = cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data;
			SparseMatrixType cupTranspose;
			transposeConjugate(cupTranspose,cup);
			SparseMatrixType cdown = cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;
			SparseMatrixType A = cupTranspose * cdown;
			SparseMatrixType Atranspose;
			transposeConjugate(Atranspose,A);

			hmatrix += modelParameters_.magneticField(0,actualIndexOfSite) * A;

			hmatrix += modelParameters_.magneticField(1,actualIndexOfSite) * Atranspose;

			SparseMatrixType nup = n(cup);
			SparseMatrixType ndown = n(cdown);

			SparseMatrixType tmp = nup;
			tmp += (-1.0)*ndown;
			hmatrix += modelParameters_.magneticField(2,actualIndexOfSite) * tmp;

		}

		void addPotentialV(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t actualIndexOfSite) const
		{
			for (size_t orb=0;orb<modelParameters_.orbitals;orb++)
				addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb);
		}

		void addPotentialV(SparseMatrixType &hmatrix,const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t actualIndexOfSite,size_t orbital) const
		{
			int dof=2*modelParameters_.orbitals;
			SparseMatrixType nup = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data);
			SparseMatrixType ndown = n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data);


			size_t linSize = geometry_.numberOfSites();

			size_t iUp = actualIndexOfSite + (orbital + 0*modelParameters_.orbitals)*linSize;
			hmatrix += modelParameters_.potentialV[iUp] * nup;
			size_t iDown = actualIndexOfSite + (orbital + 1*modelParameters_.orbitals)*linSize;
			hmatrix += modelParameters_.potentialV[iDown] * ndown;
		}

		SparseMatrixType n(const SparseMatrixType& c) const
		{
			SparseMatrixType tmpMatrix;
			SparseMatrixType cdagger;
			transposeConjugate(cdagger,c);
			multiply(tmpMatrix,c,cdagger);

			return tmpMatrix;
		}

		SparseMatrixType nBar(const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t orb1,size_t orb2,size_t spin) const
		{
			size_t dofs = 2 * modelParameters_.orbitals;
			SparseMatrixType tmpMatrix,cdagger=cm[orb1+spin*modelParameters_.orbitals+i*dofs].data;
			SparseMatrixType cbar;
			transposeConjugate(cbar,cm[orb2+(1-spin)*modelParameters_.orbitals+i*dofs].data);
			multiply(tmpMatrix,cdagger,cbar);
			return tmpMatrix;
		}

		SparseMatrixType nSummedOverSpin(const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t orbital) const
		{
			size_t dofs = 2 * modelParameters_.orbitals;
			SparseMatrixType tmpMatrix = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dofs].data);
			tmpMatrix += n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dofs].data);
			return tmpMatrix;
		}

		SparseMatrixType spinOperator(const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t orbital,size_t component) const
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

		SparseMatrixType spinOperatorAux(const typename PsimagLite::Vector<OperatorType>::Type& cm,size_t i,size_t orbital,size_t spin1,size_t spin2) const
		{
			size_t dofs = 2 * modelParameters_.orbitals;
			SparseMatrixType result,temp;
			transposeConjugate(temp,cm[orbital+spin2*modelParameters_.orbitals+i*dofs].data);
			multiply(
				result, // =
				cm[orbital+spin1*modelParameters_.orbitals+i*dofs].data, // times
				temp
			);
					
			return result;
		}

//		void setPauliMatrix()
//		{
//			PsimagLite::Matrix<RealType> matrixTmp(2,2);
//			// x component
//			matrixTmp(0,0)=matrixTmp(1,0)=matrixTmp(0,1)=matrixTmp(1,1)=0;
//			matrixTmp(1,0)=matrixTmp(0,1)=1;
//			pauliMatrix_.push_back(matrixTmp);
//			// y component
//			matrixTmp(0,0)=matrixTmp(1,0)=matrixTmp(0,1)=matrixTmp(1,1)=0;
//			matrixTmp(1,0)=matrixTmp(0,1)=1;
//			pauliMatrix_.push_back(matrixTmp);
//			// z component
//			matrixTmp(0,0)=matrixTmp(1,0)=matrixTmp(0,1)=matrixTmp(1,1)=0;
//			matrixTmp(0,0)= 1; matrixTmp(1,1)= -1;
//			pauliMatrix_.push_back(matrixTmp);
//		}

		void diagTest(const SparseMatrixType& fullm,const std::string& str) const
		{
			if (fullm.rank()!=256) return;
			PsimagLite::Matrix<SparseElementType> fullm2;
			crsMatrixToFullMatrix(fullm2,fullm);
			typename PsimagLite::Vector<SparseElementType>::Type eigs(fullm2.n_row());
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
