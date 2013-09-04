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
#include "ModelCommon.h"

namespace Dmrg {
	template<typename ModelBaseType>
	class ModelFeBasedSc : public ModelBaseType {
		
		typedef unsigned int long long WordType;
		typedef  HilbertSpaceFeAs<WordType> HilbertSpaceFeAsType;

	public:

		typedef typename ModelBaseType::ModelHelperType ModelHelperType;
		typedef typename ModelBaseType::GeometryType GeometryType;
		typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
		typedef typename ModelBaseType::LinkProductStructType LinkProductStructType;
		typedef typename ModelBaseType::LinkType LinkType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef typename HilbertSpaceFeAsType::HilbertState HilbertState;
		typedef typename ModelHelperType::BlockType BlockType;
		typedef typename ModelBaseType::SolverParamsType SolverParamsType;
		typedef typename ModelBaseType::VectorType VectorType;

		static const int FERMION_SIGN = -1;
		static const int SPIN_UP=HilbertSpaceFeAsType::SPIN_UP;
		static const int SPIN_DOWN=HilbertSpaceFeAsType::SPIN_DOWN;

		typedef typename PsimagLite::Vector<HilbertState>::Type HilbertBasisType;
		typedef LinkProductFeAs<ModelHelperType> LinkProductType;
		typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
		typedef	 typename ModelBaseType::MyBasis MyBasis;
		typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef typename ModelBaseType::InputValidatorType InputValidatorType;

		static SizeType const REINTERPRET  = 1;
		ModelFeBasedSc(const SolverParamsType& solverParams,
		               InputValidatorType& io,
		               GeometryType const &geometry)
			: ModelBaseType(solverParams,io,geometry,new ModelCommonType(geometry)),
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
				PsimagLite::String str(__FILE__);
				str += " " + ttos(__LINE__) + "\n";
				str += "potentialV length must be 2*orbitals times the number of sites\n";
				throw PsimagLite::RuntimeError(str.c_str());
			}
		}

// 		SizeType orbitals() const { return NUMBER_OF_ORBITALS; }

		SizeType hilbertSize(SizeType site) const { return (SizeType)pow(2,modelParameters_.orbitals*2); }

		void print(std::ostream& os) const { operator<<(os,modelParameters_); }

		//! find creation operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(VectorOperatorType& creationMatrix,
		                     SparseMatrixType &hamiltonian,
		                     BasisDataType &q,
		                     const BlockType& block,
		                     const RealType& time)  const
		{
			typename PsimagLite::Vector<HilbertState>::Type natBasis;
			typename PsimagLite::Vector<SizeType>::Type qvector;
			setNaturalBasis(natBasis,qvector,block);			

			setOperatorMatrices(creationMatrix,block);

			//! Set symmetry related
			setSymmetryRelated(q,natBasis,block.size());

			//! set hamiltonian
			this->calcHamiltonian(hamiltonian,creationMatrix,block,time);

			SparseMatrixType tmpMatrix2;
			tmpMatrix2.makeDiagonal(natBasis.size(),0.0);
		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(VectorOperatorType& creationMatrix,
		                         const BlockType& block) const
		{
			typename PsimagLite::Vector<HilbertState>::Type natBasis;
			SparseMatrixType tmpMatrix;
			typename PsimagLite::Vector<SizeType>::Type qvector;
			setNaturalBasis(natBasis,qvector,block);

			//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
			creationMatrix.clear();
			SizeType dofs = 2*modelParameters_.orbitals;
			for (SizeType i=0;i<block.size();i++) {
				for (SizeType sigma=0;sigma<dofs;sigma++) {
					findOperatorMatrices(tmpMatrix,i,sigma,natBasis);
					SizeType m=0;
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

		PsimagLite::Matrix<SparseElementType> naturalOperator(const PsimagLite::String& what,
		                                                      SizeType site,
		                                                      SizeType dof) const
		{
			BlockType block;
			block.resize(1);
			block[0]=site;
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			setOperatorMatrices(creationMatrix,block);
			SizeType orbital = dof % modelParameters_.orbitals;
			SizeType spin = dof/modelParameters_.orbitals;
			assert(creationMatrix.size()>0);
			SizeType nrow = creationMatrix[0].data.row();

			if (what == "i" || what=="identity") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (SizeType i = 0; i < tmp.n_row(); ++i) tmp(i,i) = 1.0;
				return tmp;
			}

			if (what=="+") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (SizeType x=0;x<modelParameters_.orbitals;x++)
					tmp += multiplyTc(creationMatrix[x].data,creationMatrix[x+modelParameters_.orbitals].data);
				return tmp;
			}
			if (what=="-") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (SizeType x=0;x<modelParameters_.orbitals;x++)
					tmp += multiplyTc(creationMatrix[x+modelParameters_.orbitals].data,creationMatrix[x].data);
				return tmp;
			}
			if (what=="z") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				PsimagLite::Matrix<SparseElementType> tmp2(nrow,nrow);
				for (SizeType x=0;x<modelParameters_.orbitals;x++) {
					tmp += multiplyTc(creationMatrix[x].data,creationMatrix[x].data);
					tmp2 += multiplyTc(creationMatrix[x+modelParameters_.orbitals].data,creationMatrix[x+modelParameters_.orbitals].data);
				}
				return tmp-tmp2;
			}
			if (what=="n") {
				PsimagLite::Matrix<SparseElementType> tmp =
				        multiplyTc(creationMatrix[dof].data,creationMatrix[dof].data);
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

			std::cerr<<"what="<<what<<"\n";
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}
		
		//! find all states in the natural basis for a block of n sites
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(typename PsimagLite::Vector<HilbertState>  ::Type&basis,
		                     typename PsimagLite::Vector<SizeType>::Type& q,
		                     const typename PsimagLite::Vector<SizeType>::Type& block) const
		{
			assert(block.size()==1);
			HilbertState a=0;
			int sitesTimesDof=2*modelParameters_.orbitals;
			HilbertState total = (1<<sitesTimesDof);

			typename PsimagLite::Vector<HilbertState>::Type  basisTmp;
			for (a=0;a<total;a++) basisTmp.push_back(a);

			// reorder the natural basis (needed for MULTIPLE BANDS)
			findQuantumNumbers(q,basisTmp,block.size());
			typename PsimagLite::Vector<SizeType>::Type iperm(q.size());
			PsimagLite::Sort<typename PsimagLite::Vector<SizeType>::Type > sort;
			sort.sort(q,iperm);
			basis.clear();
			for (a=0;a<total;a++) basis.push_back(basisTmp[iperm[a]]);
		}
		
		void findElectrons(typename PsimagLite::Vector<SizeType>::Type& electrons,
		                   const typename PsimagLite::Vector<HilbertState>::Type& basis,
		                   SizeType site) const
		{
			electrons.resize(basis.size());
			for (SizeType i=0;i<basis.size();i++) {
				// nup
				SizeType nup = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_UP);
				// ndown
				SizeType ndown = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_DOWN);
				electrons[i] = nup + ndown;
			}
		}

		void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
		                                const VectorOperatorType& cm,
		                                const BlockType& block,
		                                RealType time,
		                                RealType factorForDiagonals=1.0) const
		{
			SizeType n=block.size();

			for (SizeType i=0;i<n;i++) {
				addInteraction(hmatrix,cm,i,factorForDiagonals);
				addMagneticField(hmatrix,cm,i,block[i],factorForDiagonals);

				if (modelParameters_.potentialT.size()==0 || time==0) {
				  addPotentialV(hmatrix,cm,i,block[i],factorForDiagonals,modelParameters_.potentialV);
				} else {
				  addPotentialV(hmatrix,cm,i,block[i],factorForDiagonals,modelParameters_.potentialT);
				}
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
		RealType sign(HilbertState const &ket, int i,SizeType sigma) const
		{
			int value=0;
			SizeType dofs=2*modelParameters_.orbitals;
			for (SizeType alpha=0;alpha<dofs;alpha++)
				value += HilbertSpaceFeAsType::calcNofElectrons(ket,0,i,alpha);
			// add electron on site 0 if needed
			if (i>0) value += HilbertSpaceFeAsType::electronsAtGivenSite(ket,0);

			//order for sign is: a up, a down, b up, b down, etc
			unsigned int x = HilbertSpaceFeAsType::get(ket,i);
			int spin = sigma/modelParameters_.orbitals;
			SizeType orb = sigma % modelParameters_.orbitals;

			for (SizeType j=0;j<orb;j++) {
				for (SizeType k=0;k<2;k++) {
					SizeType ind = j + k * modelParameters_.orbitals;
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

			for (SizeType ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				
				if (HilbertSpaceFeAsType::isNonZero(ket,i,sigma)) {
					
				} else {
					HilbertSpaceFeAsType::create(bra,i,sigma);
					int jj = PsimagLite::isInVector(natBasis,bra);
					if (jj<0) throw PsimagLite::RuntimeError("findOperatorMatrices: internal error while"
							"creating.\n");
					if (ii==SizeType(jj)) {
						std::cerr<<"ii="<<i<<" ket="<<ket<<" bra="<<bra<<" sigma="<<sigma<<"\n";
						throw PsimagLite::RuntimeError("Creation operator cannot be diagonal\n");
					}
					cm(ii,jj) =sign(ket,i,sigma);
				}
			}
			if (REINTERPRET && modelParameters_.orbitals==2) reinterpret(cm,natBasis);

			SparseMatrixType temp;
			fullMatrixToCrsMatrix(temp,cm);
			transposeConjugate(creationMatrix,temp);
		}

		void findQuantumNumbers(typename PsimagLite::Vector<SizeType>::Type& q,const typename PsimagLite::Vector<HilbertState>  ::Type&basis,int n) const
		{
			BasisDataType qq;
			setSymmetryRelated(qq,basis,n);
			MyBasis::findQuantumNumbers(q,qq);
		}

		void setSymmetryRelated(BasisDataType& q,
		                        const typename PsimagLite::Vector<HilbertState>::Type& basis,
		                        int n) const
		{
			if (n!=1) PsimagLite::RuntimeError("ModelFeAs::setSymmetryRelated() implemented for n=1 only\n");

			// find j,m and flavors (do it by hand since we assume n==1)
			// note: we use 2j instead of j
			// note: we use m+j instead of m
			// This assures us that both j and m are SizeType
			typedef std::pair<SizeType,SizeType> PairType;
			typename PsimagLite::Vector<PairType>::Type jmvalues;
			typename PsimagLite::Vector<SizeType>::Type flavors; 
			PairType jmSaved = calcJmvalue<PairType>(basis[0]);
			jmSaved.first++;
			jmSaved.second++;

			typename PsimagLite::Vector<SizeType>::Type electronsUp(basis.size());
			typename PsimagLite::Vector<SizeType>::Type electronsDown(basis.size());
			for (SizeType i=0;i<basis.size();i++) {
				PairType jmpair = calcJmvalue<PairType>(basis[i]);

				jmvalues.push_back(jmpair);

				SizeType na = HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,0) +
						HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,0+2);
				SizeType nb = HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,1) +
						HilbertSpaceFeAsType::calcNofElectrons(basis[i],0,1+2); //

				SizeType flavor = na  + 3*nb;

				flavors.push_back(flavor);
				jmSaved = jmpair;

				// nup
				electronsUp[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_UP);
				// ndown
				electronsDown[i] = HilbertSpaceFeAsType::electronsWithGivenSpin(basis[i],HilbertSpaceFeAsType::SPIN_DOWN);
			}
			q.jmValues=jmvalues;
			q.flavors = flavors;
			q.electrons = electronsUp + electronsDown;
			q.szPlusConst = electronsUp;
		}

		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		// Reinterprets 6 and 9
		template<typename PairType>
		PairType calcJmvalue(const HilbertState& ket) const
		{
			PairType jm(0,0);
			if (modelParameters_.orbitals!=2) return jm;
			SizeType x=reinterpretX_,y=reinterpretY_; // these states need reinterpretation

			if (ket==x) {
				jm=std::pair<SizeType,SizeType>(2,1);
			} else if (ket==y) {
				jm=std::pair<SizeType,SizeType>(0,0);
			} else jm=calcJmValueAux<PairType>(ket);
			//std::cerr<<jm.first<<" "<<jm.second<<" |--------------\n";
			return jm; 
		}

		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		// does not work for 6 or 9
		template<typename PairType>
		PairType calcJmValueAux(const HilbertState& ket) const
		{
			SizeType site0=0;
			SizeType site1=0;

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
			if (n!=16) throw PsimagLite::RuntimeError("reinterpret (unimplemented case): "
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

		int nini(const HilbertState& state,SizeType i,SizeType gamma1,SizeType spin1,SizeType gamma2,SizeType spin2) const
		{
			int tmp1 = HilbertSpaceFeAsType::calcNofElectrons(state,i,gamma1+spin1*modelParameters_.orbitals);
			int tmp2 = HilbertSpaceFeAsType::calcNofElectrons(state,i,gamma2+spin2*modelParameters_.orbitals);
			return tmp1*tmp2;
		}

		void addInteraction(SparseMatrixType &hmatrix,
		                    const VectorOperatorType& cm,
		                    SizeType i,
		                    RealType factorForDiagonals) const
		{
			addInteractionU1(hmatrix,cm,i,factorForDiagonals);
			addInteractionU2(hmatrix,cm,i,factorForDiagonals);
			if (!modelParameters_.decay) {
				addInteractionJ1(hmatrix,cm,i,factorForDiagonals);
				addInteractionJ2(hmatrix,cm,i,factorForDiagonals);
			} else {
				addInteractionV(hmatrix,cm,i,factorForDiagonals);
			}
		}

		RealType findHubbardU(SizeType index, SizeType orb1, SizeType orb2) const
		{
			if (!modelParameters_.decay) {
				assert(index < modelParameters_.hubbardU.size());
				return modelParameters_.hubbardU[index];
			}

			assert(orb1 + orb2*modelParameters_.orbitals < modelParameters_.hubbardU.size());
			return modelParameters_.hubbardU[orb1 + orb2*modelParameters_.orbitals];
		}

		//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
		void addInteractionU1(SparseMatrixType &hmatrix,
		                      const VectorOperatorType& cm,
		                      SizeType i,
		                      RealType factorForDiagonals) const
		{
			int dof=2*modelParameters_.orbitals;
			SparseMatrixType tmpMatrix,tmpMatrix2;

			for (SizeType alpha=0;alpha<SizeType(modelParameters_.orbitals);alpha++) {
				SparseMatrixType m1=cm[alpha+SPIN_UP*modelParameters_.orbitals+i*dof].data;
				SparseMatrixType m2=cm[alpha+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;

				multiply(tmpMatrix,n(m1),n(m2));
				multiplyScalar(tmpMatrix2,
				               tmpMatrix,
				               factorForDiagonals*findHubbardU(0,alpha,alpha));
				hmatrix += tmpMatrix2;
			}
		}

		//! Term is U[1] n_{i BAND0 } n_{i BAND1}
		void addInteractionU2(SparseMatrixType &hmatrix,
		                      const VectorOperatorType& cm,
		                      SizeType i,
		                      RealType factorForDiagonals) const
		{

			for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
				for (SizeType orb2=orb1+1;orb2<modelParameters_.orbitals;orb2++) {
				SparseMatrixType tmpMatrix,tmpMatrix2;
				multiply(tmpMatrix, nSummedOverSpin(cm,i,orb1),nSummedOverSpin(cm,i,orb2));
				multiplyScalar(tmpMatrix2,
				               tmpMatrix,
				               factorForDiagonals*findHubbardU(1,orb1,orb2));
				hmatrix += tmpMatrix2;
				}
			}
		}

		//! Term is U[2] S_{i BAND0 } S_{i BAND1}
		void addInteractionJ1(SparseMatrixType &hmatrix,
		                      const VectorOperatorType& cm,
		                      SizeType i,
		                      RealType factorForDiagonals) const
		{
			RealType val=0;
			RealType val2=2.0;
			RealType val3=4.0;

			for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
				for (SizeType orb2=orb1+1;orb2<modelParameters_.orbitals;orb2++) {
					SparseMatrixType tmpMatrix2,tmpMatrix;


					multiply(tmpMatrix, spinOperator(cm,i,orb1,0),spinOperator(cm,i,orb2,1));
					val = modelParameters_.hubbardU[2]/val2;
					multiplyScalar(tmpMatrix2,tmpMatrix,factorForDiagonals*val);// this is -2*J
					hmatrix += tmpMatrix2;

					multiply(tmpMatrix, spinOperator(cm,i,orb1,1),spinOperator(cm,i,orb2,0));
					val = modelParameters_.hubbardU[2]/val2;
					multiplyScalar(tmpMatrix2,tmpMatrix,factorForDiagonals*val);// this is -2*J
					hmatrix += tmpMatrix2;

					multiply(tmpMatrix, spinOperator(cm,i,orb1,2),spinOperator(cm,i,orb2,2));
					val = modelParameters_.hubbardU[2]/val3;
					multiplyScalar(tmpMatrix2,tmpMatrix,factorForDiagonals*val);// this is -2*J
					hmatrix += tmpMatrix2;
				}
			}
		}

		//! Term is U[3] \sum_{\alpha}\bar{n}_{i\alpha UP} \bar{n}_{i\alpha DOWN} 
		//! where \bar{n}_{i\alpha \spin} = c^\dagger_{i\alpha\spin} c_{i\bar{\alpha}\bar{spin}}
		void addInteractionJ2(SparseMatrixType &hmatrix,
		                      const VectorOperatorType& cm,
		                      SizeType i,
		                      RealType factorForDiagonals) const
		{
			SparseMatrixType tmpMatrix,tmpMatrix2;

			for (SizeType orb1=0;orb1<modelParameters_.orbitals;orb1++) {
				for (SizeType orb2=0;orb2<modelParameters_.orbitals;orb2++) {
					if (orb1==orb2) continue;
					multiply(tmpMatrix,nBar(cm,i,orb1,orb2,SPIN_UP),nBar(cm,i,orb1,orb2,SPIN_DOWN));
					multiplyScalar(tmpMatrix2,
					               tmpMatrix,
					               factorForDiagonals*modelParameters_.hubbardU[3]); // this is -J
					hmatrix += tmpMatrix2;
				}
			}
		}

		void addInteractionV(SparseMatrixType& hmatrix,
		                     const VectorOperatorType& cm,
		                     SizeType site,
	                         RealType factorForDiagonals) const
		{
			assert(modelParameters_.orbitals == 3);

			SizeType dofs = 2 * modelParameters_.orbitals;

			SparseMatrixType tmpMatrix1;
			SparseMatrixType tmpMatrix2;
			SparseMatrixType tmpMatrix;

			RealType value = modelParameters_.coulombV;

			for (SizeType spin1 = 0; spin1 < 2; ++spin1) {
				for (SizeType spin2 = 0; spin2 < 2; ++spin2) {
					if (spin1 == spin2) continue;

					multiply(tmpMatrix1,
					         cm[2+spin1*modelParameters_.orbitals+site*dofs].data,
					         cm[0+spin2*modelParameters_.orbitals+site*dofs].data);


					SparseMatrixType cmTranspose1;
					transposeConjugate(cmTranspose1,
					                   cm[1+spin2*modelParameters_.orbitals+site*dofs].data);

					SparseMatrixType cmTranspose2;
					transposeConjugate(cmTranspose2,
					                   cm[1+spin1*modelParameters_.orbitals+site*dofs].data);

					multiply(tmpMatrix2,cmTranspose1,cmTranspose2);


					multiply(tmpMatrix,tmpMatrix1,tmpMatrix2);
					multiplyScalar(tmpMatrix1,
					               tmpMatrix,
					               factorForDiagonals*value);
					hmatrix += tmpMatrix1;

					transposeConjugate(tmpMatrix2,tmpMatrix1);
					hmatrix += tmpMatrix2;
				}
			}
		}

		void addMagneticField(SparseMatrixType &hmatrix,
		                      const VectorOperatorType& cm,
		                      SizeType i,
		                      SizeType actualIndexOfSite,
		                      RealType factorForDiagonals) const
		{
			if (modelParameters_.magneticField.n_row()<3) return;
			assert(actualIndexOfSite<modelParameters_.magneticField.n_col());
			for (SizeType orb=0;orb<modelParameters_.orbitals;orb++)
				addMagneticField(hmatrix,cm,i,actualIndexOfSite,orb,factorForDiagonals);
		}

		void addMagneticField(SparseMatrixType &hmatrix,
		                      const VectorOperatorType& cm,
		                      SizeType i,
		                      SizeType actualIndexOfSite,
		                      SizeType orbital,
		                      RealType factorForDiagonals) const
		{
			int dof=2*modelParameters_.orbitals;
			SparseMatrixType cup = cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data;
			SparseMatrixType cupTranspose;
			transposeConjugate(cupTranspose,cup);
			SparseMatrixType cdown = cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data;
			SparseMatrixType A = cupTranspose * cdown;
			SparseMatrixType Atranspose;
			transposeConjugate(Atranspose,A);

			hmatrix += factorForDiagonals * modelParameters_.magneticField(0,actualIndexOfSite) * A;

			hmatrix += factorForDiagonals * modelParameters_.magneticField(1,actualIndexOfSite) * Atranspose;

			SparseMatrixType nup = n(cup);
			SparseMatrixType ndown = n(cdown);

			SparseMatrixType tmp = nup;
			tmp += (-1.0)*ndown;
			hmatrix += factorForDiagonals * modelParameters_.magneticField(2,actualIndexOfSite) * tmp;

		}

		void addPotentialV(SparseMatrixType &hmatrix,
		                   const VectorOperatorType& cm,
		                   SizeType i,
		                   SizeType actualIndexOfSite,
		                   RealType factorForDiagonals,
				   const typename PsimagLite::Vector<RealType>::Type& V) const
		{
			for (SizeType orb=0;orb<modelParameters_.orbitals;orb++)
			  addPotentialV(hmatrix,cm,i,actualIndexOfSite,orb,factorForDiagonals,V);
		}

		void addPotentialV(SparseMatrixType &hmatrix,
		                   const VectorOperatorType& cm,
		                   SizeType i,
		                   SizeType actualIndexOfSite,
		                   SizeType orbital,
		                   RealType factorForDiagonals,
				   const typename PsimagLite::Vector<RealType>::Type& V) const
		{
			int dof=2*modelParameters_.orbitals;
			SparseMatrixType nup = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dof].data);
			SparseMatrixType ndown = n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dof].data);

			SizeType linSize = geometry_.numberOfSites();

			SizeType iUp = actualIndexOfSite + (orbital + 0*modelParameters_.orbitals)*linSize;
			hmatrix += factorForDiagonals * V[iUp] * nup;
			SizeType iDown = actualIndexOfSite + (orbital + 1*modelParameters_.orbitals)*linSize;
			hmatrix += factorForDiagonals * V[iDown] * ndown;
		}

		SparseMatrixType n(const SparseMatrixType& c) const
		{
			SparseMatrixType tmpMatrix;
			SparseMatrixType cdagger;
			transposeConjugate(cdagger,c);
			multiply(tmpMatrix,c,cdagger);

			return tmpMatrix;
		}

		SparseMatrixType nBar(const typename PsimagLite::Vector<OperatorType>::Type& cm,SizeType i,SizeType orb1,SizeType orb2,SizeType spin) const
		{
			SizeType dofs = 2 * modelParameters_.orbitals;
			SparseMatrixType tmpMatrix,cdagger=cm[orb1+spin*modelParameters_.orbitals+i*dofs].data;
			SparseMatrixType cbar;
			transposeConjugate(cbar,cm[orb2+(1-spin)*modelParameters_.orbitals+i*dofs].data);
			multiply(tmpMatrix,cdagger,cbar);
			return tmpMatrix;
		}

		SparseMatrixType nSummedOverSpin(const typename PsimagLite::Vector<OperatorType>::Type& cm,SizeType i,SizeType orbital) const
		{
			SizeType dofs = 2 * modelParameters_.orbitals;
			SparseMatrixType tmpMatrix = n(cm[orbital+SPIN_UP*modelParameters_.orbitals+i*dofs].data);
			tmpMatrix += n(cm[orbital+SPIN_DOWN*modelParameters_.orbitals+i*dofs].data);
			return tmpMatrix;
		}

		SparseMatrixType spinOperator(const typename PsimagLite::Vector<OperatorType>::Type& cm,SizeType i,SizeType orbital,SizeType component) const
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

		SparseMatrixType spinOperatorAux(const typename PsimagLite::Vector<OperatorType>::Type& cm,SizeType i,SizeType orbital,SizeType spin1,SizeType spin2) const
		{
			SizeType dofs = 2 * modelParameters_.orbitals;
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

		void diagTest(const SparseMatrixType& fullm,const PsimagLite::String& str) const
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

} // namespace Dmrg
/*@}*/
#endif
