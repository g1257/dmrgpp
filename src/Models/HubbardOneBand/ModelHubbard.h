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

/*! \file ModelHubbard.h
 *
 *  An implementation of the Hubbard Model to use with the DmrgSolver
 *
 */
#ifndef MODEL_HUBBARD_DMRG
#define MODEL_HUBBARD_DMRG
#include <cassert>
#include "Sort.h" // in PsimagLite
#include "ModelBase.h"
#include "ParametersModelHubbard.h"
#include "HilbertSpaceHubbard.h"
#include "LinkProductHubbardOneBand.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "ProgramGlobals.h"

namespace Dmrg {
	//! Model Hubbard for DMRG solver, inherits from ModelBase and implements its interface:
	template<typename ModelHelperType_,
	typename SparseMatrixType,
	typename DmrgGeometryType,
	template<typename> class SharedMemoryTemplate>
	class ModelHubbard : public ModelBase<ModelHelperType_,SparseMatrixType,DmrgGeometryType,
		LinkProductHubbardOneBand<ModelHelperType_>,SharedMemoryTemplate> {
	public:
		typedef ModelHelperType_ ModelHelperType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename SparseMatrixType::value_type SparseElementType;
	private:

		static int const maxNumberOfSites=ProgramGlobals::MaxNumberOfSites;;
		static const int FERMION_SIGN = -1;
		static const int DEGREES_OF_FREEDOM=2;
		static const int NUMBER_OF_ORBITALS=OperatorsType::NUMBER_OF_ORBITALS;

		typedef unsigned int long long WordType;

		typedef typename ModelHelperType::BlockType Block;
		typedef typename ModelHelperType::ReflectionSymmetryType ReflectionSymmetryType;

	public:

		typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
		typedef  HilbertSpaceHubbard<WordType> HilbertSpaceHubbardType;
		typedef typename HilbertSpaceHubbardType::HilbertState HilbertState;
		typedef LinkProductHubbardOneBand<ModelHelperType> LinkProductType;
		typedef ModelBase<ModelHelperType,SparseMatrixType,DmrgGeometryType,LinkProductType,SharedMemoryTemplate> ModelBaseType;
		typedef	typename ModelBaseType::MyBasis MyBasis;
		typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef std::vector<HilbertState> HilbertBasisType;
		
		ModelHubbard(ParametersModelHubbard<RealType> const &mp,DmrgGeometryType const &dmrgGeometry,
			    size_t offset = DEGREES_OF_FREEDOM) 
			: ModelBaseType(dmrgGeometry),modelParameters_(mp), dmrgGeometry_(dmrgGeometry),offset_(offset),
					    spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM),
					   reinterpretX_(maxNumberOfSites),reinterpretY_(maxNumberOfSites)
		{
		}

		void print(std::ostream& os) const { operator<<(os,modelParameters_); }

		size_t orbitals() const { return NUMBER_OF_ORBITALS; }

		size_t hilbertSize(size_t site) const { return (size_t)pow(2,2*NUMBER_OF_ORBITALS); }

		//! find creation operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(std::vector<OperatorType> &creationMatrix,SparseMatrixType &hamiltonian,
				BasisDataType &q,Block const &block) const
		{
			HilbertBasisType natBasis;
			std::vector<size_t> quantumNumbs;
			setNaturalBasis(natBasis,quantumNumbs,block);

			setOperatorMatrices(creationMatrix,block);

			//! Set symmetry related
			setSymmetryRelated(q,natBasis,block.size());

			//! set hamiltonian
			calcHamiltonian(hamiltonian,creationMatrix,block);
		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(std::vector<OperatorType> &creationMatrix,Block const &block) const
		{
			std::vector<typename HilbertSpaceHubbardType::HilbertState> natBasis;
			SparseMatrixType tmpMatrix;
			std::vector<size_t> quantumNumbs;
			setNaturalBasis(natBasis,quantumNumbs,block);

			//! Set the operators c^\daggger_{i\sigma} in the natural basis
			creationMatrix.clear();
			for (size_t i=0;i<block.size();i++) {
				for (int sigma=0;sigma<DEGREES_OF_FREEDOM;sigma++) {
					tmpMatrix = findOperatorMatrices(i,sigma,natBasis);
					int asign= 1;
					if (sigma>0) asign= 1;
					typename OperatorType::Su2RelatedType su2related;
					if (sigma==0) {
						su2related.source.push_back(i*offset_);
						su2related.source.push_back(i*offset_+1);	
						su2related.transpose.push_back(-1);
						su2related.transpose.push_back(-1);
						su2related.offset = NUMBER_OF_ORBITALS;
					}
					OperatorType myOp(tmpMatrix,-1,typename OperatorType::PairType(1,1-sigma),asign,su2related);
					
					creationMatrix.push_back(myOp);
				}
			}
		}

		PsimagLite::Matrix<SparseElementType> getOperator(const std::string& what,size_t gamma=0,size_t spin=0) const
		{
			Block block;
			block.resize(1);
			block[0]=0;
			std::vector<OperatorType> creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="+" or what=="i") {
				PsimagLite::Matrix<SparseElementType> tmp = multiplyTc(creationMatrix[0].data,creationMatrix[1].data);
				return tmp;
			} else if (what=="-") {
				PsimagLite::Matrix<SparseElementType> tmp = multiplyTc(creationMatrix[1].data,creationMatrix[0].data);
				return tmp;
			} else if (what=="z") {
				PsimagLite::Matrix<SparseElementType> tmp =multiplyTc(creationMatrix[0].data,creationMatrix[0].data);
				PsimagLite::Matrix<SparseElementType> tmp2 =multiplyTc(creationMatrix[1].data,creationMatrix[1].data);
				return tmp-tmp2;
			} else if (what=="n") {
				PsimagLite::Matrix<SparseElementType> tmp =  multiplyTc(creationMatrix[0].data,creationMatrix[0].data)
						+ multiplyTc(creationMatrix[1].data,creationMatrix[1].data);
				return tmp;
			} else if (what=="c") {
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[spin].data);
				return tmp;
			} else if (what=="nup") {
				PsimagLite::Matrix<SparseElementType> cup = getOperator("c",0,0);
				PsimagLite::Matrix<SparseElementType> nup = multiplyTransposeConjugate(cup,cup);
				return nup;
			} else if (what=="ndown") {
				PsimagLite::Matrix<SparseElementType> cdown = getOperator("c",0,1);
				PsimagLite::Matrix<SparseElementType> ndown = multiplyTransposeConjugate(cdown,cdown);
				return ndown;
			}
			std::cerr<<"Argument: "<<what<<" "<<__FILE__<<"\n";
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}
		
		//! find total number of electrons for each state in the basis
		void findElectrons(std::vector<size_t> &electrons,
		                   const std::vector<HilbertState>& basis,
		                   size_t site) const
		{
			int nup,ndown;
			electrons.clear();
			for (size_t i=0;i<basis.size();i++) {
				nup = HilbertSpaceHubbardType::getNofDigits(basis[i],0);
				ndown = HilbertSpaceHubbardType::getNofDigits(basis[i],1);
				electrons.push_back(nup+ndown);
			}
		}

// 		//! find all states in the natural basis for a block of n sites
// 		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(HilbertBasisType  &basis,
		                     std::vector<size_t>& q,
		                     const std::vector<size_t>& block) const
		{
			assert(block.size()==1);
			HilbertState a=0;
			int sitesTimesDof=DEGREES_OF_FREEDOM;
			HilbertState total = (1<<sitesTimesDof);

			HilbertBasisType  basisTmp;
			for (a=0;a<total;a++) basisTmp.push_back(a);

			// reorder the natural basis (needed for MULTIPLE BANDS)
			findQuantumNumbers(q,basisTmp,1);
			std::vector<size_t> iperm(q.size());

			Sort<std::vector<size_t> > sort;
			sort.sort(q,iperm);
			basis.clear();
			for (a=0;a<total;a++) basis.push_back(basisTmp[iperm[a]]);
		}
		
		const DmrgGeometryType& geometry() const { return dmrgGeometry_; }

	private:
		const ParametersModelHubbard<RealType>&  modelParameters_;
		const DmrgGeometryType &dmrgGeometry_;
		size_t offset_;
		SpinSquaredHelper<RealType,WordType> spinSquaredHelper_;
		SpinSquared<SpinSquaredHelper<RealType,WordType> > spinSquared_;
		size_t reinterpretX_,reinterpretY_;

		//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
		RealType sign(typename HilbertSpaceHubbardType::HilbertState const &ket, int i,int sigma) const
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
			if (value==0 || value%2==0) return 1.0;

			return FERMION_SIGN;
		}

		//! Compute correct value of Hubbard U parameter for  site ind
		RealType computeHubbardUValue(int type,int ind,int smax,int emin) const 
		{
			//! There are two cases:
			//! 1. (ind,jnd) in SUX --> use input hoppings
			//! 2. (ind,jnd) in YUE --> use reflected hoppings

			RealType x=0;
			switch (type) {
				case ProgramGlobals::ENVIRON_ENVIRON:
					x=modelParameters_.hubbardU[dmrgGeometry_.findReflection(ind)];
					break;

				case ProgramGlobals::SYSTEM_SYSTEM:
					x=modelParameters_.hubbardU[ind];
					break;
			}
			return x;
		}

		//! Compute correct value of onsite potential  for  site ind and internal degree of freedom sigma
		RealType computeOnsitePotential(int type,int ind,int sigma,int smax,int emin) const 
		{
			//! There are two cases:
			//! 1. (ind,jnd) in SUX --> use input hoppings
			//! 2. (ind,jnd) in YUE --> use reflected hoppings
			int totalS=modelParameters_.linSize;
			RealType x=0;
			switch (type) {
				case ProgramGlobals::ENVIRON_ENVIRON:
					x=modelParameters_.potentialV[dmrgGeometry_.findReflection(ind)+sigma*totalS];
					break;

				case ProgramGlobals::SYSTEM_SYSTEM:
					x=modelParameters_.potentialV[ind+sigma*totalS];
					break;
			}
			return x;
		}

		//! Find c^\dagger_isigma in the natural basis natBasis
		SparseMatrixType findOperatorMatrices(int i,int sigma,std::vector<typename HilbertSpaceHubbardType::HilbertState> const &natBasis) const
		{
			typename HilbertSpaceHubbardType::HilbertState bra,ket;
			int n = natBasis.size();
			PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);
			
			for (size_t ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				if (HilbertSpaceHubbardType::isNonZero(ket,i,sigma)) {
					
				} else {
					HilbertSpaceHubbardType::create(bra,i,sigma);
					int jj = PsimagLite::isInVector(natBasis,bra);
					if (jj<0) throw std::runtime_error("findOperatorMatrices: internal error while"
							"creating.\n");
					cm(ii,jj) =sign(ket,i,sigma);
				}
			}

			SparseMatrixType creationMatrix(cm);
			return creationMatrix;
		}

		//! find quantum numbers for each state of this basis, 
		//! considered symmetries for this model are: n_up and n_down
// 		void findQuantumNumbers(std::vector<int> &q,std::vector<typename HilbertSpaceHubbardType::HilbertState>  const &basis) const
// 		{
// 			int nup,ndown;
// 			q.clear();
// 			for (size_t i=0;i<basis.size();i++) {
// 				nup = HilbertSpaceHubbardType::getNofDigits(basis[i],0);
// 				ndown = HilbertSpaceHubbardType::getNofDigits(basis[i],1);
// 				q.push_back(nup +maxNumberOfSites*ndown);
// 			}
// 		}

		void findQuantumNumbers(std::vector<size_t>& q,const HilbertBasisType  &basis,int n) const
		{
			BasisDataType qq;
			setSymmetryRelated(qq,basis,n);
			MyBasis::findQuantumNumbers(q,qq);
		}

		void setSymmetryRelated(BasisDataType& q,HilbertBasisType  const &basis,int n) const
		{
			if (n!=1) std::runtime_error("ModelHubbard::setSymmetryRelated() implemented for n=1 only\n");

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
				electronsUp[i] = HilbertSpaceHubbardType::getNofDigits(basis[i],HilbertSpaceHubbardType::SPIN_UP);
				// ndown
				electronsDown[i] = HilbertSpaceHubbardType::getNofDigits(basis[i],HilbertSpaceHubbardType::SPIN_DOWN);
				
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
		PairType calcJmvalue(const HilbertState& ket) const
		{
			PairType jm(0,0);
			
			size_t x=reinterpretX_,y=reinterpretY_; // these states need reinterpretation
			if (ket==x) {
				jm=std::pair<size_t,size_t>(2,1);
			} else if (ket==y) {
				jm=std::pair<size_t,size_t>(0,0);
			} else jm=calcJmValueAux<PairType>(ket);
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

		//! Full hamiltonian from creation matrices cm
		void calcHamiltonian(SparseMatrixType &hmatrix,std::vector<OperatorType> const &cm,Block const &block) const
		{
			size_t n=block.size();
			//int type,sigma;
			SparseMatrixType tmpMatrix,tmpMatrix2,niup,nidown;
			
			hmatrix.makeDiagonal(cm[0].data.rank());
			size_t linSize = dmrgGeometry_.numberOfSites();
			assert(block.size()==1);
			
			for (size_t i=0;i<n;i++) {
				//! hopping part
				for (size_t j=0;j<n;j++) {
					for (size_t term=0;term<dmrgGeometry_.terms();term++) {
						for (size_t dofs=0;dofs<LinkProductType::dofs(term);dofs++) {
							std::pair<size_t,size_t> edofs = LinkProductType::connectorDofs(term,dofs);
							RealType tmp = dmrgGeometry_(block[i],edofs.first,block[j],edofs.second,term);
						
							if (i==j || tmp==0.0) continue;
				
							size_t sigma = dofs;
							transposeConjugate(tmpMatrix2,cm[sigma+j*offset_].data);
							multiply(tmpMatrix,cm[sigma+i*offset_].data,tmpMatrix2);
							multiplyScalar(tmpMatrix2,tmpMatrix,static_cast<SparseElementType>(tmp));
							hmatrix += tmpMatrix2;
						}
					}
				}
				// onsite U hubbard 
				//n_i up
				size_t sigma =0; // up sector
				transposeConjugate(tmpMatrix,cm[sigma+i*offset_].data);
				multiply(niup,tmpMatrix,cm[sigma+i*offset_].data);
				//n_i down
				sigma =1; // down sector
				transposeConjugate(tmpMatrix,cm[sigma+i*offset_].data);
				multiply(nidown,tmpMatrix,cm[sigma+i*offset_].data);
				 
				multiply(tmpMatrix,niup,nidown);
				//type = dmrgGeometry_.calcConnectorType(block[i],block[i]);
				RealType tmp = modelParameters_.hubbardU[block[i]]; //computeHubbardUValue(type,block[i],smax,emin);
				multiplyScalar(tmpMatrix2,tmpMatrix,static_cast<SparseElementType>(tmp));

				hmatrix += tmpMatrix2;

				// V_iup term
				tmp = modelParameters_.potentialV[block[i]+0*linSize]; //computeOnsitePotential(type,block[i],0,smax,emin);
				multiplyScalar(tmpMatrix,niup,static_cast<SparseElementType>(tmp));
				hmatrix += tmpMatrix;

				// V_idown term
				tmp = modelParameters_.potentialV[block[i]+1*linSize]; //computeOnsitePotential(type,block[i],1,smax,emin);
				multiplyScalar(tmpMatrix,nidown,static_cast<SparseElementType>(tmp));
				hmatrix += tmpMatrix;
			}
		}
	};	//class ModelHubbard

	template<typename ModelHelperType,
	typename SparseMatrixType,
	typename DmrgGeometryType,
	template<typename> class SharedMemoryTemplate>
	std::ostream &operator<<(std::ostream &os,const ModelHubbard<ModelHelperType,SparseMatrixType,DmrgGeometryType,SharedMemoryTemplate>& model)
	{
		model.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif
