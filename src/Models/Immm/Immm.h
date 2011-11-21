// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009-2011, UT-Battelle, LLC
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

/*! \file Immm.h
 *
 *  DOC TBW FIXME
 *
 */
#ifndef IMMM_HEADER_H
#define IMMM_HEADER_H
#include "ModelBase.h"
#include "ParametersImmm.h"
#include "HilbertSpaceImmm.h"
#include "CrsMatrix.h"
#include "SpinSquaredHelper.h"
#include "SpinSquared.h"
#include "VerySparseMatrix.h"
#include "LinkProductImmm.h"
#include "ProgramGlobals.h"
#include <cassert>

namespace Dmrg {
	template<
		typename ModelHelperType_,
		typename SparseMatrixType,
		typename GeometryType,
  		template<typename> class SharedMemoryTemplate>
	class Immm : public ModelBase<ModelHelperType_,SparseMatrixType,GeometryType,
 	LinkProductImmm<ModelHelperType_>,SharedMemoryTemplate> {
		
		typedef unsigned int long long WordType;

	public:
		typedef ModelHelperType_ ModelHelperType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef PsimagLite::Matrix<SparseElementType> MatrixType;
		typedef typename ModelHelperType::BlockType Block;
		typedef typename ModelHelperType::ReflectionSymmetryType ReflectionSymmetryType;
		typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
		typedef LinkProductImmm<ModelHelperType> LinkProductType;
		typedef  HilbertSpaceImmm<WordType> HilbertSpaceImmmType;
		typedef typename HilbertSpaceImmmType::HilbertState HilbertState;
		typedef std::vector<HilbertState> HilbertBasisType;
		typedef   ModelBase<ModelHelperType,SparseMatrixType,GeometryType,LinkProductType,SharedMemoryTemplate> ModelBaseType;
		typedef	 typename ModelBaseType::MyBasis MyBasis;
		typedef	 typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;

		static int const maxNumberOfSites=ProgramGlobals::MaxNumberOfSites;;
		static const int FERMION_SIGN = -1;
		static const int SPIN_UP=HilbertSpaceImmmType::SPIN_UP;
		static const int SPIN_DOWN=HilbertSpaceImmmType::SPIN_DOWN;
		static const int NUMBER_OF_SPINS=HilbertSpaceImmmType::NUMBER_OF_SPINS;

		Immm(ParametersImmm<RealType> const &mp,GeometryType const &geometry) 
		: ModelBaseType(geometry),
		  modelParameters_(mp),
		  geometry_(geometry),
		  degreesOfFreedom_(geometry_.numberOfSites()),
		  index2Op_(geometry_.numberOfSites()),
		  hilbertSpace_(degreesOfFreedom_)
		{
			size_t cooperEach = 4;
			buildDofs(cooperEach);
			buildIndex2Op();
		}
		
		size_t dOf(size_t site) const { return degreesOfFreedom_[site]; }

		size_t hilbertSize(size_t site) const
		{
			return (1<<dOf(site));
		} 

		void print(std::ostream& os) const { operator<<(os,modelParameters_); }

		//! find creation operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		void setNaturalBasis(std::vector<OperatorType>& creationMatrix,
		                     SparseMatrixType& hamiltonian,
		                     BasisDataType& q,
		                     const Block& block)  const
		{
			assert(block.size()==1);
			std::vector<HilbertState> natBasis;
			std::vector<size_t> qvector;
			setNaturalBasis(natBasis,qvector,block);			

			setOperatorMatrices(creationMatrix,block);

			//! Set symmetry related
			setSymmetryRelated(q,natBasis,block.size());

			//! set hamiltonian
			calcHamiltonian(hamiltonian,creationMatrix,block);
		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(std::vector<OperatorType>& creationMatrix,
		                         const Block& block) const
		{
			assert(block.size()==1);
			std::vector<HilbertState> natBasis;
			SparseMatrixType tmpMatrix;
			std::vector<size_t> qvector;
			setNaturalBasis(natBasis,qvector,block);

			//! Set the operators c^\daggger_{i\gamma\sigma} in the natural basis
			creationMatrix.clear();
			size_t dof = dOf(block[0]);

			for (size_t sigma=0;sigma<dof;sigma++) {
				findOperatorMatrices(tmpMatrix,block[0],sigma,natBasis);
				typename OperatorType::Su2RelatedType su2related;
					
				OperatorType myOp(tmpMatrix,-1,typename OperatorType::PairType(0,0),1,su2related);
				creationMatrix.push_back(myOp);
			}
		}

		MatrixType getOperator(const std::string& what,size_t orbital=0,size_t spin=0) const
		{
			throw std::runtime_error("Need site\n");
			Block block;
			block.resize(1);
			block[0]=0;
			size_t norb = dOf(block[0])/HilbertSpaceImmmType::NUMBER_OF_SPINS;
			std::vector<OperatorType> creationMatrix;
			setOperatorMatrices(creationMatrix,block);

			if (what=="+" or what=="i") {
				return cDaggerCi(block,SPIN_UP,SPIN_DOWN);
			}
			if (what=="-") {
				return cDaggerCi(block,SPIN_DOWN,SPIN_UP);
			}
			if (what=="z") {
				return nUpOrDown(block,SPIN_UP)-nUpOrDown(block,SPIN_DOWN);
			}
			if (what=="n") {
				return nUpOrDown(block,SPIN_UP)+nUpOrDown(block,SPIN_DOWN);
			} 
			if (what=="c") {
				MatrixType tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[index2Op_[block[0]] + orbital + spin*norb].data);
				return tmp;
			}
			if (what=="d") { // delta = c^\dagger * c^dagger
				throw std::runtime_error("delta unimplemented for this model\n");
			}
			std::cerr<<"what="<<what<<"\n";
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}

		//! find all states in the natural basis for a block of n sites
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(std::vector<HilbertState>  &basis,
		                     std::vector<size_t>& q,
		                     const std::vector<size_t>& block) const
		{
			assert(block.size()==1);
			size_t dof =  dOf(block[0]);
			HilbertState total = (1<<dof);

			std::vector<HilbertState>  basisTmp;
			for (HilbertState a=0;a<total;a++) basisTmp.push_back(a);

			// reorder the natural basis (needed for MULTIPLE BANDS)
			findQuantumNumbers(q,basisTmp,block[0]);
			std::vector<size_t> iperm(q.size());
			Sort<std::vector<size_t> > sort;
			sort.sort(q,iperm);
			basis.clear();
			for (HilbertState a=0;a<total;a++)
				basis.push_back(basisTmp[iperm[a]]);
		}

		void findElectrons(std::vector<size_t>& electrons,
		                   const std::vector<HilbertState>& basis,
		                   size_t site) const
		{
			electrons.resize(basis.size());
			for (size_t i=0;i<basis.size();i++) {
				// nup
				size_t nup = hilbertSpace_.electronsWithGivenSpin(basis[i],site,HilbertSpaceImmmType::SPIN_UP);
				// ndown
				size_t ndown = hilbertSpace_.electronsWithGivenSpin(basis[i],site,HilbertSpaceImmmType::SPIN_DOWN);
				electrons[i] = nup + ndown;
			}
		}

	private:


		//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS 
		RealType sign(HilbertState const &ket, size_t site,size_t sigma) const
		{
			int value=0;
			size_t dof = dOf(site);
			for (size_t alpha=0;alpha<dof;alpha++) 
				value += hilbertSpace_.calcNofElectrons(ket,0,site,alpha);
			// add electron on site 0 if needed
			if (site>0) value += hilbertSpace_.electronsAtGivenSite(ket,0);

			//order for sign is: up a (sigma==0), down a (sigma==2), up b (sigma==1), down b(sigma==3)
			unsigned int x = hilbertSpace_.get(ket,site);	
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
		void findOperatorMatrices(SparseMatrixType& creationMatrix,
		                          size_t site,
		                          size_t sigma,
		                          const std::vector<HilbertState>& natBasis) const
		{
			HilbertState bra,ket;
			size_t n = natBasis.size();
			MatrixType cm(n,n);

			for (size_t ii=0;ii<n;ii++) {
				bra=ket=natBasis[ii];
				
				if (hilbertSpace_.isNonZero(ket,0,sigma)) {
					
				} else {
					hilbertSpace_.create(bra,0,sigma);
					int jj = PsimagLite::isInVector(natBasis,bra);
					assert(jj>=0);
					assert(ii!=size_t(jj));
					cm(ii,jj) =sign(ket,site,sigma);
				}
			}
			// here reinterpret for SU(2) if needed

			SparseMatrixType temp;
			fullMatrixToCrsMatrix(temp,cm);
			transposeConjugate(creationMatrix,temp);
		}

		void findQuantumNumbers(std::vector<size_t>& q,
		                        const std::vector<HilbertState>& basis,
		                        size_t site) const
		{
			BasisDataType qq;
			setSymmetryRelated(qq,basis,site);
			MyBasis::findQuantumNumbers(q,qq);
		}

		void setSymmetryRelated(BasisDataType& q,
		                        const std::vector<HilbertState>& basis,
		                        size_t site) const
		{

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

// 				size_t na = hilbertSpace_.calcNofElectrons(basis[i],0,0) +
// 						hilbertSpace_.calcNofElectrons(basis[i],0,0+2);
// 				size_t nb = hilbertSpace_.calcNofElectrons(basis[i],0,1) +
// 						hilbertSpace_.calcNofElectrons(basis[i],0,1+2);

				size_t flavor = 0; // na  + 3*nb;

				flavors.push_back(flavor);
				jmSaved = jmpair;

				// nup
				electronsUp[i] = hilbertSpace_.electronsWithGivenSpin(basis[i],site,HilbertSpaceImmmType::SPIN_UP);
				// ndown
				electronsDown[i] = hilbertSpace_.electronsWithGivenSpin(basis[i],site,HilbertSpaceImmmType::SPIN_DOWN);
			}
			q.jmValues=jmvalues;
			q.flavors = flavors;
			q.electronsUp = electronsUp;
			q.electronsDown = electronsDown;
		}

		//! Not implemented, su(2) symmetry won't work
		template<typename PairType>
		PairType calcJmvalue(const HilbertState& ket) const
		{
			PairType jm(0,0);
			return jm; 
		}

		//! Full hamiltonian from creation matrices cm
		void calcHamiltonian(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,Block const &block) const
		{
			assert(block.size()==1);
			SparseMatrixType tmpMatrix,tmpMatrix2;

			hmatrix.makeDiagonal(cm[0].data.rank());
			//addInteraction(hmatrix,cm,block);
		}

// 		void addInteraction(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
// 		{
// 			addInteractionU1(hmatrix,cm,i);
// 			addInteractionU2(hmatrix,cm,i);
// 			addInteractionJ1(hmatrix,cm,i);
// 			addInteractionJ2(hmatrix,cm,i);
// 		}
// 
// 		//! Term is U[0]\sum_{\alpha}n_{i\alpha UP} n_{i\alpha DOWN}
// 		void addInteractionU1(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
// 		{
// 			int dof=DEGREES_OF_FREEDOM;
// 			SparseMatrixType tmpMatrix,tmpMatrix2;
// 
// 			for (size_t alpha=0;alpha<size_t(NUMBER_OF_ORBITALS);alpha++) {
// 				SparseMatrixType m1=cm[alpha+SPIN_UP*NUMBER_OF_ORBITALS+i*dof].data;
// 				SparseMatrixType m2=cm[alpha+SPIN_DOWN*NUMBER_OF_ORBITALS+i*dof].data;
// 
// 				multiply(tmpMatrix,n(m1),n(m2));
// 				multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[0]); // this is U
// 				hmatrix += tmpMatrix2;
// 			}
// 		}
// 
// 		//! Term is U[1] n_{i BAND0 } n_{i BAND1}
// 		void addInteractionU2(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
// 		{
// 			size_t orbital0=0,orbital1=1;
// 			SparseMatrixType tmpMatrix,tmpMatrix2;
// 
// 			multiply(tmpMatrix, nSummedOverSpin(cm,i,orbital0),nSummedOverSpin(cm,i,orbital1));
// 			multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[1]);// this is U'-J/2
// 			hmatrix += tmpMatrix2;
// 		}
// 
// 		//! Term is U[2] S_{i BAND0 } S_{i BAND1}
// 		void addInteractionJ1(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
// 		{
// 			size_t orbital0=0,orbital1=1;
// 			SparseMatrixType tmpMatrix2,tmpMatrix;
// 			RealType val=0;
// 			RealType val2=2.0;
// 			RealType val3=4.0;
// 
// 			multiply(tmpMatrix, spinOperator(cm,i,orbital0,0),spinOperator(cm,i,orbital1,1));
// 			val = modelParameters_.hubbardU[2]/val2;
// 			multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
// 			hmatrix += tmpMatrix2;
// 
// 			multiply(tmpMatrix, spinOperator(cm,i,orbital0,1),spinOperator(cm,i,orbital1,0));
// 			val = modelParameters_.hubbardU[2]/val2;
// 			multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
// 			hmatrix += tmpMatrix2;
// 
// 			multiply(tmpMatrix, spinOperator(cm,i,orbital0,2),spinOperator(cm,i,orbital1,2));
// 			val = modelParameters_.hubbardU[2]/val3;
// 			multiplyScalar(tmpMatrix2,tmpMatrix,val);// this is -2*J
// 			hmatrix += tmpMatrix2;
// 			return;
// 		}
// 
// 		//! Term is U[3] \sum_{\alpha}\bar{n}_{i\alpha UP} \bar{n}_{i\alpha DOWN} 
// 		//! where \bar{n}_{i\alpha \spin} = c^\dagger_{i\alpha\spin} c_{i\bar{\alpha}\bar{spin}}
// 		void addInteractionJ2(SparseMatrixType &hmatrix,const std::vector<OperatorType>& cm,size_t i) const
// 		{
// 			SparseMatrixType tmpMatrix,tmpMatrix2;
// 
// 			for (size_t alpha=0;alpha<size_t(NUMBER_OF_ORBITALS);alpha++) {
// 				multiply(tmpMatrix,nBar(cm,i,alpha,SPIN_UP),nBar(cm,i,alpha,SPIN_DOWN));
// 				multiplyScalar(tmpMatrix2,tmpMatrix,modelParameters_.hubbardU[3]); // this is -J
// 				hmatrix += tmpMatrix2;
// 			}
// 		}
// 
// 		SparseMatrixType n(const SparseMatrixType& c) const
// 		{
// 			SparseMatrixType tmpMatrix;
// 			SparseMatrixType cdagger;
// 			transposeConjugate(cdagger,c);
// 			multiply(tmpMatrix,c,cdagger);
// 
// 			return tmpMatrix;
// 		}
// 
// 		SparseMatrixType nBar(const std::vector<OperatorType>& cm,size_t i,size_t orbital,size_t spin) const
// 		{
// 			SparseMatrixType tmpMatrix,cdagger=cm[orbital+spin*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data;
// 			SparseMatrixType cbar;
// 			transposeConjugate(cbar,cm[1-orbital+(1-spin)*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
// 			multiply(tmpMatrix,cdagger,cbar);
// 			return tmpMatrix;
// 		}
// 
// 		SparseMatrixType nSummedOverSpin(const std::vector<OperatorType>& cm,size_t i,size_t orbital) const
// 		{
// 			SparseMatrixType tmpMatrix = n(cm[orbital+SPIN_UP*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
// 			tmpMatrix += n(cm[orbital+SPIN_DOWN*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
// 			return tmpMatrix;
// 		}
// 
// 		SparseMatrixType spinOperator(const std::vector<OperatorType>& cm,size_t i,size_t orbital,size_t component) const
// 		{
// 			switch (component) {
// 				case 0: // S^+
// 					return spinOperatorAux(cm,i,orbital,SPIN_UP,SPIN_DOWN);
// 					break;
// 				case 1: // S^-
// 					return spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_UP);
// 					break;
// 			}
// 			SparseMatrixType tmpMatrix=spinOperatorAux(cm,i,orbital,SPIN_UP,SPIN_UP);
// 			SparseMatrixType tmpMatrix2=spinOperatorAux(cm,i,orbital,SPIN_DOWN,SPIN_DOWN);
// 			SparseMatrixType tmpMatrix3;
// 			multiplyScalar(tmpMatrix3,tmpMatrix2,-1.0);
// 			tmpMatrix += tmpMatrix3;
// 			return tmpMatrix;
// 		}
// 
// 		SparseMatrixType spinOperatorAux(const std::vector<OperatorType>& cm,size_t i,size_t orbital,size_t spin1,size_t spin2) const
// 		{
// 			SparseMatrixType result,temp;
// 			transposeConjugate(temp,cm[orbital+spin2*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data);
// 			multiply(
// 				result, // =
// 				cm[orbital+spin1*NUMBER_OF_ORBITALS+i*DEGREES_OF_FREEDOM].data, // times
// 				temp
// 			);
// 					
// 			return result;
// 		}

		void diagTest(const SparseMatrixType& fullm,const std::string& str) const
		{
			if (fullm.rank()!=256) return;
			MatrixType fullm2;
			crsMatrixToFullMatrix(fullm2,fullm);
			std::vector<SparseElementType> eigs(fullm2.n_row());
			PsimagLite::diag(fullm2,eigs,'V');
			std::cout<<str<<" diagTest size="<<fullm.rank()<<" eigs[0]="<<eigs[0]<<"\n";
			std::cout<<fullm;
		}

		//! This is a mapping from site to the start of the operators for 
		//! that site. Note that there are more than one operator per site,
		//! and the number of operators per site varies.
		//! So, you have for example operators 0 1 2 3 belonging to site 0;
		//! 4, 5, 6, 7 to site 1; 8, 9 belonging to site 2, etc.
		//! That was just an example and might not reflect the actual numbers
		//! for this model
		void buildIndex2Op()
		{
			size_t counter = 0;
			for (size_t i=0;i<geometry_.numberOfSites();i++) {
				index2Op_[i] = counter;
				counter += dOf(i);
			}
		}
		
		MatrixType cDaggerCi(const std::vector<size_t>& block,
		                            size_t spin1,
		                            size_t spin2) const
		{
			assert(block.size()==1);
			size_t site = block[0];
			MatrixType tmp;
			size_t norb = dOf(site)/NUMBER_OF_SPINS;
			std::vector<OperatorType> creationMatrix;
			setOperatorMatrices(creationMatrix,block);
			for (size_t orb=0;orb<norb;orb++)
				tmp += multiplyTc(creationMatrix[orb+SPIN_UP*norb].data,
				                  creationMatrix[orb+SPIN_DOWN*norb].data);
			return tmp;
		}
		
		MatrixType nUpOrDown(const std::vector<size_t>& block,size_t spin) const
		{
			return cDaggerCi(block,spin,spin);
		}

		//! If there's only spin  at site i degreesOfFreedom_[i]=2
		//! If there's spin an norb orbitals then degreesOfFreedom_[i]=2*norb
		//! etc.
		void buildDofs(size_t copperEach)
		{
			size_t counter = 5;
			for (size_t i=0;i<degreesOfFreedom_.size();i++) {
				if (counter%copperEach==0) degreesOfFreedom_[i]=2;
				else degreesOfFreedom_[i]=4;
				counter++;
			}
		}

		const ParametersImmm<RealType>&  modelParameters_;
		GeometryType const &geometry_;
		std::vector<size_t> degreesOfFreedom_;
		std::vector<size_t> index2Op_;
		HilbertSpaceImmmType hilbertSpace_;
	};     //class Immm

	template<typename ModelHelperType,
	         typename SparseMatrixType,
	         typename GeometryType,
	         template<typename> class SharedMemoryTemplate>
	std::ostream &operator<<(std::ostream &os,
	                         const Immm<ModelHelperType,
	                                    SparseMatrixType,
	                                    GeometryType,
	                                    SharedMemoryTemplate>& model)
	{
		model.print(os);
		return os;
	}
} // namespace Dmrg
/*@}*/
#endif // IMMM_HEADER_H
