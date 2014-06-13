/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file Tj1Orb.h
 *
 *  Hubbard + Heisenberg
 *
 */
#ifndef TJ_1ORB_H
#define TJ_1ORB_H
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "../Models/HeisenbergSpinOneHalf/ModelHeisenberg.h"
#include "../Models/Tj1Orb/LinkProductTj1Orb.h"
#include "../Models/Tj1Orb/ParametersModelTj1Orb.h"
#include "ModelCommon.h"

namespace Dmrg {
	//! t-J model for DMRG solver, uses ModelHubbard and ModelHeisenberg by containment
	template<typename ModelBaseType>
	class Tj1Orb : public ModelBaseType {

	public:

		typedef ModelHubbard<ModelBaseType> ModelHubbardType;
		typedef ModelHeisenberg<ModelBaseType> ModelHeisenbergType;
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
		typedef LinkProductTj1Orb<ModelHelperType> LinkProductType;
		typedef ModelCommon<ModelBaseType,LinkProductType> ModelCommonType;
		typedef	typename ModelBaseType::MyBasis MyBasis;
		typedef	typename ModelBaseType::BasisWithOperatorsType MyBasisWithOperators;
		typedef typename MyBasis::BasisDataType BasisDataType;
		typedef typename ModelHubbardType::HilbertState HilbertStateType;
		typedef typename ModelHubbardType::HilbertBasisType HilbertBasisType;
		typedef typename ModelHelperType::BlockType BlockType;
		typedef typename ModelBaseType::SolverParamsType SolverParamsType;
		typedef typename ModelBaseType::VectorType VectorType;
		typedef typename ModelHubbardType::HilbertSpaceHubbardType HilbertSpaceHubbardType;
		typedef typename ModelBaseType::InputValidatorType InputValidatorType;

		static const int DEGREES_OF_FREEDOM=2;
		static const int NUMBER_OF_ORBITALS = 1;
		static const int FERMION_SIGN = -1;

		enum {SPIN_UP, SPIN_DOWN};

		Tj1Orb(const SolverParamsType& solverParams,
		       InputValidatorType& io,
		       GeometryType const &geometry)
		    : ModelBaseType(io,new ModelCommonType(solverParams,geometry)),
		      modelParameters_(io),
		      geometry_(geometry),
		      offset_(DEGREES_OF_FREEDOM+3), // c^\dagger_up, c^\dagger_down, S+, Sz, n
		      spinSquared_(spinSquaredHelper_,NUMBER_OF_ORBITALS,DEGREES_OF_FREEDOM)
		{}

		SizeType memResolv(PsimagLite::MemResolv& mres,
		                   SizeType x,
		                   PsimagLite::String msg = "") const
		{
			return 0;
		}

		SizeType hilbertSize(SizeType site) const
		{
			return 3;
		}

		//! find creation operator matrices for (i,sigma) in the natural basis, find quantum numbers and number of electrons
		//! for each state in the basis
		virtual void setNaturalBasis(VectorOperatorType& creationMatrix,
		                             SparseMatrixType &hamiltonian,
		                             BasisDataType& q,
		                             const BlockType& block,
		                             const RealType& time) const
		{

			HilbertBasisType natBasis;
			typename PsimagLite::Vector<SizeType>::Type quantumNumbs;
			setNaturalBasis(natBasis,quantumNumbs,block);

			setOperatorMatrices(creationMatrix,block);

			//! Set symmetry related
			setSymmetryRelated(q,natBasis,block.size());

			//! set hamiltonian
			this->calcHamiltonian(hamiltonian,creationMatrix,block,time);
		}

		//! set creation matrices for sites in block
		void setOperatorMatrices(typename PsimagLite::Vector<OperatorType> ::Type&creationMatrix,
		                         const BlockType& block) const
		{
			typename PsimagLite::Vector<HilbertStateType>::Type natBasis;
			SparseMatrixType tmpMatrix;
			typename PsimagLite::Vector<SizeType>::Type quantumNumbs;
			setNaturalBasis(natBasis,quantumNumbs,block);

			// Set the operators c^\daggger_{i\sigma} in the natural basis
			creationMatrix.clear();
			for (SizeType i=0;i<block.size();i++) {
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
				creationMatrix.push_back(myOp);

				// Set the operators S^z_i in the natural basis
				tmpMatrix = findSzMatrices(i,natBasis);
				typename OperatorType::Su2RelatedType su2related2;
				OperatorType myOp2(tmpMatrix,1,typename OperatorType::PairType(2,1),1.0/sqrt(2.0),su2related2);
				creationMatrix.push_back(myOp2);

				// Set ni matrices:
				SparseMatrixType tmpMatrix = findNiMatrices(0,natBasis);
				RealType angularFactor= 1;
				typename OperatorType::Su2RelatedType su2related3;
				su2related3.offset = 1; //check FIXME
				OperatorType myOp3(tmpMatrix,1,typename OperatorType::PairType(0,0),angularFactor,su2related3);

				creationMatrix.push_back(myOp3);
			}
		}

		/** \cppFunction{!PTEX_THISFUNCTION} returns the operator in the unmangled (natural) basis of one-site */
		PsimagLite::Matrix<SparseElementType> naturalOperator(const PsimagLite::String& what,
										  SizeType site,
										  SizeType dof) const
		{
			BlockType block;
			block.resize(1);
			block[0]=site;
			typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
			setOperatorMatrices(creationMatrix,block);
			assert(creationMatrix.size()>0);
			SizeType nrow = creationMatrix[0].data.row();

			if (what == "i" || what=="identity") {
				PsimagLite::Matrix<SparseElementType> tmp(nrow,nrow);
				for (SizeType i = 0; i < tmp.n_row(); ++i) tmp(i,i) = 1.0;
				return tmp;
			}

			if (what=="+") {
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[2].data);
				return tmp;
			}

			if (what=="-") {
				SparseMatrixType tmp2;
				transposeConjugate(tmp2,creationMatrix[2].data);
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,tmp2);
				return tmp;
			}

			if (what=="z") {
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[3].data);
				return tmp;
			}

			if (what=="c") {
				PsimagLite::Matrix<SparseElementType> tmp;
				assert(dof<2);
				crsMatrixToFullMatrix(tmp,creationMatrix[dof].data);
				return tmp;
			}

			if (what=="n") {
				PsimagLite::Matrix<SparseElementType> tmp;
				crsMatrixToFullMatrix(tmp,creationMatrix[4].data);
				return tmp;
			}

			if (what=="nup") {
				PsimagLite::Matrix<SparseElementType> cup = naturalOperator("c",site,SPIN_UP);
				PsimagLite::Matrix<SparseElementType> nup = multiplyTransposeConjugate(cup,cup);
				return nup;
			}

			if (what=="ndown") {
				PsimagLite::Matrix<SparseElementType> cdown = naturalOperator("c",site,SPIN_DOWN);
				PsimagLite::Matrix<SparseElementType> ndown = multiplyTransposeConjugate(cdown,cdown);
				return ndown;
			}

			std::cerr<<"Argument: "<<what<<" "<<__FILE__<<"\n";
			throw std::logic_error("DmrgObserve::spinOperator(): invalid argument\n");
		}

		//! find total number of electrons for each state in the basis
		void findElectrons(typename PsimagLite::Vector<SizeType> ::Type&electrons,
					   const typename PsimagLite::Vector<HilbertStateType>::Type& basis,
					   SizeType site) const
		{
			int nup,ndown;
			electrons.clear();
			for (SizeType i=0;i<basis.size();i++) {
				nup = HilbertSpaceHubbardType::getNofDigits(basis[i],0);
				ndown = HilbertSpaceHubbardType::getNofDigits(basis[i],1);
				electrons.push_back(nup+ndown);
			}
		}

		//! find all states in the natural basis for a block of n sites
		//! N.B.: HAS BEEN CHANGED TO ACCOMODATE FOR MULTIPLE BANDS
		void setNaturalBasis(HilbertBasisType  &basis,
					 typename PsimagLite::Vector<SizeType>::Type& q,
					 const typename PsimagLite::Vector<SizeType>::Type& block) const
		{
			assert(block.size()==1);
			HilbertStateType a=0;
			int sitesTimesDof=DEGREES_OF_FREEDOM;
			HilbertStateType total = (1<<sitesTimesDof);
			total--;

			HilbertBasisType  basisTmp;
			for (a=0;a<total;a++) basisTmp.push_back(a);

			// reorder the natural basis (needed for MULTIPLE BANDS)
			findQuantumNumbers(q,basisTmp,1);
			typename PsimagLite::Vector<SizeType>::Type iperm(q.size());

			PsimagLite::Sort<typename PsimagLite::Vector<SizeType>::Type > sort;
			sort.sort(q,iperm);
			basis.clear();
			for (a=0;a<total;a++) basis.push_back(basisTmp[iperm[a]]);
		}

		void print(std::ostream& os) const
		{
			os<<modelParameters_;
		}

	private:

		//! Calculate fermionic sign when applying operator c^\dagger_{i\sigma} to basis state ket
		RealType sign(HilbertStateType const &ket, int i,int sigma) const
		{
			return 1;
		}

		//! Find c^\dagger_isigma in the natural basis natBasis
		SparseMatrixType findOperatorMatrices(int i,
		                                      int sigma,
		                                      const typename PsimagLite::Vector<HilbertStateType>::Type& natBasis) const
		{
			HilbertStateType bra,ket;
			int n = natBasis.size();
			PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

			for (SizeType ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				if (HilbertSpaceHubbardType::isNonZero(ket,i,sigma)) {

				} else {
					HilbertSpaceHubbardType::create(bra,i,sigma);
					int jj = PsimagLite::isInVector(natBasis,bra);
					if (jj<0) continue;
					cm(ii,jj) =sign(ket,i,sigma);
				}
			}
			//std::cout<<"Cm\n";
			//std::cout<<cm;
			SparseMatrixType creationMatrix(cm);
			return creationMatrix;
		}

		//! Find S^+_i in the natural basis natBasis
		SparseMatrixType findSplusMatrices(int i,
		                                   const typename PsimagLite::Vector<HilbertStateType>::Type& natBasis) const
		{
			HilbertStateType bra,ket;
			int n = natBasis.size();
			PsimagLite::Matrix<SparseElementType> cm(n,n);

			for (SizeType ii=0;ii<natBasis.size();ii++) {
				bra=ket=natBasis[ii];
				if (HilbertSpaceHubbardType::get(ket,i)==2) {
					// it is a down electron, then flip it:
					HilbertSpaceHubbardType::destroy(bra,i,1);
					HilbertSpaceHubbardType::create(bra,i,0);
					int jj = PsimagLite::isInVector(natBasis,bra);
					assert(jj>=0);
					cm(ii,jj)=1.0;
				}
			}
			SparseMatrixType operatorMatrix(cm);
			return operatorMatrix;
		}

		//! Find S^z_i in the natural basis natBasis
		SparseMatrixType findSzMatrices(int i,
		                                const typename PsimagLite::Vector<HilbertStateType>::Type& natBasis) const
		{
			HilbertStateType ket;
			int n = natBasis.size();
			PsimagLite::Matrix<SparseElementType> cm(n,n);

			for (SizeType ii=0;ii<natBasis.size();ii++) {
				ket=natBasis[ii];
				SizeType value = HilbertSpaceHubbardType::get(ket,i);
				switch (value) {
					case 1:
						cm(ii,ii)=0.5;
						break;
					case 2:
						cm(ii,ii)= -0.5;
						break;
				}
			}
			SparseMatrixType operatorMatrix(cm);
			return operatorMatrix;
		}

		//! Find n_i in the natural basis natBasis
		SparseMatrixType findNiMatrices(int i,
		                                const typename PsimagLite::Vector<HilbertStateType>::Type& natBasis) const
		{
			SizeType n = natBasis.size();
			PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n,n);

			for (SizeType ii=0;ii<natBasis.size();ii++) {
				HilbertStateType ket=natBasis[ii];
				cm(ii,ii) = 0.0;
				for (SizeType sigma=0;sigma<2;sigma++)
					if (HilbertSpaceHubbardType::isNonZero(ket,i,sigma))
						cm(ii,ii) += 1.0;
			}
			SparseMatrixType creationMatrix(cm);
			return creationMatrix;
		}

		void addDiagonalsInNaturalBasis(SparseMatrixType &hmatrix,
		                                const VectorOperatorType& cm,
		                                const BlockType& block,
		                                RealType time,
		                                RealType factorForDiagonals=1.0) const
		{
			SizeType n=block.size();

			SizeType linSize = geometry_.numberOfSites();
			for (SizeType i=0;i<n;i++) {
				// potentialV
				SparseMatrixType nup(naturalOperator("nup",i,0));
				SparseMatrixType ndown(naturalOperator("ndown",i,0));
				SparseMatrixType m = nup;
				assert(block[i]+linSize<modelParameters_.potentialV.size());
				m *= modelParameters_.potentialV[block[i]];
				m += modelParameters_.potentialV[block[i]+linSize]*ndown;
				hmatrix += factorForDiagonals * m;
			}
		}

		void findQuantumNumbers(typename PsimagLite::Vector<SizeType>::Type& q,const HilbertBasisType  &basis,int n) const
		{
			BasisDataType qq;
			setSymmetryRelated(qq,basis,n);
			MyBasis::findQuantumNumbers(q,qq);
		}

		void setSymmetryRelated(BasisDataType& q,HilbertBasisType  const &basis,int n) const
		{
			assert(n==1);

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
				// nup
				electronsUp[i] = HilbertSpaceHubbardType::getNofDigits(basis[i],HilbertSpaceHubbardType::SPIN_UP);
				// ndown
				electronsDown[i] = HilbertSpaceHubbardType::getNofDigits(basis[i],HilbertSpaceHubbardType::SPIN_DOWN);

				flavors.push_back(electronsUp[i]+electronsDown[i]);
				jmSaved = jmpair;
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
		PairType calcJmvalue(const HilbertStateType& ket) const
		{
			return calcJmValueAux<PairType>(ket);
		}

		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		// does not work for 6 or 9
		template<typename PairType>
		PairType calcJmValueAux(const HilbertStateType& ket) const
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

		//serializr start class Tj1Orb
		//serializr vptr
		//serializr normal modelParameters_
		ParametersModelTj1Orb<RealType>  modelParameters_;
		//serializr ref geometry_ end
		const GeometryType &geometry_;
		//serializr normal offset_
		SizeType offset_;
		//serializr normal spinSquaredHelper_
		SpinSquaredHelper<RealType,HilbertStateType> spinSquaredHelper_;
		//serializr normal spinSquared_
		SpinSquared<SpinSquaredHelper<RealType,HilbertStateType> > spinSquared_;

	};	//class Tj1Orb

} // namespace Dmrg
/*@}*/
#endif // TJ_1ORB_H

