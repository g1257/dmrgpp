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


#ifndef BASISWITHOPERATORS_HEADER_H
#define BASISWITHOPERATORS_HEADER_H

#include "ApplyFactors.h"
#include "Basis.h"

namespace Dmrg {

	/*! A class to represent a Hilbert Space for a strongly correlated electron model
	    Derives from Basis

	 C++ class \\cppClass{Basis} (and \\cppClass{BasisImplementation}) implement only 
	 certain functionality associated with a Hilbert space basis, as mentioned in 
	 the previous section. However, more capabilities related to a Hilbert space basis are needed.

	 C++ class \\cppClass{!PTEX_THISCLASS} inherits from \\cppClass{Basis}, and contains 
	 certain local operators for the basis in question, as well as the Hamiltonian matrix.
	 The operators that need to be considered here are operators necessary to compute 
	 the Hamiltonian across the system and environment, and to compute observables. 
	 Therefore, the specific operators vary from model to model.
	 For example, for the Hubbard model, we consider $c_{i\\sigma}$ operators, 
	 that destroy an electron with spin $\\sigma$ on site $i$.
	 For the Heisenberg model, we consider operators $S^+_i$ and $S^z_i$ for each site $i$. 
	 In each case these operators are calculated by the model class (see Section~\\ref{subsec:models}) 
	 on the ``natural'' basis,  and added to the basis in question with a call to 
	 \\cppFunction{setOperators()}.  
	 These local operators are stored as sparse matrices to save memory, 
	 although the matrix type is templated and could be anything.
	 For details on the implementation of these operators, see \\cppClass{OperatorsBase},
	 its common implementation \\cppClass{OperatorsImplementation}, and the two examples provided
	 \\cppClass{OperatorsHubbard} and \\cppClass{OperatorsHeisenberg} for the Hubbard 
	 and Heisenberg models, respectively.
	 Additionally, \\cppClass{!PTEX_THISCLASS} has a number of member functions to 
	 handle operations that the DMRG method performs on 
	 local operators in a Hilbert space basis. These include functions to create 
	 an outer product of two given Hilbert spaces, to transform a basis, to truncate a basis, etc. 
	 */
	template<typename OperatorsType_,typename ConcurrencyType_>
	class	BasisWithOperators : public  OperatorsType_::BasisType {

		typedef std::pair<size_t,size_t> PairType;


	public:

		typedef typename OperatorsType_::BasisType::RealType RealType;
		typedef OperatorsType_ OperatorsType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename OperatorsType::BasisType BasisType;
		typedef typename BasisType::BlockType BlockType;
		typedef typename OperatorType::SparseMatrixType SparseMatrixType;
		typedef BasisWithOperators<OperatorsType,ConcurrencyType> 	ThisType;
		typedef typename BasisType::BasisDataType BasisDataType;
		typedef typename BasisType::FactorsType FactorsType;

		enum {GROW_RIGHT,GROW_LEFT};

		BasisWithOperators(const std::string& s) :BasisType(s),operators_(this) {}

		template<typename IoInputter>
		BasisWithOperators(IoInputter& io,
		                   const std::string& ss,
		                   size_t counter=0)
		: BasisType(io,ss,counter),operators_(io,0,this)
		{
			io.read(operatorsPerSite_,"#OPERATORSPERSITE");
		}

		template<typename IoInputter>
		void load(IoInputter& io)
		{
			BasisType::load(io); // parent loads
			operators_.load(io);
			io.read(operatorsPerSite_,"#OPERATORSPERSITE");
		}

		//! set this basis to the outer product of   basis2 and basis3 or basis3 and basis2  depending on dir
		void setToProduct(const ThisType& basis2,const ThisType& basis3,int dir)
		{
			if (dir==GROW_RIGHT) setToProduct(basis2,basis3);
			else setToProduct(basis3,basis2);
		}

		//! set this basis to the outer product of   basis2 and basis3 
		//!PTEX_LABEL{setToProductOps}
		void setToProduct(const ThisType& basis2,const ThisType& basis3)
		{
			BasisType &parent = *this;
			// reorder the basis
			parent.setToProduct(basis2,basis3);

			std::vector<double> fermionicSigns;
			size_t x = basis2.numberOfOperators()+basis3.numberOfOperators();

			if (this->useSu2Symmetry()) setMomentumOfOperators(basis2);
			operators_.setToProduct(basis2,basis3,x,this);
			ApplyFactors<FactorsType> apply(this->getFactors(),this->useSu2Symmetry());
			int savedSign = 0;

			for (size_t i=0;i<this->numberOfOperators();i++) {
				if (i<basis2.numberOfOperators()) {
					if (!this->useSu2Symmetry()) {
						const OperatorType& myOp =  basis2.getOperatorByIndex(i);
						if (savedSign != myOp.fermionSign) {
							fillFermionicSigns(fermionicSigns,basis2.electronsVector(),myOp.fermionSign);
							savedSign = myOp.fermionSign;
						}
						operators_.externalProduct(i,myOp,basis3.size(),fermionicSigns,true,apply);
					} else {
						operators_.externalProductReduced(i,basis2,basis3,true,basis2.getReducedOperatorByIndex(i));
					}
				} else {
					if (!this->useSu2Symmetry()) {
						const OperatorType& myOp = basis3.getOperatorByIndex(i-basis2.numberOfOperators());

						if (savedSign != myOp.fermionSign) {
							fillFermionicSigns(fermionicSigns,basis2.electronsVector(),myOp.fermionSign);
							savedSign = myOp.fermionSign;
						}
						operators_.externalProduct(i,myOp,basis2.size(),fermionicSigns,false,apply);
					} else {
						operators_.externalProductReduced(i,basis2,basis3,false,
						   basis3.getReducedOperatorByIndex(i-basis2.numberOfOperators()));
					}
				}
			}

			//! Calc. hamiltonian
			operators_.outerProductHamiltonian(basis2.hamiltonian(),basis3.hamiltonian(),apply);
			operators_.outerProductHamiltonianReduced(basis2,basis3,basis2.reducedHamiltonian(),basis3.reducedHamiltonian());
			//! re-order operators and hamiltonian 
			operators_.reorder(this->permutationVector());
			
			size_t offset1 = basis2.operatorsPerSite_.size();
			operatorsPerSite_.resize(offset1+basis3.operatorsPerSite_.size());
			for (size_t i=0;i<offset1;i++)
				operatorsPerSite_[i] =  basis2.operatorsPerSite_[i];
			
			for (size_t i=0;i<basis3.operatorsPerSite_.size();i++)
				operatorsPerSite_[i+offset1] =  basis3.operatorsPerSite_[i];
		}

//		template<typename SolverParametersType>
//		void changeBasis(std::vector<size_t>& removedIndices,
//				 std::vector<RealType>& eigs,
//				 size_t kept,
//				 const SolverParametersType& solverParams)
//		{
//			this->changeBasis(removedIndices,eigs,kept,solverParams);
//		}

		//! transform this basis by transform 
		//! note: basis change must conserve total number of electrons and all quantum numbers
		template<typename BlockMatrixType>
		RealType truncateBasis(SparseMatrixType& ftransform,
				       const BlockMatrixType& transform,
				       const std::vector<RealType>& eigs,
				       const std::vector<size_t>& removedIndices,
				       ConcurrencyType& concurrency)
		{
			BasisType &parent = *this;
			RealType error = parent.truncateBasis(ftransform,transform,eigs,removedIndices);
			size_t newSize = transform.rank() - removedIndices.size();
			changeBasisDirect(ftransform,concurrency,newSize);

			return error;
		}

		void changeBasisDirect(const SparseMatrixType& ftransform, ConcurrencyType& concurrency,size_t newSize)
		{
			operators_.changeBasis(ftransform,this,concurrency,newSize); //startEnd);
		}

		void setHamiltonian(SparseMatrixType const &h) { operators_.setHamiltonian(h); }

		//! Returns the Hamiltonian as stored in this basis
		SparseMatrixType hamiltonian() const { return operators_.hamiltonian(); }

		SparseMatrixType reducedHamiltonian() const { return operators_.reducedHamiltonian(); }

		void setVarious(BlockType const &block,
		                SparseMatrixType const &h,
		                BasisDataType const &qm,
		                const std::vector<OperatorType>& ops)
		{
			set(block);
			this->setSymmetryRelated(qm);
			setHamiltonian(h);
			operators_.setOperators(ops);
			operatorsPerSite_.push_back(ops.size());
		}

// 		void getOperatorByIndex(std::string& s,int i) const
// 		{
// 			operators_.getOperatorByIndex(i);
// 		}

		PairType getOperatorIndices(size_t i,size_t sigma) const
		{
			size_t sum = 0;
			for (size_t j=0;j<i;j++)
				sum += operatorsPerSite_[j];

			return PairType(sum + sigma,operatorsPerSite_[i]);
		}

		const OperatorType& getOperatorByIndex(int i) const
		{
			return operators_.getOperatorByIndex(i);
		}

		const OperatorType& getReducedOperatorByIndex(int i) const
		{
			return operators_.getReducedOperatorByIndex(i);
		}
		
		const OperatorType& getReducedOperatorByIndex(char modifier,const PairType& p) const
		{
			return operators_.getReducedOperatorByIndex(modifier,p);
		}
// 
// 		const OperatorType& getOperator(int i,int sigma) const
// 		{
// 			return operators_.getOperator(i,sigma);
// 		}

// 		const OperatorType& getOperator(const std::string& s,int i) const
// 		{
// 			return operators_.getOperator(s,i);
// 		}

		size_t numberOfOperators() const { return operators_.numberOfOperators(); }

// 		static size_t numberOfOperatorsPerSite()
// 		{
// 			return 	OperatorsType::numberOfOperatorsPerSite();
// 		}

		int fermionicSign(size_t i,int fsign) const
		{
			const BasisType &parent = *this;
			return parent.fermionicSign(i,fsign);
		}

		template<typename IoOutputter>
		void save(IoOutputter& io,const std::string& s) const
		{
			BasisType::save(io,s); // parent saves
			operators_.save(io,s);
			io.printVector(operatorsPerSite_,"#OPERATORSPERSITE");
		}

		template<typename IoOutputter>
		void save(IoOutputter& io) const
		{
			BasisType::save(io); // parent saves
			operators_.save(io,this->name());
			io.printVector(operatorsPerSite_,"#OPERATORSPERSITE");
		}

	private:
		OperatorsType operators_;
		std::vector<size_t> operatorsPerSite_;

		template<typename SomeType>
		void fillFermionicSigns(std::vector<SomeType>& fermionicSigns,const std::vector<size_t>& electrons,int f)
		{
			fermionicSigns.resize(electrons.size());
			for (size_t i=0;i<fermionicSigns.size();i++)
				fermionicSigns[i]= (electrons[i]%2==0) ? 1.0 : static_cast<double>(f);
		}

		void setMomentumOfOperators(const ThisType& basis)
		{
			std::vector<size_t> momentum;
			for (size_t i=0;i<basis.numberOfOperators();i++) {
				int x = PsimagLite::isInVector(
				     momentum,basis.getReducedOperatorByIndex(i).jm.first);
				if (x<0)
					momentum.push_back(basis.getReducedOperatorByIndex(i).jm.first);
			}
			operators_.setMomentumOfOperators(momentum);
		}
	}; // class BasisWithOperators

	template<typename OperatorsType,typename ConcurrencyType>
	std::ostream& operator<<(std::ostream& os,
	                         const BasisWithOperators<OperatorsType,ConcurrencyType>& bwo)
	{
		throw std::runtime_error("Unimplemented <<");
		return os;
	}

	template<typename OperatorsType,typename ConcurrencyType>
	std::istream& operator>>(std::istream& is,
	                         BasisWithOperators<OperatorsType,ConcurrencyType>& bwo)
	{
		throw std::runtime_error("Unimplemented >>");
		return is;
	}
} // namespace Dmrg

/*@}*/
#endif

