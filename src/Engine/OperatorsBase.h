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

/*! \file OperatorsBase.h
 *
 *  
 *
 */
#ifndef OPERATORS_BASE_H
#define OPERATORS_BASE_H

#include "OperatorsImplementation.h"
#include "ProgressIndicator.h"

namespace Dmrg {
	//! This class must be inherited
	template<typename OperatorType_,typename BasisType_>
	class OperatorsBase {
	public:	
		typedef BasisType_ BasisType;
		typedef OperatorType_ OperatorType;
		typedef typename OperatorType::SparseMatrixType SparseMatrixType;

		OperatorsBase(const BasisType* thisBasis,size_t dof,size_t nOrbitals) :
			operatorsImpl_(thisBasis,dof,nOrbitals),
					progress_("Operators",0) 
		{ }
					
		template<typename IoInputter>
		OperatorsBase(IoInputter& io,size_t level,const BasisType* thisBasis,
			      size_t dof,size_t orbitals) 
			: operatorsImpl_(io,level,thisBasis,dof,orbitals),progress_("Operators",0) 
		{
		}

		void setOperators(const std::vector<OperatorType>& ops)
		{
			operatorsImpl_.setOperators(ops);
		}

		const OperatorType& getOperatorByIndex(int i) const 
		{
			return operatorsImpl_.getOperatorByIndex(i);			
		}

		const OperatorType& getReducedOperatorByIndex(int i) const 
		{
			return operatorsImpl_.getReducedOperatorByIndex(i);
		}

		const OperatorType& getReducedOperatorByIndex(char modifier,int i) const 
		{
			return operatorsImpl_.getReducedOperatorByIndex(modifier,i);
		}

		void getOperatorByIndex(std::string& s,int i) const
		{
			operatorsImpl_.getOperatorByIndex(i);
		}

		size_t numberOfOperators() const { return operatorsImpl_.size(); }

		template<typename TransformElementType,typename ConcurrencyType>
		void changeBasis(psimag::Matrix<TransformElementType> const &ftransform,const BasisType* thisBasis,
					ConcurrencyType &concurrency)
		{
			return operatorsImpl_.changeBasis(ftransform,thisBasis,concurrency);
			std::ostringstream msg;
			msg<<"Done with changeBasis";
			progress_.printline(msg,std::cerr);
		}

		template<typename TransformElementType>
		void changeBasis(SparseMatrixType &v,psimag::Matrix<TransformElementType> const &ftransform)
		{
			operatorsImpl_.changeBasis(v,ftransform);
		}

		void reorder(const std::vector<size_t>& permutation)
		{
			operatorsImpl_.reorder(permutation);
		}

		void reorder(SparseMatrixType &v,const std::vector<size_t>& permutation)
		{
			operatorsImpl_.reorder(v,permutation);
		}

		template<typename Field>
		void reorder(std::vector<Field> &data,const std::vector<size_t>& permutation)
		{
			operatorsImpl_.reorder(data,permutation);
		}

		void setMomentumOfOperators(const std::vector<size_t>& momentum)
		{
			operatorsImpl_.setMomentumOfOperators(momentum);
		}

		void setToProduct(const BasisType& basis2,const BasisType& basis3,size_t x,const BasisType* thisBasis)
		{
			operatorsImpl_.setToProduct(basis2,basis3,x,thisBasis);
		}

		template<typename ApplyFactorsType>
		void externalProduct(
				size_t i,
				const OperatorType& m,
    				int x,
				const std::vector<double>& fermionicSigns,
    				bool option,
				ApplyFactorsType& apply)
		{
			operatorsImpl_.externalProduct(i,m,x,fermionicSigns,option,apply);
		}

		void externalProductReduced(
				size_t i,
    				const BasisType& basis2,
				const BasisType& basis3,
				bool option,
				const OperatorType& A)
		{
			operatorsImpl_.externalProductReduced(i,basis2,basis3,option,A);
		}

		template<typename ApplyFactorsType>
		void outerProductHamiltonian(const SparseMatrixType& h2,const SparseMatrixType& h3,ApplyFactorsType& apply)
		{
			operatorsImpl_.outerProductHamiltonian(h2,h3,apply);
		}

		void outerProductHamiltonianReduced(const BasisType& basis2,const BasisType& basis3,
				const SparseMatrixType& h2,const SparseMatrixType& h3)
		{
			operatorsImpl_.outerProductHamiltonianReduced(basis2,basis3,h2,h3);
		}

		void setHamiltonian(SparseMatrixType const &h) { operatorsImpl_.setHamiltonian(h); }

		const SparseMatrixType& hamiltonian() const { return operatorsImpl_.hamiltonian(); }	

		const SparseMatrixType& reducedHamiltonian() const { return operatorsImpl_.reducedHamiltonian(); }	

		const std::vector<size_t>& electrons() const {return operatorsImpl_.electrons(); }

		void print(int i= -1) const { operatorsImpl_.print(i); }

		template<typename IoOutputter>
		void save(IoOutputter& io,const std::string& s) const
		{
			operatorsImpl_.save(io,s);
		}

	private:
		OperatorsImplementation<OperatorType,BasisType> operatorsImpl_;
		ProgressIndicator progress_;
	}; //class OperatorsBase 
} // namespace Dmrg

/*@}*/
#endif
