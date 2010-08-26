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

/*! \file OperatorImplementation.h
 *
 *  
 *
 */
#ifndef OPERATOR_IMPL_H
#define OPERATOR_IMPL_H

#include "ReducedOperators.h"

namespace Dmrg {
	//! 
	template<typename OperatorType,typename DmrgBasisType>
	class OperatorsImplementation {
	public:	
		typedef typename OperatorType::SparseMatrixType SparseMatrixType;

		OperatorsImplementation(const DmrgBasisType* thisBasis,
				       size_t dof,size_t nOrbitals) :
			useSu2Symmetry_(DmrgBasisType::useSu2Symmetry()),reducedOpImpl_(thisBasis,dof,nOrbitals) 
		{
		}
		
		template<typename IoInputter>
		OperatorsImplementation(IoInputter& io,size_t level,const DmrgBasisType* thisBasis,
					size_t dof,size_t orbitals) 
			: useSu2Symmetry_(DmrgBasisType::useSu2Symmetry()),
			  reducedOpImpl_(io,level,thisBasis,dof,orbitals) 
		{
			if (!useSu2Symmetry_)
				io.read(operators_,"#OPERATORS");
			//else reducedOpImpl_.load(io,level);

			io.readMatrix(hamiltonian_,"#HAMILTONIAN");
			reducedOpImpl_.setHamiltonian(hamiltonian_);
		}

		void setOperators(const std::vector<OperatorType>& ops)
		{
			if (!useSu2Symmetry_) operators_=ops;
			else reducedOpImpl_.setOperators(ops);
		}

		const OperatorType& getOperatorByIndex(int i) const 
		{
			if (useSu2Symmetry_) throw std::runtime_error("EERRRRRRRRRRRRRRRORRRRRRRRRRRR\n");
			return operators_[i];
		}

		const OperatorType& getReducedOperatorByIndex(int i) const 
		{
			if (!useSu2Symmetry_) throw std::runtime_error("EERRRRRRRRRRRRRRRORRRRRRRRRRRR\n");
			return reducedOpImpl_.getReducedOperatorByIndex(i);
		}

		const OperatorType& getReducedOperatorByIndex(char modifier,int i) const 
		{
			if (!useSu2Symmetry_) throw std::runtime_error("EERRRRRRRRRRRRRRRORRRRRRRRRRRR\n");
			return reducedOpImpl_.getReducedOperatorByIndex(modifier,i);		
		}

		size_t size() const 
		{
			if (useSu2Symmetry_) return reducedOpImpl_.size();
			return operators_.size(); 
		}

		template<typename TransformElementType,typename ConcurrencyType>
		void changeBasis(psimag::Matrix<TransformElementType> const &ftransform,const DmrgBasisType* thisBasis,
					ConcurrencyType &concurrency)
		{
			size_t total = size();
			size_t k;

			concurrency.loopCreate(total);

			reducedOpImpl_.prepareTransform(ftransform,thisBasis);
			size_t dof = total / thisBasis->block().size();	
			while(concurrency.loop(k)) {
				if (isExcluded(k,thisBasis,dof)) {
					operators_[k].data.resize(ftransform.n_col(),ftransform.n_col());
					continue;
				}
				if (!useSu2Symmetry_) changeBasis(operators_[k].data,ftransform);
				reducedOpImpl_.changeBasis(k);
			}

			if (!useSu2Symmetry_) {
				gather(operators_,concurrency);
				broadcast(operators_,concurrency);
			} else {	
				reducedOpImpl_.gather(concurrency);
				reducedOpImpl_.broadcast(concurrency);
			}
			
			changeBasis(hamiltonian_,ftransform);
			reducedOpImpl_.changeBasisHamiltonian();
		}
		
		bool isExcluded(size_t k,const DmrgBasisType* thisBasis,size_t dof)
		{
			return false; // disabled for now
		}

		void changeBasis(SparseMatrixType &v,psimag::Matrix<typename SparseMatrixType::value_type> const &ftransform)
		{
			fullMatrixToCrsMatrix(v,transformFullFast(v,ftransform));
		}

		void reorder(const std::vector<size_t>& permutation)
		{
			for (size_t k=0;k<size();k++) {
				if (!useSu2Symmetry_) reorder(operators_[k].data,permutation);
				reducedOpImpl_.reorder(k,permutation);
			}
			reorder(hamiltonian_,permutation);
			reducedOpImpl_.reorderHamiltonian(permutation);
		}

		void setToProduct(const DmrgBasisType& basis2, const DmrgBasisType& basis3,size_t x,const DmrgBasisType* thisBasis)
		{
			if (!useSu2Symmetry_) operators_.resize(x);
			reducedOpImpl_.setToProduct(basis2,basis3,x,thisBasis); 
		}

		void setMomentumOfOperators(const std::vector<size_t>& momentum)
		{
			reducedOpImpl_.setMomentumOfOperators(momentum);
		}

		template<typename ApplyFactorsType>
		void externalProduct(size_t i,const OperatorType& m,int x,const std::vector<double>& fermionicSigns,
				     bool option,ApplyFactorsType& apply)
		{
			if (useSu2Symmetry_) std::cerr<<"#######EEEEEEEEERRRRRRRRRRRRRROOOOOOOOOOORRRRRRRRRR\n"; 
			Dmrg::externalProduct(operators_[i].data,m.data,x,fermionicSigns,option);
			// don't forget to set fermion sign and j:
			operators_[i].fermionSign=m.fermionSign;
			operators_[i].jm=m.jm;
			operators_[i].angularFactor=m.angularFactor;
			apply(operators_[i].data);
		}

		void externalProductReduced(size_t i,const DmrgBasisType& basis2,const DmrgBasisType& basis3,bool option,
					    const OperatorType& A)
		{
			reducedOpImpl_.externalProduct(i,basis2,basis3,option,A);
		}

		template<typename ApplyFactorsType>
		void outerProductHamiltonian(const SparseMatrixType& h2,const SparseMatrixType& h3,ApplyFactorsType& apply)
		{
			SparseMatrixType tmpMatrix;
			std::vector<double> ones(h2.rank(),1.0);
			Dmrg::externalProduct(hamiltonian_,h2,h3.rank(),ones,true);

			Dmrg::externalProduct(tmpMatrix,h3,h2.rank(),ones,false);

			hamiltonian_ += tmpMatrix;

			apply(hamiltonian_);
		}

		void outerProductHamiltonianReduced(const DmrgBasisType& basis2,const DmrgBasisType& basis3,
				const SparseMatrixType& h2,const SparseMatrixType& h3)
		{
			reducedOpImpl_.outerProductHamiltonian(basis2,basis3,h2,h3);
		}

		void setHamiltonian(SparseMatrixType const &h) 
		{ 
			hamiltonian_=h; 
			reducedOpImpl_.setHamiltonian(h);
		}

		const SparseMatrixType& hamiltonian() const { return hamiltonian_; }

		const SparseMatrixType& reducedHamiltonian() const { return reducedOpImpl_.hamiltonian(); }

		void print(int ind) const
		{
			if (!useSu2Symmetry_) {
				if (ind<0) for (size_t i=0;i<operators_.size();i++) std::cerr<<operators_[i];
				else std::cerr<<operators_[ind];
			} else {
				reducedOpImpl_.print(ind);
			}
		}

		template<typename IoOutputter>
		void save(IoOutputter& io,const std::string& s) const
		{
			if (!useSu2Symmetry_) io.printVector(operators_,"#OPERATORS");
			else reducedOpImpl_.save(io,s);
			io.printMatrix(hamiltonian_,"#HAMILTONIAN");
		}

	private:
		bool useSu2Symmetry_;
		ReducedOperators<OperatorType,DmrgBasisType> reducedOpImpl_;
		std::vector<OperatorType> operators_;
		SparseMatrixType hamiltonian_;
		
		void reorder(SparseMatrixType &v,const std::vector<size_t>& permutation)
		{
			SparseMatrixType matrixTmp;

			permute(matrixTmp,v,permutation);
			permuteInverse(v,matrixTmp,permutation);
		}
	}; //class OperatorsImplementation

} // namespace Dmrg

/*@}*/
#endif
