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

/*! \file OperatorsBase.h
 *
 *  Documentation needed FIXME 
 *
 */
#ifndef OPERATORS_BASE_H
#define OPERATORS_BASE_H

#include "ReducedOperators.h"
#include "Range.h"
#include <cassert>
#include "ProgressIndicator.h"

namespace Dmrg {
/**
The \\cppClass{!PTEX_THISCLASS} class stores the local operators for this basis.
Only the local operators corresponding to the most recently added sites
will be meaningful. Indeed, if we  apply transformation $W$ (possibly
truncating the basis, see Eq.~(\\ref{eq:transformation})) then
\\begin{equation}
(W^\\dagger A W)  (W^\\dagger BW) \\neq W^\\dagger  (AB)  W,
\\end{equation}
since $WW^\\dagger\\neq 1$ because the DMRG truncation does not assure us
that $W^\\dagger$ will be the right inverse of $W$ (but $W^\\dagger W=1$
always holds). Because of this reason we cannot construct the Hamiltonian
simply from transformed local operators, even if we store them for all sites,
but we need to store also the Hamiltonian in the most recently transformed
basis. The fact that \\cppClass{!PTEX_THISCLASS} stores local operators in
the most recently transformed basis for \\emph{all sites} does not increase
memory usage too much, and simplifies the writing of code for complicated
geometries or connections, because all local opeators are availabel at all
times. Each SCE model class is responsible for determining whether a
transformed operator can be used (or not because of the reason limitation above).
*/
	template<typename OperatorType_,typename BasisType_>
	class OperatorsBase {
		
		typedef std::pair<size_t,size_t> PairType;
		static const bool EXCLUDE = false;

	public:
		
		typedef BasisType_ BasisType;
		typedef OperatorType_ OperatorType;
		typedef typename OperatorType::SparseMatrixType SparseMatrixType;

		OperatorsBase(const BasisType* thisBasis)
		: useSu2Symmetry_(BasisType::useSu2Symmetry()),
		  reducedOpImpl_(thisBasis),
		  progress_("Operators",0)
		{}

		template<typename IoInputter>
		OperatorsBase(IoInputter& io,
		              size_t level,
		              const BasisType* thisBasis)
		: useSu2Symmetry_(BasisType::useSu2Symmetry()),
		  reducedOpImpl_(io,level,thisBasis),
		  progress_("Operators",0)
		{
			if (!useSu2Symmetry_) io.read(operators_,"#OPERATORS");

			io.readMatrix(hamiltonian_,"#HAMILTONIAN");
			reducedOpImpl_.setHamiltonian(hamiltonian_);
		}

		template<typename IoInputter>
		void load(IoInputter& io)
		{
			if (!useSu2Symmetry_)
				io.read(operators_,"#OPERATORS");
			else reducedOpImpl_.load(io);

			io.readMatrix(hamiltonian_,"#HAMILTONIAN");
			reducedOpImpl_.setHamiltonian(hamiltonian_);
		}

		void setOperators(const std::vector<OperatorType>& ops)
		{
			if (!useSu2Symmetry_) operators_=ops;
			else reducedOpImpl_.setOperators(ops);
		}
		
		const OperatorType& getReducedOperatorByIndex(char modifier,const PairType& p) const
		{
			assert(useSu2Symmetry_);
			return reducedOpImpl_.getReducedOperatorByIndex(modifier,p);
		}

		const OperatorType& getOperatorByIndex(int i) const
		{
			assert(!useSu2Symmetry_);
			return operators_[i];
		}

		const OperatorType& getReducedOperatorByIndex(int i) const
		{
			assert(useSu2Symmetry_);
			return reducedOpImpl_.getReducedOperatorByIndex(i);
		}

		size_t numberOfOperators() const
		{
			if (useSu2Symmetry_) return reducedOpImpl_.size();
			return operators_.size();
		}

		template<typename ConcurrencyType>
		void changeBasis(const SparseMatrixType& ftransform,
		                 const BasisType* thisBasis,
				 ConcurrencyType &concurrency)
// 		                 const std::pair<size_t,size_t>& startEnd)
		{
			reducedOpImpl_.prepareTransform(ftransform,thisBasis);
			size_t total = numberOfOperators();

			PsimagLite::Range<ConcurrencyType> range(0,total,concurrency);
			for (;!range.end();range.next()) {
				size_t k = range.index();
				if (isExcluded(k,thisBasis)) {
					operators_[k].data.clear(); //resize(ftransform.n_col(),ftransform.n_col());
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

		void changeBasis(SparseMatrixType &v,const SparseMatrixType& ftransform)
		{
			SparseMatrixType transformConj;
			transposeConjugate(transformConj,ftransform);
			SparseMatrixType tmp = v*ftransform;
			multiply(v,transformConj,tmp);
		}

		void reorder(const std::vector<size_t>& permutation)
		{
			for (size_t k=0;k<numberOfOperators();k++) {
				if (!useSu2Symmetry_) reorder(operators_[k].data,permutation);
				reducedOpImpl_.reorder(k,permutation);
			}
			reorder(hamiltonian_,permutation);
			reducedOpImpl_.reorderHamiltonian(permutation);
		}

		void setMomentumOfOperators(const std::vector<size_t>& momentum)
		{
			reducedOpImpl_.setMomentumOfOperators(momentum);
		}

		void setToProduct(const BasisType& basis2,const BasisType& basis3,size_t x,const BasisType* thisBasis)
		{
			if (!useSu2Symmetry_) operators_.resize(x);
			reducedOpImpl_.setToProduct(basis2,basis3,x,thisBasis);
		}

		/**
		I will know explain how the full outer product between two operators
		is implemented. If local operator $A$ lives in Hilbert space
		$\\mathcal{A}$ and local operator $B$ lives in Hilbert space
		$\\mathcal{B}$, then $C=AB$ lives in Hilbert space
		$\\mathcal{C}=\\mathcal{A}\\otimes\\mathcal{B}$. Let $\\alpha_1$ and
		$\\alpha_2$ represent states of $\\mathcal{A}$, and let $\\beta_1$ and
		$\\beta_2$ represent states of   $\\mathcal{B}$. Then, in the product
		basis, $C_{\\alpha_1,\\beta_1;\\alpha_2,\\beta_2}=A_{\\alpha_1,\\alpha_2}
		B_{\\beta_1,\\beta_2}$. Additionally,  $\\mathcal{C}$ is reordered
		such that each state of this outer product basis is labeled in
		increasing effective quantum number (see
		Section~\\ref{sec:dmrgbasis}). In the previous example, if the Hilbert
		spaces  $\\mathcal{A}$ and $\\mathcal{B}$ had sizes $a$ and $b$,
		respectively, then their outer product would have size $ab$.
		When we add sites to the system (or the environment) the memory
		usage remains bounded by the truncation, and it is usually not a
		problem to store full product matrices, as long as we do it in a
		sparse way (DMRG++ uses compressed row storage). In short, local
		operators are always stored in the most recently transformed basis
		for \\emph{all sites} and, if applicable, \\emph{all values} of the
		internal degree of freedom $\\sigma$. See !PTEX\\_REF{setToProductOps}
		and !PTEX\\_REF{HERE}.
		*/
		template<typename ApplyFactorsType>
		void externalProduct(size_t i,
		                     const OperatorType& m,
		                     int x,
		                     const std::vector<double>& fermionicSigns,
		                     bool option,
		                     ApplyFactorsType& apply)
		{
			assert(!useSu2Symmetry_);
			PsimagLite::externalProduct(operators_[i].data,m.data,x,fermionicSigns,option);
			// don't forget to set fermion sign and j:
			operators_[i].fermionSign=m.fermionSign;
			operators_[i].jm=m.jm;
			operators_[i].angularFactor=m.angularFactor;
			apply(operators_[i].data);
		}

		void externalProductReduced(size_t i,
		                            const BasisType& basis2,
		                            const BasisType& basis3,
		                            bool option,
		                            const OperatorType& A)
		{
			reducedOpImpl_.externalProduct(i,basis2,basis3,option,A);
		}

		template<typename ApplyFactorsType>
		void outerProductHamiltonian(const SparseMatrixType& h2,const SparseMatrixType& h3,ApplyFactorsType& apply)
		{
			SparseMatrixType tmpMatrix;
			assert(h2.row()==h2.col());
			std::vector<double> ones(h2.row(),1.0);
			PsimagLite::externalProduct(hamiltonian_,h2,h3.row(),ones,true);

			PsimagLite::externalProduct(tmpMatrix,h3,h2.row(),ones,false);

			hamiltonian_ += tmpMatrix;

			apply(hamiltonian_);
		}

		void outerProductHamiltonianReduced(const BasisType& basis2,
		                                    const BasisType& basis3,
		                                    const SparseMatrixType& h2,
		                                    const SparseMatrixType& h3)
		{
			reducedOpImpl_.outerProductHamiltonian(basis2,basis3,h2,h3);
		}

		void setHamiltonian(SparseMatrixType const &h)
		{
			hamiltonian_=h;
			reducedOpImpl_.setHamiltonian(h);
		}

		const SparseMatrixType& hamiltonian() const { return hamiltonian_; }

		const SparseMatrixType& reducedHamiltonian() const
		{
			return reducedOpImpl_.hamiltonian();
		}

		//const std::vector<size_t>& electrons() const {return operatorsImpl_.electrons(); }

		void print(int ind= -1) const
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

		void reorder(SparseMatrixType &v,const std::vector<size_t>& permutation)
		{
			SparseMatrixType matrixTmp;

			permute(matrixTmp,v,permutation);
			permuteInverse(v,matrixTmp,permutation);
		}

		bool isExcluded(size_t k,
			       const BasisType* thisBasis)
// 		               const std::pair<size_t,size_t>& startEnd)
		{
			if (!EXCLUDE) return false; // <-- this is the safest answer
// 			if (k<startEnd.first || k>=startEnd.second) return true;
			return false;
		}

		bool useSu2Symmetry_;
		ReducedOperators<OperatorType,BasisType> reducedOpImpl_;
		std::vector<OperatorType> operators_;
		SparseMatrixType hamiltonian_;
		PsimagLite::ProgressIndicator progress_;
	}; //class OperatorsBase 
} // namespace Dmrg

/*@}*/
#endif

