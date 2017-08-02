/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

/*! \file Operators.h
 *
 *  Documentation needed FIXME
 *
 */
#ifndef OPERATORS_H
#define OPERATORS_H

#include "ReducedOperators.h"
#include <cassert>
#include "ProgressIndicator.h"
#include "Complex.h"
#include "Concurrency.h"
#include "Parallelizer.h"

namespace Dmrg {
/* PSIDOC Operators
The \cppClass{Operators} class stores the local operators for this basis.
Only the local operators corresponding to the most recently added sites
will be meaningful. Indeed, if we  apply transformation $W$ (possibly
truncating the basis, then
\begin{equation}
(W^\dagger A W)  (W^\dagger BW) \neq W^\dagger  (AB)  W,
\end{equation}
since $WW^\dagger\neq 1$ because the DMRG truncation does not assure us
that $W^\dagger$ will be the right inverse of $W$ (but $W^\dagger W=1$
always holds). Because of this reason we cannot construct the Hamiltonian
simply from transformed local operators, even if we store them for all sites,
but we need to store also the Hamiltonian in the most recently transformed
basis. The fact that \cppClass{Operators} stores local operators in
the most recently transformed basis for \emph{all sites} does not increase
memory usage too much, and simplifies the writing of code for complicated
geometries or connections, because all local opeators are availabel at all
times. Each SCE model class is responsible for determining whether a
transformed operator can be used (or not because of the reason limitation above).
*/
template<typename BasisType_>
class Operators {

	typedef std::pair<SizeType,SizeType> PairType;

public:

	typedef BasisType_ BasisType;
	typedef ReducedOperators<BasisType> ReducedOperatorsType;
	typedef typename ReducedOperatorsType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename ReducedOperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType,SizeType> PairSizeSizeType;

	class MyLoop {

	public:

		MyLoop(bool useSu2Symmetry,
		       ReducedOperatorsType& reducedOpImpl,
		       typename PsimagLite::Vector<OperatorType>::Type& operators,
		       const BlockDiagonalMatrixType& ftransform1,
		       const BasisType* thisBasis1,
		       const PairSizeSizeType& startEnd)
		    : useSu2Symmetry_(useSu2Symmetry),
		      reducedOpImpl_(reducedOpImpl),
		      operators_(operators),
		      ftransform(ftransform1),
		      thisBasis(thisBasis1),
		      hasMpi_(ConcurrencyType::hasMpi()),
		      startEnd_(startEnd)
		{
			reducedOpImpl_.prepareTransform(ftransform,thisBasis);
		}

		void doTask(SizeType taskNumber , SizeType threadNum)
		{
			SizeType k = taskNumber;
			if (isExcluded(k) && k < operators_.size()) {
				operators_[k].data.clear();
				return;
			}

			if (!useSu2Symmetry_)
				reducedOpImpl_.changeBasis(operators_[k].data);
			else
				reducedOpImpl_.changeBasis(k);
		}

		SizeType tasks() const
		{
			if (useSu2Symmetry_) return reducedOpImpl_.size();
			return operators_.size();
		}

		void gather()
		{
			if (ConcurrencyType::isMpiDisabled("Operators")) return;

			if (!useSu2Symmetry_) {
				gatherOperators();
				bcastOperators();
			} else {
				reducedOpImpl_.gather();
				reducedOpImpl_.bcast();
			}
		}

	private:

		bool isExcluded(SizeType k) const
		{
#ifdef OPERATORS_CHANGE_ALL
			return false; // <-- this is the safest answer
#endif
			if (k < startEnd_.first || k >= startEnd_.second) return true;
			return false;
		}

		void gatherOperators()
		{
			if (!hasMpi_) return;
			PsimagLite::MPI::pointByPointGather(operators_);
		}

		void bcastOperators()
		{
			if (!hasMpi_) return;
			for (SizeType i = 0; i < operators_.size(); i++)
				Dmrg::bcast(operators_[i]);
		}

		bool useSu2Symmetry_;
		ReducedOperatorsType& reducedOpImpl_;
		typename PsimagLite::Vector<OperatorType>::Type& operators_;
		const BlockDiagonalMatrixType& ftransform;
		const BasisType* thisBasis;
		bool hasMpi_;
		const PairSizeSizeType& startEnd_;
	};

	Operators(const BasisType* thisBasis)
	    : useSu2Symmetry_(BasisType::useSu2Symmetry()),
	      reducedOpImpl_(thisBasis),
	      progress_("Operators")
	{
		announceChangeAll();
	}

	template<typename IoInputter>
	Operators(IoInputter& io,
	          SizeType level,
	          const BasisType* thisBasis,
	          bool isObserveCode)
	    : useSu2Symmetry_(BasisType::useSu2Symmetry()),
	      reducedOpImpl_(io,level,thisBasis),
	      progress_("Operators")
	{
		if (isObserveCode) return;

		announceChangeAll();

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

	void setOperators(const typename PsimagLite::Vector<OperatorType>::Type& ops)
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
		assert(i>=0 && SizeType(i)<operators_.size());
		return operators_[i];
	}

	const OperatorType& getReducedOperatorByIndex(int i) const
	{
		assert(useSu2Symmetry_);
		return reducedOpImpl_.getReducedOperatorByIndex(i);
	}

	SizeType numberOfOperators() const
	{
		if (useSu2Symmetry_) return reducedOpImpl_.size();
		return operators_.size();
	}

	void changeBasis(const BlockDiagonalMatrixType& ftransform,
	                 const BasisType* thisBasis,
	                 const PairSizeSizeType& startEnd)
	{
		typedef PsimagLite::Parallelizer<MyLoop> ParallelizerType;
		ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
		                              PsimagLite::MPI::COMM_WORLD);

		MyLoop helper(useSu2Symmetry_,reducedOpImpl_,operators_,ftransform,thisBasis,startEnd);

		threadObject.loopCreate(helper); // FIXME: needs weights

		helper.gather();

		reducedOpImpl_.changeBasisHamiltonian(hamiltonian_,ftransform);
	}

	void reorder(const   VectorSizeType& permutation)
	{
		for (SizeType k=0;k<numberOfOperators();k++) {
			if (!useSu2Symmetry_) reorder(operators_[k].data,permutation);
			reducedOpImpl_.reorder(k,permutation);
		}
		reorder(hamiltonian_,permutation);
		reducedOpImpl_.reorderHamiltonian(permutation);
	}

	void setMomentumOfOperators(const   VectorSizeType& momentum)
	{
		reducedOpImpl_.setMomentumOfOperators(momentum);
	}

	void setToProduct(const BasisType& basis2,
	                  const BasisType& basis3,
	                  SizeType x,
	                  const BasisType* thisBasis)
	{
		if (!useSu2Symmetry_) operators_.resize(x);
		reducedOpImpl_.setToProduct(basis2,basis3,x,thisBasis);
	}

	/* PSIDOC OperatorsExternalProduct
		I will know explain how the full outer product between two operators
		is implemented. If local operator $A$ lives in Hilbert space
		$\mathcal{A}$ and local operator $B$ lives in Hilbert space
		$\mathcal{B}$, then $C=AB$ lives in Hilbert space
		$\mathcal{C}=\mathcal{A}\otimes\mathcal{B}$. Let $\alpha_1$ and
		$\alpha_2$ represent states of $\mathcal{A}$, and let $\beta_1$ and
		$\beta_2$ represent states of   $\mathcal{B}$. Then, in the product
		basis, $C_{\alpha_1,\beta_1;\alpha_2,\beta_2}=A_{\alpha_1,\alpha_2}
		B_{\beta_1,\beta_2}$. Additionally,  $\mathcal{C}$ is reordered
		such that each state of this outer product basis is labeled in
		increasing effective quantum number (see
		Section~\ref{sec:dmrgbasis}). In the previous example, if the Hilbert
		spaces  $\mathcal{A}$ and $\mathcal{B}$ had sizes $a$ and $b$,
		respectively, then their outer product would have size $ab$.
		When we add sites to the system (or the environment) the memory
		usage remains bounded by the truncation, and it is usually not a
		problem to store full product matrices, as long as we do it in a
		sparse way (DMRG++ uses compressed row storage). In short, local
		operators are always stored in the most recently transformed basis
		for \emph{all sites} and, if applicable, \emph{all values} of the
		internal degree of freedom $\sigma$. See PTEXREF{setToProductOps}
		and PTEXREF{HERE}.
		*/
	template<typename ApplyFactorsType>
	void externalProduct(SizeType i,
	                     const OperatorType& m,
	                     int x,
	                     const   VectorRealType& fermionicSigns,
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

	void externalProductReduced(SizeType i,
	                            const BasisType& basis2,
	                            const BasisType& basis3,
	                            bool option,
	                            const OperatorType& A)
	{
		reducedOpImpl_.externalProduct(i,basis2,basis3,option,A);
	}

	template<typename ApplyFactorsType>
	void outerProductHamiltonian(const SparseMatrixType& h2,
	                             const SparseMatrixType& h3,
	                             ApplyFactorsType& apply)
	{
		SparseMatrixType tmpMatrix;
		assert(h2.row()==h2.col());
		VectorRealType ones(h2.row(),1.0);
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

	//const   VectorSizeType& electrons() const {return operatorsImpl_.electrons(); }

	void print(int ind= -1) const
	{
		if (!useSu2Symmetry_) {
			if (ind<0)
				for (SizeType i=0;i<operators_.size();i++) std::cerr<<operators_[i];
			else std::cerr<<operators_[ind];
		} else {
			reducedOpImpl_.print(ind);
		}
	}

	template<typename IoOutputter>
	void save(IoOutputter& io,const PsimagLite::String& s) const
	{
		if (!useSu2Symmetry_) io.printVector(operators_,"#OPERATORS");
		else reducedOpImpl_.save(io,s);
		io.printMatrix(hamiltonian_,"#HAMILTONIAN");
	}

	template<typename IoOutputter>
	void saveEmpty(IoOutputter& io,const PsimagLite::String& s) const
	{
		PsimagLite::Vector<SizeType>::Type tmp;
		if (!useSu2Symmetry_) io.printVector(tmp,"#OPERATORS");
		else reducedOpImpl_.saveEmpty(io,s);
		PsimagLite::Matrix<SizeType> tmp2(0,0);
		io.printMatrix(tmp2,"#HAMILTONIAN");
	}

	SizeType size() const { return operators_.size(); }

private:

	void reorder(SparseMatrixType &v,const   VectorSizeType& permutation)
	{
		if (v.row() == 0 || v.col() == 0) {
			assert(v.row() == 0 && v.col() == 0);
			return;
		}

		SparseMatrixType matrixTmp;

		permute(matrixTmp,v,permutation);
		permuteInverse(v,matrixTmp,permutation);
	}

	void announceChangeAll() const
	{
#ifdef OPERATORS_CHANGE_ALL
		static bool flag = false;
		if (flag) return;
		PsimagLite::String msg("OPERATORS_CHANGE_ALL in use: ");
		msg += "GeometryMaxConnections value might not be used\n";
		std::cerr<<msg;
		std::cout<<msg;
		flag = true;
#endif
	}

	bool useSu2Symmetry_;
	ReducedOperatorsType reducedOpImpl_;
	typename PsimagLite::Vector<OperatorType>::Type operators_;
	SparseMatrixType hamiltonian_;
	PsimagLite::ProgressIndicator progress_;
}; //class Operators
} // namespace Dmrg

/*@}*/
#endif

