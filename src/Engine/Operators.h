/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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
#include "ApplyFactors.h"

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
	typedef typename OperatorType::StorageType StorageType;
	typedef typename StorageType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType,SizeType> PairSizeSizeType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef typename BasisType::FactorsType FactorsType;

	// law of the excluded middle went out the window here:
	enum class ChangeAllEnum { UNSET, TRUE_SET, FALSE_SET};

	class MyLoop {

	public:

		MyLoop(ReducedOperatorsType& reducedOpImpl,
		       typename PsimagLite::Vector<OperatorType>::Type& operators,
		       const BlockDiagonalMatrixType& ftransform1,
		       const BasisType* thisBasis1,
		       const PairSizeSizeType& startEnd)
		    : reducedOpImpl_(reducedOpImpl),
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

			if (!BasisType::useSu2Symmetry())
				reducedOpImpl_.changeBasis(operators_[k].data);
			else
				reducedOpImpl_.changeBasis(k);
		}

		SizeType tasks() const
		{
			if (BasisType::useSu2Symmetry()) return reducedOpImpl_.size();
			return operators_.size();
		}

		void gather()
		{
			if (ConcurrencyType::isMpiDisabled("Operators")) return;

			if (!BasisType::useSu2Symmetry()) {
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
			if (changeAll_ == ChangeAllEnum::TRUE_SET)
				return false; // <-- this is the safest answer

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
				bcast(operators_[i]);
		}

		ReducedOperatorsType& reducedOpImpl_;
		typename PsimagLite::Vector<OperatorType>::Type& operators_;
		const BlockDiagonalMatrixType& ftransform;
		const BasisType* thisBasis;
		bool hasMpi_;
		const PairSizeSizeType& startEnd_;
	};

	Operators(const BasisType* thisBasis)
	    : reducedOpImpl_(thisBasis),
	      progress_("Operators")
	{
		if (changeAll_ == ChangeAllEnum::UNSET)
			changeAll_ = ChangeAllEnum::FALSE_SET;
	}

	template<typename IoInputter>
	Operators(IoInputter& io,
	          PsimagLite::String prefix,
	          SizeType level,
	          const BasisType* thisBasis,
	          bool isObserveCode)
	    : reducedOpImpl_(io,level,thisBasis),
	      progress_("Operators")
	{
		if (changeAll_ == ChangeAllEnum::UNSET)
			changeAll_ = ChangeAllEnum::FALSE_SET;

		if (isObserveCode) return;

		read(io, prefix, false);
	}

	template<typename IoInputter>
	void read(IoInputter& io,
	          PsimagLite::String prefix,
	          bool roi = true, // it is false only when called from constructor
	          typename PsimagLite::EnableIf<
	          PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	{
		prefix += "/";

		if (!BasisType::useSu2Symmetry()) {
			io.read(operators_, prefix + "Operators");
		} else {
			if (roi) reducedOpImpl_.read(io);
		}

		io.read(hamiltonian_, prefix + "Hamiltonian");
		reducedOpImpl_.setHamiltonian(hamiltonian_);
	}

	static void setChangeAll(bool flag)
	{
		if (!flag) {
			changeAll_ = ChangeAllEnum::FALSE_SET;
			printChangeAll();
			return;
		}

		assert(flag);

		if (changeAll_ == ChangeAllEnum::UNSET) {
			changeAll_ = ChangeAllEnum::TRUE_SET;
			printChangeAll();
			return;
		}

		err("Operators::setChangeAll(true) called to late\n");
	}

	void setOperators(const typename PsimagLite::Vector<OperatorType>::Type& ops)
	{
		if (!BasisType::useSu2Symmetry()) operators_=ops;
		else reducedOpImpl_.setOperators(ops);
	}

	const OperatorType& getReducedOperatorByIndex(char modifier,const PairType& p) const
	{
		assert(BasisType::useSu2Symmetry());
		return reducedOpImpl_.getReducedOperatorByIndex(modifier,p);
	}

	const OperatorType& getOperatorByIndex(int i) const
	{
		assert(!BasisType::useSu2Symmetry());
		assert(i>=0 && SizeType(i)<operators_.size());
		return operators_[i];
	}

	const OperatorType& getReducedOperatorByIndex(int i) const
	{
		assert(BasisType::useSu2Symmetry());
		return reducedOpImpl_.getReducedOperatorByIndex(i);
	}

	SizeType numberOfOperators() const
	{
		if (BasisType::useSu2Symmetry()) return reducedOpImpl_.size();
		return operators_.size();
	}

	void changeBasis(const BlockDiagonalMatrixType& ftransform,
	                 const BasisType* thisBasis,
	                 const PairSizeSizeType& startEnd)
	{
		typedef PsimagLite::Parallelizer<MyLoop> ParallelizerType;
		ParallelizerType threadObject(PsimagLite::Concurrency::codeSectionParams);

		MyLoop helper(reducedOpImpl_,operators_,ftransform,thisBasis,startEnd);

		threadObject.loopCreate(helper); // FIXME: needs weights

		helper.gather();

		reducedOpImpl_.changeBasisHamiltonian(hamiltonian_, ftransform);
	}

	void reorder(const VectorSizeType& permutation)
	{
		for (SizeType k=0;k<numberOfOperators();k++) {
			if (!BasisType::useSu2Symmetry())
				reorder2(operators_[k].data, permutation);
			reducedOpImpl_.reorder(k,permutation);
		}

		reorder2(hamiltonian_,permutation);
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
		if (!BasisType::useSu2Symmetry())
			operators_.resize(x);
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
	void externalProduct(SizeType i,
	                     const OperatorType& m,
	                     int x,
	                     const VectorRealType& fermionicSigns,
	                     bool option,
	                     ApplyFactors<FactorsType>& apply)
	{
		assert(!BasisType::useSu2Symmetry());
		externalProduct2(operators_[i].data,m.data,x,fermionicSigns,option);
		// don't forget to set fermion sign and j:
		operators_[i].fermionOrBoson=m.fermionOrBoson;
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
	void outerProductHamiltonian(const StorageType& h2,
	                             const StorageType& h3,
	                             ApplyFactorsType& apply)
	{
		StorageType tmpMatrix;
		assert(h2.rows()==h2.cols());
		VectorRealType ones(h2.rows(),1.0);
		externalProduct2(hamiltonian_,h2,h3.rows(),ones,true);

		externalProduct2(tmpMatrix,h3,h2.rows(),ones,false);

		hamiltonian_ += tmpMatrix;

		apply(hamiltonian_);
	}

	void outerProductHamiltonianReduced(const BasisType& basis2,
	                                    const BasisType& basis3,
	                                    const StorageType& h2,
	                                    const StorageType& h3)
	{
		reducedOpImpl_.outerProductHamiltonian(basis2,basis3,h2,h3);
	}

	void setHamiltonian(StorageType const &h)
	{
		hamiltonian_ = h;
		reducedOpImpl_.setHamiltonian(h);
	}

	void setHamiltonian(const SparseMatrixType& h)
	{
		fromCRS(hamiltonian_, h);
		reducedOpImpl_.setHamiltonian(hamiltonian_);
	}

	const StorageType& hamiltonian() const { return hamiltonian_; }

	const StorageType& reducedHamiltonian() const
	{
		return reducedOpImpl_.hamiltonian();
	}

	//const   VectorSizeType& electrons() const {return operatorsImpl_.electrons(); }

	void print(int ind= -1) const
	{
		if (!BasisType::useSu2Symmetry()) {
			if (ind<0)
				for (SizeType i=0;i<operators_.size();i++) std::cerr<<operators_[i];
			else std::cerr<<operators_[ind];
		} else {
			reducedOpImpl_.print(ind);
		}
	}

	template<typename SomeIoOutType>
	void overwrite(SomeIoOutType& io,
	               const PsimagLite::String& s,
	               typename PsimagLite::EnableIf<
	               PsimagLite::IsOutputLike<SomeIoOutType>::True, int*>::Type = 0) const
	{
		if (!BasisType::useSu2Symmetry())
			io.overwrite(operators_, s + "/Operators");
		else
			reducedOpImpl_.overwrite(io,s);

		io.overwrite(hamiltonian_, s + "/Hamiltonian");
	}

	void write(PsimagLite::IoNg::Out& io,
	           const PsimagLite::String& s,
	           PsimagLite::IoNgSerializer::WriteMode mode) const
	{
		if (!BasisType::useSu2Symmetry()) {
			if (mode == PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
				io.overwrite(operators_, s + "/Operators");
			else
				io.write(operators_, s + "/Operators");
		} else {
			reducedOpImpl_.write(io, s, mode);
		}

		if (mode == PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
			io.overwrite(hamiltonian_, s + "/Hamiltonian");
		else
			io.write(hamiltonian_, s + "/Hamiltonian");
	}

	SizeType size() const { return operators_.size(); }

	void clear()
	{
		reducedOpImpl_.clear();
		operators_.clear();
		hamiltonian_.clear();
	}

private:

	static void printChangeAll()
	{
		PsimagLite::String msg("INFO: Operators::changeAll_=");
		msg += toString(changeAll_) + "\n";
		if (changeAll_ == ChangeAllEnum::TRUE_SET)
			msg += "GeometryMaxConnections value might not be used\n";
		std::cerr<<msg;
		std::cout<<msg;
	}

	static PsimagLite::String toString(ChangeAllEnum value)
	{
		switch (value) {
		case ChangeAllEnum::UNSET:
			return "UNSET";
			break;
		case ChangeAllEnum::FALSE_SET:
			return "FALSE";
			break;
		case ChangeAllEnum::TRUE_SET:
		default:
			return "TRUE";
			break;
		}
	}

	static ChangeAllEnum changeAll_;
	ReducedOperatorsType reducedOpImpl_;
	typename PsimagLite::Vector<OperatorType>::Type operators_;
	StorageType hamiltonian_;
	PsimagLite::ProgressIndicator progress_;
}; //class Operators

template<typename T>
typename Operators<T>::ChangeAllEnum Operators<T>::changeAll_ =
        Operators<T>::ChangeAllEnum::UNSET;

} // namespace Dmrg

/*@}*/
#endif

