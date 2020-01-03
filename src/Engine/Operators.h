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

#include <cassert>
#include "ProgressIndicator.h"
#include "Complex.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "BlockOffDiagMatrix.h"
#include "ChangeOfBasis.h"
#include "Operator.h"
#include "Matrix.h"

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

public:

	typedef std::pair<SizeType,SizeType> PairType;
	typedef BasisType_ BasisType;
	typedef Operators<BasisType_> ThisType;
	typedef typename BasisType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef Operator<SparseElementType> OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef PsimagLite::Matrix<SparseElementType> DenseMatrixType;
	typedef ChangeOfBasis<OperatorStorageType, DenseMatrixType> ChangeOfBasisType;
	typedef typename OperatorType::StorageType StorageType;
	typedef typename StorageType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType,SizeType> PairSizeSizeType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef typename ChangeOfBasisType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;

	// law of the excluded middle went out the window here:
	enum class ChangeAllEnum { UNSET, TRUE_SET, FALSE_SET};

	class MyLoop {

	public:

		MyLoop(VectorOperatorType& operators,
		       VectorOperatorType& superOps,
		       ChangeOfBasisType& changeOfBasis,
		       const BlockDiagonalMatrixType& ftransform1,
		       const PairSizeSizeType& startEnd)
		    : operators_(operators),
		      superOps_(superOps),
		      changeOfBasis_(changeOfBasis),
		      ftransform(ftransform1),
		      startEnd_(startEnd)
		{
			changeOfBasis.update(ftransform);
		}

		void doTask(SizeType taskNumber, SizeType)
		{
			const SizeType nLocals =  operators_.size();
			if (taskNumber < nLocals)
				doTaskForLocal(taskNumber);
			else
				doTaskForSuper(taskNumber - nLocals);
		}

		SizeType tasks() const
		{
			return operators_.size() + superOps_.size();
		}

	private:

		void doTaskForLocal(SizeType k)
		{
			if (isLocalExcluded(k) && k < operators_.size()) {
				operators_[k].clear();
				return;
			}

			changeOfBasis_(operators_[k].getStorageNonConst());
		}

		void doTaskForSuper(SizeType k)
		{
			if (isSuperExcluded(k) && k < superOps_.size()) {
				superOps_[k].clear();
				return;
			}

			changeOfBasis_(superOps_[k].getStorageNonConst());
		}

		bool isLocalExcluded(SizeType k) const
		{
			if (changeAll_ == ChangeAllEnum::TRUE_SET)
				return false; // <-- this is the safest answer

			if (k < startEnd_.first || k >= startEnd_.second) return true;
			return false;
		}

		bool isSuperExcluded(SizeType) const
		{
			throw PsimagLite::RuntimeError("Operators.h: isSuperExcluded not written yet\n");
		}

		VectorOperatorType& operators_;
		VectorOperatorType& superOps_;
		ChangeOfBasisType& changeOfBasis_;
		const BlockDiagonalMatrixType& ftransform;
		const PairSizeSizeType& startEnd_;
	};

	Operators() : progress_("Operators")
	{
		if (changeAll_ == ChangeAllEnum::UNSET)
			changeAll_ = ChangeAllEnum::FALSE_SET;
	}

	template<typename IoInputter>
	Operators(IoInputter& io, PsimagLite::String prefix, bool isObserveCode)
	    : progress_("Operators")
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
		const SizeType last = prefix.length() - 1;
		if (last >= prefix.length())
			err("Operators.h: read\n");

		if (prefix[last] != '/') prefix += "/";

		io.read(operators_, prefix + "Operators");
		//io.read(superOps_, prefix + "SuperOperators");
		io.read(hamiltonian_, prefix + "Hamiltonian");
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

	void setLocal(const VectorOperatorType& ops)
	{
		operators_ = ops;
	}

	const OperatorType& getLocalByIndex(SizeType i) const
	{
		assert(i < operators_.size());
		return operators_[i];
	}

	SizeType sizeOfLocal() const
	{
		return operators_.size();
	}

	void changeBasis(const BlockDiagonalMatrixType& ftransform,
	                 const PairSizeSizeType& startEnd,
	                 bool blasIsThreadSafe)
	{
		typedef PsimagLite::Parallelizer<MyLoop> ParallelizerType;
		SizeType threads = (blasIsThreadSafe) ? PsimagLite::Concurrency::
		                                        codeSectionParams.npthreads : 1;
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerType threadObject(codeSectionParams);

		MyLoop helper(operators_, superOps_, changeOfBasis_, ftransform, startEnd);

		threadObject.loopCreate(helper); // FIXME: needs weights

		hamiltonian_.checkValidity();
		ChangeOfBasisType::changeBasis(hamiltonian_, ftransform);
	}

	template<typename SomeSuperOperatorHelperType>
	void setToProduct(const BasisType& basis1,
	                  const ThisType& ops1,
	                  const BasisType& basis2,
	                  const ThisType& ops2,
	                  const VectorSizeType& permutationInverse,
	                  const SomeSuperOperatorHelperType& someSuperOpHelper)
	{
		setToProductLocal(basis1, ops1, basis2, ops2, permutationInverse);
		setToProductSuper(basis1, ops1, basis2, ops2, permutationInverse, someSuperOpHelper);
	}

	void outerProductHamiltonian(const StorageType& h2,
	                             const StorageType& h3,
	                             const VectorSizeType& permutationFull)
	{
		StorageType tmpMatrix;
		assert(h2.rows()==h2.cols());
		VectorRealType ones(h2.rows(),1.0);
		externalProduct2(hamiltonian_,h2,h3.rows(),ones,true, permutationFull);

		externalProduct2(tmpMatrix,h3,h2.rows(),ones,false, permutationFull);

		hamiltonian_ += tmpMatrix;
	}

	void setHamiltonian(StorageType const &h)
	{
		hamiltonian_ = h;
	}

	void setHamiltonian(const SparseMatrixType& h)
	{
		fromCRS(hamiltonian_, h);
	}

	const StorageType& hamiltonian() const { return hamiltonian_; }

	void print(int ind= -1) const
	{
		if (ind<0)
			for (SizeType i=0;i<operators_.size();i++) std::cerr<<operators_[i];
		else std::cerr<<operators_[ind];
	}

	template<typename SomeIoOutType>
	void overwrite(SomeIoOutType& io,
	               const PsimagLite::String& s,
	               typename PsimagLite::EnableIf<
	               PsimagLite::IsOutputLike<SomeIoOutType>::True, int*>::Type = 0) const
	{
		io.overwrite(operators_, s + "/Operators");
		//		 io.overwrite(superOps_, s + "/SuperOperators");
		io.overwrite(hamiltonian_, s + "/Hamiltonian");
	}

	void write(PsimagLite::IoNg::Out& io,
	           const PsimagLite::String& s,
	           PsimagLite::IoNgSerializer::WriteMode mode) const
	{
		if (mode == PsimagLite::IoNgSerializer::ALLOW_OVERWRITE) {
			io.overwrite(operators_, s + "/Operators");
			//			io.overwrite(superOps_, s + "/SuperOperators");
			io.overwrite(hamiltonian_, s + "/Hamiltonian");
		} else {
			io.write(operators_, s + "/Operators");
			//			io.write(superOps_, s + "/SuperOperators");
			io.write(hamiltonian_, s + "/Hamiltonian");
		}
	}

	void clear()
	{
		operators_.clear();
		superOps_.clear();
		hamiltonian_.clear();
	}

	const OperatorType& getSuperByIndex(SizeType ind) const
	{
		assert(ind < superOps_.size());
		return superOps_[ind];
	}

	SizeType superIndices(const VectorSizeType&, SizeType) const
	{
		PsimagLite::String msg(__FILE__);
		throw PsimagLite::RuntimeError(msg + "::superOperatorIndices() not implemented yet\n");
	}

private:

	void setToProductLocal(const BasisType& basis2,
	                       const ThisType& ops2,
		                   const BasisType& basis3,
		                   const ThisType& ops3,
	                       const VectorSizeType& permutationInverse)
	{
		typename PsimagLite::Vector<RealType>::Type fermionicSigns;
		SizeType nlocalOps = ops2.sizeOfLocal() + ops3.sizeOfLocal();
		operators_.resize(nlocalOps);
		ProgramGlobals::FermionOrBosonEnum savedSign = ProgramGlobals::FermionOrBosonEnum::BOSON;

		for (SizeType i = 0; i < nlocalOps; ++i) {
			if (i < ops2.sizeOfLocal()) {
				const OperatorType& myOp = ops2.getLocalByIndex(i);
				bool isFermion = (myOp.fermionOrBoson() ==
				                  ProgramGlobals::FermionOrBosonEnum::FERMION);
				if (savedSign != myOp.fermionOrBoson() || fermionicSigns.size() == 0) {
					utils::fillFermionicSigns(fermionicSigns,
					                          basis2.signs(),
					                          (isFermion) ? -1 : 1);
					savedSign = myOp.fermionOrBoson();
				}

				crossProductForLocal(i,
				                     myOp,
				                     basis3.size(),
				                     fermionicSigns,
				                     true,
				                     permutationInverse);

			} else {
				const OperatorType& myOp = ops3.getLocalByIndex(i - ops2.sizeOfLocal());

				bool isFermion = (myOp.fermionOrBoson() ==
				                  ProgramGlobals::FermionOrBosonEnum::FERMION);

				if (savedSign != myOp.fermionOrBoson() || fermionicSigns.size() == 0) {
					utils::fillFermionicSigns(fermionicSigns,
					                          basis2.signs(),
					                          (isFermion) ? -1 : 1);
					savedSign = myOp.fermionOrBoson();
				}

				crossProductForLocal(i,
				                     myOp,
				                     basis2.size(),
				                     fermionicSigns,
				                     false,
				                     permutationInverse);
			}
		}
	}

	template<typename SomeSuperOperatorHelperType>
	void setToProductSuper(const BasisType& basis2,
	                       const ThisType& ops2,
		                   const BasisType& basis3,
		                   const ThisType& ops3,
	                       const VectorSizeType& permutationInverse,
	                       const SomeSuperOperatorHelperType& someSuperOpHelper)
	{
		typename PsimagLite::Vector<RealType>::Type fermionicSigns;
		SizeType nSuperOps = someSuperOpHelper.size();
		superOps_.resize(nSuperOps);
		ProgramGlobals::FermionOrBosonEnum savedSign = ProgramGlobals::FermionOrBosonEnum::BOSON;
		typedef typename SomeSuperOperatorHelperType::PairBoolSizeType PairBoolSizeType;
		const bool option = (basis3.block().size() == 1);
		for (SizeType i = 0; i < nSuperOps; ++i) {
			const PairBoolSizeType op2Index  = someSuperOpHelper.leftOperatorIndex(i);
			const PairBoolSizeType op3Index = someSuperOpHelper.rightOperatorIndex(i);
			const OperatorType& op1 = (!op2Index.first) ? ops2.getLocalByIndex(op2Index.second)
			                                                   : ops2.getSuperByIndex(op2Index.
			                                                                          second);
			const OperatorType& op3 = (!op3Index.first) ? ops3.getLocalByIndex(op3Index.second)
			                                                   : ops3.getSuperByIndex(op3Index.
			                                                                          second);
			bool isFermion = (op3.fermionOrBoson() == ProgramGlobals::FermionOrBosonEnum::FERMION);

			if (savedSign != op3.fermionOrBoson() || fermionicSigns.size() == 0) {
				utils::fillFermionicSigns(fermionicSigns, basis2.signs(), (isFermion) ? -1 : 1);
				savedSign = op3.fermionOrBoson();
			}


			superOps_[i].outerProduct(op1, op3, fermionicSigns, option, permutationInverse);

		}
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
	void crossProductForLocal(SizeType i,
	                          const OperatorType& m,
	                          int x,
	                          const VectorRealType& fermionicSigns,
	                          bool option,
	                          const VectorSizeType& permutationFull)
	{
		assert(!BasisType::useSu2Symmetry());
		operators_[i].outerProduct(m,
		                           x,
		                           fermionicSigns,
		                           option,
		                           permutationFull);
		// don't forget to set fermion sign and j:
		operators_[i].set(m.fermionOrBoson(), m.jm(), m.angularFactor());
		// apply(operators_[i]);
	}

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
	ChangeOfBasisType changeOfBasis_;
	VectorOperatorType operators_;
	VectorOperatorType superOps_;
	StorageType hamiltonian_;
	PsimagLite::ProgressIndicator progress_;
}; //class Operators

template<typename T>
typename Operators<T>::ChangeAllEnum Operators<T>::changeAll_ =
        Operators<T>::ChangeAllEnum::UNSET;

} // namespace Dmrg

/*@}*/
#endif

