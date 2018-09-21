/*
Copyright (c) 2009-2018, UT-Battelle, LLC
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

/*! \file Basis.h
 *
 *   A class to represent in a light way a Dmrg basis (used only to implement symmetries).
 *
 */
#ifndef BASIS_HEADER_H
#define BASIS_HEADER_H
#include "Utils.h"
#include "Sort.h" // in PsimagLite
#include "HamiltonianSymmetryLocal.h"
#include "HamiltonianSymmetrySu2.h"
#include "ProgressIndicator.h"
#include "Qn.h"
#include "NotReallySort.h"

namespace Dmrg {
// A class to represent in a light way a Dmrg basis (used only to implement symmetries).
// (See corresponding section in paper)
template<typename SparseMatrixType_>
class Basis {

	typedef typename SparseMatrixType_::value_type SparseElementType;
	typedef typename PsimagLite::Real<SparseElementType>::Type RealType_;
	typedef Basis<SparseMatrixType_> ThisType;
	typedef HamiltonianSymmetryLocal<SparseMatrixType_>  HamiltonianSymmetryLocalType;

public:

	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef VectorSizeType BlockType;
	typedef SparseMatrixType_ SparseMatrixType;
	typedef RealType_ RealType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef Qn QnType;
	typedef HamiltonianSymmetrySu2<SparseMatrixType_, QnType> HamiltonianSymmetrySu2Type;
	typedef typename HamiltonianSymmetrySu2Type::FactorsType FactorsType;
	typedef typename HamiltonianSymmetrySu2Type::PairType PairType;
	typedef typename QnType::VectorQnType VectorQnType;

	enum {SAVE_ALL, SAVE_PARTIAL};

	//! Constructor, s=name of this basis
	Basis(const PsimagLite::String& s)
	    : dmrgTransformed_(false), name_(s), progress_(s)
	{}

	//! Loads this basis from memory or disk
	template<typename IoInputter>
	Basis(IoInputter& io,
	      const PsimagLite::String& ss,
	      bool minimizeRead)
	    : dmrgTransformed_(false), name_(ss), progress_(ss)
	{
		correctNameIfNeeded();
		PsimagLite::String prefix =  ss + "/";
		loadInternal(io, prefix, minimizeRead);
	}

	//! Loads this basis from memory or disk
	template<typename IoInputter>
	void read(IoInputter& io, PsimagLite::String prefix)
	{
		name_ = "";
		loadInternal(io, prefix, false);
	}

	//! Returns the name of this basis
	const PsimagLite::String& name() const { return name_; }

	const VectorBoolType& signs() const { return signs_; }

	//! Sets the block of sites for this basis
	void set(BlockType const &B) { block_ = B; }

	/* PSIDOC BasisSetToProduct
		The quantum numbers of the original (untransformed) real-space basis
		are set by the model class (to be described in Section~\ref{subsec:models}),
		whereas the quantum numbers of outer products are handled
		by the class \cppClass{Basis} and \cppClass{BasisImplementation},
		function \cppFunction{setToProduct}. This can be done because if $|a\rangle$
		has quantum number $q_a$ and $|b\rangle$ has quantum number
		$q_b$, then $|a\rangle\otimes|b\rangle$ has quantum number
		$q_a+q_b$.  \cppClass{Basis} knows how quantum
		numbers change when we change the basis: they
		do not change since the DMRG transformation
		preserves quantum numbers; and  \cppClass{Basis} also
		knows what happens to quantum numbers when we truncate the basis:
		quantum numbers of discarded states are discarded.
		In this way, symmetries are implemented efficiently,
		with minimal dependencies and in a model-independent way.
		*/
	void setToProduct(const ThisType& basis1,
	                  const ThisType& basis2,
	                  const QnType* pseudoQn = 0)
	{
		block_.clear();
		utils::blockUnion(block_,basis1.block_,basis2.block_);
		VectorQnType qns;

		if (useSu2Symmetry_) {
			std::cout<<"Basis: SU(2) Symmetry is in use\n";

			symmSu2_.setToProduct(basis1.symmSu2_,
			                      basis2.symmSu2_,
			                      pseudoQn,
			                      basis1.qns_,
			                      basis2.qns_,
			                      qns_,
			                      qns);
		} else {
			SizeType ns = basis2.size();
			SizeType ne = basis1.size();

			unsigned long int check = ns*ne;
			unsigned int shift = 8*sizeof(SizeType)-1;
			unsigned long int max = 1;
			max <<= shift;
			if (check >= max) {
				PsimagLite::String msg("Basis::setToProduct: Basis too large. ");
				msg += "Current= "+ ttos(check) + " max " + ttos(max) + " ";
				msg += "Please recompile with -DUSE_LONG\n";
				throw PsimagLite::RuntimeError(msg);
			}

			SizeType npe = basis2.partition_.size();
			if (npe > 0) --npe;
			SizeType nps = basis1.partition_.size();
			if (nps > 0) --nps;

			SizeType total = basis1.size() * basis2.size();
			qns.reserve(total);
			signs_.clear(); // reserve isn't affected
			signs_.reserve(total);
			for (SizeType pe = 0; pe < npe; ++pe) {
				for (SizeType i = basis2.partition_[pe]; i < basis2.partition_[pe + 1]; ++i) {
					for (SizeType ps = 0; ps < nps; ++ps) {
						for (SizeType j = basis1.partition_[ps];
						     j < basis1.partition_[ps + 1];
						     ++j) {
							qns.push_back(QnType(basis2.qns_[pe], basis1.qns_[ps]));
							signs_.push_back(basis1.signs_[j] ^ basis2.signs_[i]);
						}
					}
				}
			}
		}

		bool notSuper = (basis1.block().size() == 1 || basis2.block().size() == 1);
		ProgramGlobals::VerboseEnum verbose = (notSuper) ? ProgramGlobals::VERBOSE_NO :
		                                                   ProgramGlobals::VERBOSE_YES;

		// order quantum numbers of combined basis:
		findPermutationAndPartitionAndQns(qns, true, false, verbose);
		reorder();
		signsOld_ = signs_;
	}

	//! returns the effective quantum number of basis state i
	const QnType& qnEx(SizeType i) const
	{
		assert(i < qns_.size());
		return qns_[i];
	}

	//! returns the start of basis partition i (see paper)
	SizeType partition(SizeType i) const
	{
		assert(i < partition_.size());
		return partition_[i];
	}

	//! returns number of partitions for this basis (see paper)
	SizeType partition() const { return partition_.size(); }

	//! returns the permutation of i
	SizeType permutation(SizeType i) const
	{
		assert(i<permutationVector_.size());
		return  permutationVector_[i];
	}

	//! Return the permutation vector
	const VectorSizeType& permutationVector() const
	{
		return  permutationVector_;
	}

	//! returns the inverse permutation of i
	int permutationInverse(SizeType i) const
	{
		assert(i<permInverse_.size());
		return permInverse_[i];
	}

	//! returns the inverse permutation vector
	const VectorSizeType& permutationInverse() const
	{
		return permInverse_;
	}

	//! returns the block of sites over which this basis is built
	const BlockType& block() const { return block_; }

	//! returns the size of this basis
	SizeType size() const
	{
		assert(partition_.size() > 0);
		return partition_[partition_.size()-1];
	}

	//! finds the partition that contains basis state i
	SizeType findPartitionNumber(SizeType i) const
	{
		for (SizeType j=0;j<partition_.size()-1;j++)
			if (i>=partition_[j] && i<partition_[j+1]) return j;

		throw PsimagLite::RuntimeError("Basis: No partition found for this state\n");
	}

	//! removes the indices contained in removedIndices and
	//! transforms this basis by transform
	//! eigs cannot be const because su(2) needs to sort them due to being obtained in blocks
	template<typename SolverParametersType>
	void changeBasis(VectorSizeType& removedIndices,
	                 typename PsimagLite::Vector<RealType>::Type& eigs,
	                 SizeType kept,
	                 const SolverParametersType& solverParams)
	{
		removedIndices.clear();
		if (useSu2Symmetry_)
			symmSu2_.calcRemovedIndices(removedIndices,eigs,kept,solverParams);
		else
			symmLocal_.calcRemovedIndices(removedIndices,eigs,kept,solverParams);

		if (removedIndices.size()==0) return;

		VectorSizeType perm(removedIndices.size());
		PsimagLite::Sort<VectorSizeType > sort;
		sort.sort(removedIndices,perm);
	}

	RealType truncateBasis(const typename PsimagLite::Vector<RealType>::Type& eigs,
	                       const VectorSizeType& removedIndices)
	{
		dmrgTransformed_=true;

		if (removedIndices.size()==0) return 0;

		PsimagLite::OstringStream msg2;
		msg2<<"Truncating indices...";
		progress_.printline(msg2,std::cout);

		VectorQnType qns;
		unShrinkVector(qns, qns_, partition_);
		truncate(qns, removedIndices);

		// N.B.: false below means that we don't truncate the permutation vectors
		//	because they're needed for the WFT
		findPermutationAndPartitionAndQns(qns, false, false, ProgramGlobals::VERBOSE_NO);

		PsimagLite::OstringStream msg;
		msg<<"Done with changeBasis";
		progress_.printline(msg,std::cout);
		return calcError(eigs,removedIndices);
	}

	//! Returns the factors that mix this basis
	//! If not using SU(2) this is trivial
	const FactorsType* getFactors() const
	{
		return (useSu2Symmetry_) ? &(symmSu2_.getFactors()) : 0;
	}

	//! returns the flavor of state i of this basis
	SizeType getFlavor(SizeType i) const
	{
		if (useSu2Symmetry_) return symmSu2_.getFlavor(i);
		else return  symmLocal_.getFlavor(i);
	}

	SizeType flavor2Index(SizeType f1,
	                      SizeType f2,
	                      SizeType ne1,
	                      SizeType ne2,
	                      SizeType j1,
	                      SizeType j2) const
	{
		assert(useSu2Symmetry_);
		return symmSu2_.flavor2Index(f1,f2,ne1,ne2,j1,j2);
	}

	void flavor2Index(PsimagLite::Map<SizeType,SizeType>::Type& mm,
	                  const PairType& jm) const
	{
		assert(useSu2Symmetry_);
		symmSu2_.flavor2Index(mm,jm);
	}

	const VectorSizeType& flavorsOld() const
	{
		assert(useSu2Symmetry_);
		return symmSu2_.flavorsOld();
	}

	const VectorBoolType& oldSigns() const { return signsOld_; }

	//! Returns the fermionic sign for state i
	int fermionicSign(SizeType i, int f) const
	{
		assert(i < signs_.size());
		return (signs_[i]) ? f : 1;
	}

	//! Returns the (j,m) for state i of this basis
	PairType jmValue(SizeType i) const
	{
		if (!useSu2Symmetry_) return PairType(0,0);
		return symmSu2_.jmValue(i);
	}

	//! Returns true if using SU(2) symmetry or false otherwise
	static bool useSu2Symmetry()  { return useSu2Symmetry_; }

	//! Tells this basis to use SU(2) symmetry or not
	static void useSu2Symmetry(bool flag)  { useSu2Symmetry_=flag; }

	//! Returns true if this basis has been DMRG transformed, or false if it hasn't
	bool dmrgTransformed() const { return dmrgTransformed_; }

	//! Returns the reduced (by the Wigner Eckart theorem) index corresponding to state i
	SizeType reducedIndex(SizeType i) const
	{
		assert(useSu2Symmetry_);
		return symmSu2_.reducedIndex(i);
	}

	//! Returns the size of this basis when reduced (by the Wigner-Eckart theorem)
	SizeType reducedSize() const
	{
		assert(useSu2Symmetry_);
		return symmSu2_.reducedSize();
	}

	//! Returns the i-th distinct j value for this basis
	SizeType jVals(SizeType i) const
	{
		assert(useSu2Symmetry_);
		return symmSu2_.jVals(i);
	}

	//! Returns the number of distinct j values for this basis
	SizeType jVals() const
	{
		assert(useSu2Symmetry_);
		return symmSu2_.jVals();
	}

	//! Returns the maximum j value in this basis
	SizeType jMax() const
	{
		//assert(useSu2Symmetry_);
		return symmSu2_.jMax();
	}

	PsimagLite::String pseudoQnToString(SizeType i) const
	{
		return "unimplemented";
	}

	QnType pseudoQn(SizeType i) const
	{
		if (useSu2Symmetry_) {
			assert(i < partition_.size());
			SizeType ind = partition_[i];
			PairType jmPair(symmSu2_.jmValue(ind).first, 0);
			assert(ind < signs_.size());
			QnType q(signs_[ind], VectorSizeType(), jmPair, 0);
			return q;
		} else {
			return qnEx(i);
		}
	}

	void su2ElectronsBridge(VectorSizeType &v) const
	{
		QnType::su2ElectronsBridge(v, qns_);
	}

	//! saves this basis to disk
	template<typename SomeIoType>
	void write(SomeIoType& io,
	           const PsimagLite::String& ss,
	           typename SomeIoType::Serializer::WriteMode mode,
	           bool minimizeWrite) const
	{
		checkSigns();
		PsimagLite::String label = ss + "/";
		io.createGroup(ss);
		io.write(useSu2Symmetry_, label + "useSu2Symmetry");
		io.write(block_, label + "BLOCK");

		if (!minimizeWrite) {
			io.write(signs_, label + "signs_");
			io.write(signsOld_, label + "SignsOld");
		}

		io.write(partition_, label + "PARTITION");
		io.write(permInverse_, label + "PERMUTATIONINVERSE");
		io.write(qns_, label + "QNShrink");

		if (useSu2Symmetry_) symmSu2_.write(io, label);
		else symmLocal_.write(io, label);
	}

	//! saves this basis to disk
	template<typename SomeIoType>
	void write(SomeIoType& io,
	           typename SomeIoType::Serializer::WriteMode mode,
	           PsimagLite::String prefix,
	           bool minimizeWrite,
	           typename PsimagLite::EnableIf<
	           PsimagLite::IsOutputLike<SomeIoType>::True, int>::Type = 0) const
	{
		write(io, prefix + "/" + name_, mode, minimizeWrite);
	}

	//! The operator<< is a friend
	friend std::ostream& operator<<(std::ostream& os,
	                                const Basis<SparseMatrixType>& x)
	{
		os<<"dmrgTransformed="<<x.dmrgTransformed_<<"\n";
		os<<"name="<<x.name_<<"\n";
		os<<"quantumNumbers\n";
		os<<x.quantumNumbers_;
		os<<"electrons\n";
		os<<x.electrons_;
		os<<"partition\n";
		os<<x.partition_;
		os<<"permutation\n";
		os<<x.permutationVector_;
		os<<"block\n";
		os<<x.block_;
		return os;
	}

	friend std::istream& operator>>(std::istream& is,
	                                Basis<SparseMatrixType>&)
	{
		throw PsimagLite::RuntimeError("Unimplemented >>");
	}

protected:

	//! Sets symmetry information for this basis, see SymmetryElectronsSz.h for more
	void setSymmetryRelated(const VectorQnType& basisData)
	{
		if (useSu2Symmetry_) symmSu2_.set(basisData);

		SizeType n = basisData.size();
		signs_.resize(n);
		for (SizeType i = 0; i < n; ++i)
			signs_[i] = basisData[i].oddElectrons;

		VectorQnType basisData2 = basisData;
		if (!useSu2Symmetry()) flattenQns(basisData2);

		findPermutationAndPartitionAndQns(basisData2,
		                                  true,
		                                  true,
		                                  ProgramGlobals::VERBOSE_NO);
		reorder();
		signsOld_ = signs_;
	}

private:

	template<typename IoInputter>
	void loadInternal(IoInputter& io,
	                  PsimagLite::String prefix,
	                  bool minimizeRead,
	                  typename PsimagLite::EnableIf<
	                  PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	{
		useSu2Symmetry_=false;
		prefix += "/";
		io.read(useSu2Symmetry_, prefix + "useSu2Symmetry");
		io.read(block_, prefix + "BLOCK");

		if (!minimizeRead) {
			io.read(signs_, prefix + "signs_");
			io.read(signsOld_, prefix + "SignsOld");
		}

		io.read(partition_, prefix + "PARTITION");
		io.read(permInverse_, prefix + "PERMUTATIONINVERSE");
		permutationVector_.resize(permInverse_.size());
		for (SizeType i=0;i<permInverse_.size();i++)
			permutationVector_[permInverse_[i]]=i;

		QnType::readVector(qns_, prefix + "QNShrink", io);
		if (!minimizeRead) checkSigns();
		dmrgTransformed_=false;
		if (useSu2Symmetry_)
			symmSu2_.read(io,prefix, minimizeRead);
		else
			symmLocal_.read(io,prefix, minimizeRead);
	}


	void checkSigns() const
	{
#ifndef NDEBUG
		if (!Qn::ifPresentOther0IsElectrons || Qn::modalStruct.size() == 0)
			return;

		SizeType n = partition_.size();
		assert(n > 0);
		--n;
		assert(partition_[n] == signs_.size());
		for (SizeType p = 0; p < n; ++p) {
			SizeType start = partition_[p];
			SizeType end = partition_[p + 1];
			SizeType expected = (qns_[p].other[0] & 1);
			for (SizeType i = start; i < end; ++i) {
				if (signs_[i] != expected)
					err("Unexpected sign\n");
			}
		}
#endif
	}

	static void flattenQns(VectorQnType& qns)
	{
		SizeType n = qns.size();
		for (SizeType i = 0; i < n; ++i) {
			qns[i].flavors = 0;
			qns[i].jmPair.first = qns[i].jmPair.second = 0;
		}
	}

	void unShrinkVector(VectorQnType& dest,
	                    const VectorQnType& src,
	                    const VectorSizeType& partition) const
	{
		SizeType n = partition.size();
		assert(n > 0);
		assert(src.size() == n -1);

		dest.resize(partition[n - 1], QnType::zero());
		for (SizeType i = 0; i < n - 1; ++i) {
			SizeType start = partition[i];
			SizeType end = partition[i + 1];
			assert(end < 1 + dest.size());
			for (SizeType j = start; j < end; ++j)
				dest[j] = src[i];
		}
	}

	RealType calcError(const typename PsimagLite::Vector<RealType>::Type& eigs,
	                   const VectorSizeType& removedIndices) const
	{
		RealType sum=static_cast<RealType>(0.0);
		for (SizeType i=0;i<eigs.size();i++)
			if (PsimagLite::indexOrMinusOne(removedIndices,i) < 0) sum += eigs[i];
		return 1.0 - sum;
	}

	void truncate(VectorQnType& qns, const VectorSizeType& removedIndices)
	{
		utils::truncateVector(qns, removedIndices);
		utils::truncateVector(signs_,removedIndices);
		if (useSu2Symmetry_) {
			VectorSizeType electrons;
			QnType::su2ElectronsBridge(electrons, qns);
			symmSu2_.truncate(removedIndices, electrons);
		}
	}

	void reorder()
	{
		utils::reorder(signs_,permutationVector_);
		if (useSu2Symmetry_) symmSu2_.reorder(permutationVector_);
	}

	void findPermutationAndPartitionAndQns(const VectorQnType& qns,
	                                       bool changePermutation,
	                                       bool doNotSort,
	                                       ProgramGlobals::VerboseEnum verbose)
	{
		SizeType n = qns.size();

		VectorSizeType numbers(n);
		for (SizeType i = 0; i < n; ++i) numbers[i] = i;
		VectorSizeType permutationVector;
		NotReallySort notReallySort;
		notReallySort(permutationVector,
		              qns_,
		              partition_,
		              numbers,
		              qns,
		              doNotSort,
		              verbose);

		if (changePermutation) {
			permutationVector_ = (useSu2Symmetry_) ? numbers : permutationVector;
			permInverse_.resize(permutationVector_.size());
			for (SizeType i=0;i<permInverse_.size();i++)
				permInverse_[permutationVector_[i]]=i;
		}
	}

	void correctNameIfNeeded()
	{
		if (name_.find("/") == PsimagLite::String::npos)
			return;

		if (name_.find("system") != PsimagLite::String::npos)
			name_ = "system";
		else
			name_ = "environ";
	}

	/* PSIDOC BasisQuantumNumbers
		Symmetries will allow the solver to block the Hamiltonian matrix in blocks,
using less memory, speeding up
		the computation and allowing the code to parallelize matrix blocks related by symmetry.
		Let us assume that our particular model has $N_s$ symmetries labeled by $0\le \alpha < N_s$.
		Therefore, each element $k$  of the basis has $N_s$ associated ``good'' quantum numbers
		 $\tilde{q}_{k,\alpha}$. These quantum numbers can refer to practically anything,
		 for example, to number of particles with a given spin or orbital or to the $z$
component of the spin.
		We do not need to know the details to block the matrix. We know, however, that
these numbers are
		finite, and let $Q$ be an integer such that $\tilde{q}_{k,\alpha}< Q$ $\forall k,\alpha$.
		We can then combine all these quantum numbers into a single one,
		like this: $q_k = \sum_\alpha \tilde{q}_{k,\alpha} Q^\alpha$,
		and this mapping is bijective. In essence, we combined all ``good''
		quantum numbers into a single one and from now on we
		will consider that we have only one Hamiltonian symmetry called the
		``effective'' symmetry, and only one corresponding number $q_k$, the
		``effective'' quantum number. These numbers are stored in the  member
		{\it quantumNumbers} of C++ class \cppClass{Basis}.
		(Note that if one has 100 sites or less,\footnote{This is probably a
		maximum for systems of correlated electrons such as the Hubbard model
		or the t-J model.} then the number $Q$ defined above is probably of the
		order of hundreds for usual symmetries, making this implementation very practical for
		systems of correlated electrons.)
		*/
	VectorQnType qns_;
	VectorBoolType signs_;
	VectorBoolType signsOld_;

	/* PSIDOC BasisPartition
		What remains to be done is to find a partition of the basis which
		labels where the quantum number changes. Let us say that the
		quantum numbers of the reordered basis states are
		\[
		\{3,3,3,3,8,8,9,9,9,15,\cdots\}.
		\]
		Then we define a vector named ``partition'', such that partition[0]=0,
		partition[1]=4, because the quantum number changes in the 4th position
		(from 3 to 8), and then partition[2]=6, because the quantum number
		changes again (from 8 to 9) in the 6th position, etc.
		Now we know that our Hamiltonian matrix will be composed first of a
		block of 4x4, then of a block of 2x2, etc.
		*/
	VectorSizeType partition_;

	/* PSIDOC BasisPermutationVector
		We then reorder our basis such that its elements are given in
		increasing $q$ number. There will be a permutation vector associated
		with this reordering, that will be stored in the member
		\verb!permutationVector! of class \cppClass{Basis}.
		For ease of coding we also store its inverse in \verb!permInverse!.
		*/
	VectorSizeType permutationVector_;
	VectorSizeType permInverse_;
	HamiltonianSymmetryLocalType symmLocal_;
	HamiltonianSymmetrySu2Type symmSu2_;
	/* PSIDOC BasisBlock
		The variable block of a \cppClass{DmrgBasis} object indicates over
		which sites the basis represented by this object is being built.
		*/
	BlockType block_;
	bool dmrgTransformed_;
	PsimagLite::String name_;
	PsimagLite::ProgressIndicator progress_;
	static bool useSu2Symmetry_;

}; // class Basis

template<typename SparseMatrixType>
bool Basis<SparseMatrixType>::useSu2Symmetry_=false;

template<typename SparseMatrixType_>
struct IsBasisType<Basis<SparseMatrixType_> > {
	enum {True = true};
};

} // namespace Dmrg

/*@}*/
#endif

