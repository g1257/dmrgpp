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

namespace Dmrg {
	//! A class to represent in a light way a Dmrg basis (used only to implement symmetries).
	//! (See corresponding section in paper)
	template<typename RealType_,typename SparseMatrixType>
	class	Basis {

		typedef  Basis<RealType_,SparseMatrixType> ThisType;
		typedef HamiltonianSymmetryLocal<RealType_,SparseMatrixType>  HamiltonianSymmetryLocalType;
		typedef HamiltonianSymmetrySu2<RealType_,SparseMatrixType>  HamiltonianSymmetrySu2Type;
//		typedef Reflection ReflectionSymmetryType;

	public:
		typedef typename HamiltonianSymmetrySu2Type::FactorsType FactorsType;
		typedef typename HamiltonianSymmetrySu2Type::PairType PairType;
		typedef  BasisData<PairType> BasisDataType;
		typedef typename PsimagLite::Vector<size_t>::Type BlockType;

		enum {BEFORE_TRANSFORM,AFTER_TRANSFORM};

		typedef RealType_ RealType;

		//! Constructor, s=name of this basis 
		Basis(const PsimagLite::String& s)
		: dmrgTransformed_(false), name_(s), progress_(s,0)
		{
			symmLocal_.createDummyFactors(1,1);
		}

		//! Loads this basis from memory or disk
		template<typename IoInputter>
		Basis(IoInputter& io,const PsimagLite::String& ss,size_t counter=0,bool bogus = false)
		: dmrgTransformed_(false), name_(ss), progress_(ss,0)
		{
			io.advance("#NAME="+ss,counter);
			loadInternal(io);
		}

		//! Loads this basis from memory or disk
		template<typename IoInputter>
		void load(IoInputter& io)
		{
			PsimagLite::String nn="#NAME=";
			std::pair<PsimagLite::String,size_t> sc = io.advance(nn);
			name_ = sc.first.substr(nn.size(),sc.first.size());
			loadInternal(io);
		}

		//! Returns the name of this basis
		const PsimagLite::String& name() const { return name_; }

		//! Sets the block of sites for this basis
		void set(BlockType const &B) { block_ = B; }

		//! Sets symmetry information for this basis, see BasisData.h for more
		void setSymmetryRelated(const BasisDataType& basisData)
		{
			if (useSu2Symmetry_) symmSu2_.set(basisData);
			electrons_.resize(basisData.electronsUp.size());
			for (size_t i=0;i<basisData.electronsUp.size();i++)
				electrons_[i]=basisData.electronsUp[i]+basisData.electronsDown[i];
			findQuantumNumbers(quantumNumbers_,basisData);
			findPermutationAndPartition();
			electronsOld_=electrons_;
		}


		/**
		The quantum numbers of the original (untransformed) real-space basis
		are set by the model class (to be described in Section~\\ref{subsec:models}),
		whereas the quantum numbers of outer products are handled
		by the class \\cppClass{Basis} and \\cppClass{BasisImplementation},
		function \cppFunction{setToProduct}. This can be done because if $|a\\rangle$
		has quantum number $q_a$ and $|b\\rangle$ has quantum number
		$q_b$, then $|a\\rangle\\otimes|b\\rangle$ has quantum number
		$q_a+q_b$.  \\cppClass{!PTEX_THISCLASS} knows how quantum
		numbers change when we change the basis: they
		do not change since the DMRG transformation
		preserves quantum numbers; and  \\cppClass{!PTEX_THISCLASS} also
		knows what happens to quantum numbers when we truncate the basis:
		quantum numbers of discarded states are discarded.
		In this way, symmetries are implemented efficiently,
		with minimal dependencies and in a model-independent way.
		*/
		void setToProduct(const ThisType& su2Symmetry2,const ThisType& su2Symmetry3,int pseudoQn = -1)
		{
			block_.clear();
			utils::blockUnion(block_,su2Symmetry2.block_,su2Symmetry3.block_); //! B= pS.block Union X

			if (useSu2Symmetry_) {
				symmSu2_.setToProduct(su2Symmetry2.symmSu2_,su2Symmetry3.symmSu2_,pseudoQn,
						su2Symmetry2.electrons_,su2Symmetry3.electrons_,electrons_,quantumNumbers_);
			} else {
				size_t ns = su2Symmetry2.size();
				size_t ne = su2Symmetry3.size();

				quantumNumbers_.clear();
				electrons_.clear();

				for (size_t j=0;j<ne;j++) for (size_t i=0;i<ns;i++) {
					quantumNumbers_.push_back(su2Symmetry2.quantumNumbers_[i]+su2Symmetry3.quantumNumbers_[j]);
					electrons_.push_back(su2Symmetry2.electrons(i)+su2Symmetry3.electrons(j));
				}

				symmLocal_.createDummyFactors(ns,ne);
			}
			// order quantum numbers of combined basis:
			findPermutationAndPartition();

			reorder();
			electronsOld_ = electrons_;
		}

		//! returns the effective quantum number of basis state i
		int qn(int i,size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			if (beforeOrAfterTransform==AFTER_TRANSFORM)
				return quantumNumbers_[i];
			return quantumNumbersOld_[i];
		}

		//! Returns the partition that corresponds to quantum number qn
		int partitionFromQn(size_t qn,size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			const typename PsimagLite::Vector<size_t>::Type *quantumNumbers, *partition;

			if (beforeOrAfterTransform==AFTER_TRANSFORM) {
				quantumNumbers = &quantumNumbers_;
				partition = &partition_;
			} else {
				quantumNumbers = &quantumNumbersOld_;
				partition = &partitionOld_;
			}

			for (size_t i=0;i<partition->size();i++) {
				size_t state = (*partition)[i];
				if ((*quantumNumbers)[state]==qn) return i;
			}
			return -1;
		}

		//! returns the start of basis partition i (see paper)
		size_t partition(size_t i,size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			if (beforeOrAfterTransform==AFTER_TRANSFORM)
				return partition_[i];
			return partitionOld_[i];
		}

		//! returns number of partitions for this basis (see paper)
		size_t partition() const { return partition_.size(); }

		//! returns the permutation of i 
		size_t permutation(size_t i) const
		{
			assert(i<permutationVector_.size());
			return  permutationVector_[i];
		}

		//! Return the permutation vector
		const typename PsimagLite::Vector<size_t>::Type& permutationVector() const
		{
			return  permutationVector_;
		}

		//! returns the inverse permutation of i 
		int permutationInverse(size_t i) const
		{
			assert(i<permInverse_.size());
			return permInverse_[i];
		}

		//! returns the inverse permutation vector
		const typename PsimagLite::Vector<size_t>::Type& permutationInverse() const { return permInverse_; }

		//! returns the block of sites over which this basis is built
		const BlockType& block() const { return block_; }

		//! returns the size of this basis
		size_t size() const { return quantumNumbers_.size(); }

		//! finds the partition that contains basis state i
		size_t findPartitionNumber(size_t i) const
		{
			for (size_t j=0;j<partition_.size()-1;j++)
				if (i>=partition_[j] && i<partition_[j+1]) return j;
			throw PsimagLite::RuntimeError("BasisImplementation:: No partition found for this state\n");
		}

		//! Encodes the quantum numbers into a single unique size_t and returns it
		static size_t encodeQuantumNumber(const typename PsimagLite::Vector<size_t>::Type& quantumNumbers)
		{
			if (useSu2Symmetry_) return HamiltonianSymmetrySu2Type::encodeQuantumNumber(quantumNumbers);
			else return HamiltonianSymmetryLocalType::encodeQuantumNumber(quantumNumbers);
		}

		//! Inverse for encodeQuantumNumber
		static typename PsimagLite::Vector<size_t>::Type decodeQuantumNumber(int q)
		{
			if (useSu2Symmetry_) return HamiltonianSymmetrySu2Type::decodeQuantumNumber(q);
			else return HamiltonianSymmetryLocalType::decodeQuantumNumber(q);
		}

		//! Encodes (flavor,jvalue,density) into a unique number and returns it
		static size_t pseudoQuantumNumber(const typename PsimagLite::Vector<size_t>::Type& targets)
		{
			if (useSu2Symmetry_)
				return HamiltonianSymmetrySu2Type::pseudoQuantumNumber(targets);
			else
				return HamiltonianSymmetryLocalType::pseudoQuantumNumber(targets);
		}

		//! Inverse of pseudoQuantumNumber
		size_t pseudoEffectiveNumber(size_t i) const
		{
			if (useSu2Symmetry_) return symmSu2_.pseudoEffectiveNumber(electrons_[i],symmSu2_.jmValue(i).first);
			else return quantumNumbers_[i];
		}

		//! Given the information in the structure bdt, calculates the quantum numbers in q
		static void findQuantumNumbers(typename PsimagLite::Vector<size_t>::Type& qn,const BasisDataType& basisData)
		{
			if (useSu2Symmetry_) HamiltonianSymmetrySu2Type::findQuantumNumbers(qn,basisData);
			else HamiltonianSymmetryLocalType::findQuantumNumbers(qn,basisData);
		}

		//! removes the indices contained in removedIndices and
		//! transforms this basis by transform
		//! eigs cannot be const because su(2) needs to sort them due to being obtained in blocks
		template<typename SolverParametersType>
		void changeBasis(typename PsimagLite::Vector<size_t>::Type& removedIndices,
				 typename PsimagLite::Vector<RealType>::Type& eigs,
				 size_t kept,
				 const SolverParametersType& solverParams)
		{
			removedIndices.clear();
			if (useSu2Symmetry_) symmSu2_.calcRemovedIndices(removedIndices,eigs,kept,solverParams);
			else symmLocal_.calcRemovedIndices(removedIndices,eigs,kept,solverParams);

			if (removedIndices.size()==0) return;

			typename PsimagLite::Vector<size_t>::Type perm(removedIndices.size());
			PsimagLite::Sort<typename PsimagLite::Vector<size_t>::Type > sort;
			sort.sort(removedIndices,perm);
		}

		template<typename BlockMatrixType>
		RealType truncateBasis(SparseMatrixType& ftransform,
				       const BlockMatrixType& transform,
				       const typename PsimagLite::Vector<RealType>::Type& eigs,
				       const typename PsimagLite::Vector<size_t>::Type& removedIndices)
		{
			quantumNumbersOld_ = quantumNumbers_;
			partitionOld_ = partition_;
			dmrgTransformed_=true;

			blockMatrixToSparseMatrix(ftransform,transform);

			if (removedIndices.size()==0) return 0;

			PsimagLite::OstringStream msg0;
			msg0<<"Truncating transform...";
			utils::truncate(ftransform,removedIndices,false);
			progress_.printline(msg0,std::cerr);

			PsimagLite::OstringStream msg2;
			msg2<<"Truncating indices...";
			progress_.printline(msg2,std::cerr);
			truncate(removedIndices);

			// N.B.: false below means that we don't truncate the permutation vectors
			//	because they're needed for the WFT
			findPermutationAndPartition(false);
			PsimagLite::OstringStream msg;
			msg<<"Done with changeBasis";
			progress_.printline(msg,std::cerr);
			return calcError(eigs,removedIndices);
		}

		//! Finds a partition of the basis given the effecitve quantum numbers
		//! Find a partition of the basis given the effecitve quantum numbers
		//! (see section about Symmetries in paper)
		void findPartition()
		{
			size_t qtmp = quantumNumbers_[0]+1;
			partition_.clear();
			for (size_t i=0;i<size();i++) {
				if (quantumNumbers_[i]!=qtmp) {
					partition_.push_back(i);
					qtmp = quantumNumbers_[i];
				}
			}
			partition_.push_back(size());
		}

		//! Returns the factors that mix this basis 
		//! If not using SU(2) this is trivial
		const FactorsType& getFactors() const
		{
			if (useSu2Symmetry_) return symmSu2_.getFactors();
			else return symmLocal_.getFactors();
		}

		//! returns the number of electrons for state i of this basis
		size_t electrons(size_t i) const { return electrons_[i];  }

		//! returns the flavor of state i of this basis 
		size_t getFlavor(size_t i) const
		{
			if (useSu2Symmetry_) return symmSu2_.getFlavor(i);
			else return  symmLocal_.getFlavor(i);
		}

		size_t flavor2Index(size_t f1,size_t f2,size_t ne1,size_t ne2,size_t j1,size_t j2) const
		{
			assert(useSu2Symmetry_);
			return symmSu2_.flavor2Index(f1,f2,ne1,ne2,j1,j2);
		}

		void flavor2Index(PsimagLite::Map<size_t,size_t>::Type& mm,const PairType& jm) const
		{
			assert(useSu2Symmetry_);
			symmSu2_.flavor2Index(mm,jm);
		}

		const typename PsimagLite::Vector<size_t>::Type& flavorsOld() const
		{
			assert(useSu2Symmetry_);
			return symmSu2_.flavorsOld();
		}

		//! Returns the vector of electrons for this basis
		const typename PsimagLite::Vector<size_t>::Type& electronsVector(size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			if (beforeOrAfterTransform == AFTER_TRANSFORM) return electrons_;
			return electronsOld_;
		}

		//! Returns the fermionic sign for state i
		int fermionicSign(size_t i,int f) const { return (electrons_[i]&1) ? f : 1; }

		//! Returns the (j,m) for state i of this basis
		PairType jmValue(size_t i) const
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

		//! Returns the reduced (by the Wigner Eckart theorem) index corresponding to state i of this basis
		size_t reducedIndex(size_t i) const
		{
			assert(useSu2Symmetry_);
			return symmSu2_.reducedIndex(i);
		}

		//! Returns the size of this basis when reduced (by the Wigner-Eckart theorem)
		size_t reducedSize() const
		{
			assert(useSu2Symmetry_);
			return symmSu2_.reducedSize();
		}

		//! Returns the i-th distinct j value for this basis
		size_t jVals(size_t i) const
		{
			assert(useSu2Symmetry_);
			return symmSu2_.jVals(i);
		}

		//! Returns the number of distinct j values for this basis
		size_t jVals() const
		{
			assert(useSu2Symmetry_);
			return symmSu2_.jVals();
		}

		//! Returns the maximum j value in this basis
		size_t jMax() const
		{
			//assert(useSu2Symmetry_);
			return symmSu2_.jMax();
		}

		//! saves this basis to disk
		template<typename IoOutputter>
		void save(IoOutputter& io,const PsimagLite::String& ss) const
		{
			io.printline("#NAME="+ss);
			saveInternal(io);
		}

		//! saves this basis to disk
		template<typename IoOutputter>
		void save(IoOutputter& io) const
		{
			//PsimagLite::OstringStream msg;
			//msg<<"Now saving to disk";
			//progress_.printline(msg,std::cerr);
			io.printline("#NAME="+name_);
			saveInternal(io);

		}

		//! The operator<< is a friend
		template<typename RealType2,typename SparseMatrixType2>
		friend std::ostream& operator<<(std::ostream& os,const Basis<RealType2,SparseMatrixType2>& x);

		template<typename RealType2,typename SparseMatrixType2>
		friend std::istream& operator>>(std::istream& is,Basis<RealType2,SparseMatrixType2>& x);

	private:

		template<typename IoInputter>
		void loadInternal(IoInputter& io)
		{
			int x=0;
			useSu2Symmetry_=false;
			io.readline(x,"#useSu2Symmetry=");
			if (x>0) useSu2Symmetry_=true;
			io.read(block_,"#BLOCK");
			io.read(quantumNumbers_,"#QN");
			io.read(electrons_,"#ELECTRONS");
			io.read(electronsOld_,"#0OLDELECTRONS");
			io.read(partition_,"#PARTITION");
			io.read(permInverse_,"#PERMUTATIONINVERSE");
			permutationVector_.resize(permInverse_.size());
			for (size_t i=0;i<permInverse_.size();i++) permutationVector_[permInverse_[i]]=i;
			dmrgTransformed_=false;
			if (useSu2Symmetry_)
				symmSu2_.load(io);
			else
				symmLocal_.load(io);
		}

		template<typename IoOutputter>
		void saveInternal(IoOutputter& io) const
		{
			PsimagLite::String s="#useSu2Symmetry="+ttos(useSu2Symmetry_);
			io.printline(s);
			io.printVector(block_,"#BLOCK");
			io.printVector(quantumNumbers_,"#QN");
			io.printVector(electrons_,"#ELECTRONS");
			io.printVector(electronsOld_,"#0OLDELECTRONS");
			io.printVector(partition_,"#PARTITION");
			io.printVector(permInverse_,"#PERMUTATIONINVERSE");

			if (useSu2Symmetry_) symmSu2_.save(io);
			else symmLocal_.save(io);
		}

		RealType calcError(const typename PsimagLite::Vector<RealType>::Type& eigs,
		                   const typename PsimagLite::Vector<size_t>::Type& removedIndices) const
		{
			RealType sum=static_cast<RealType>(0.0);
			for (size_t i=0;i<eigs.size();i++)
				if (PsimagLite::isInVector(removedIndices,i)<0) sum+=eigs[i];
			return 1.0-sum;
		}

		void truncate(const typename PsimagLite::Vector<size_t>::Type& removedIndices)
		{
			utils::truncateVector(quantumNumbers_,removedIndices);
			utils::truncateVector(electrons_,removedIndices);
			if (useSu2Symmetry_) symmSu2_.truncate(removedIndices,electrons_);
		}

		void reorder()
		{
			utils::reorder(electrons_,permutationVector_);
			if (useSu2Symmetry_) symmSu2_.reorder(permutationVector_);
		}

		void findPermutationAndPartition(bool changePermutation=true)
		{
			if (changePermutation) {
				permutationVector_.resize(size());
				if (useSu2Symmetry_) 	{
					//symmSu2_.orderFlavors(permutationVector_,partition_);
					for (size_t i=0;i<permutationVector_.size();i++)
						permutationVector_[i]=i;
				} else 	{
					PsimagLite::Sort<typename PsimagLite::Vector<size_t>::Type > sort;
					sort.sort(quantumNumbers_,permutationVector_);
				}
			}

			findPartition();

			if (changePermutation) {
				permInverse_.resize(permutationVector_.size());
				for (size_t i=0;i<permInverse_.size();i++)
					permInverse_[permutationVector_[i]]=i;

			}
		}

		/**
		Symmetries will allow the solver to block the Hamiltonian matrix in blocks, using less memory, speeding up
		the computation and allowing the code to parallelize matrix blocks related by symmetry.
		Let us assume that our particular model has $N_s$ symmetries labeled by $0\\le \\alpha < N_s$.
		Therefore, each element $k$  of the basis has $N_s$ associated ``good'' quantum numbers
		 $\\tilde{q}_{k,\\alpha}$. These quantum numbers can refer to practically anything,
		 for example, to number of particles with a given spin or orbital or to the $z$ component of the spin.
		We do not need to know the details to block the matrix. We know, however, that these numbers are
		finite, and let $Q$ be an integer such that $\\tilde{q}_{k,\\alpha}< Q$ $\\forall k,\\alpha$.
		We can then combine all these quantum numbers into a single one,
		like this: $q_k = \\sum_\\alpha \\tilde{q}_{k,\\alpha} Q^\\alpha$,
		and this mapping is bijective. In essence, we combined all ``good''
		quantum numbers into a single one and from now on we
		will consider that we have only one Hamiltonian symmetry called the
		``effective'' symmetry, and only one corresponding number $q_k$, the
		``effective'' quantum number. These numbers are stored in the  member
		{\\it quantumNumbers} of C++ class \\cppClass{!PTEX_THISCLASS}.
		(Note that if one has 100 sites or less,\\footnote{This is probably a
		maximum for systems of correlated electrons such as the Hubbard model
		or the t-J model.} then the number $Q$ defined above is probably of the
		order of hundreds for usual symmetries, making this implementation very practical for
		systems of correlated electrons.)
		*/
		typename PsimagLite::Vector<size_t>::Type quantumNumbers_;
		typename PsimagLite::Vector<size_t>::Type quantumNumbersOld_;
		typename PsimagLite::Vector<size_t>::Type electrons_;
		typename PsimagLite::Vector<size_t>::Type electronsOld_;

		/**
		What remains to be done is to find a partition of the basis which
		labels where the quantum number changes. Let us say that the
		quantum numbers of the reordered basis states are
		\\[
		\\{3,3,3,3,8,8,9,9,9,15,\\cdots\\}.
		\\]
		Then we define a vector named ``partition'', such that partition[0]=0,
		partition[1]=4, because the quantum number changes in the 4th position
		(from 3 to 8), and then partition[2]=6, because the quantum number
		changes again (from 8 to 9) in the 6th position, etc.
		Now we know that our Hamiltonian matrix will be composed first of a
		block of 4x4, then of a block of 2x2, etc.
		*/
		typename PsimagLite::Vector<size_t>::Type partition_;
		typename PsimagLite::Vector<size_t>::Type partitionOld_;

		/**
		We then reorder our basis such that its elements are given in
		increasing $q$ number. There will be a permutation vector associated
		with this reordering, that will be stored in the member
		\\verb!permutationVector! of class \\cppClass{!PTEX_THISCLASS}.
		For ease of coding we also store its inverse in \\verb!permInverse!.
		*/
		typename PsimagLite::Vector<size_t>::Type permutationVector_;
		typename PsimagLite::Vector<size_t>::Type permInverse_;
		HamiltonianSymmetryLocalType symmLocal_;
		HamiltonianSymmetrySu2Type symmSu2_;
		/**
		The variable block\_ of a \cppClass{DmrgBasis} object indicates over
		which sites the basis represented by this object is being built.
		*/
		BlockType block_;
		bool dmrgTransformed_;
		PsimagLite::String name_;
		PsimagLite::ProgressIndicator progress_;
//		ReflectionSymmetryType reflection_;
		static bool useSu2Symmetry_;

	}; // class Basis

	template<typename RealType,typename SparseMatrixType>
	std::ostream& operator<<(std::ostream& os,const Basis<RealType,SparseMatrixType>& x)
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

	template<typename RealType,typename SparseMatrixType>
	std::istream& operator>>(std::istream& is,Basis<RealType,SparseMatrixType>& x)
	{
		throw PsimagLite::RuntimeError("Unimplemented >>");
		return is;
	}

	template<typename RealType,typename SparseMatrixType>
	bool Basis<RealType,SparseMatrixType>::useSu2Symmetry_=false;

} // namespace Dmrg

/*@}*/
#endif

