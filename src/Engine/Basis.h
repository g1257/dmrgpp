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

/*! \file Basis.h
 *
 *   A class to represent in a light way a Dmrg basis (used only to implement symmetries).
 * 
 */
#ifndef BASIS_HEADER_H
#define BASIS_HEADER_H

#include "BasisImplementation.h"

namespace Dmrg {
	//! A class to represent in a light way a Dmrg basis (used only to implement symmetries).
	//! (See corresponding section in paper)
	template<typename RealType_,typename SparseMatrixType>
	class	Basis {
		typedef  Basis<RealType_,SparseMatrixType> ThisType;
		typedef BasisImplementation<RealType_,SparseMatrixType> BasisImplementationType;
	public:
		typedef  typename BasisImplementationType::FactorsType FactorsType;
		typedef RealType_ RealType;

		//! Structure that contains information about the basis 
		//! see BasisData.h for more info
		typedef typename BasisImplementationType::BasisDataType BasisDataType;

		//! A vector of sites
		typedef typename BasisImplementationType::BlockType BlockType;

		//! Constant used by the WaveFunctionTransformation
		static size_t const BEFORE_TRANSFORM =  BasisImplementationType::BEFORE_TRANSFORM;

		//! Contant used by the WaveFunctionTransformation
		static size_t const AFTER_TRANSFORM =  BasisImplementationType::AFTER_TRANSFORM;

		//! Constructor, s=name of this basis 
		Basis(const std::string& s) : basisImplementation_(s) {  }

		//! Loads this basis from memory or disk
		template<typename IoInputter>
		Basis(IoInputter& io,const std::string& ss,size_t counter=0,bool bogus = false)
		: basisImplementation_(io,ss,counter,bogus)
		{}

		//! Loads this basis from memory or disk
		template<typename IoInputter>
		void load(IoInputter& io)
		{
			basisImplementation_.load(io);
		}

		//! Returns the name of this basis
		const std::string& name() const { return  basisImplementation_.name(); }

		//! Sets the block of sites for this basis
		void set(BlockType const &B) { basisImplementation_.set(B); }

		//! Sets symmetry information for this basis, see BasisData.h for more
		void setSymmetryRelated(BasisDataType const &q) { basisImplementation_.setSymmetryRelated(q); }

		//! sets this basis to the outer product of basis2 basis3
		void setToProduct(ThisType const &basis2,ThisType const &basis3,int q =  -1)
		{
			basisImplementation_.setToProduct(basis2.basisImplementation_,basis3.basisImplementation_,q);
		}

		//! returns the effective quantum number of basis state i
		int qn(int i,size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			return basisImplementation_.qn(i,beforeOrAfterTransform);
		}

		//! Returns the partition that corresponds to quantum number qn
		int partitionFromQn(size_t qn,size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			return basisImplementation_.partitionFromQn(qn,beforeOrAfterTransform);
		}

		//! returns the start of basis partition i (see paper)
		size_t partition(size_t i,size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			return basisImplementation_.partition(i,beforeOrAfterTransform);
		}

		//! returns number of partitions for this basis (see paper)
		size_t partition() const { return basisImplementation_.partition(); }

		//! returns the permutation of i 
		size_t permutation(size_t i) const { return  basisImplementation_.permutation(i); }

		//! Return the permutation vector
		const std::vector<size_t>& permutationVector() const
		{
			return  basisImplementation_.permutationVector();
		}

		//! returns the inverse permutation of i 
		int permutationInverse(int i) const { return basisImplementation_.permutationInverse(i); }

		//! returns the inverse permutation vector
		const std::vector<size_t>& permutationInverse() const { return basisImplementation_.permutationInverse(); }

		//! returns the block of sites over which this basis is built
		const BlockType& block() const { return basisImplementation_.block(); }

		//! returns the size of this basis
		size_t size() const { return basisImplementation_.size(); }

		//! finds the partition that contains basis state i
		size_t findPartitionNumber(size_t i) const { return basisImplementation_.findPartitionNumber(i); }

		//! Encodes the quantum numbers into a single unique size_t and returns it
		static size_t encodeQuantumNumber(const std::vector<size_t>& quantumNumbers)
		{
			return BasisImplementationType::encodeQuantumNumber(quantumNumbers);
		}

		//! Inverse for encodeQuantumNumber
		static std::vector<size_t> decodeQuantumNumber(int q)
		{
			return BasisImplementationType::decodeQuantumNumber(q);
		}

		//! Encodes (flavor,jvalue,density) into a unique number and returns it
		static size_t pseudoQuantumNumber(const std::vector<size_t>& targets)
		{
			return BasisImplementationType::pseudoQuantumNumber(targets);
		}

		//! Inverse of pseudoQuantumNumber
		size_t pseudoEffectiveNumber(size_t i) const
		{
			return basisImplementation_.pseudoEffectiveNumber(i);
		}

		//! Given the information in the structure bdt, calculates the quantum numbers in q
		static void findQuantumNumbers(std::vector<size_t>& q,const BasisDataType& bdt)
		{
			BasisImplementationType::findQuantumNumbers(q,bdt);
		}

		//! removes the indices contained in removedIndices and
		//! transforms this basis by transform 
		template<typename BlockMatrixType,typename SolverParametersType>
		RealType changeBasis(typename BlockMatrixType::BuildingBlockType& ftransform,
		                     BlockMatrixType &transform,
		                     std::vector<RealType>& eigs,
		                     size_t kept,
		                     const SolverParametersType& solverParams)
		{
			return basisImplementation_.changeBasis(ftransform,transform,eigs,kept,solverParams);
		}

		//! Finds a partition of the basis given the effecitve quantum numbers
		void findPartition() { return basisImplementation_.findPartition(); }

		//! Returns the factors that mix this basis 
		//! If not using SU(2) this is trivial
		const FactorsType& getFactors() const { return basisImplementation_.getFactors(); }

		//! returns the number of electrons for state i of this basis
		size_t electrons(size_t i) const { return basisImplementation_.getNe(i); }

		//! returns the flavor of state i of this basis 
		size_t getFlavor(size_t i) const { return basisImplementation_.getFlavor(i); }

		//! Given the double triplet (f1f2,ne1ne2,j1j2) returns a unique size_t
		size_t flavor2Index(size_t f1,size_t f2,size_t ne1,size_t ne2,size_t j1,size_t j2) const
		{
			return basisImplementation_.flavor2Index(f1,f2,ne1,ne2,j1,j2);
		}

		//! 
		void flavor2Index(std::map<size_t,size_t>& mm,const std::pair<size_t,size_t>& jm) const
		{
			basisImplementation_.flavor2Index(mm,jm);
		}

		const std::vector<size_t>& flavorsOld() const
		{
			return basisImplementation_.flavorsOld();
		}

		//! Returns the vector of electrons for this basis
		const std::vector<size_t>& electronsVector(size_t beforeOrAfterTransform=AFTER_TRANSFORM) const
		{
			return basisImplementation_.electronsVector(beforeOrAfterTransform);
		}

		//! Returns the fermionic sign for state i
		int fermionicSign(size_t i,int f) const { return basisImplementation_.fermionicSign(i,f); }

		//! Returns the (j,m) for state i of this basis
		typename BasisImplementationType::PairType jmValue(size_t i) const { return basisImplementation_.jmValue(i); }

		//! Returns true if using SU(2) symmetry or false otherwise
		static bool useSu2Symmetry() { return BasisImplementationType::useSu2Symmetry(); }

		//! Tells this basis to use SU(2) symmetry or not
		static void useSu2Symmetry(bool tf) {  BasisImplementationType::useSu2Symmetry(tf); }

		//! Returns true if this basis has been DMRG transformed, or false if it hasn't
		bool dmrgTransformed() const { return basisImplementation_.dmrgTransformed(); }

		//! Returns the reduced (by the Wigner Eckart theorem) index corresponding to state i of this basis
		size_t reducedIndex(size_t i) const
		{
			return basisImplementation_.reducedIndex(i);
		}

		//! Returns the size of this basis when reduced (by the Wigner-Eckart theorem)
		size_t reducedSize() const
		{
			return basisImplementation_.reducedSize();
		}

		//! Returns the i-th distinct j value for this basis
		size_t jVals(size_t i) const
		{
			return basisImplementation_.jVals(i);
		}

		//! Returns the number of distinct j values for this basis
		size_t jVals() const
		{
			return basisImplementation_.jVals();
		}

		//! Returns the maximum j value in this basis
		size_t jMax() const
		{
			return basisImplementation_.jMax();
		}

		//! saves this basis to disk
		template<typename IoOutputter>
		void save(IoOutputter& io,const std::string& s) const
		{
			basisImplementation_.save(io,s);
		}

		//! saves this basis to disk
		template<typename IoOutputter>
		void save(IoOutputter& io) const
		{
			basisImplementation_.save(io);
		}

		//! The operator<< is a friend
		template<typename RealType2,typename SparseMatrixType2>
		friend std::ostream& operator<<(std::ostream& os,const Basis<RealType2,SparseMatrixType2>& x);

		template<typename RealType2,typename SparseMatrixType2>
		friend std::istream& operator>>(std::istream& is,Basis<RealType2,SparseMatrixType2>& x);

	private:
		BasisImplementationType basisImplementation_;
	}; // class Basis

	template<typename RealType,typename SparseMatrixType>
	std::ostream& operator<<(std::ostream& os,const Basis<RealType,SparseMatrixType>& x)
	{
		os<<x.basisImplementation_;
		return os;
	}

	template<typename RealType,typename SparseMatrixType>
	std::istream& operator>>(std::istream& is,Basis<RealType,SparseMatrixType>& x)
	{
		is>>x.basisImplementation_;
		return is;
	}

} // namespace Dmrg

/*@}*/
#endif

