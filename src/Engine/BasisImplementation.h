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

/*! \file BasisImplementation.h
 *
 *  
 *
 */
#ifndef BASIS_IMPL_H
#define BASIS_IMPL_H

#include "HamiltonianSymmetryLocal.h"
#include "HamiltonianSymmetrySu2.h"
#include "ProgressIndicator.h"

namespace Dmrg {
	
	
	//!
	template<typename RealType,typename SparseMatrixType>
	class BasisImplementation {
		typedef BasisImplementation<RealType,SparseMatrixType> ThisType;
		typedef HamiltonianSymmetryLocal<RealType,SparseMatrixType>  HamiltonianSymmetryLocalType;
		typedef HamiltonianSymmetrySu2<RealType,SparseMatrixType>  HamiltonianSymmetrySu2Type;

	public:
		typedef typename HamiltonianSymmetrySu2Type::FactorsType FactorsType;
		typedef typename HamiltonianSymmetrySu2Type::PairType PairType;
		typedef  BasisData<PairType> BasisDataType;
		typedef std::vector<int> BlockType;
		
		enum {BEFORE_TRANSFORM,AFTER_TRANSFORM};
		
		BasisImplementation(const std::string& s) : dmrgTransformed_(false), name_(s), progress_(s,0) 
		{
			symmLocal_.createDummyFactors(1,1);
		}
		
		// use this if you know the name
		template<typename IoInputter>
		BasisImplementation(IoInputter& io,const std::string& ss,size_t counter=0,bool bogus = false)
				: dmrgTransformed_(false), name_(ss), progress_(ss,0) 
		{
			io.advance("#NAME="+ss,counter);
			loadInternal(io);
		}

		// use this if you don't know the name
		/*template<typename IoInputter>
		BasisImplementation(IoInputter& io,size_t counter=0,bool bogus = false)
				: dmrgTransformed_(false), name_("#NAME"), progress_("#NAME",0)
		{
			std::string nn="#NAME=";
			std::pair<std::string,size_t> sc = io.advance(nn,counter);
			name_ = sc.first.substr(nn.size(),sc.first.size());
			loadInternal(io);
			
		}*/

		const std::string& name() const { return name_; }

		static bool useSu2Symmetry()  { return useSu2Symmetry_; }

		static void useSu2Symmetry(bool flag)  { useSu2Symmetry_=flag; }

		static size_t pseudoQuantumNumber(const std::vector<size_t>& targets)
		{
			if (useSu2Symmetry_) return HamiltonianSymmetrySu2Type::pseudoQuantumNumber(targets);
			else return HamiltonianSymmetryLocalType::pseudoQuantumNumber(targets);
		}

		size_t pseudoEffectiveNumber(size_t i) const
		{
			if (useSu2Symmetry_) return symmSu2_.pseudoEffectiveNumber(electrons_[i],symmSu2_.jmValue(i).first);
			else return quantumNumbers_[i];
		}

		size_t size() const { return quantumNumbers_.size(); }

		int qn(size_t i,size_t beforeOrAfterTransform) const 
		{
			if (beforeOrAfterTransform==AFTER_TRANSFORM)
				return quantumNumbers_[i];
			return quantumNumbersOld_[i];
		}

		int partitionFromQn(size_t qn,size_t beforeOrAfterTransform) const 
		{
			const std::vector<size_t> *quantumNumbers, *partition;
			
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

		static void findQuantumNumbers(std::vector<size_t>& qn,const BasisDataType& basisData) 
		{
			if (useSu2Symmetry_) HamiltonianSymmetrySu2Type::findQuantumNumbers(qn,basisData);
			else HamiltonianSymmetryLocalType::findQuantumNumbers(qn,basisData);
		}
		
		size_t findPartitionNumber(size_t i) const
		{
			for (size_t j=0;j<partition_.size()-1;j++) 
				if (i>=partition_[j] && i<partition_[j+1]) return j;
			throw std::runtime_error("BasisImplementation:: No partition found for this state\n");
		}

		const BlockType& block() const { return block_; }

		void setToProduct(const ThisType& su2Symmetry2,const ThisType& su2Symmetry3,int pseudoQn)
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
					electrons_.push_back(su2Symmetry2.getNe(i)+su2Symmetry3.getNe(j));
				}
				
				symmLocal_.createDummyFactors(ns,ne);
			}
			// order quantum numbers of combined basis:
			findPermutationAndPartition();

			reorder();
			electronsOld_ = electrons_;
		}

		size_t getFlavor(size_t i) const
		{
			if (useSu2Symmetry_) return symmSu2_.getFlavor(i);
			else return  symmLocal_.getFlavor(i);
		}

		size_t flavor2Index(size_t f1,size_t f2,size_t ne1,size_t ne2,size_t j1,size_t j2) const
		{
			if (!useSu2Symmetry_) std::cerr<<"flavor2Index: ERRRRRRRRRRRRRRROOOOOOOOOOOOOOR\n";
			return symmSu2_.flavor2Index(f1,f2,ne1,ne2,j1,j2);
		}

		void flavor2Index(std::map<size_t,size_t>& mm,const PairType& jm) const
		{
			if (!useSu2Symmetry_) std::cerr<<"flavor2Index: ERRRRRRRRRRRRRRROOOOOOOOOOOOOOR\n";
			symmSu2_.flavor2Index(mm,jm);
		}

		const std::vector<size_t>& flavorsOld() const
		{
			if (!useSu2Symmetry_) std::cerr<<"flavor2Index: ERRRRRRRRRRRRRRROOOOOOOOOOOOOOR\n";
			return symmSu2_.flavorsOld();
		}

		size_t getNe(size_t i) const { return electrons_[i]; }

		const std::vector<size_t>& electronsVector(size_t beforeOrAfterTransform) const 
		{
			if (beforeOrAfterTransform == AFTER_TRANSFORM) return electrons_;
			return electronsOld_;
		}

		int fermionicSign(size_t i,int f) const { return ((electrons_[i]%2)==0) ? 1 : f; }

		PairType jmValue(size_t i) const 
		{
			if (!useSu2Symmetry_) return PairType(0,0); 
			return symmSu2_.jmValue(i); 
		}

		template<typename BlockMatrixType,typename SolverParametersType>
		RealType changeBasis(typename BlockMatrixType::BuildingBlockType &ftransform,BlockMatrixType& transform,
				      std::vector<RealType>& eigs,size_t kept,const SolverParametersType& solverParams)
			{
			/* if (!isUnitary(transform)) {
				std::cerr<<"------------------------------------\n";
				operator<<(transform,std::cerr);
				throw std::runtime_error("BasisWithOperators::changeBasis(): transform is not unitary.\n");
		
			}*/
			quantumNumbersOld_ = quantumNumbers_;
			partitionOld_ = partition_;
			dmrgTransformed_=true;
			
			blockMatrixToFullMatrix(ftransform,transform);

			std::vector<size_t> removedIndices;
			if (useSu2Symmetry_) symmSu2_.calcRemovedIndices(removedIndices,eigs,kept,solverParams);
			else symmLocal_.calcRemovedIndices(removedIndices,eigs,kept,solverParams);

			if (removedIndices.size()>0) {
				std::vector<size_t> perm(removedIndices.size());
				utils::sort(removedIndices,perm);
				std::ostringstream msg;
				msg<<"Truncating transform...";
				psimag::truncate(ftransform,removedIndices,false);
				progress_.printline(msg,std::cerr);
				/*if (!isUnitary(ftransform)) { // only used for debugging
					//std::cerr<<"------------------------------------\n";
					//operator<<(transform,std::cerr);
					std::cerr<<"------------------------------------\n";
					std::cerr<<ftransform;
					throw std::runtime_error("BasisWithOperators::changeBasis(): transform is not unitary.\n");
		
				}*/
				std::ostringstream msg2;	
				msg2<<"Truncating indices...";
				progress_.printline(msg2,std::cerr);
				truncate(removedIndices);
			}
			
			// N.B.: false below means that we don't truncate the permutation vectors
			//	because they're needed for the WFT
			findPermutationAndPartition(false);
			std::ostringstream msg;
			msg<<"Done with changeBasis";
			progress_.printline(msg,std::cerr);
			return calcError(eigs,removedIndices);
		}

		//! return number of partitions for this basis (see paper)
		size_t partition(size_t i,size_t beforeOrAfterTransform) const 
		{
			if (beforeOrAfterTransform==AFTER_TRANSFORM)
				return partition_[i];
			return partitionOld_[i];
		}

		//! return number of partitions for this basis (see paper)
		size_t partition() const { return partition_.size(); }

		//! return the permutation of i 
		size_t permutation(size_t i) const { return   permutationVector_[i]; }

		const std::vector<size_t>& permutationVector() const 
		{
			return  permutationVector_;
		}

		//! return the inverse permutation of i 
		size_t permutationInverse(size_t i) const { return permInverse_[i]; }

		//! return the inverse permutation vector
		const std::vector<size_t>& permutationInverse() const { return permInverse_; }

		void set(BlockType const &B) { block_=B; }

		static size_t encodeQuantumNumber(const std::vector<size_t>& quantumNumbers)
		{
			if (useSu2Symmetry_) return HamiltonianSymmetrySu2Type::encodeQuantumNumber(quantumNumbers);
			else return HamiltonianSymmetryLocalType::encodeQuantumNumber(quantumNumbers);
		}

		static std::vector<size_t> decodeQuantumNumber(int q)
		{
			if (useSu2Symmetry_) return HamiltonianSymmetrySu2Type::decodeQuantumNumber(q);
			else return HamiltonianSymmetryLocalType::decodeQuantumNumber(q);
		}	
	
		const FactorsType& getFactors() const 
		{
			if (useSu2Symmetry_) return symmSu2_.getFactors();
			else return symmLocal_.getFactors();
		}

		// reduced:
		size_t reducedIndex(size_t i) const 
		{ 
			if (!useSu2Symmetry_) std::cerr<<"########EEEEEEEEERRRRRRRRRRRRRORRRRRRRRRRRR\n";
			return symmSu2_.reducedIndex(i); 
		}

		size_t reducedSize() const 
		{ 
			if (!useSu2Symmetry_) throw std::runtime_error("########EEEEEEEEERRRRRRRRRRRRRORRRRRRRRRRRR\n");
			return symmSu2_.reducedSize(); 
		}

		size_t jVals(size_t i) const 
		{ 
			if (!useSu2Symmetry_)throw std::runtime_error("########EEEEEEEEERRRRRRRRRRRRRORRRRRRRRRRRR\n");
			return symmSu2_.jVals(i); 
		}

		size_t jVals() const 
		{ 
			if (!useSu2Symmetry_) throw std::runtime_error("########EEEEEEEEERRRRRRRRRRRRRORRRRRRRRRRRR\n");
			return symmSu2_.jVals(); 
		}

		size_t jMax() const 
		{ 
			//if (!useSu2Symmetry_) throw std::runtime_error("########EEEEEEEEERRRRRRRRRRRRRORRRRRRRRRRRR\n");
			return symmSu2_.jMax(); 
		}
		
		bool dmrgTransformed() const { return dmrgTransformed_; }

		template<typename IoOutputter>
		void save(IoOutputter& io,const std::string& ss) const
		{
			io.printline("#NAME="+ss);
			saveInternal(io);
		}

		template<typename IoOutputter>
		void save(IoOutputter& io) const
		{
			io.printline("#NAME="+name_);
			saveInternal(io);
		}
		
		template<typename RealType2,typename SparseMatrixType2>
		friend std::ostream& operator<<(std::ostream& os,const BasisImplementation<RealType2,SparseMatrixType2>& x);

	private:
		bool dmrgTransformed_;
		std::string name_;
		ProgressIndicator progress_;
		static bool useSu2Symmetry_;
		std::vector<size_t> quantumNumbers_;
		std::vector<size_t> quantumNumbersOld_;
		std::vector<size_t> electrons_;
		std::vector<size_t> electronsOld_;
		std::vector<size_t> partition_;
		std::vector<size_t> partitionOld_;
		std::vector<size_t> permutationVector_;
		std::vector<size_t> permInverse_;
		HamiltonianSymmetryLocalType symmLocal_;
		HamiltonianSymmetrySu2Type symmSu2_;
		BlockType block_;
		
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
			
			if (useSu2Symmetry_)
				symmSu2_.load(io);
			else
				symmLocal_.load(io);
		}

		template<typename IoOutputter>
		void saveInternal(IoOutputter& io) const
		{
			std::string s="#useSu2Symmetry="+utils::ttos(useSu2Symmetry_);
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

		RealType calcError(std::vector<RealType> const &eigs,std::vector<size_t> const &removedIndices)
		{
			RealType sum=static_cast<RealType>(0.0);
			for (size_t i=0;i<eigs.size();i++) if (utils::isInVector(removedIndices,i)<0) sum+=eigs[i];
			return 1.0-sum;
		}

		void truncate(std::vector<size_t> const &removedIndices)
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
					utils::sort(quantumNumbers_,permutationVector_);
				}
			}
			
			findPartition();
			
			if (changePermutation) {
				permInverse_.resize(permutationVector_.size());
				for (size_t i=0;i<permInverse_.size();i++) 
					permInverse_[permutationVector_[i]]=i;
				
			}
		}

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
	}; // class BasisImplementation

	template<typename RealType,typename SparseMatrixType>
	std::ostream& operator<<(std::ostream& os,const BasisImplementation<RealType,SparseMatrixType>& x)
	{
		os<<"dmrgTransformed="<<x.dmrgTransformed_<<"\n";
		os<<"name="<<x.name_<<"\n";
		utils::vectorPrint(x.quantumNumbers_,"quantumNumbers",os);
		//std::vector<size_t> quantumNumbersOld_;
		utils::vectorPrint(x.electrons_,"electrons",os);
		utils::vectorPrint(x.partition_,"partition",os);
		//std::vector<size_t> partitionOld_;
		utils::vectorPrint(x.permutationVector_,"permutation",os);
		//std::vector<size_t> permInverse_;
		//HamiltonianSymmetryLocalType symmLocal_;
		//HamiltonianSymmetrySu2Type symmSu2_;
		utils::vectorPrint(x.block_,"block",os);
		//for (size_t i=0;i<x.size();i++) 
		//	os<<"i="<<i<<" jm="<<x.jmValue(i)<<" f="<<x.getFlavor(i)<<" e="<<x.getNe(i)<<"\n";
		return os;
	}

	template<typename RealType,typename SparseMatrixType>
	bool BasisImplementation<RealType,SparseMatrixType>::useSu2Symmetry_=false;

} // namespace Dmrg

/*@}*/
#endif


