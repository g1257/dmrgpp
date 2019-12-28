/*
Copyright (c) 2009-2012-2019, UT-Battelle, LLC
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

/*! \file BasisWithOperators.h
*/

#ifndef BASISWITHOPERATORS_HEADER_H
#define BASISWITHOPERATORS_HEADER_H

#include "ApplyFactors.h"
#include "Basis.h"
#include "Operators.h"

namespace Dmrg {

/* PSIDOC BasisWithOperators
	A class to represent a Hilbert Space for a strongly correlated electron model
		Derives from Basis

	 C++ class \cppClass{Basis} (and \cppClass{BasisImplementation}) implement only
	 certain functionality associated with a Hilbert space basis, as mentioned in
	 the previous section. However, more capabilities related to a Hilbert space basis are needed.

	 C++ class \cppClass{BasisWithOperators} inherits from \cppClass{Basis}, and contains
	 certain local operators for the basis in question, as well as the Hamiltonian matrix.
	 The operators that need to be considered here are operators necessary to compute
	 the Hamiltonian across the system and environment, and to compute observables.
	 Therefore, the specific operators vary from model to model.
	 For example, for the Hubbard model, we consider $c_{i\sigma}$ operators,
	 that destroy an electron with spin $\sigma$ on site $i$.
	 For the Heisenberg model, we consider operators $S^+_i$ and $S^z_i$ for each site $i$.
	 In each case these operators are calculated by the model class (see Section~\ref{subsec:models})
	 on the ``natural'' basis,  and added to the basis in question with a call to
	 \cppFunction{setOperators()}.
	 These local operators are stored as sparse matrices to save memory,
	 although the matrix type is templated and could be anything.
	 For details on the implementation of these operators, see \cppClass{OperatorsBase},
	 its common implementation \cppClass{OperatorsImplementation}, and the two examples provided
	 \cppClass{OperatorsHubbard} and \cppClass{OperatorsHeisenberg} for the Hubbard
	 and Heisenberg models, respectively.
	 Additionally, \cppClass{BasisWithOperators} has a number of member functions to
	 handle operations that the DMRG method performs on
	 local operators in a Hilbert space basis. These include functions to create
	 an outer product of two given Hilbert spaces, to transform a basis, to truncate a basis, etc.
	 */
template<typename BasisType_>
class BasisWithOperators : public BasisType_ {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorIntegerType;
	typedef BasisType_ BaseType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename BaseType::RealType RealType;
	typedef Operators<BasisType_> OperatorsType;
	typedef typename OperatorsType::PairSizeSizeType PairSizeSizeType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType VectorSizeType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef BasisWithOperators<BasisType_> ThisType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef typename OperatorStorageType::value_type ComplexOrRealType;
	typedef typename PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;

	enum class SaveEnum {ALL, PARTIAL};

	BasisWithOperators(const PsimagLite::String& s)
	    : BasisType(s)
	{}

	template<typename IoInputter>
	BasisWithOperators(IoInputter& io,
	                   const PsimagLite::String& ss,
	                   bool isObserveCode)
	    : BasisType(io, ss, false),
	      operators_(io, ss + "/", isObserveCode)
	{
		const PsimagLite::String prefix = ss + "/";
		io.read(operatorsPerSite_, prefix + "OperatorPerSite");
	}

	template<typename IoInputter>
	void read(IoInputter& io,
	          PsimagLite::String prefix,
	          typename PsimagLite::EnableIf<
	          PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	{
		BasisType::read(io, prefix); // parent loads
		operators_.read(io, prefix);
		io.read(operatorsPerSite_, prefix + "/OperatorPerSite");
	}

	void dontCopyOperators(const BasisWithOperators& b)
	{
		BaseType& base = *this;
		const BaseType& b1 = static_cast<BaseType>(b);
		base = b1;
		operatorsPerSite_ = b.operatorsPerSite_;
		operators_.clear();
	}

	// set this basis to the outer product of
	// basis2 and basis3 or basis3 and basis2  depending on dir
	void setToProduct(const ThisType& basis2,
	                  const ThisType& basis3,
	                  ProgramGlobals::DirectionEnum dir)
	{
		if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			setToProduct(basis2, basis3);
		else
			setToProduct(basis3, basis2);
	}

	//! transform this basis by transform
	//! note: basis change must conserve total number of electrons and all quantum numbers
	RealType truncateBasis(const BlockDiagonalMatrixType& ftransform,
	                       const typename PsimagLite::Vector<RealType>::Type& eigs,
	                       const VectorSizeType& removedIndices,
	                       const PairSizeSizeType& startEnd,
	                       bool blasIsThreadSafe)
	{
		RealType error = BasisType::truncateBasis(eigs, removedIndices);

		operators_.changeBasis(ftransform, this, startEnd, blasIsThreadSafe);

		return error;
	}

	void setHamiltonian(const SparseMatrixType& h)
	{
		operators_.setHamiltonian(h);
	}

	const OperatorStorageType& hamiltonian() const
	{
		return operators_.hamiltonian();
	}

	const OperatorStorageType& reducedHamiltonian() const
	{
		return operators_.reducedHamiltonian();
	}

	template<typename SomeModelType>
	void setOneSite(const VectorSizeType& block,
	                const SomeModelType& model,
	                RealType time)
	{
		typename BaseType::VectorQnType qm;

		BaseType::set(block);
		typename PsimagLite::Vector<OperatorType>::Type ops;
		SparseMatrixType h;

		model.setOperatorMatrices(ops, qm, block);
		BaseType::setSymmetryRelated(qm);

		model.calcHamiltonian(h, ops, block, time);

		OperatorStorageType hOp(h);
		operators_.setHamiltonian(hOp);
		operators_.setLocal(ops);

		// one site basis is assumed ordered

		operatorsPerSite_.clear();
		for (SizeType i = 0; i < block.size(); ++i)
			operatorsPerSite_.push_back(ops.size()/block.size());

		assert(operatorsPerSite_.size() > 0);
	}

	PairType getOperatorIndices(SizeType i,SizeType sigma) const
	{
		SizeType sum = 0;
		for (SizeType j=0;j<i;j++) {
			assert(j<operatorsPerSite_.size());
			sum += operatorsPerSite_[j];
		}

		assert(i<operatorsPerSite_.size());
		return PairType(sum + sigma,operatorsPerSite_[i]);
	}

	const OperatorType& getOperatorByIndex(int i) const
	{
		return operators_.getLocalByIndex(i);
	}

	SizeType numberOfLocalOperators() const { return operators_.sizeOfLocal(); }

	SizeType superOperatorIndices(const VectorSizeType& sites, SizeType sigma) const
	{
		return operators_.superIndices(sites, sigma);
	}

	SizeType operatorsPerSite(SizeType i) const
	{
		assert(i < operatorsPerSite_.size());
		return operatorsPerSite_[i];
	}

	int fermionicSign(SizeType i,int fsign) const
	{
		const BasisType &parent = *this;
		return parent.fermionicSign(i,fsign);
	}

	template<typename SomeOutputType>
	void write(SomeOutputType& io,
	           typename SomeOutputType::Serializer::WriteMode mode,
	           PsimagLite::String prefix,
	           SaveEnum option,
	           typename PsimagLite::EnableIf<
	           PsimagLite::IsOutputLike<SomeOutputType>::True, int*>::Type = 0) const
	{
		write(io, prefix + "/" + BasisType::name(), mode, option);
	}

	template<typename SomeOutputType>
	void write(SomeOutputType& io,
	           const PsimagLite::String& s,
	           typename SomeOutputType::Serializer::WriteMode mode,
	           SaveEnum option,
	           typename PsimagLite::EnableIf<
	           PsimagLite::IsOutputLike<SomeOutputType>::True, int*>::Type = 0) const
	{
		BasisType::write(io, s, mode, false); // parent saves
		if (option == SaveEnum::ALL)
			operators_.write(io, s, mode);

		assert(operatorsPerSite_.size() > 0);
		io.write(operatorsPerSite_, s + "/OperatorPerSite", mode);
	}

private:

	//! set this basis to the outer product of   basis2 and basis3
	//!PTEX_LABEL{setToProductOps}
	void setToProduct(const ThisType& basis2,
	                  const ThisType& basis3)
	{
		// reorder the basis
		BasisType::setToProduct(basis2, basis3);

		typename PsimagLite::Vector<RealType>::Type fermionicSigns;
		SizeType x = basis2.numberOfLocalOperators()+basis3.numberOfLocalOperators();

		operators_.setToProduct(basis2, basis3, x);
		ApplyFactors<FactorsType> apply(BasisType::getFactors(), BasisType::useSu2Symmetry());
		ProgramGlobals::FermionOrBosonEnum savedSign = ProgramGlobals::FermionOrBosonEnum::BOSON;

		const SizeType nlocalOps = numberOfLocalOperators();
		for (SizeType i = 0; i < nlocalOps; ++i) {
			if (i<basis2.numberOfLocalOperators()) {
				const OperatorType& myOp =  basis2.getOperatorByIndex(i);
				bool isFermion = (myOp.fermionOrBoson() ==
				                  ProgramGlobals::FermionOrBosonEnum::FERMION);
				if (savedSign != myOp.fermionOrBoson() || fermionicSigns.size() == 0) {
					utils::fillFermionicSigns(fermionicSigns,
					                          basis2.signs(),
					                          (isFermion) ? -1 : 1);
					savedSign = myOp.fermionOrBoson();
				}

				operators_.crossProductForLocal(i,
				                                myOp,
				                                basis3.size(),
				                                fermionicSigns,
				                                true,
				                                BaseType::permutationInverse());

			} else {
				const OperatorType& myOp = basis3.
				        getOperatorByIndex(i - basis2.numberOfLocalOperators());

				bool isFermion = (myOp.fermionOrBoson() ==
				                  ProgramGlobals::FermionOrBosonEnum::FERMION);

				if (savedSign != myOp.fermionOrBoson() || fermionicSigns.size() == 0) {
					utils::fillFermionicSigns(fermionicSigns,
					                          basis2.signs(),
					                          (isFermion) ? -1 : 1);
					savedSign = myOp.fermionOrBoson();
				}

				operators_.crossProductForLocal(i,
				                                myOp,
				                                basis2.size(),
				                                fermionicSigns,
				                                false,
				                                BaseType::permutationInverse());
			}
		}

		//! Calc. hamiltonian
		operators_.outerProductHamiltonian(basis2.hamiltonian(),
		                                   basis3.hamiltonian(),
		                                   apply,
		                                   BaseType::permutationInverse());

		SizeType offset1 = basis2.operatorsPerSite_.size();
		operatorsPerSite_.resize(offset1+basis3.operatorsPerSite_.size());
		for (SizeType i=0;i<offset1;i++)
			operatorsPerSite_[i] =  basis2.operatorsPerSite_[i];

		for (SizeType i=0;i<basis3.operatorsPerSite_.size();i++)
			operatorsPerSite_[i+offset1] =  basis3.operatorsPerSite_[i];
		assert(operatorsPerSite_.size() > 0);
	}

	// BasisWithOperators(const BasisWithOperators&);

	// BasisWithOperators& operator=(const BasisWithOperators&);

	OperatorsType operators_;
	PsimagLite::Vector<SizeType>::Type operatorsPerSite_;
}; // class BasisWithOperators

template<typename T>
struct IsBasisType<BasisWithOperators<T> > {
	enum {True = true};
};
} // namespace Dmrg

/*@}*/
#endif

