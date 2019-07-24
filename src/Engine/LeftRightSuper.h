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

/*! \file LeftRightSuper.h
 *
 *  A class that contains the left block or system, the
 *  right block or environ, and the superblock
 */
#ifndef LEFT_RIGHT_SUPER_H
#define LEFT_RIGHT_SUPER_H

#include "ProgressIndicator.h"
#include "KroneckerDumper.h"
#include "Io/IoNg.h"

namespace Dmrg {

template<typename BasisWithOperatorsType_,typename SuperBlockType>
class LeftRightSuper {

public:

	typedef typename SuperBlockType::RealType RealType;
	typedef BasisWithOperatorsType_ BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef typename BasisType::BlockType BlockType;
	typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
	typedef  LeftRightSuper<BasisWithOperatorsType_,SuperBlockType> ThisType;
	typedef KroneckerDumper<ThisType> KroneckerDumperType;
	typedef typename KroneckerDumperType::ParamsForKroneckerDumper ParamsForKroneckerDumperType;
	typedef typename BasisType::QnType QnType;

	template<typename IoInputter>
	LeftRightSuper(IoInputter& io,
	               PsimagLite::String prefix,
	               bool isObserveCode,
	               typename PsimagLite::EnableIf<
	               PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	    : progress_("LeftRightSuper"),
	      left_(0),right_(0),super_(0),refCounter_(0)
	{
		bool minimizeRead = isObserveCode;

		prefix += "/LRS";

		PsimagLite::String nameSuper;
		io.read(nameSuper, prefix + "/NameSuper");

		PsimagLite::String nameLeft;
		io.read(nameLeft, prefix + "/NameSystem");

		PsimagLite::String nameRight;
		io.read(nameRight, prefix + "/NameEnviron");

		super_ = new SuperBlockType(io, prefix + "/" + nameSuper, minimizeRead);
		left_ = new BasisWithOperatorsType(io, prefix + "/" + nameLeft, true);
		right_ = new BasisWithOperatorsType(io, prefix + "/" + nameRight, true);
	}

	LeftRightSuper(const PsimagLite::String& slabel,
	               const PsimagLite::String& elabel,
	               const PsimagLite::String& selabel)
	    : progress_("LeftRightSuper"),
	      left_(0),right_(0),super_(0),refCounter_(0)
	{
		left_ = new BasisWithOperatorsType(slabel);
		right_ = new BasisWithOperatorsType(elabel);
		super_ = new SuperBlockType(selabel);
	}

	~LeftRightSuper()
	{
		if (refCounter_ > 0) {
			--refCounter_;
			return;
		}

		delete left_;
		left_ = 0;
		delete right_;
		right_ = 0;
		delete super_;
		super_ = 0;
	}

	LeftRightSuper(BasisWithOperatorsType& left,
	               BasisWithOperatorsType& right,
	               SuperBlockType& super)
	    : progress_("LeftRightSuper"),
	      left_(&left),right_(&right),super_(&super),refCounter_(1)
	{}

	LeftRightSuper(const ThisType& rls)
	    : progress_("LeftRightSuper"),refCounter_(1)
	{
		left_ = rls.left_;
		right_ = rls.right_;
		super_ = rls.super_;
	}

	void dontCopyOperators(const ThisType& rls)
	{
		assert(left_);
		left_->dontCopyOperators(rls.left());
		assert(right_);
		right_->dontCopyOperators(rls.right());
		assert(super_);
		assert(rls.super_);
		*super_ = *rls.super_;
		if (refCounter_ > 0) --refCounter_;
	}

	template<typename SomeModelType>
	void growLeftBlock(const SomeModelType& model,
	                   BasisWithOperatorsType &pS,
	                   BlockType const &X,
	                   RealType time)
	{
		assert(left_);
		grow(*left_,
		     model,
		     pS,
		     X,
		     ProgramGlobals::DirectionEnum::EXPAND_SYSTEM,
		     time);
	}

	template<typename SomeModelType>
	void growRightBlock(const SomeModelType& model,
	                    BasisWithOperatorsType &pE,
	                    BlockType const &X,
	                    RealType time)
	{
		assert(right_);
		grow(*right_,
		     model,
		     pE,
		     X,
		     ProgramGlobals::DirectionEnum::EXPAND_ENVIRON,
		     time);
	}

	void printSizes(const PsimagLite::String& label,std::ostream& os) const
	{
		assert(left_);
		assert(right_);

		PsimagLite::OstringStream msg;
		msg<<label<<": left-block basis="<<left_->size();
		msg<<", right-block basis="<<right_->size();
		msg<<" sites="<<left_->block().size()<<"+";
		msg<<right_->block().size();
		progress_.printline(msg,os);
	}

	SizeType sites() const
	{
		assert(left_);
		assert(right_);
		return left_->block().size() + right_->block().size();
	}

	/*!PTEX_LABEL{setToProductLrs} */
	void setToProduct(QnType quantumSector, SizeType initialSizeOfHashTable)
	{
		assert(left_);
		assert(right_);
		assert(super_);
		super_->setToProduct(*left_, *right_, &quantumSector, initialSizeOfHashTable);
	}

	void write(PsimagLite::IoNg::Out& io,
	           PsimagLite::String prefix,
	           typename BasisWithOperatorsType::SaveEnum option,
	           bool minimizeWrite) const
	{
		prefix += "/LRS";
		io.createGroup(prefix);

		io.write(super_->name(), prefix + "/NameSuper");
		io.write(left_->name(), prefix + "/NameSystem");
		io.write(right_->name(), prefix + "/NameEnviron");

		assert(left_);
		assert(right_);
		assert(super_);

		super_->write(io, PsimagLite::IoSerializer::NO_OVERWRITE, prefix, minimizeWrite);
		left_->write(io, PsimagLite::IoSerializer::NO_OVERWRITE, prefix, option);
		right_->write(io, PsimagLite::IoSerializer::NO_OVERWRITE, prefix, option);
	}

	const BasisWithOperatorsType& left()  const
	{
		assert(left_);
		return *left_;
	}

	const BasisWithOperatorsType& right() const
	{
		assert(right_);
		return *right_;
	}

	BasisWithOperatorsType& leftNonConst()
	{
		assert(left_);
		return *left_;
	}

	BasisWithOperatorsType& rightNonConst()
	{
		assert(right_);
		return *right_;
	}

	const SuperBlockType& super() const
	{
		assert(super_);
		return *super_;
	}

	void left(const BasisWithOperatorsType& left)
	{
		if (refCounter_ > 0)
			err("LeftRightSuper::left(...): not the owner\n");
		assert(left_);
		*left_=left; // deep copy
	}

	void right(const BasisWithOperatorsType& right)
	{
		if (refCounter_ > 0)
			err("LeftRightSuper::right(...): not the owner\n");
		assert(right_);
		*right_=right; // deep copy
	}

	template<typename IoInputType>
	void read(IoInputType& io,
	          PsimagLite::String prefix,
	          typename PsimagLite::EnableIf<
	          PsimagLite::IsInputLike<IoInputType>::True, int>::Type = 0)
	{
		prefix += "/LRS";

		PsimagLite::String nameSuper;
		io.read(nameSuper, prefix + "/NameSuper");
		PsimagLite::String nameSys;
		io.read(nameSys, prefix + "/NameSystem");
		PsimagLite::String nameEnviron;
		io.read(nameEnviron, prefix + "/NameEnviron");

		super_->read(io,  prefix + "/" + nameSuper);
		left_->read(io, prefix + "/" + nameSys);
		right_->read(io, prefix + "/" + nameEnviron);
	}

private:

	//! add block X to basis pS and put the result in left_:
	/* PSIDOC LeftRightSuperGrow
		Local operators are set for the basis in question with a call to
		\cppClass{BasisWithOperators}'s member function \cppFunction{setOperators()}.
		When adding sites to the system or environment the program does a
		full outer product, i.e., it increases the size of all local operators.
		This is performed by the call to \verb!setToProduct!
		\verb!(pSprime,pS,Xbasis,dir,option)!
		in the grow function, which actually calls \verb!pSprime.setToProduct!
		\verb!(pS,xBasis,dir)!
		This function also recalculates the Hamiltonian in the outer product
		of (i) the previous system basis $pS$, and (ii) the basis $Xbasis$
		corresponding to the site(s) that is (are) being added.
		To do this, the Hamiltonian connection between the two parts
		needs to be calculated and added, and this is done in the call to
		\cppFunction{addHamiltonianConnection}.
		Finally, the resulting dmrgBasis object for the outer product,
		pSprime, is set to contain this full Hamiltonian with the call
		to  \cppFunction{pSprime.setHamiltonian(matrix)}.
		*/
	template<typename SomeModelType>
	void grow(BasisWithOperatorsType& leftOrRight,
	          const SomeModelType& model,
	          BasisWithOperatorsType &pS,
	          const BlockType& X,
	          ProgramGlobals::DirectionEnum dir,
	          RealType time)
	{
		BasisWithOperatorsType Xbasis("Xbasis");
		typedef LeftRightSuper<BasisWithOperatorsType, BasisType> LeftRightSuper2Type;
		Xbasis.setVarious(X, model, time);
		leftOrRight.setToProduct(pS, Xbasis, dir);

		SparseMatrixType matrix = leftOrRight.hamiltonian().getCRS();

		LeftRightSuper2Type* lrs;
		BasisType* leftOrRightL =  &leftOrRight;
		if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			lrs = new LeftRightSuper2Type(pS,Xbasis,*leftOrRightL);
		} else {
			lrs = new  LeftRightSuper2Type(Xbasis,pS,*leftOrRightL);
		}
		//!PTEX_LABEL{295}
		model.addHamiltonianConnection(matrix, *lrs, time);
		delete lrs;
		leftOrRight.setHamiltonian(matrix);
	}

	LeftRightSuper(LeftRightSuper&);

	LeftRightSuper& operator=(const LeftRightSuper&);

	ProgressIndicatorType progress_;
	BasisWithOperatorsType* left_;
	BasisWithOperatorsType* right_;
	SuperBlockType* super_;
	SizeType refCounter_;
}; // class LeftRightSuper

} // namespace Dmrg

/*@}*/
#endif // LEFT_RIGHT_SUPER_H

