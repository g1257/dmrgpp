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

/*! \file ModelBase.h
 *
 *
 */

#ifndef MODEL_BASE_H
#define MODEL_BASE_H

#include "Vector.h"
#include "Sort.h"
#include "MemResolv.h"
#include "TargetQuantumElectrons.h"
#include "Io/IoSerializerStub.h"
#include "ModelCommon.h"
#include "QnHash.h"
#include "ParallelHamiltonianConnection.h"

namespace Dmrg {

template<typename ModelHelperType_,
         typename ParametersType_,
         typename InputValidatorType_,
         typename GeometryType_>
class ModelBase  {

public:

	struct OpaqueOp {

		OpaqueOp(PsimagLite::String name_, SizeType dof_ = 0, SizeType edof_ = 0)
		    : name(name_), dof(dof_), edof(edof_)
		{}

		PsimagLite::String name;
		SizeType dof;
		SizeType edof;
	};

	typedef ParametersType_ ParametersType;
	typedef InputValidatorType_ InputValidatorType;
	typedef ModelHelperType_ ModelHelperType;
	typedef GeometryType_ GeometryType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::BlockType BlockType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::BasisType MyBasis;
	typedef typename ModelHelperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename MyBasis::QnType QnType;
	typedef TargetQuantumElectrons<RealType, QnType> TargetQuantumElectronsType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::SparseElementType ComplexOrRealType;
	typedef ModelCommon<ParametersType, GeometryType, ModelHelperType> ModelCommonType;
	typedef typename ModelCommonType::HamiltonianConnectionType HamiltonianConnectionType;
	typedef typename ModelCommonType::VectorLinkType VectorLinkType;
	typedef typename ModelCommonType::VectorType VectorType;
	typedef ParametersType SolverParamsType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef VectorSizeType HilbertBasisType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename QnType::PairSizeType PairSizeType;
	typedef typename BasisWithOperatorsType::VectorBoolType VectorBoolType;
	typedef typename ModelCommonType::LabeledOperatorsType LabeledOperatorsType;
	typedef typename ModelCommonType::ModelLinksType ModelLinksType;
	typedef typename LabeledOperatorsType::LabelType OpsLabelType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename ModelCommonType::VerySparseMatrixType VerySparseMatrixType;
	typedef ParallelHamiltonianConnection<HamiltonianConnectionType> ParallelHamConnectionType;
	typedef typename ModelLinksType::TermType ModelTermType;
	typedef OpaqueOp OpForLinkType;

	ModelBase(const ParametersType& params,
	          const GeometryType_& geometry,
	          InputValidatorType& io)
	    : modelCommon_(params, geometry),
	      targetQuantum_(io),
	      ioIn_(io)
	{
		labeledOperators_.setModelName(params.model);
	}

	void postCtor()
	{
		fillLabeledOperators(qns_); // fills qns_ and labeledOperators_
		modelLinks_.postCtor1(labeledOperators_, modelCommon_.geometry().terms());
		fillModelLinks(); // fills modelLinks_
		modelLinks_.postCtor2();

		ProgramGlobals::init(maxElectronsOneSpin());
	}

	virtual ~ModelBase() {}

	/* PSIDOC ModelInterface
	These are the functions that each model must implement.

		PSIDOCCOPY ModelBaseWrite

\noindent\dotfill\\

		PSIDOCCOPY addDiagonalsInNaturalBasis

\noindent\dotfill\\

		PSIDOCCOPY fillLabeledOperators

\noindent\dotfill\\

		PSIDOCCOPY fillModelLinks
	*/

	/* PSIDOC ModelBaseWrite FirstProtoBelow
	PSIDOCCOPY $FirstProtoBelow
	For information purposes only.
	String is the first argument and contains the group.
	Serializer object is second argument.
	As an example, we describe the write() virtual memeber function
	for the Hubbard model.

	PSIDOCCOPY Hubbard::write
	*/
	virtual void write(PsimagLite::String,
	                   PsimagLite::IoNg::Out::Serializer&) const = 0;


	/* PSIDOC addDiagonalsInNaturalBasis
	PSIDOCCOPY $FirstProtoBelow
	Takes 3 arguments: the first one is output, the last two are input.
	The first argument is a CRS matrix that you need to
	fill with the on-site Hamiltonian terms in the on-site basis.
	The second argument is a vector of sites containing only one site; this
	argument is an input. The last argument is
	the current physical time; it's an input, and is given in case your
	model has time-dependent on-site terms. We briefly discuss the
	addition of a magnetic field Zeeman term to the Heisenberg model.

	PSIDOCCOPY Heisenberg::addDiagonalsInNaturalBasis
	*/
	virtual void addDiagonalsInNaturalBasis(SparseMatrixType&,
	                                        const BlockType& block,
	                                        RealType)  const = 0;

	/* PSIDOC fillLabeledOperators
	PSIDOCCOPY $FirstProtoBelow
	  The only argument given to this function is the qns that you must fill.
	  This is a vector of \texttt{Qn} objects.
	  You need to loop over the one-site Hilbert space, something like
\begin{lstlisting}
for (SizeType i = 0; i < numberOfStates; ++i) {
  ...
  qns[i] = QnType(sign, other, jmpair, flavor);
}
	   \end{lstlisting}
	   where sign is true is state i has odd number of fermions and false
	   otherwise, other is a vector of conserved quantities evaluated on state i,
	   jmpair is the 2j and j+m of this state, and flavor is the flavor of state i.
	   You may set jmpair to PairSizeType(0,0) and flavor to 0, in which case
	   your model won't support SU(2) symmetry.

	   After setting qns, you need to create all one-site labeled operators for this
	   model. These are the operators that the engine must keep track, and, optionally,
	   other labeled operators to be used in operator specifications.
	   You create labeled operators with the following syntax.
\begin{lstlisting}
OpsLabelType& myoperator = this->createOperatorLabel("myoperator");
myoperator.makeTrackable(); // only if needs tracking
\end{lstlisting}
	   then you fill it with
\begin{lstlisting}
for (SizeType dof = 0; dof < numberOfDofs; ++dof) {
  ...
  OperatorType someOperator(sparseMatrix, fermionSign, ...);
  myoperator.push(someOperator);
}
\end{lstlisting}
		And you must mark the operators that need tracking as shown above.
	*/
	virtual void fillLabeledOperators(VectorQnType&) = 0;

	/* PSIDOC fillModelLinks
	PSIDOCCOPY $FirstProtoBelow
	Give the Hamiltonian connections that this model has.
	These are the terms in the Hamiltonian that connect \emph{different} sites.
	We'll go first though the simpler example of the Hubbard model.

	PSIDOCCOPY Hubbard::fillModelLinks
	\vspace{1em}

	OK, let's discuss a more complicated example: the case of the FeAs model.
	PSIDOCCOPY FeAs::fillModelLinks
	 */
	virtual void fillModelLinks() = 0;

	// END ^^^^^^^^^^^Functions that each model needs to implement

	virtual void findOddElectronsOfOneSite(VectorBoolType& oddElectrons,
	                                       SizeType site) const
	{
		typename PsimagLite::Vector<SizeType>::Type block(1, site);
		typename PsimagLite::Vector<OperatorType>::Type cm;
		VectorQnType qq;
		setOperatorMatrices(cm, qq, block);
		SizeType n = qq.size();
		oddElectrons.resize(n);
		for (SizeType i = 0; i < n; ++i)
			oddElectrons[i] = qq[i].oddElectrons;
	}

	//! Full hamiltonian from creation matrices cm
	virtual void calcHamiltonian(SparseMatrixType &hmatrix,
	                             const VectorOperatorType& cm,
	                             const BlockType& block,
	                             RealType time)  const
	{
		hmatrix.makeDiagonal(cm[0].getStorage().rows());

		modelCommon_.addConnectionsInNaturalBasis(hmatrix,cm,block,time);

		addDiagonalsInNaturalBasis(hmatrix, block, time);
	}

	virtual SizeType maxElectronsOneSpin() const
	{
		SizeType tmp = hilbertSize(0);
		tmp = static_cast<SizeType>(log(tmp)/log(2.0));
		SizeType maxElectrons = static_cast<SizeType>(tmp/2);
		if (tmp & 1) maxElectrons++;

		return maxElectrons*modelCommon_.geometry().numberOfSites() + 1;
	}

	virtual SizeType siteToAtomKind(SizeType) const { return 0; }

	virtual SizeType kindsOfAtoms() const { return 1; }

	virtual PsimagLite::String oracle() const { return ""; }

	/**
		The function \cppFunction{addHamiltonianConnection} implements
		the Hamiltonian connection (e.g. tight-binding links in the case of the Hubbard Model
		or products $S_i\cdot S_j$ in the case of the Heisenberg model) between
		two basis, $basis2$ and $basis3$, in the order of the outer product,
		$basis1={\rm SymmetryOrdering}(basis2\otimes basis3)$. This was
		explained before in Section~\ref{subsec:dmrgBasisWithOperators}.
		This function has a default implementation.
		*/
	void addHamiltonianConnection(SparseMatrixType& matrix,
	                              const LeftRightSuperType& lrs,
	                              RealType currentTime) const
	{
		PsimagLite::Profiling profiling("addHamiltonianConnection",
		                                "",
		                                std::cout);

		assert(lrs.super().partition() > 0);
		SizeType total = lrs.super().partition()-1;

		typename PsimagLite::Vector<VerySparseMatrixType*>::Type vvsm(total, 0);
		VectorSizeType nzs(total, 0);

		for (SizeType m = 0; m < total; ++m) {
			SizeType offset = lrs.super().partition(m);
			assert(lrs.super().partition(m + 1) >= offset);
			SizeType bs = lrs.super().partition(m + 1) - offset;

			vvsm[m] = new VerySparseMatrixType(bs, bs);
			VerySparseMatrixType& vsm = *(vvsm[m]);
			HamiltonianConnectionType hc(m,
			                             lrs,
			                             modelCommon_.geometry(),
			                             modelLinks_,
			                             currentTime,
			                             0);

			hc.matrixBond(vsm);
			nzs[m] = vsm.nonZeros();
			if (nzs[m] > 0) continue;
			delete vvsm[m];
			vvsm[m] = 0;
		}

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType permutation(total, 0);
		sort.sort(nzs, permutation);

		typename PsimagLite::Vector<const SparseMatrixType*>::Type vectorOfCrs;

		assert(total == permutation.size());
		for (SizeType i = 0; i < total; ++i) { // loop over new order

			SizeType m = permutation[i]; // get old index from new index

			if (vvsm[m] == 0) continue;

			const VerySparseMatrixType& vsm = *(vvsm[m]);
			SparseMatrixType matrixBlock2;
			matrixBlock2 = vsm;
			delete vvsm[m];
			vvsm[m] = 0;

			SizeType offset = lrs.super().partition(m);
			SparseMatrixType* full = new SparseMatrixType(matrix.rows(),
			                                              matrix.cols(),
			                                              matrixBlock2.nonZeros());
			fromBlockToFull(*full, matrixBlock2, offset);
			vectorOfCrs.push_back(full);
		}

		if (vectorOfCrs.size() == 0) return;

		vectorOfCrs.push_back(&matrix);
		SizeType effectiveTotal = vectorOfCrs.size();

		VectorType ones(effectiveTotal, 1.0);
		SparseMatrixType sumCrs;
		sum(sumCrs, vectorOfCrs, ones);
		vectorOfCrs.pop_back();
		effectiveTotal = vectorOfCrs.size();
		for (SizeType i = 0; i < effectiveTotal; ++i) {
			delete vectorOfCrs[i];
			vectorOfCrs[i] = 0;
		}

		matrix.swap(sumCrs);
	}

	/** Let H be the hamiltonian of the  model for basis1 and partition m
	 * consisting of the external product
		 * of basis2 \otimes basis3
		 * This function does x += H*y
		 * The \cppFunction{matrixVectorProduct} function implements the operation $x+=Hy$.
		 * This function
		 * has a default implementation.
		 */
	void matrixVectorProduct(VectorType& x,
	                         const VectorType& y,
	                         const HamiltonianConnectionType& hc) const
	{
		typedef PsimagLite::Parallelizer<ParallelHamConnectionType> ParallelizerType;
		ParallelizerType parallelConnections(PsimagLite::Concurrency::codeSectionParams);

		ParallelHamConnectionType phc(x, y, hc);
		parallelConnections.loopCreate(phc);

		phc.sync();
	}

	void fullHamiltonian(SparseMatrixType& matrix,
	                     const HamiltonianConnectionType& hc) const
	{
		return modelCommon_.fullHamiltonian(matrix, hc);
	}

	// Return the size of the one-site Hilbert space basis for this model
	// site MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	// should be static
	SizeType hilbertSize(SizeType actualSite) const
	{
		const SizeType kindOfSite = siteToAtomKind(actualSite);
		return modelLinks_.hilbertSize(kindOfSite, labeledOperators_);
	}

	// Fill the VectorOperatorType with operators that need to be kept
	// track by the DMRG++ Engine.
	// Fill VectorQnType with the qns of the one site basis in the order
	// you chose to give the operators
	// You can check that block.size() == 1 or throw otherwise
	// The contents of block MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	// should be static
	void setOperatorMatrices(VectorOperatorType& cm,
	                         VectorQnType& qns,
	                         const BlockType& block) const
	{
		assert(block.size() == 1);

		const SizeType kindOfSite = siteKind(block[0]);
		modelLinks_.setOperatorMatrices(cm, labeledOperators_, kindOfSite);

		const SizeType k = kindsOfAtoms();
		SizeType start = 0;
		for (SizeType i = 0; i < k; ++i) {
			if (i == kindOfSite)
				break;
			start += modelLinks_.hilbertSize(i, labeledOperators_);
		}

		const SizeType end = start + modelLinks_.hilbertSize(kindOfSite, labeledOperators_);

		qns.resize(end - start);
		std::copy(qns_.begin() + start, qns_.end() + end, qns.begin());
	}

	static const ModelLinksType& modelLinks()
	{
		return modelLinks_;
	}

	// should be static
	const OperatorType& naturalOperator(const PsimagLite::String& what,
	                                    SizeType kindOfSite, // ignore, legacy
	                                    SizeType dof) const
	{
		return labeledOperators_(what, dof);
	}

	// should be static
	bool instrospect() const
	{
		labeledOperators_.instrospect();
		return true;
	}

	// should be static
	void printBasis(SizeType site) const
	{
		BlockType block(1, site);
		typename PsimagLite::Vector<OperatorType>::Type cm;
		VectorQnType qq;
		setOperatorMatrices(cm, qq, block);
		std::cout<<"block="<<block;
		std::cout<<"qq="<<qq;

		SizeType n = cm.size();
		for (SizeType i = 0; i < n; ++i) {
			std::cout<<"Matrix "<<i<<"\n";
			cm[i].write(std::cout);
		}
	}

	const GeometryType& geometry() const { return modelCommon_.geometry(); }

	const ParametersType& params() const { return modelCommon_.params(); }

	const TargetQuantumElectronsType& targetQuantum() const
	{
		return targetQuantum_;
	}

	static void orderByQuantum(VectorSizeType& basis, VectorQnType& qns)
	{
		const SizeType initialSizeOfHashTable = 10;
		const SizeType n = qns.size();
		std::hash<QnType> myhash(true);
		std::unordered_map<QnType, SizeType> qnSizes(initialSizeOfHashTable, myhash);
		std::unordered_map<QnType, SizeType> seenThisQns(initialSizeOfHashTable, myhash);
		VectorQnType uniqueQns;

		for (SizeType i = 0; i < n; ++i) {
			const QnType& qn = qns[i];
			++qnSizes[qn];
			if (seenThisQns[qn] == 1) continue;
			seenThisQns[qn] = 1;
			uniqueQns.push_back(qn);
		}

		std::unordered_map<QnType, SizeType> offsets(initialSizeOfHashTable, myhash);
		offsetsFromSizes(offsets, qnSizes, uniqueQns);

		std::unordered_map<QnType, SizeType> extraOffsets(initialSizeOfHashTable, myhash);
		VectorSizeType basisNew(basis.size());
		assert(0 < qns.size());
		VectorQnType qnNew(qns.size(), qns[0]);
		for (SizeType i = 0; i < n; ++i) {
			const QnType& thisQn = qns[i];
			SizeType sum = extraOffsets[thisQn];
			const SizeType offset = offsets[thisQn];
			const SizeType ipos = offset + sum;
			++extraOffsets[thisQn];
			basisNew[ipos] = i;
			qnNew[ipos] = thisQn;
		}

		basis = basisNew;
		qns = qnNew;
	}

	InputValidatorType_& ioIn() const { return ioIn_; }

protected:

	PsimagLite::String oracle(const RealType& energy,
	                          const PsimagLite::String formula) const
	{
		if (modelCommon_.params().options.isSet("TargetingAncilla"))
			return "";

		const PsimagLite::String s = " Oracle: " + ttos(energy) + " : " + formula;
		return s;
	}

	static OpsLabelType& createOpsLabel(PsimagLite::String name,
	                                    SizeType kindOfSite = 0)
	{
		return labeledOperators_.createLabel(name, kindOfSite);
	}

	static void makeTrackable(PsimagLite::String name)
	{
		modelLinks_.makeTrackable(name);
	}

//	static void makeTrackableOrderMatters(VectorStringType vname, SizeType site = 0)
//	{
//		SizeType n = vname.size();
//		for (SizeType i = 0; i < n; ++i)
//			labeledOperators_.makeTrackableOrderMatters(vname[i], site);
//	}

	static ModelTermType& createTerm(PsimagLite::String name)
	{
		return modelLinks_.createTerm(name);
	}

private:

	static void offsetsFromSizes(std::unordered_map<QnType, SizeType>& offsets,
	                             std::unordered_map<QnType, SizeType>& sizes,
	                             const VectorQnType& qns)
	{
		const SizeType total = sizes.size();

		SizeType offset = 0;
		for (SizeType i = 0; i < total; ++i) {
			const QnType& qn = qns[i];
			const SizeType thisSize = sizes[qn];
			offsets[qn] = offset;
			offset  += thisSize;
		}
	}

	ModelCommonType modelCommon_;
	TargetQuantumElectronsType targetQuantum_;
	InputValidatorType_& ioIn_;
	static LabeledOperatorsType labeledOperators_;
	static ModelLinksType modelLinks_;
	static VectorQnType qns_;
}; //class ModelBase

template<typename T1, typename T2, typename T3, typename T4>
typename ModelBase<T1, T2, T3, T4>::VectorQnType ModelBase<T1, T2, T3, T4>::qns_;

template<typename T1, typename T2, typename T3, typename T4>
typename ModelBase<T1, T2, T3, T4>::ModelLinksType ModelBase<T1, T2, T3, T4>::modelLinks_;

} // namespace Dmrg
/*@}*/
#endif
