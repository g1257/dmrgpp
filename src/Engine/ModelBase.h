/*
Copyright (c) 2009-2018-2021, UT-Battelle, LLC
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
#include "Braket.h"
#include "SuperOpHelperBase.h"
#include "OutputFileOrNot.h"
#include "CanonicalExpression.h"
#include "OperatorSpec.h"

namespace Dmrg {

template<typename ModelHelperType_,
         typename ParametersType_,
         typename InputValidatorType_,
         typename SuperGeometryType_>
class ModelBase  {

public:

	struct OpaqueOp {

		OpaqueOp(PsimagLite::String name_, SizeType dof_ = 0, SizeType edof_ = 0)
		    : name(name_), dof(dof_), edof(edof_)
		{
			fermionOrBoson = labeledOperators_(name_, dof_).fermionOrBoson();
			kindOfSite = labeledOperators_.findLabel(name_).kindOfSite();
		}

		PsimagLite::String name;
		SizeType dof;
		SizeType edof;
		SizeType kindOfSite;
		ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
	};

	typedef ParametersType_ ParametersType;
	typedef InputValidatorType_ InputValidatorType;
	typedef ModelHelperType_ ModelHelperType;
	typedef SuperGeometryType_ SuperGeometryType;
	typedef ModelBase<ModelHelperType_, ParametersType_, InputValidatorType_, SuperGeometryType>
	ThisType;
	typedef Braket<ThisType> BraketType;
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
	typedef ModelCommon<ParametersType, SuperGeometryType, ModelHelperType> ModelCommonType;
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
	typedef typename ModelLinksType::AtomKindBase AtomKindBaseType;
	typedef SuperOpHelperBase<SuperGeometryType, ParametersType> SuperOpHelperBaseType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename ModelLinksType::LabelType LabelType;
	typedef OperatorSpec<ThisType, OperatorType> OperatorSpecType;
	typedef PsimagLite::CanonicalExpression<OperatorSpecType> CanonicalExpressionType;

	ModelBase(const ParametersType& params,
	          const SuperGeometryType& superGeometry,
	          InputValidatorType& io)
	    : modelCommon_(params, superGeometry),
	      targetQuantum_(io),
	      ioIn_(io),
	      atomKind_(nullptr),
	      superOpHelper_(nullptr)
	{
		labeledOperators_.setModelName(params.model);
		modelLinks_.clear();
		qns_.clear();
		readOnSiteH(io);
	}

	void postCtor()
	{
		modelLinks_.setAtomKind(&getAtomKind());
		fillLabeledOperators(qns_); // fills qns_ and labeledOperators_
		createIdentity();
		if (modelLinks_.kindsOfAtoms() == 1)
			checkThatQnsAreNotReallySorted();
		modelLinks_.postCtor1(labeledOperators_, modelCommon_.superGeometry().terms());
		fillModelLinks(); // fills modelLinks_
		customOperators();
		modelLinks_.postCtor2();

		superOpHelper_ = setSuperOpHelper();
		assert(superOpHelper_);

		ProgramGlobals::init(maxElectronsOneSpin());
	}

	// START OF VIRTUAL FUNCTIONS

	virtual ~ModelBase()
	{
		delete atomKind_;
		atomKind_ = nullptr;
		delete superOpHelper_;
		superOpHelper_ = nullptr;
	}

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

	virtual void findOddElectronsOfOneSite(VectorBoolType& oddElectrons, SizeType site) const
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

		return maxElectrons*modelCommon_.superGeometry().numberOfSites() + 1;
	}

	virtual const AtomKindBaseType& getAtomKind()
	{
		if (!atomKind_)
			atomKind_ = new AtomKindBaseType();
		return *atomKind_;
	}

	// No need to override unless doing OneSiteTruncation
	// Fill the VectorOperatorType with operators that need to be kept
	// track by the DMRG++ Engine.
	// Fill VectorQnType with the qns of the one site basis in the order
	// you chose to give the operators
	// You can check that block.size() == 1 or throw otherwise
	// The contents of block MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	// should be static
	virtual SizeType setOperatorMatrices(VectorOperatorType& cm,
	                                     VectorQnType& qns,
	                                     const BlockType& block) const
	{
		assert(block.size() == 1);

		const SizeType kindOfSite = modelLinks_.siteToAtomKind(block[0]);
		modelLinks_.setOperatorMatrices(cm, labeledOperators_, kindOfSite);

		const SizeType k = modelLinks_.kindsOfAtoms();
		SizeType start = 0;
		for (SizeType i = 0; i < k; ++i) {
			if (i == kindOfSite)
				break;
			start += modelLinks_.hilbertSize(i);
		}

		const SizeType end = start + modelLinks_.hilbertSize(kindOfSite);

		assert(end >= start);
		assert(start < qns_.size());
		qns.resize(end - start, qns_[start]);
		std::copy(qns_.begin() + start, qns_.begin() + end, qns.begin());
		return 0;
	}

	// Models may ignore announcements from the engine
	virtual void announce(PsimagLite::String) const
	{}

	virtual PsimagLite::String oracle() const { return ""; }

	virtual void oneSiteTruncationUpdate(OutputFileOrNot&, const MatrixType&, SizeType)
	{
		qns_.clear();
		labeledOperators_.clear();
		modelLinks_.clear();
		// we could also clear atomKind and superOpHelper here if needed (?)

		postCtor();
	}

	// for models with entanglers only
	virtual bool isCorrectlyPaired(SizeType) const
	{
		throw PsimagLite::RuntimeError("This model does not support entanglers\n");
	}

	// END OF VIRTUAL FUNCTIONS

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

		HamiltonianConnectionType hc(lrs,
		                             modelLinks_,
		                             currentTime,
		                             superOpHelper());

		for (SizeType m = 0; m < total; ++m) {
			SizeType offset = lrs.super().partition(m);
			assert(lrs.super().partition(m + 1) >= offset);
			SizeType bs = lrs.super().partition(m + 1) - offset;

			vvsm[m] = new VerySparseMatrixType(bs, bs);
			VerySparseMatrixType& vsm = *(vvsm[m]);

			typename HamiltonianConnectionType::AuxType aux(m, lrs);

			hc.matrixBond(vsm, aux);
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
	                         const HamiltonianConnectionType& hc,
	                         const typename ModelHelperType::Aux& aux) const
	{
		typedef PsimagLite::Parallelizer<ParallelHamConnectionType> ParallelizerType;
		ParallelizerType parallelConnections(PsimagLite::Concurrency::codeSectionParams);

		ParallelHamConnectionType phc(x, y, hc, aux);
		parallelConnections.loopCreate(phc);

		phc.sync();
	}

	void fullHamiltonian(SparseMatrixType& matrix,
	                     const HamiltonianConnectionType& hc,
	                     const typename ModelHelperType::Aux& aux) const
	{
		return modelCommon_.fullHamiltonian(matrix, hc, aux);
	}

	// Return the size of the one-site Hilbert space basis for this model
	// site MUST be ignored unless your model has a site-dependent
	// Hilbert space (SDHS)
	// should be static
	SizeType hilbertSize(SizeType actualSite) const
	{
		const SizeType kindOfSite = modelLinks_.siteToAtomKind(actualSite);
		return modelLinks_.hilbertSize(kindOfSite);
	}

	static const ModelLinksType& modelLinks()
	{
		return modelLinks_;
	}

	static OperatorType naturalOperator(const PsimagLite::String& what,
	                                    SizeType site,
	                                    SizeType dof)
	{
		static const PsimagLite::String expipi = "exp_i_pi_";
		static const SizeType l = expipi.length();
		PsimagLite::String what2 = what;
		OperatorType op;

		if (what.substr(0, l) == expipi) {
			what2 = what.substr(l, what.length() - l);
			op = labeledOperators_(what2, dof);
			if (op.fermionOrBoson() == ProgramGlobals::FermionOrBosonEnum::FERMION)
				err("Don't know how to exponentiate a fermionic operator\n");

			MatrixType m2 = op.getCRS().toDense();
			expIpi(m2);
			op.fromStorage(m2);
		} else {
			op = labeledOperators_(what, dof);
		}

		invalidateIfNeeded(op, site, what2);
		return op;
	}

	static bool introspect()
	{
		labeledOperators_.introspect();
		return true;
	}

	void printBasis(SizeType site) const
	{
		BlockType block(1, site);
		typename PsimagLite::Vector<OperatorType>::Type cm;
		VectorQnType qq;
		setOperatorMatrices(cm, qq, block);
		std::cout<<"block="<<block;
		std::cout<<"qq="<<qq;

		const SizeType kindOfSite = modelLinks_.siteToAtomKind(site);
		const SizeType n = labeledOperators_.size();
		for (SizeType i = 0; i < n; ++i) {
			const LabelType& ll = labeledOperators_[i];

			if (ll.kindOfSite() != kindOfSite)
				continue;

			const SizeType dofs = ll.dofs();
			std::cout<<"Operator name="<<ll.name()<<" has "<<dofs<<" dofs";
			if (ll.isTrackable())
				std::cout<<" Trackable: YES\n";
			else
				std::cout<<" Trackable: NO\n";

			for (SizeType j = 0; j < dofs; ++j) {
				std::cout<<"Operator name="<<ll.name()<<" dof="<<j;
				PsimagLite::String desc = ll.description(j);
				std::cout<<" "<<desc<<"\n";
				ll(j).write(std::cout);
			}
		}
	}

	void printTerms() const
	{
		modelLinks_.printTerms(std::cout, labeledOperators_);
	}

	const SuperGeometryType& superGeometry() const
	{
		return modelCommon_.superGeometry();
	}

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

	SuperOpHelperBaseType& superOpHelper() const
	{
		return *superOpHelper_;
	}

	// protected:

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
		labeledOperators_.makeTrackable(name);
	}

	static ModelTermType& createTerm(PsimagLite::String name,
	                                 bool wantsHermitian = true,
	                                 PsimagLite::String geometryFrom = "")
	{
		return modelLinks_.createTerm(name, wantsHermitian, geometryFrom);
	}

	virtual SuperOpHelperBaseType* setSuperOpHelper()
	{
		return (superOpHelper_) ? superOpHelper_
		                        : new SuperOpHelperBaseType(modelCommon_.superGeometry());
	}

	static void notReallySort(VectorSizeType& basis, VectorQnType& qns)
	{
		const SizeType n = qns.size();
		if (n == 0) return;
		VectorQnType qunique;
		findQunique(qunique, qns);
		assert(qunique.size() > 0);

		VectorQnType qns2;
		VectorSizeType basis2;
		VectorBoolType done(n, false);
		const SizeType m = qunique.size();
		for (SizeType i = 0; i < m; ++i) {
			const QnType& thisq = qunique[i];
			for (SizeType j = 0; j < n; ++j) {
				if (done[j]) continue;

				// add all qns with thisq
				if (qns[j] == thisq) {
					qns2.push_back(thisq);
					basis2.push_back(basis[j]);
					done[j] = true;
				}
			}
		}

		basis = basis2;
		qns = qns2;
	}

protected:

	void onSiteHLegacyFix(const VectorStringType& legacyPotentialT)
	{
		if (onSiteHadd_.size() > 0) {
			if (legacyPotentialT.size() > 0)
				err("AddOnSiteHamiltonian: You cannot give both legacy and standard entries\n");
		} else {
			onSiteHadd_ = legacyPotentialT;
		}
	}

	void additionalOnSiteHamiltonian(SparseMatrixType& hmatrix,
	                                 const BlockType& block,
	                                 RealType time) const
	{
		if (onSiteHadd_.size() != superGeometry().numberOfSites())
			return;

		assert(block.size() == 1);
		const SizeType site = block[0];
		assert(onSiteHadd_.size() > site);
		const PsimagLite::String hOnSite = onSiteHadd_[site];
		if (hOnSite == "") return;

		OperatorSpecType opSpec(*this);
		CanonicalExpressionType canonicalExpression(opSpec);
		OperatorType hOft;
		OperatorType opEmpty;
		PsimagLite::String expression = CanonicalExpressionType::replaceAll(hOnSite,
		                                                                    "%t",
		                                                                    time).second;
		int bogus = 0;
		canonicalExpression(hOft, expression, opEmpty, bogus);
		hmatrix += hOft.getStorage().getCRS();
	}

private:

	void readOnSiteH(InputValidatorType& io)
	{
		PsimagLite::String tmp2;

		try {
			io.readline(tmp2, "AddOnSiteHamiltonian=");
			onSiteHadd_.resize(superGeometry().numberOfSites());
			stringToVectorOfStrings(onSiteHadd_, tmp2);
		} catch (std::exception&) {}
	}

	static void stringToVectorOfStrings(VectorStringType& vec, PsimagLite::String str)
	{
		if (str[0] == '[') {
			stringToVectorOfStringsCommaMode(vec, str);
		} else {
			stringToVectorOfStringsPlusMode(vec, str);
		}
	}

	static void stringToVectorOfStringsPlusMode(VectorStringType& vec, PsimagLite::String str)
	{
		str = ProgramGlobals::killSpaces(str);

		// break on plus
		const SizeType nsites = vec.size();
		VectorStringType tokens;
		PsimagLite::split(tokens, str, "+");
		const SizeType n = tokens.size();
		for (SizeType i = 0; i < n; ++i) {
			std::pair<PsimagLite::String, SizeType> oneSummand = getSiteAndContent(tokens[i]);
			const SizeType site = oneSummand.second;
			if (site >= nsites)
				err("You provided a site " + ttos(site) + " >= " + ttos(nsites) + "\n");
			vec[site] = oneSummand.first;
		}
	}

	static std::pair<PsimagLite::String, SizeType> getSiteAndContent(PsimagLite::String str)
	{
		const SizeType n = str.length();
		SizeType status = 0; // 0 = closed, 1 = open
		SizeType site = 0;
		bool foundSite = false;
		PsimagLite::String buffer;
		PsimagLite::String content;
		for (SizeType i = 0; i < n; ++i) {
			const char c = str[i];
			if (c == '[') {
				if (status == 1)
					err("Nested brakets found\n");
				status = 1; // open
				continue;
			} else if (c == ']') {
				if (status != 1)
					err("Closing braket without opening one\n");
				site = PsimagLite::atoi(buffer);
				buffer = "";
				foundSite = true;
				status = 0; // closing
				continue;
			}

			if (status == 1) buffer += c;
			else content += c;
		}

		if (!foundSite)
			err("A term for AddOnSiteHamiltonian was given without a site\n");

		return std::pair<PsimagLite::String, SizeType>(content, site);
	}

	static void stringToVectorOfStringsCommaMode(VectorStringType& vec, PsimagLite::String str)
	{
		const SizeType last = str.length() - 1;
		if (str.length() < 3 || str[0] != '[' || str[last] != ']')
			err("Expected [...] in comma mode\n");

		str = str.substr(1, str.length() - 2); // remove [ and ]

		// break on ,
		const SizeType nsites = vec.size();
		VectorStringType tokens;
		PsimagLite::split(tokens, str, ",");
		const SizeType n = tokens.size();
		if (n != nsites)
			err("Expected " + ttos(nsites) + " entries but got " + ttos(n) + "\n");
		vec.swap(tokens);
	}

	static void invalidateIfNeeded(OperatorType& op, SizeType site, PsimagLite::String what)
	{
		SizeType siteKind = modelLinks_.siteToAtomKind(site);
		SizeType opKind = labeledOperators_.findLabel(what).kindOfSite();

		if (siteKind == opKind) return;

		op.clear();
	}

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

	virtual void customOperators()
	{
		PsimagLite::String ops;
		try {
			ioIn_.readline(ops, "DefineOperators=");
		} catch (std::exception&) {
			return;
		}

		VectorStringType tokens;
		PsimagLite::split(tokens, ops, ",");
		const SizeType n = tokens.size();
		for (SizeType i = 0; i < n; ++i) {
			customOperator(tokens[i]);
		}
	}

	virtual void customOperator(PsimagLite::String opStr)
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, opStr, ":");
		if (tokens.size() != 2)
			err("Custom Operator " + opStr + " must have exactly one colon\n");

		if (tokens[0].length() == 0) return;
		const char firstChar = tokens[0][0];
		OpsLabelType& c = this->createOpsLabel(tokens[0]);
		BraketType braket(*this, "<gs|" + tokens[1] + "|gs>");

		c.push(braket.op(0));
		if (firstChar == '_')
			this->makeTrackable(tokens[0]);
	}

	static void findQunique(VectorQnType& qunique, const VectorQnType& qns)
	{
		qunique.clear();
		const SizeType n = qns.size();
		if (n == 0) return;
		QnType qprev = qns[0];
		qunique.push_back(qprev);
		for (SizeType i = 1; i < n; ++i) {
			QnType thisq = qns[i];
			if (thisq == qprev) continue;
			qunique.push_back(thisq);
			qprev = thisq;
		}
	}

	void checkThatQnsAreNotReallySorted() const
	{
#ifdef NDEBUG
		return;
#endif

		VectorQnType qunique;
		findQunique(qunique, qns_);
		if (qunique.size() == 0) return;

		const SizeType m = qunique.size();
		for (SizeType i = 0; i < m; ++i) {
			const QnType thisq = qunique[i];
			for (SizeType j = i + 1; j < m; ++j) {
				if (thisq != qunique[j]) continue;
				err("QNS of one site: not ordered: Model must order QNS\n");
			}
		}
	}

	static void createIdentity()
	{
		if (labeledOperators_.size() == 0)
			err("createIdentity: INTERNAL ERROR\n");

		const SizeType n = labeledOperators_[0].rows();
		MatrixType m(n, n);
		for (SizeType i = 0; i < n; ++i) m(i, i) = 1;
		typename OperatorType::Su2RelatedType su2related;
		OperatorType myOp(SparseMatrixType(m),
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  PairSizeType(0, 0),
		                  1,
		                  su2related);
		createOpsLabel("identity").push(myOp);
	}

	ModelCommonType modelCommon_;
	TargetQuantumElectronsType targetQuantum_;
	InputValidatorType_& ioIn_;
	VectorStringType onSiteHadd_;
	AtomKindBaseType* atomKind_;
	mutable SuperOpHelperBaseType* superOpHelper_;
	static LabeledOperatorsType labeledOperators_;
	static ModelLinksType modelLinks_;
	static VectorQnType qns_;
}; //class ModelBase

template<typename T1, typename T2, typename T3, typename T4>
typename ModelBase<T1, T2, T3, T4>::LabeledOperatorsType
ModelBase<T1, T2, T3, T4>::labeledOperators_;

template<typename T1, typename T2, typename T3, typename T4>
typename ModelBase<T1, T2, T3, T4>::ModelLinksType ModelBase<T1, T2, T3, T4>::modelLinks_;

template<typename T1, typename T2, typename T3, typename T4>
typename ModelBase<T1, T2, T3, T4>::VectorQnType ModelBase<T1, T2, T3, T4>::qns_;
} // namespace Dmrg
/*@}*/
#endif
