#ifndef MODEL_LINKS_H
#define MODEL_LINKS_H
#include "OneLink.hh"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include "Vector.h"
#include <map>

namespace Dmrg {

template <typename LabeledOperatorsType, typename SuperGeometryType_> class ModelLinks {

public:

	typedef std::pair<SizeType, SizeType>                         PairSizeType;
	typedef std::pair<char, char>                                 PairCharType;
	typedef std::pair<PsimagLite::String, PsimagLite::String>     PairStringType;
	typedef std::pair<PsimagLite::String, SizeType>               PairStringSizeType;
	typedef typename LabeledOperatorsType::OperatorType::RealType RealType_;
	typedef PsimagLite::Vector<SizeType>::Type                    VectorSizeType;
	typedef typename LabeledOperatorsType::LabelType              LabelType;
	typedef typename LabeledOperatorsType::ComplexOrRealType      ComplexOrRealType;
	typedef typename LabeledOperatorsType::OperatorType           OperatorType;
	typedef typename OperatorType::RealType                       RealType;
	typedef typename OperatorType::StorageType                    OperatorStorageType;
	typedef typename PsimagLite::Vector<OperatorType>::Type       VectorOperatorType;
	typedef typename PsimagLite::Vector<typename LabelType::PairStringSizeType>::Type
	    VectorPairStringSizeType;
	using OneLinkType = OneLink<ComplexOrRealType>;

	class AtomKindBase {

	public:

		virtual ~AtomKindBase() { }

		virtual SizeType siteToAtomKind(SizeType) const { return 0; }

		virtual SizeType kindsOfAtoms() const { return 1; }
	};

	class Term {

		using VectorOneLinkType = std::vector<OneLinkType>;

	public:

		struct Su2Properties {
			Su2Properties(SizeType a = 0, RealType_ f = 1.0, SizeType c = 0)
			    : angularMomentum(a)
			    , angularFactor(f)
			    , category(c)
			{ }

			SizeType  angularMomentum;
			RealType_ angularFactor;
			SizeType  category;
		};

		using LambdaType    = typename OneLinkType::LambdaType;
		using OldLambdaType = typename OneLinkType::OldLambdaType;

		// pair of sites should actually be pair of kinds of sites
		Term(PsimagLite::String name, bool wantsHermitian = true) // name of term,
		                                                          // not name of operator
		    : name_(name)
		    , wantsHermitian_(wantsHermitian)
		{ }

		bool wantsHermitian() const { return wantsHermitian_; }

		bool areSitesCompatible(const VectorSizeType& actualSites) const
		{
			const SizeType n = actualSites.size();
			assert(n == vectorKind_.size());
			for (SizeType i = 0; i < n; ++i) {
				if (vectorKind_[i] != atomKind_->siteToAtomKind(actualSites[i]))
					return false;
			}

			return true;
		}

		bool areSitesCompatible2(const VectorSizeType& kinds) const
		{
			const SizeType n = kinds.size();
			assert(n == vectorKind_.size());
			for (SizeType i = 0; i < n; ++i) {
				if (vectorKind_[i] != kinds[i])
					return false;
			}

			return true;
		}

		// give only su2properties
		template <typename OpaqueOp>
		void push(const OpaqueOp& op1,
		          char            mod1,
		          const OpaqueOp& op2,
		          char            mod2,
		          Su2Properties   su2properties)
		{
			push(
			    op1,
			    mod1,
			    op2,
			    mod2,
			    [](ComplexOrRealType&, RealType, SizeType) { },
			    su2properties);
		}

		// give only lambda (new)
		template <typename OpaqueOp>
		void push(const OpaqueOp& op1,
		          char            mod1,
		          const OpaqueOp& op2,
		          char            mod2,
		          LambdaType      modifier)
		{
			push(op1, mod1, op2, mod2, modifier, Su2Properties());
		}

		// give only lambda (old)
		template <typename OpaqueOp>
		void push(const OpaqueOp& op1,
		          char            mod1,
		          const OpaqueOp& op2,
		          char            mod2,
		          OldLambdaType   modifier)
		{
			LambdaType newModif
			    = [modifier](ComplexOrRealType& value, RealType, SizeType)
			{ modifier(value); };
			push(op1, mod1, op2, mod2, newModif, Su2Properties());
		}

		// give nothing
		template <typename OpaqueOp>
		void push(const OpaqueOp& op1, char mod1, const OpaqueOp& op2, char mod2)
		{
			return push(
			    op1,
			    mod1,
			    op2,
			    mod2,
			    [](ComplexOrRealType&, RealType, SizeType) { },
			    Su2Properties());
		}

		// give all
		template <typename OpaqueOp>
		void push(const OpaqueOp& op1,
		          char            mod1,
		          const OpaqueOp& op2,
		          char            mod2,
		          LambdaType      vModifier,
		          Su2Properties   su2properties)
		{
			if (links_.size() > 0) {
				if (!areSitesCompatible2(
				        VectorSizeType { op1.kindOfSite, op2.kindOfSite }))
					err("Term " + name_ + " incompatible atom kinds at push\n");
			} else {
				vectorKind_ = VectorSizeType { op1.kindOfSite, op2.kindOfSite };
			}

			SizeType index1 = findIndexOfOp(op1.name, op1.dof);
			SizeType index2 = findIndexOfOp(op2.name, op2.dof);

			ProgramGlobals::FermionOrBosonEnum fermionOrBoson
			    = ProgramGlobals::FermionOrBosonEnum::BOSON;
			if (op1.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION
			    && op2.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION)
				fermionOrBoson = ProgramGlobals::FermionOrBosonEnum::FERMION;
			// can we also infer angularMomentum, angularFactor, and category? FIXME
			// TODO

			PsimagLite::String modStr("NN");
			modStr[0] = mod1;
			modStr[1] = mod2;
			links_.push_back(OneLink(VectorSizeType { index1, index2 },
			                         VectorSizeType { op1.edof, op2.edof },
			                         fermionOrBoson,
			                         modStr,
			                         su2properties.angularMomentum,
			                         su2properties.angularFactor,
			                         su2properties.category,
			                         vModifier));
		}

		template <typename OpaqueOp>
		void push4(
		    const OpaqueOp& op1,
		    char            mod1,
		    const OpaqueOp& op2,
		    char            mod2,
		    const OpaqueOp& op3,
		    char            mod3,
		    const OpaqueOp& op4,
		    char            mod4,
		    OldLambdaType   vModifier     = [](ComplexOrRealType&) { },
		    Su2Properties   su2properties = Su2Properties())
		{
			if (links_.size() > 0) {
				if (!areSitesCompatible2(VectorSizeType { op1.kindOfSite,
				                                          op2.kindOfSite,
				                                          op3.kindOfSite,
				                                          op4.kindOfSite }))
					err("Term " + name_ + " incompatible atom kinds at push\n");
			} else {
				vectorKind_ = VectorSizeType { op1.kindOfSite,
					                       op2.kindOfSite,
					                       op3.kindOfSite,
					                       op4.kindOfSite };
			}

			SizeType index1 = findIndexOfOp(op1.name, op1.dof);
			SizeType index2 = findIndexOfOp(op2.name, op2.dof);
			SizeType index3 = findIndexOfOp(op3.name, op2.dof);
			SizeType index4 = findIndexOfOp(op3.name, op2.dof);

			ProgramGlobals::FermionOrBosonEnum fermionOrBoson
			    = ProgramGlobals::FermionOrBosonEnum::BOSON;
			// FIXME:
			if (op1.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION
			    || op2.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION
			    || op3.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION
			    || op4.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION)
				err(PsimagLite::String(__FILE__)
				    + "::push4(): Unsupported fermionic ops\n");

			// can we also infer angularMomentum, angularFactor, and category? FIXME
			// TODO

			PsimagLite::String modStr("NNNN");
			modStr[0] = mod1;
			modStr[1] = mod2;
			modStr[2] = mod3;
			modStr[3] = mod4;
			links_.push_back(
			    OneLink(VectorSizeType { index1, index2, index3, index4 },
			            VectorSizeType { op1.edof, op2.edof, op3.edof, op4.edof },
			            fermionOrBoson,
			            modStr,
			            su2properties.angularMomentum,
			            su2properties.angularFactor,
			            su2properties.category,
			            vModifier));
		}

		SizeType size() const { return links_.size(); }

		const OneLinkType& operator()(SizeType dof) const
		{
			assert(dof < links_.size());
			return links_[dof];
		}

		const PsimagLite::String& name() const { return name_; }

		void print(std::ostream& os, const LabeledOperatorsType& labeledOps) const
		{
			os << "Term Name=" << name_ << " Dofs=" << links_.size() << "\n";
			const SizeType n = links_.size();
			for (SizeType i = 0; i < n; ++i) {
				const OneLinkType&    onelink = links_[i];
				const VectorSizeType& indices = onelink.indices;
				const SizeType        m       = indices.size();
				PsimagLite::String    fermOrBos
				    = (onelink.fermionOrBoson
				       == ProgramGlobals::FermionOrBosonEnum::FERMION)
				    ? "[Fermionic]"
				    : "[Bosonic]";
				os << "\t" << fermOrBos << "    ";

				for (SizeType j = 0; j < m; ++j) {
					const SizeType index = indices[j];

					assert(j < onelink.mods.length());
					const PairSizeType lPair
					    = findOperatorIndex(index, labeledOps);
					const PsimagLite::String opName
					    = labeledOps[lPair.first].name();
					const PsimagLite::String dof
					    = (lPair.second == 0) ? "" : "?" + ttos(lPair.second);
					const char               modChar = onelink.mods[j];
					const PsimagLite::String mod = (modChar == 'N') ? "" : "'";
					const SizeType           orbital = onelink.orbs[j];
					const PsimagLite::String orbitalStr
					    = (orbital == 0) ? "" : "!" + ttos(orbital);
					const SizeType kind = labeledOps[lPair.first].kindOfSite();
					const PsimagLite::String kindStr
					    = (kind == 0) ? "" : "[@" + ttos(kind) + "]";
					os << "\t" << opName << dof << orbitalStr << mod << kindStr
					   << " ";
				}

				os << "\n";
			}
		}

	private:

		SizeType findIndexOfOp(PsimagLite::String name, SizeType dof) const
		{
			return offsets_[name] + dof;
		}

		static PairSizeType findOperatorIndex(SizeType                    index,
		                                      const LabeledOperatorsType& labeledOps)
		{
			const SizeType n = labeledOps.size();
			SizeType       k = 0;
			for (SizeType i = 0; i < n; ++i) {
				if (!labeledOps[i].isTrackable())
					continue;

				const SizeType dofs = labeledOps[i].dofs();
				for (SizeType j = 0; j < dofs; ++j)
					if (k++ == index)
						return PairSizeType(i, j);
			}

			throw PsimagLite::RuntimeError("findOperatorName: Not found for index "
			                               + ttos(index) + "\n");
		}

		Term(const Term&);

		Term& operator=(const Term&);

		PsimagLite::String name_; // name of term, not name of operator
		const bool         wantsHermitian_;
		VectorOneLinkType  links_;
		VectorSizeType     vectorKind_;
	};

	class IsValue {

	public:

		IsValue(PsimagLite::String name)
		    : name_(name)
		{ }

		bool operator()(const Term* term) const { return (term->name() == name_); }

	private:

		PsimagLite::String name_;
	};

	enum HermitianEnum
	{
		HERMIT_NEITHER,
		HERMIT_PLUS,
		HERMIT_MINUS
	};

	typedef SuperGeometryType_                               SuperGeometryType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type     VectorStringType;
	typedef typename PsimagLite::Vector<HermitianEnum>::Type VectorHermitianEnum;
	typedef typename PsimagLite::Vector<Term*>::Type         VectorTermType;
	typedef Term                                             TermType;

	ModelLinks()
	    : maxDofs_(0)
	{ }

	~ModelLinks() { clear(); }

	void clear()
	{
		// offsets_.clear();
		SizeType n = terms_.size();
		for (SizeType i = 0; i < n; ++i)
			delete terms_[i];

		terms_.clear();
		hermit_.clear();
		termGeomReplacement_.clear();
		hilbert_.clear();
	}

	void setAtomKind(const AtomKindBase* ptr) { atomKind_ = ptr; }

	void postCtor1(const LabeledOperatorsType& labeledOps, SizeType geometryTerms)
	{
		if (terms_.size() > geometryTerms) {
			PsimagLite::String str("ModelBase: NumberOfTerms must be ");
			err(str + ttos(terms_.size()) + " in input file for this model\n");
		}

		VectorOperatorType cm; // only for hermit
		SizeType           n = labeledOps.size();
		VectorSizeType     dofsByKind(n);
		hilbert_.resize(atomKind_->kindsOfAtoms());
		for (SizeType i = 0; i < n; ++i) {
			const LabelType& ll = labeledOps[i];

			if (!ll.isTrackable())
				continue;

			const SizeType kindOfSite = ll.kindOfSite();

			offsets_[ll.name()] = (i == 0) ? 0 : dofsByKind[kindOfSite];

			dofsByKind[kindOfSite] += ll.dofs();

			for (SizeType j = 0; j < ll.dofs(); ++j)
				cm.push_back(ll(j));

			hilbert_[kindOfSite] = ll(0).getCRS().rows();
		}

		SizeType m = cm.size();
		if (m == 0)
			return;
		hermit_.resize(m);
		for (SizeType i = 0; i < m; ++i)
			hermit_[i] = getHermitianProperty(cm[i].getStorage());
	}

	void postCtor2()
	{
		const SizeType   fromModel = terms_.size();
		VectorStringType input;
		for (SizeType termIndex = 0; termIndex < fromModel; ++termIndex) {
			const SizeType termIndexForGeom = termIndexForGeometry(termIndex);
			if (termIndex == termIndexForGeom) {
				input.push_back(terms_[termIndex]->name());
			}
		}

		VectorSizeType termGeomReplacement(fromModel);
		for (SizeType termIndex = 0; termIndex < fromModel; ++termIndex) {
			const SizeType     termIndexForGeom = termIndexForGeometry(termIndex);
			PsimagLite::String name             = terms_[termIndexForGeom]->name();
			typename VectorStringType::const_iterator x
			    = std::find(input.begin(), input.end(), name);
			if (x == input.end())
				err("ModelLinks: INTERNAL ERROR term " + name
				    + " termIndex= " + ttos(termIndex) + "\n");

			termGeomReplacement[termIndex] = x - input.begin();
		}

		termGeomReplacement_.swap(termGeomReplacement);

		SizeType t = terms_.size();
		for (SizeType i = 0; i < t; ++i) {
			SizeType dof = terms_[i]->size();
			if (dof > maxDofs_)
				maxDofs_ = dof;
		}
	}

	Term&
	createTerm(PsimagLite::String name, bool wantsHermitian, PsimagLite::String geometryFrom)
	{
		typename VectorTermType::const_iterator x
		    = std::find_if(terms_.begin(), terms_.end(), IsValue(name));

		if (x != terms_.end())
			err("Repeated term " + name + "\n");

		termGeomReplacement_.push_back(terms_.size());
		if (geometryFrom != "") {
			typename VectorTermType::const_iterator y
			    = std::find_if(terms_.begin(), terms_.end(), IsValue(geometryFrom));
			termGeomReplacement_[termGeomReplacement_.size() - 1] = y - terms_.begin();
		}

		Term* term = new Term(name, wantsHermitian);
		terms_.push_back(term);

		assert(termGeomReplacement_.size() == terms_.size());
		return *term;
	}

	// FIXME: For Immm and SDHS
	HermitianEnum getHermitianProperty(SizeType opsIndex) const
	{
		assert(opsIndex < hermit_.size());
		return hermit_[opsIndex];
	}

	SizeType dofsAllocationSize() const { return maxDofs_; }

	SizeType termIndexForGeometry(SizeType termIndex) const
	{
		assert(termIndex < termGeomReplacement_.size());
		return termGeomReplacement_[termIndex];
	}

	SizeType numberOfTerms() const { return terms_.size(); }

	const TermType& term(SizeType term) const
	{
		assert(term < terms_.size());
		return *(terms_[term]);
	}

	static void setOperatorMatrices(VectorOperatorType&         cm,
	                                const LabeledOperatorsType& labeledOps,
	                                SizeType                    kindOfSite)
	{
		cm.clear();
		SizeType n = labeledOps.size();
		for (SizeType i = 0; i < n; ++i) {
			const LabelType& ll = labeledOps[i];

			if (!ll.isTrackable())
				continue;

			if (ll.kindOfSite() != kindOfSite)
				continue;

			const SizeType dofs = ll.dofs();
			for (SizeType j = 0; j < dofs; ++j)
				cm.push_back(ll(j));
		}
	}

	SizeType hilbertSize(SizeType kindOfSite) const
	{
		assert(kindOfSite < hilbert_.size());
		return hilbert_[kindOfSite];
	}

	bool areSitesCompatibleForThisTerm(SizeType              termIndex,
	                                   const VectorSizeType& actualSites) const
	{
		assert(atomKind_);
		assert(termIndex < terms_.size());
		return terms_[termIndex]->areSitesCompatible(actualSites);
	}

	SizeType siteToAtomKind(SizeType actualSite) const
	{
		assert(atomKind_);
		return atomKind_->siteToAtomKind(actualSite);
	}

	SizeType kindsOfAtoms() const
	{
		assert(atomKind_);
		return atomKind_->kindsOfAtoms();
	}

	void printTerms(std::ostream& os, const LabeledOperatorsType& labeledOps) const
	{
		const SizeType n = terms_.size();
		os << "Model " << labeledOps.modelName() << " has " << n << " Hamiltonian terms\n";
		for (SizeType i = 0; i < n; ++i) {
			os << "Term " << i << " ";
			terms_[i]->print(os, labeledOps);
		}
	}

private:

	static HermitianEnum getHermitianProperty(const OperatorStorageType& m)
	{
		if (isHermitian(m))
			return HERMIT_PLUS;
		return (isAntiHermitian(m)) ? HERMIT_MINUS : HERMIT_NEITHER;
	}

	VectorTermType                                terms_;
	VectorHermitianEnum                           hermit_;
	SizeType                                      maxDofs_;
	VectorSizeType                                termGeomReplacement_;
	const static AtomKindBase*                    atomKind_;
	static std::map<PsimagLite::String, SizeType> offsets_;
	VectorSizeType                                hilbert_;
};

template <typename T1, typename T2>
std::map<PsimagLite::String, SizeType> ModelLinks<T1, T2>::offsets_;

template <typename T1, typename T2>
const typename ModelLinks<T1, T2>::AtomKindBase* ModelLinks<T1, T2>::atomKind_ = 0;
}
#endif // MODEL_LINKS_H
