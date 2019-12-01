#ifndef MODEL_LINKS_H
#define MODEL_LINKS_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include <functional>
#include <map>

namespace Dmrg {

template<typename LabeledOperatorsType, typename GeometryType_>
class ModelLinks {

	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef std::pair<char, char> PairCharType;
	typedef std::pair<PsimagLite::String, PsimagLite::String> PairStringType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;
	typedef typename LabeledOperatorsType::OperatorType::RealType RealType_;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename LabeledOperatorsType::LabelType LabelType;
	typedef typename LabeledOperatorsType::ComplexOrRealType ComplexOrRealType;
	typedef typename LabeledOperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::RealType RealType;
	typedef typename OperatorType::StorageType OperatorStorageType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<typename LabelType::PairStringSizeType>::Type
	VectorPairStringSizeType;

public:

	class OneLink {

	public:

		typedef std::function<void(ComplexOrRealType&)> LambdaType;

		OneLink(PairSizeType indices_,
		        PairSizeType orbs_,
		        ProgramGlobals::FermionOrBosonEnum fermionOrBoson_,
		        PairCharType mods_,
		        SizeType angularMomentum_,
		        RealType_ angularFactor_,
		        SizeType category_,
		        LambdaType vModifier_)
		    : indices(indices_),
		      orbs(orbs_),
		      fermionOrBoson(fermionOrBoson_),
		      mods(mods_),
		      angularMomentum(angularMomentum_),
		      angularFactor(angularFactor_),
		      category(category_),
		      modifier(vModifier_)
		{}

		PairSizeType indices;
		PairSizeType orbs;
		ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
		PairCharType mods;
		SizeType angularMomentum;
		RealType_ angularFactor;
		SizeType category;
		LambdaType modifier;
	}; // OneLink

	class AtomKindBase {

	public:

		virtual ~AtomKindBase() {}

		virtual SizeType siteToAtomKind(SizeType) const  { return 0; }

		virtual SizeType kindsOfAtoms() const { return 1; }
	};

	class Term {

		typedef typename PsimagLite::Vector<OneLink>::Type VectorOneLinkType;

	public:

		typedef OneLink OneLinkType;
		typedef typename OneLinkType::LambdaType LambdaType;

		// pair of sites should actually be pair of kinds of sites
		Term(PsimagLite::String name) // name of term, not name of operator
		    : name_(name), pairKind_(0, 0)
		{}

		bool areSitesCompatible(SizeType kind0, SizeType kind1) const
		{
			return (pairKind_.first == kind0 && pairKind_.second == kind1);
		}

		bool areSiteKindsEqual() const
		{
			return (pairKind_.first == pairKind_.second);
		}

		template<typename OpaqueOp>
		void push(const OpaqueOp& op1,
		          char mod1,
		          const OpaqueOp& op2,
		          char mod2,
		          SizeType angularMomentum = 0,
		          RealType_ angularFactor = 1.0,
		          SizeType category = 0
		        ,LambdaType vModifier = [](ComplexOrRealType&) {})
		{
			if (links_.size() > 0) {
				if (!areSitesCompatible(op1.kindOfSite, op2.kindOfSite))
					err("Term " + name_ + " incompatible atom kinds at push\n");
			} else {
				pairKind_ = PairSizeType(op1.kindOfSite, op2.kindOfSite);
			}

			SizeType index1 = findIndexOfOp(op1.name, op1.dof);
			SizeType index2 = findIndexOfOp(op2.name, op2.dof);

			ProgramGlobals::FermionOrBosonEnum fermionOrBoson =
			        ProgramGlobals::FermionOrBosonEnum::BOSON;
			if (op1.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION &&
			        op2.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION)
				fermionOrBoson = ProgramGlobals::FermionOrBosonEnum::FERMION;
			// can we also infer angularMomentum, angularFactor, and category? FIXME TODO

			links_.push_back(OneLink(PairSizeType(index1, index2),
			                         PairSizeType(op1.edof, op2.edof),
			                         fermionOrBoson,
			                         PairCharType(mod1, mod2),
			                         angularMomentum,
			                         angularFactor,
			                         category,
			                         vModifier));
		}

		SizeType size() const { return links_.size(); }

		const OneLinkType& operator()(SizeType dof) const
		{
			assert(dof < links_.size());
			return links_[dof];
		}

		const PsimagLite::String& name() const { return name_; }

	private:

		SizeType findIndexOfOp(PsimagLite::String name, SizeType dof) const
		{
			return offsets_[name] + dof;
		}

		Term(const Term&);

		Term& operator=(const Term&);

		PsimagLite::String name_; // name of term, not name of operator
		VectorOneLinkType links_;
		PairSizeType pairKind_;
	};

	class IsValue {

	public:

		IsValue(PsimagLite::String name)
		    : name_(name) {}

		bool operator()(const Term* term) const
		{
			return (term->name() == name_);
		}

	private:

		PsimagLite::String name_;
	};

	enum HermitianEnum { HERMIT_NEITHER, HERMIT_PLUS, HERMIT_MINUS};

	typedef GeometryType_ GeometryType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<HermitianEnum>::Type VectorHermitianEnum;
	typedef typename PsimagLite::Vector<Term*>::Type VectorTermType;
	typedef Term TermType;

	ModelLinks() : maxDofs_(0), atomKind_(0) {}

	~ModelLinks()
	{
		SizeType n = terms_.size();
		for (SizeType i = 0; i < n; ++i)
			delete terms_[i];

		terms_.clear();
	}

	void postCtor1(const LabeledOperatorsType& labeledOps,
	               const AtomKindBase* ptr,
	               SizeType geometryTerms)
	{
		if (terms_.size() > geometryTerms) {
			PsimagLite::String str("ModelBase: NumberOfTerms must be ");
			err(str + ttos(terms_.size()) + " in input file for this model\n");
		}

		atomKind_ = ptr;
		VectorOperatorType cm; // only for hermit
		SizeType n = trackables_.size();
		VectorSizeType dofsByKind(n);
		for (SizeType i = 0; i < n; ++i) {
			const LabelType& ll = labeledOps.findLabel(trackables_[i]);

			const SizeType kindOfSite = ll.kindOfSite();

			offsets_[ll.name()] = (i == 0) ? 0 : dofsByKind[kindOfSite];

			dofsByKind[kindOfSite] += ll.dofs();

			for (SizeType j = 0; j < ll.dofs(); ++j)
				cm.push_back(ll(j));
		}

		SizeType m = cm.size();
		if (m == 0) return;
		hermit_.resize(m);
		for (SizeType i = 0; i < m; ++i)
			hermit_[i] = getHermitianProperty(cm[i].getStorage());
	}

	void postCtor2()
	{
		SizeType t = terms_.size();
		for (SizeType i = 0; i < t; ++i) {
			SizeType dof = terms_[i]->size();
			if (dof > maxDofs_) maxDofs_ = dof;
		}
	}

	void makeTrackable(PsimagLite::String what)
	{
		if (std::find(trackables_.begin(), trackables_.end(), what) != trackables_.end())
			err("makeTrackable: cannot find operator " + what + "\n");
		trackables_.push_back(what);
	}

	Term& createTerm(PsimagLite::String name)
	{
		typename VectorTermType::const_iterator x = std::find_if(terms_.begin(),
		                                                         terms_.end(),
		                                                         IsValue(name));

		if (x != terms_.end())
			err("Repeated term " + name + "\n");

		Term* term = new Term(name);
		terms_.push_back(term);
		return *term;
	}

	// FIXME: For Immm and SDHS
	HermitianEnum getHermitianProperty(SizeType opsIndex, SizeType) const
	{
		assert(opsIndex < hermit_.size());
		return hermit_[opsIndex];
	}

	SizeType dofsAllocationSize() const { return maxDofs_; }

	const TermType& term(SizeType term) const
	{
		assert(term < terms_.size());
		return *(terms_[term]);
	}

	void setOperatorMatrices(VectorOperatorType& cm,
	                         const LabeledOperatorsType& labeledOps,
	                         SizeType kindOfSite) const
	{
		cm.clear();
		SizeType n = trackables_.size();
		for (SizeType i = 0; i < n; ++i) {
			const LabelType& ll = labeledOps.findLabel(trackables_[i]);

			if (ll.kindOfSite() != kindOfSite)
				continue;

			const SizeType dofs = ll.dofs();
			for (SizeType j = 0; j < dofs; ++j)
				cm.push_back(ll(j));
		}
	}

	SizeType hilbertSize(SizeType kindOfSite,
	                     const LabeledOperatorsType& labeledOps) const
	{
		SizeType n = trackables_.size();
		for (SizeType i = 0; i < n; ++i) {
			const LabelType& ll = labeledOps.findLabel(trackables_[i]);

			if (ll.kindOfSite() != kindOfSite)
				continue;

			return ll(0).getStorage().rows();
		}

		throw PsimagLite::RuntimeError("hilbertSize FATAL: " + ttos(kindOfSite) + "\n");
	}

	bool areSiteKindsEqual(SizeType termIndex) const
	{
		assert(termIndex < terms_.size());
		return terms_[termIndex]->areSiteKindsEqual();
	}

	bool areSitesCompatibleForThisTerm(SizeType termIndex,
	                                   SizeType actualSite0,
	                                   SizeType actualSite1) const
	{
		assert(atomKind_);
		const SizeType kind0 = atomKind_->siteToAtomKind(actualSite0);
		const SizeType kind1 = atomKind_->siteToAtomKind(actualSite1);

		assert(termIndex < terms_.size());
		return terms_[termIndex]->areSitesCompatible(kind0, kind1);
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

private:

	static HermitianEnum getHermitianProperty(const OperatorStorageType& m)
	{
		if (isHermitian(m)) return HERMIT_PLUS;
		return (isAntiHermitian(m)) ? HERMIT_MINUS : HERMIT_NEITHER;
	}

	VectorTermType terms_;
	VectorHermitianEnum hermit_;
	SizeType maxDofs_;
	const AtomKindBase* atomKind_;
	static std::map<PsimagLite::String, SizeType> offsets_;
	static VectorStringType trackables_;
};

template<typename T1, typename T2>
std::map<PsimagLite::String, SizeType> ModelLinks<T1, T2>::offsets_;

template<typename T1, typename T2>
typename ModelLinks<T1, T2>::VectorStringType ModelLinks<T1, T2>::trackables_;
}
#endif // MODEL_LINKS_H
