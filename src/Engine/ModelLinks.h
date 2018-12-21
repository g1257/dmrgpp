#ifndef MODEL_LINKS_H
#define MODEL_LINKS_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"

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
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;

	class OneLink {

	public:

		static void nullModifier(ComplexOrRealType&, void*) {}

		typedef void (*ModifierType)(ComplexOrRealType&, void*);

		OneLink(PairSizeType indices_,
		        PairSizeType orbs_,
		        ProgramGlobals::FermionOrBosonEnum fermionOrBoson_,
		        PairCharType mods_,
		        SizeType angularMomentum_,
		        RealType_ angularFactor_,
		        SizeType category_,
		        ModifierType vModifier_,
		        void* modifierData_)
		    : indices(indices_),
		      orbs(orbs_),
		      fermionOrBoson(fermionOrBoson_),
		      mods(mods_),
		      angularMomentum(angularMomentum_),
		      angularFactor(angularFactor_),
		      category(category_),
		      modifier(vModifier_),
		      modifierData(modifierData_)
		{}

		PairSizeType indices;
		PairSizeType orbs;
		ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
		PairCharType mods;
		SizeType angularMomentum;
		RealType_ angularFactor;
		SizeType category;
		ModifierType modifier;
		void* modifierData;
	}; // OneLink

	class Term {

		typedef typename PsimagLite::Vector<OneLink>::Type VectorOneLinkType;

	public:

		typedef OneLink OneLinkType;

		// pair of sites should actually be pair of kinds of sites
		Term(PsimagLite::String name, // name of term, not name of operator
		     const VectorSizeType& sites)
		    : name_(name), sites_(sites)
		{}

		static void init(const LabeledOperatorsType& l,
		                 const VectorOperatorType& cm,
		                 const VectorSizeType& offsets)
		{
			labeledOps_ = &l;
			cm_ = &cm;
			offsets_ = &offsets;
		}

		void push(PsimagLite::String name1,
		          PsimagLite::String name2,
		          SizeType dof1 = 0,
		          SizeType dof2 = 0,
		          char mod1 = 'N',
		          char mod2 = 'C',
		          SizeType angularMomentum = 0,
		          RealType_ angularFactor = 1.0,
		          SizeType category = 0,
		          typename OneLinkType::ModifierType vModifier = OneLink::nullModifier,
		          void* modifierData = 0)
		{
			assert(sites_.size() == 2);
			SizeType index1 = findIndexOfOp(name1, sites_[0], dof1);
			SizeType index2 = findIndexOfOp(name2, sites_[1], dof2);

			assert(cm_);
			const VectorOperatorType& cm = *cm_;
			assert(index1 < cm.size() && index2 < cm.size());
			ProgramGlobals::FermionOrBosonEnum fermionOrBoson = ProgramGlobals::BOSON;
			if (cm[index1].fermionSign < 0 && cm[index2].fermionSign < 0)
				fermionOrBoson = ProgramGlobals::FERMION;

			links_.push_back(OneLink(PairSizeType(index1, index2),
			                         PairSizeType(dof1, dof2),
			                         fermionOrBoson,
			                         PairCharType(mod1, mod2),
			                         angularMomentum,
			                         angularFactor,
			                         category,
			                         vModifier,
			                         modifierData));
		}

		SizeType size() const { return links_.size(); }

		const OneLinkType& operator()(SizeType dof) const
		{
			assert(dof < links_.size());
			return links_[dof];
		}

		const VectorSizeType& sites() const { return sites_; }

		const PsimagLite::String& name() const { return name_; }

	private:

		SizeType findIndexOfOp(PsimagLite::String name,
		                       SizeType site,
		                       SizeType dof) const
		{
			assert(labeledOps_);
			assert(offsets_);
			const VectorSizeType& offsets = *offsets_;
			SizeType n = offsets.size();
			for (SizeType i = 0; i < n; ++i) {
				const typename LabelType::PairStringSizeType& nameAndSite =
				        labeledOps_->trackables(i);

				if (nameAndSite.first == name && nameAndSite.second == site)
					return offsets[i] + dof;
			}

			throw PsimagLite::RuntimeError("Cannot find TRACKABLE operator " +
			                               name + " with dof " + ttos(dof) +
			                               " and site " + ttos(site) + "\n");
		}

		Term(const Term&);

		Term& operator=(const Term&);

		PsimagLite::String name_; // name of term, not name of operator
		VectorSizeType sites_;
		VectorOneLinkType links_;
		static const LabeledOperatorsType* labeledOps_;
		static const VectorSizeType* offsets_;
		static const VectorOperatorType* cm_;
	};

	class IsValue {

	public:

		IsValue(PsimagLite::String name, VectorSizeType sites)
		    : name_(name), sites_(sites)
		{}

		bool operator()(const Term* term) const
		{
			return (term->name() == name_ && sites_ == term->sites());
		}

	private:

		PsimagLite::String name_;
		VectorSizeType sites_;
	};

public:

	enum HermitianEnum { HERMIT_NEITHER, HERMIT_PLUS, HERMIT_MINUS};

	typedef GeometryType_ GeometryType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<HermitianEnum>::Type VectorHermitianEnum;
	typedef typename PsimagLite::Vector<Term*>::Type VectorTermType;
	typedef Term TermType;

	ModelLinks() : maxDofs_(0) {}

	~ModelLinks()
	{
		SizeType n = terms_.size();
		for (SizeType i = 0; i < n; ++i)
			delete terms_[i];

		terms_.clear();
	}

	void postCtor1(const LabeledOperatorsType& labeledOps,
	               SizeType geometryTerms)
	{
		if (terms_.size() > geometryTerms) {
			PsimagLite::String str("ModelBase: NumberOfTerms must be ");
			err(str + ttos(terms_.size()) + " in input file for this model\n");
		}

		SizeType n = labeledOps.trackables();
		offsets_.resize(n, 0);
		for (SizeType i = 0; i < n; ++i) {
			const typename LabelType::PairStringSizeType& nameAndSite = labeledOps.trackables(i);
			const LabelType& ll = labeledOps.findLabel(nameAndSite.first,
			                                           nameAndSite.second);
			SizeType dofs = ll.dofs();

			offsets_[i] = (i == 0) ? 0 : offsets_[i - 1] + ll.dofs();

			for (SizeType j = 0; j < dofs; ++j)
				cm_.push_back(ll(j));
		}

		SizeType m = cm_.size();
		if (m == 0) return;
		hermit_.resize(m);
		for (SizeType i = 0; i < m; ++i)
			hermit_[i] = getHermitianProperty(cm_[i].data);

		Term::init(labeledOps, cm_, offsets_);
	}

	void postCtor2()
	{
		SizeType t = terms_.size();
		for (SizeType i = 0; i < t; ++i) {
			SizeType dof = terms_[i]->size();
			if (dof > maxDofs_) maxDofs_ = dof;
		}
	}

	Term& createTerm(PsimagLite::String name,
	                 VectorSizeType sites)
	{
		typename VectorTermType::const_iterator x = std::find_if(terms_.begin(),
		                                                         terms_.end(),
		                                                         IsValue(name, sites));

		if (x != terms_.end())
			err("Repeated term " + name + " sites=NOT DISPLAYED\n");

		Term* term = new Term(name, sites);
		terms_.push_back(term);
		return *term;
	}

	// FIXME: For Immm and SDHS
	HermitianEnum getHermitianProperty(SizeType opsIndex, SizeType) const
	{
		assert(opsIndex < hermit_.size());
		return hermit_[opsIndex];
	}

	SizeType hilbertSize(SizeType) const
	{
		assert(cm_.size() > 0);
		return cm_[0].data.rows();
	}

	void setOperatorMatrices(VectorOperatorType& cm) const
	{
		cm = cm_;
	}

	SizeType dofsAllocationSize() const { return maxDofs_; }

	const TermType& term(SizeType term) const
	{
		assert(term < terms_.size());
		return *(terms_[term]);
	}

	const VectorOperatorType& trackableOps(SizeType) const
	{
		assert(cm_.size() > 0);
		return cm_;
	}

private:

	static HermitianEnum getHermitianProperty(const SparseMatrixType& m)
	{
		if (isHermitian(m)) return HERMIT_PLUS;
		return (isAntiHermitian(m)) ? HERMIT_MINUS : HERMIT_NEITHER;
	}

	VectorTermType terms_;
	VectorOperatorType cm_;
	VectorSizeType offsets_;
	VectorHermitianEnum hermit_;
	SizeType maxDofs_;
};

template<typename T1, typename T2>
const T1* ModelLinks<T1, T2>::Term::labeledOps_ = 0;

template<typename T1, typename T2>
const typename ModelLinks<T1, T2>::VectorSizeType* ModelLinks<T1, T2>::Term::offsets_ = 0;

template<typename T1, typename T2>
const typename ModelLinks<T1, T2>::VectorOperatorType* ModelLinks<T1, T2>::Term::cm_ = 0;
}
#endif // MODEL_LINKS_H
