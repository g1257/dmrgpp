#ifndef MODEL_LINKS_H
#define MODEL_LINKS_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"

namespace Dmrg {

template<typename LabeledOperatorsType, typename GeometryType_>
class ModelLinks {

	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef std::pair<PsimagLite::String, PsimagLite::String> PairStringType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;
	typedef typename LabeledOperatorsType::OperatorType::RealType RealType_;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename LabeledOperatorsType::LabelType LabelType;

	class OneLink {

	public:

		typedef std::pair<char, char> PairCharType;

		OneLink(SizeType index1,
		        SizeType index2,
		        ProgramGlobals::FermionOrBosonEnum fermionOrBoson_,
		        char mod1,
		        char mod2,
		        SizeType angularMomentum_,
		        RealType_ angularFactor_,
		        SizeType category_,
		        RealType_ vModifier_)
		    : indices(PairSizeType(index1, index2)),
		      fermionOrBoson(fermionOrBoson_),
		      mods(PairCharType(mod1, mod2)),
		      angularMomentum(angularMomentum_),
		      angularFactor(angularFactor_),
		      category(category_),
		      vModifier(vModifier_)
		{}

		PairSizeType indices;
		ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
		PairCharType mods;
		SizeType angularMomentum;
		RealType_ angularFactor;
		SizeType category;
		RealType_ vModifier;
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

		void push(PsimagLite::String name1,
		          PsimagLite::String name2,
		          ProgramGlobals::FermionOrBosonEnum fermionOrBoson,
		          char mod1 = 'N',
		          char mod2 = 'C',
		          SizeType angularMomentum = 1,
		          RealType_ angularFactor = 1.0,
		          SizeType category = 0,
		          RealType_ vModifier = 1.0)
		{
			SizeType index1 = findIndexOfOp(name1);
			SizeType index2 = findIndexOfOp(name2);

			links_.push_back(OneLink(index1,
			                         index2,
			                         fermionOrBoson,
			                         mod1,
			                         mod2,
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

		const VectorSizeType& sites() const { return sites_; }

		const PsimagLite::String& name() const { return name_; }

	private:

		SizeType findIndexOfOp(PsimagLite::String) const
		{
			err("Wrong\n");
			return 0;
		}

		Term(const Term&);

		Term& operator=(const Term&);

		PsimagLite::String name_; // name of term, not name of operator
		VectorSizeType sites_;
		VectorOneLinkType links_;
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
	typedef typename LabeledOperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::RealType RealType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
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

	void postCtor(const LabeledOperatorsType& labeledOps,
	              SizeType geometryTerms)
	{
		if (terms_.size() > geometryTerms) {
			PsimagLite::String str("ModelBase: NumberOfTerms must be ");
			err(str + ttos(terms_.size()) + " in input file for this model\n");
		}

		SizeType n = labeledOps.trackables();
		for (SizeType i = 0; i < n; ++i) {
			const typename LabelType::PairStringSizeType& nameAndSite = labeledOps.trackables(i);
			const LabelType& ll = labeledOps.findLabel(nameAndSite.first,
			                                           nameAndSite.second);
			SizeType dofs = ll.dofs();

			for (SizeType j = 0; j < dofs; ++j)
				cm_.push_back(ll(j));
		}

		SizeType m = cm_.size();
		if (m == 0) return;
		hermit_.resize(m);
		for (SizeType i = 0; i < m; ++i)
			hermit_[i] = getHermitianProperty(cm_[i].data);

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
	VectorHermitianEnum hermit_;
	SizeType maxDofs_;
};
}
#endif // MODEL_LINKS_H
