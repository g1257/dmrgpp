#ifndef LINKPRODUCTBASE_H
#define LINKPRODUCTBASE_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"

namespace Dmrg {

template<typename LabeledOperatorsType, typename GeometryType_>
class LinkProductBase {

	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef std::pair<PsimagLite::String, PsimagLite::String> PairStringType;
	typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;
	typedef typename LabeledOperatorsType::RealType RealType_;
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
		        SizeType category_)
		    : indices(PairSizeType(index1, index2)),
		      fermionOrBoson(fermionOrBoson_),
		      mods(PairCharType(mod1, mod2)),
		      angularMomentum(angularMomentum_),
		      angularFactor(angularFactor_),
		      category(category_)
		{}

		PairSizeType indices;
		ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
		PairCharType mods;
		SizeType angularMomentum;
		RealType_ angularFactor;
		SizeType category;
	}; // OneLink

	class Term {

		typedef typename PsimagLite::Vector<OneLink>::Type VectorOneLinkType;

	public:

		// pair of sites should actually be pair of kinds of sites
		Term(PsimagLite::String name, // name of term, not name of operator
		     PairSizeType sites = PairSizeType(0,0))
		    : name_(name), sites_(sites)
		{}

		void push(PsimagLite::String name1,
		          PsimagLite::String name2,
		          ProgramGlobals::FermionOrBosonEnum fermionOrBoson,
		          char mod1 = 'N',
		          char mod2 = 'C',
		          SizeType angularMomentum = 1,
		          RealType_ angularFactor = 1.0,
		          SizeType category = 0)
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
			                         category));
		}

	private:

		//		SizeType findIndexOfOp(PsimagLite::String) const
		//		{

		//		}

		Term(const Term&);

		Term& operator=(const Term&);

		PsimagLite::String name_; // name of term, not name of operator
		PairSizeType sites_;
		VectorOneLinkType links_;
	};

	class IsValue {

	public:

		IsValue(PsimagLite::String name, VectorSizeType sites)
		    : name_(name), sites_(sites)
		{}

		bool operator()(const Term& term) const
		{
			return (term.name() == name_ && sites_ == term.sites());
		}

	private:

		PsimagLite::String name_;
		VectorSizeType sites_;
	};

public:

	enum HermitianEnum { HERMIT_NEITHER, HERMIT_PLUS, HERMIT_MINUS};

	typedef GeometryType_ GeometryType;
	typedef typename LabeledOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename LabeledOperatorsType::RealType RealType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename LabeledOperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<HermitianEnum>::Type VectorHermitianEnum;
	typedef typename PsimagLite::Vector<Term>::Type VectorTermType;

	void postCtor(const LabeledOperatorsType& labeledOps,
	              SizeType geometryTerms)
	{
		if (terms_.size() > geometryTerms) {
			PsimagLite::String str("ModelBase: NumberOfTerms must be ");
			err(str + ttos(terms_.size()) + " in input file for this model\n");
		}

		SizeType n = labeledOps.trackables();
		for (SizeType i = 0; i < n; ++i) {
			typename LabelType::PairStringSizeType nameAndSite = labeledOps.trackables(i);
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
	}

	Term& createTerm(PsimagLite::String name,
	                 VectorSizeType sites)
	{
		typename VectorTermType::const_iterator x = std::find_if(terms_.begin(),
		                                                         terms_.end(),
		                                                         IsValue(name, sites));

		if (x != terms_.end())
			err("Repeated term " + name + " sites=NOT DISPLAYED\n");

		terms_.push_back(Term(name, sites));
		return terms_[terms_.size() - 1];
	}

	// FIXME: For Immm and SDHS
	HermitianEnum getHermitianProperty(SizeType opsIndex, SizeType) const
	{
		assert(opsIndex < hermit_.size());
		return hermit_[opsIndex];
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
};
}
#endif // LINKPRODUCTBASE_H
