#ifndef MODEL_LINKS_H
#define MODEL_LINKS_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "PsimagLite.h"
#include <functional>

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

	struct OpaqueOp {

		OpaqueOp(PsimagLite::String name_, SizeType dof_ = 0, SizeType edof_ = 0)
		    : name(name_), dof(dof_), edof(edof_)
		{}

		PsimagLite::String name;
		SizeType dof;
		SizeType edof;
	};

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

	class Term {

		typedef typename PsimagLite::Vector<OneLink>::Type VectorOneLinkType;

	public:

		typedef OneLink OneLinkType;
		typedef typename OneLinkType::LambdaType LambdaType;

		// pair of sites should actually be pair of kinds of sites
		Term(PsimagLite::String name) // name of term, not name of operator
		    : name_(name)
		{}

		void push(const OpaqueOp& op1,
		          char mod1,
		          const OpaqueOp& op2,
		          char mod2,
		          SizeType angularMomentum = 0,
		          RealType_ angularFactor = 1.0,
		          SizeType category = 0,
		          LambdaType vModifier = [](ComplexOrRealType&) {})
		{
			SizeType index1 = findIndexOfOp(op1.name, op1.dof);
			SizeType index2 = findIndexOfOp(op2.name, op2.dof);


			assert(index1 < cm_.size() && index2 < cm_.size());
			ProgramGlobals::FermionOrBosonEnum fermionOrBoson =
			        ProgramGlobals::FermionOrBosonEnum::BOSON;
			if (cm_[index1].fermionOrBoson() == ProgramGlobals::FermionOrBosonEnum::FERMION &&
			cm_[index2].fermionOrBoson() == ProgramGlobals::FermionOrBosonEnum::FERMION)
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
			SizeType n = offsets_.size();
			for (SizeType i = 0; i < n; ++i) {
				if (trackables_[i] == name)
					return offsets_[i] + dof;
			}

			throw PsimagLite::RuntimeError("Cannot find TRACKABLE operator " +
			                               name + " with dof " + ttos(dof) + "\n");
		}

		Term(const Term&);

		Term& operator=(const Term&);

		PsimagLite::String name_; // name of term, not name of operator
		VectorOneLinkType links_;
	};

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

		SizeType n = trackables_.size();
		offsets_.resize(n, 0);
		SizeType dofs = 0;
		for (SizeType i = 0; i < n; ++i) {
			const LabelType& ll = labeledOps.findLabel(trackables_[i]);

			offsets_[i] = (i == 0) ? 0 : offsets_[i - 1] + dofs;

			dofs = ll.dofs();

			for (SizeType j = 0; j < dofs; ++j) {
				cm_.push_back(ll(j));
				cmNameDof_.push_back(typename LabelType::PairStringSizeType(trackables_[i],
				                                                            j));
			}
		}

		SizeType m = cm_.size();
		if (m == 0) return;
		hermit_.resize(m);
		for (SizeType i = 0; i < m; ++i)
			hermit_[i] = getHermitianProperty(cm_[i].getStorage());
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
		                                                         name);

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

	SizeType hilbertSize(SizeType) const
	{
		assert(cm_.size() > 0);
		return cm_[0].getStorage().rows();
	}

	SizeType dofsAllocationSize() const { return maxDofs_; }

	const TermType& term(SizeType term) const
	{
		assert(term < terms_.size());
		return *(terms_[term]);
	}

	// FIXME TODO SDHS Immm: add type of site to arguments here
	SizeType nameDofToIndex(PsimagLite::String name,
	                        SizeType dof) const
	{
		PairStringSizeType p(name, dof);
		auto result = std::find(cmNameDof_.begin(), cmNameDof_.end(), p);
		if (result == cmNameDof_.end())
			err("siteNameDofToIndex: not found for name=" + name + " dof= " +
			    ttos(dof) + "\n");
		return (result - cmNameDof_.begin());
	}

	const VectorOperatorType& cm() const { return cm_; }

private:

	static HermitianEnum getHermitianProperty(const OperatorStorageType& m)
	{
		if (isHermitian(m)) return HERMIT_PLUS;
		return (isAntiHermitian(m)) ? HERMIT_MINUS : HERMIT_NEITHER;
	}

	VectorTermType terms_;
	VectorHermitianEnum hermit_;
	SizeType maxDofs_;
	static VectorOperatorType cm_;
	static VectorSizeType offsets_;
	static VectorStringType trackables_;
	static VectorPairStringSizeType cmNameDof_;
};

template<typename T1, typename T2>
typename ModelLinks<T1, T2>::VectorSizeType ModelLinks<T1, T2>::offsets_;

template<typename T1, typename T2>
typename ModelLinks<T1, T2>::VectorOperatorType ModelLinks<T1, T2>::cm_;

template<typename T1, typename T2>
typename ModelLinks<T1, T2>::VectorStringType ModelLinks<T1, T2>::trackables_;

template<typename T1, typename T2>
typename ModelLinks<T1, T2>::VectorPairStringSizeType ModelLinks<T1, T2>::cmNameDof_;
}
#endif // MODEL_LINKS_H
