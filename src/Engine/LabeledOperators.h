#ifndef LABELEDOPERATORS_H
#define LABELEDOPERATORS_H
#include "Vector.h"
#include "TypeToString.h"

namespace Dmrg {

template<typename OperatorType_>
class LabeledOperators {

	class Label {

		typedef typename PsimagLite::Vector<OperatorType_>::Type VectorOperatorType;

	public:

		typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;

		// site should actually be kind of site
		Label(PsimagLite::String name, SizeType site)
		    : name_(name), site_(site) {}

		const OperatorType_& operator()(SizeType dof) const
		{
			if (dof < ops_.size())
				return ops_[dof];

			PsimagLite::String msg("FATAL: LabeledOperators:");
			msg += " dof=" + ttos(dof) + " OUT OF RANGE, for label = " + name_ + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		bool operator==(PairStringSizeType nameAndSite) const
		{
			return (nameAndSite.first == name_ && nameAndSite.second == site_);
		}

		SizeType rows() const
		{
			if (ops_.size() == 0)
				err("FATAL: LabeledOperators::Label::rows(): Internal Error\n");
			return ops_[0].data.rows();
		}

		void push(const OperatorType_& op)
		{
			ops_.push_back(op);
		}

		void instrospect() const
		{
			std::cout<<"Label "<<name_<<" site="<<site_<<" with "<<ops_.size()<<" dofs.\n";
		}

		SizeType dofs() const { return ops_.size(); }

		SizeType site() const { return site_; }

	private:

		Label(const Label&);

		Label& operator=(const Label&);

		PsimagLite::String name_;
		SizeType site_;
		VectorOperatorType ops_;
	};

	typedef typename PsimagLite::Vector<Label*>::Type VectorLabelType;

	class IsValue {

	public:

		IsValue(PsimagLite::String value, SizeType site)
		    : value_(typename Label::PairStringSizeType(value, site)) {}

		bool operator()(Label const* label) const
		{
			return (*label == value_);
		}

	private:

		typename Label::PairStringSizeType value_;

	};

public:

	typedef OperatorType_ OperatorType;
	typedef Label LabelType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<typename LabelType::PairStringSizeType>::Type
	VectorPairStringSizeType;
	typedef typename OperatorType::value_type ComplexOrRealType;

	LabeledOperators(PsimagLite::String model = "") : model_(model), sites_(0)
	{}

	~LabeledOperators()
	{
		SizeType n = labels_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete labels_[i];
			labels_[i] = 0;
		}
	}

	void setModelName(PsimagLite::String model)
	{
		model_ = model;
	}

	void postCtor()
	{
		VectorSizeType h(sites_);
		SizeType n = trackables_.size();
		for (SizeType i = 0; i < n; ++i) {
			typename Label::PairStringSizeType nameAndSite = trackables_[i];
			const LabelType& ll = findLabel(nameAndSite.first, nameAndSite.second);
			SizeType dofs = ll.dofs();

			for (SizeType j = 0; j < dofs; ++j) {
				assert(nameAndSite.second < sites_);
				h[nameAndSite.second] = ll(j).data.rows();
			}
		}

		for (SizeType site = 0; site < sites_; ++site) {
			Label* labeli = new Label("i", site);
			labels_.push_back(labeli);
			pushIdentity(*labeli, h[site]);
		}
	}

	Label& createLabel(PsimagLite::String name,
	                   SizeType site)
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		                                                          labels_.end(),
		                                                          IsValue(name, site));

		if (x != labels_.end())
			err("Repeated label " + name + " site=" + ttos(site) + "\n");

		Label* label = new Label(name, site);
		labels_.push_back(label);
		if (site + 1 > sites_) sites_ = site + 1;
		return *label;
	}

	void makeTrackableOrderMatters(PsimagLite::String what, SizeType site)
	{
		typename Label::PairStringSizeType p(what, site);
		if (std::find(trackables_.begin(), trackables_.end(), p) != trackables_.end())
			err("makeTrackableOrderMatters: repeated entry " + what + "\n");
		trackables_.push_back(p);
	}

	const OperatorType& operator()(PsimagLite::String what,
	                               SizeType site,
	                               SizeType dof) const
	{
		return findLabel(what, site)(dof);
	}

	const LabelType& findLabel(PsimagLite::String what, SizeType site) const
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		                                                          labels_.end(),
		                                                          IsValue(what, site));
		if (x != labels_.end())
			return *(labels_[x - labels_.begin()]);

		PsimagLite::String str("LabeledOperators: model=" + model_);
		str += " label=" + what + " not found\n";
		throw PsimagLite::RuntimeError(str);
	}

	SizeType trackables() const  { return trackables_.size(); }

	const typename LabelType::PairStringSizeType& trackables(SizeType ind) const
	{
		assert(ind < trackables_.size());
		return trackables_[ind];
	}

	void instrospect() const
	{
		SizeType n = labels_.size();
		std::cout<<"There are "<<n<<" labels available for the "<<model_<<" model\n";
		for (SizeType i = 0; i < n; ++i)
			labels_[i]->instrospect();
	}

	void instrospect(PsimagLite::String what, SizeType site) const
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		                                                          labels_.end(),
		                                                          IsValue(what, site));
		if (x != labels_.end())
			return labels_[x - labels_.begin()]->instrospect();

		PsimagLite::String str("LabeledOperators: model=" + model_);
		str += " label=" + what + " or site=" + ttos(site) + " not found\n";
		throw PsimagLite::RuntimeError(str);
	}

private:

	void pushIdentity(LabelType& label, SizeType nrow)
	{
		typename OperatorType::StorageType tmp(nrow, nrow);
		tmp.makeDiagonal(nrow, 1.0);
		typename OperatorType::Su2RelatedType su2Related;
		label.push(OperatorType(tmp,
		                        1.0,
		                        typename OperatorType::PairType(0,0),
		                        1.0,
		                        su2Related));

	}

	LabeledOperators(const LabeledOperators&);

	LabeledOperators& operator=(const LabeledOperators&);

	PsimagLite::String model_;
	SizeType sites_;
	VectorLabelType labels_;
	VectorPairStringSizeType trackables_;
};
}
#endif // LABELEDOPERATORS_H
