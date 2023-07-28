#ifndef LABELEDOPERATORS_H
#define LABELEDOPERATORS_H
#include "TypeToString.h"
#include "Vector.h"

namespace Dmrg
{

template <typename OperatorType_>
class LabeledOperators
{

	class Label
	{

	public:

		typedef typename PsimagLite::Vector<OperatorType_>::Type VectorOperatorType;
		typedef typename PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
		typedef std::pair<PsimagLite::String, SizeType> PairStringSizeType;

		Label(PsimagLite::String name, SizeType kindOfSite)
		    : name_(name)
		    , kindOfSite_(kindOfSite)
		    , isTrackable_(false)
		{
		}

		const OperatorType_& operator()(SizeType dof) const
		{
			if (dof < ops_.size())
				return ops_[dof];

			PsimagLite::String msg("FATAL: LabeledOperators:");
			msg += " dof=" + ttos(dof) + " OUT OF RANGE, for label = " + name_ + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		bool operator==(const PsimagLite::String& name) const
		{
			return (name == name_);
		}

		SizeType rows() const
		{
			if (ops_.size() == 0)
				err("FATAL: LabeledOperators::Label::rows(): Internal Error\n");
			return ops_[0].getStorage().rows();
		}

		void push(const OperatorType_& op, PsimagLite::String desc = "")
		{
			const SizeType n = ops_.size();
			if (n > 0 && ops_[n - 1].getStorage().rows() != op.getStorage().rows())
				err("LabeledOperators::Label::push: FATAL\n");

			ops_.push_back(op);
			descriptions_.push_back(desc);
		}

		void introspect() const
		{
			std::cout << "Label " << name_ << " kindOfSite=" << kindOfSite_;
			std::cout << " with " << ops_.size() << " dofs.\n";
		}

		PsimagLite::String description(SizeType j) const
		{
			assert(j < descriptions_.size());
			return descriptions_[j];
		}

		SizeType dofs() const { return ops_.size(); }

		PsimagLite::String name() const { return name_; }

		SizeType kindOfSite() const { return kindOfSite_; }

		bool isTrackable() const { return isTrackable_; }

		void makeTrackable()
		{
			isTrackable_ = true;
		}

	private:

		Label(const Label&);

		Label& operator=(const Label&);

		PsimagLite::String name_;
		SizeType kindOfSite_;
		bool isTrackable_;
		VectorOperatorType ops_;
		VectorStringType descriptions_;
	};

	typedef typename PsimagLite::Vector<Label*>::Type VectorLabelType;
	class IsValue
	{

	public:

		IsValue(PsimagLite::String value)
		    : value_(value)
		{
		}

		bool operator()(Label const* label) const
		{
			return (*label == value_);
		}

	private:

		PsimagLite::String value_;
	};

public:

	typedef OperatorType_ OperatorType;
	typedef Label LabelType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename OperatorType::value_type ComplexOrRealType;
	typedef std::pair<SizeType, SizeType> PairSizeType;

	LabeledOperators(PsimagLite::String model = "")
	    : model_(model)
	{
	}

	~LabeledOperators()
	{
		SizeType n = labels_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete labels_[i];
			labels_[i] = 0;
		}
	}

	void clear()
	{
		const SizeType n = labels_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete labels_[i];
			labels_[i] = nullptr;
		}

		labels_.clear();
	}

	void setModelName(PsimagLite::String model)
	{
		model_ = model;

		SizeType n = labels_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete labels_[i];
			labels_[i] = 0;
		}

		labels_.clear();
	}

	Label& createLabel(PsimagLite::String name,
	    SizeType kindOfSite)
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		    labels_.end(),
		    IsValue(name));

		if (x != labels_.end())
			err("Repeated label " + name + "\n");

		Label* label = new Label(name, kindOfSite);
		labels_.push_back(label);
		return *label;
	}

	const OperatorType& operator()(PsimagLite::String what,
	    SizeType dof) const
	{
		return findLabel(what)(dof);
	}

	const LabelType& findLabel(PsimagLite::String what) const
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		    labels_.end(),
		    IsValue(what));
		if (x != labels_.end())
			return *(labels_[x - labels_.begin()]);

		PsimagLite::String str("LabeledOperators: model=" + model_);
		str += " label=" + what + " not found\n";
		throw PsimagLite::RuntimeError(str);
	}

	void introspect() const
	{
		SizeType n = labels_.size();
		std::cout << "There are " << n << " labels available for the " << model_ << " model\n";
		for (SizeType i = 0; i < n; ++i)
			labels_[i]->introspect();
	}

	void introspect(PsimagLite::String what) const
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		    labels_.end(),
		    what);
		if (x != labels_.end())
			return labels_[x - labels_.begin()]->introspect();

		PsimagLite::String str("LabeledOperators: model=" + model_);
		str += " label=" + what + " not found\n";
		throw PsimagLite::RuntimeError(str);
	}

	PsimagLite::String modelName() const { return model_; }

	SizeType size() const { return labels_.size(); }

	const LabelType& operator[](SizeType ind) const
	{
		assert(ind < labels_.size());
		assert(labels_[ind]);
		return *(labels_[ind]);
	}

	void makeTrackable(PsimagLite::String name)
	{
		const LabelType& label = findLabel(name);
		LabelType& labelNonConst = const_cast<LabelType&>(label);
		labelNonConst.makeTrackable();
	}

private:

	void pushIdentity(LabelType& label, SizeType nrow)
	{
		typename OperatorType::StorageType tmp(nrow, nrow);
		tmp.makeDiagonal(nrow, 1.0);
		typename OperatorType::Su2RelatedType su2Related;
		label.push(OperatorType(tmp,
		    1.0,
		    typename OperatorType::PairType(0, 0),
		    1.0,
		    su2Related));
	}

	LabeledOperators(const LabeledOperators&);

	LabeledOperators& operator=(const LabeledOperators&);

	PsimagLite::String model_;
	VectorLabelType labels_;
};
}
#endif // LABELEDOPERATORS_H
