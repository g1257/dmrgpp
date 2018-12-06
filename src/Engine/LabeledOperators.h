#ifndef LABELEDOPERATORS_H
#define LABELEDOPERATORS_H
#include "Vector.h"
#include "TypeToString.h"

namespace Dmrg {

template<typename OperatorType>
class LabeledOperators {

	class Label {

		typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;

		Label(PsimagLite::String name) : name_(name) {}

		const OperatorType& operator()(SizeType site,
		                               SizeType dof) const
		{
			if (site != 0)
				std::cerr<<"WARNING: LabeledOperators::Label does not support site ! = 0\n";
			if (dof < ops_.size())
				return ops_[dof];

			PsimagLite::String msg("FATAL: LabeledOperators:");
			msg += " dof=" + ttos(dof) + " OUT OF RANGE, for label = " + name_ + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		bool operator==(PsimagLite::String name) const
		{
			return (name == name_);
		}

		Label(const Label&);

		Label& operator=(const Label&);

	public:

		void push(const OperatorType& op)
		{
			ops_.push_back(op);
		}

		friend class LabeledOperators;

	private:

		PsimagLite::String name_;
		VectorOperatorType ops_;
	};

	typedef typename PsimagLite::Vector<Label*>::Type VectorLabelType;

	class IsValue {

	public:

		IsValue(PsimagLite::String value) : value_(value) {}

		bool operator()(Label const* label) const
		{
			return (*label == value_);
		}

	private:

		PsimagLite::String value_;

	};
public:

	typedef Label LabelType;

	LabeledOperators(PsimagLite::String model) : model_(model)
	{}

	Label& createLabel(PsimagLite::String name)
	{
		Label* label = new Label(name);
		labels_.push_back(label);
		return *label;
	}

	const OperatorType& operator()(PsimagLite::String what,
	                               SizeType site,
	                               SizeType dof) const
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		                                                          labels_.end(),
		                                                          IsValue(what));
		if (x != labels_.end())
			return labels_[x - labels_.begin()]->operator()(site, dof);

		PsimagLite::String str("LabeledOperators: model=" + model_);
		str += " label=" + what + " not found\n";
		throw PsimagLite::RuntimeError(str);
	}

private:

	PsimagLite::String model_;
	VectorLabelType labels_;
};
}
#endif // LABELEDOPERATORS_H
