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

		const OperatorType& operator()(SizeType,
		                               SizeType dof) const
		{
			//if (site != 0)
			//	std::cerr<<"WARNING: LabeledOperators::Label does not support site ! = 0\n";
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

		SizeType rows() const
		{
			if (ops_.size() == 0)
				err("FATAL: LabeledOperators::Label::rows(): Internal Error\n");
			return ops_[0].data.rows();
		}

		Label(const Label&);

		Label& operator=(const Label&);

	public:

		void push(const OperatorType& op)
		{
			ops_.push_back(op);
		}

		void instrospect() const
		{
			std::cout<<"Label "<<name_<<" with "<<ops_.size()<<" dofs.\n";
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

	~LabeledOperators()
	{
		SizeType n = labels_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete labels_[i];
			labels_[i] = 0;
		}
	}

	void postCtor(SizeType h)
	{
		Label* labeli = new Label("i");
		labels_.push_back(labeli);
		pushIdentity(*labeli, h);
	}

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

	void instrospect() const
	{
		SizeType n = labels_.size();
		std::cout<<"There are "<<n<<" labels available for this model\n";
		for (SizeType i = 0; i < n; ++i)
			labels_[i]->instrospect();
	}

	void instrospect(PsimagLite::String what) const
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		                                                          labels_.end(),
		                                                          IsValue(what));
		if (x != labels_.end())
			return labels_[x - labels_.begin()]->instrospect();

		PsimagLite::String str("LabeledOperators: model=" + model_);
		str += " label=" + what + " not found\n";
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
	VectorLabelType labels_;
};
}
#endif // LABELEDOPERATORS_H
