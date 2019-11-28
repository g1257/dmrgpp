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

		Label(PsimagLite::String name, SizeType kindOfSite)
		    : name_(name), kindOfSite_(kindOfSite) {}

		const OperatorType_& operator()(SizeType dof) const
		{
			if (dof < ops_.size())
				return ops_[dof];

			PsimagLite::String msg("FATAL: LabeledOperators:");
			msg += " dof=" + ttos(dof) + " OUT OF RANGE, for label = " + name_ + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

//		bool operator==(PairStringSizeType nameAndSite) const
//		{
//			return (nameAndSite.first == name_ && nameAndSite.second == site_);
//		}

		SizeType rows() const
		{
			if (ops_.size() == 0)
				err("FATAL: LabeledOperators::Label::rows(): Internal Error\n");
			return ops_[0].data.rows();
		}

		void push(const OperatorType_& op)
		{
			const SizeType n = ops_.size();
			if (n > 0 && ops_[n - 1].getCrs().rows() != op.getCrs().rows())
				err("LabeledOperators::Label::push: FATAL\n");

			ops_.push_back(op);
		}

		void instrospect() const
		{
			std::cout<<"Label "<<name_<<" kindOfSite="<<kindOfSite_;
			std::cout<<" with "<<ops_.size()<<" dofs.\n";
		}

		SizeType dofs() const { return ops_.size(); }

		SizeType kindOfSite() const { return kindOfSite_; }

	private:

		Label(const Label&);

		Label& operator=(const Label&);

		PsimagLite::String name_;
		SizeType kindOfSite_;
		VectorOperatorType ops_;
	};

	typedef typename PsimagLite::Vector<Label*>::Type VectorLabelType;

public:

	typedef OperatorType_ OperatorType;
	typedef Label LabelType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename OperatorType::value_type ComplexOrRealType;
	typedef std::pair<SizeType, SizeType> PairSizeType;

	LabeledOperators(PsimagLite::String model = "") : model_(model)
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

//	void postCtor()
//	{
//		assert(0 < labels_.size());
//		SizeType h = labels_[0].data.rows();
//		Label* labeli = new Label("i", site);
//		labels_.push_back(labeli);
//		pushIdentity(*labeli, h);
//	}

	Label& createLabel(PsimagLite::String name,
	                   SizeType site)
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		                                                          labels_.end(),
		                                                          name);

		if (x != labels_.end())
			err("Repeated label " + name + "\n");

		Label* label = new Label(name, site);
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
		                                                          what);
		if (x != labels_.end())
			return *(labels_[x - labels_.begin()]);

		PsimagLite::String str("LabeledOperators: model=" + model_);
		str += " label=" + what + " not found\n";
		throw PsimagLite::RuntimeError(str);
	}

	void instrospect() const
	{
		SizeType n = labels_.size();
		std::cout<<"There are "<<n<<" labels available for the "<<model_<<" model\n";
		for (SizeType i = 0; i < n; ++i)
			labels_[i]->instrospect();
	}

	void instrospect(PsimagLite::String what) const
	{
		typename VectorLabelType::const_iterator x = std::find_if(labels_.begin(),
		                                                          labels_.end(),
		                                                          what);
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
