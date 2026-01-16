#ifndef LABELEDOPERATOR_H
#define LABELEDOPERATOR_H
#include "PsimagLite.h"

namespace Dmft {

class LabeledOperator {

public:

	enum class Label {OPERATOR_NIL,
		              OPERATOR_C,
		              OPERATOR_SZ,
		              OPERATOR_CDAGGER,
		              OPERATOR_N,
		              OPERATOR_SPLUS,
		              OPERATOR_SMINUS};

	LabeledOperator(const Label& l) : what_(l), x_(0) {}

	/* PSIDOC OperatorIds
	  These are \verb!c! \verb!cdagger! \verb!sz!
	  \verb!nil! \verb!n! \verb!splus! and \verb!sminus!
	  They are case sensitive and model depedent.
	 */
	explicit LabeledOperator(const PsimagLite::String& s) : what_(toId(s)), x_(0) {}

	explicit LabeledOperator(SizeType x) : what_(Label::OPERATOR_NIL), x_(x) {}

	SizeType toUint() const
	{
		assert(what_ == Label::OPERATOR_NIL);
		return x_;
	}

	static Label toId(PsimagLite::String s)
	{
		if (s == "c") {
			return Label::OPERATOR_C;
		} else if (s == "cdagger") {
			return Label::OPERATOR_CDAGGER;
		} else if (s == "sz") {
			return Label::OPERATOR_SZ;
		} else if (s == "nil") {
			return Label::OPERATOR_NIL;
		} else if (s == "n") {
			return Label::OPERATOR_N;
		} else if (s == "splus") {
			return Label::OPERATOR_SPLUS;
		} else if (s == "sminus") {
			return Label::OPERATOR_SMINUS;
		} else {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) +  "\n";
			err(str + "operatorWithType: unsupported operator " + s + "\n");
		}

		return Label::OPERATOR_NIL;
	}

	Label id() const { return what_; }

	PsimagLite::String toString() const
	{
		PsimagLite::Vector<PsimagLite::String>::Type labels = {"cdagger",
		                                                       "c",
		                                                       "n",
		                                                       "sz",
		                                                       "splus",
		                                                       "sminus",
		                                                       "nil"};
		for (SizeType i = 0; i < labels.size(); ++i)
			if (toId(labels[i]) == what_) return labels[i];

		return "UNKNOWN";
	}

	SizeType numberOfTypes() const
	{
		return 4;
	}

	bool needsNewBasis() const
	{
		if (what_ == Label::OPERATOR_C || what_ == Label::OPERATOR_CDAGGER)
			return true;
		if (what_ == Label::OPERATOR_SPLUS || what_ == Label::OPERATOR_SMINUS)
			return true;
		return false;
	}

	PsimagLite::String unknownOperator() const
	{
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Unknown operator " + toString() + "\n";
		return str;
	}

	bool isFermionic() const
	{
		if (what_ == Label::OPERATOR_C) return true;
		if (what_ == Label::OPERATOR_CDAGGER) return true;
		return false;
	}

	Label transposeConjugate() const
	{
		Label l = what_;
		if (what_ == Label::OPERATOR_C)
			l = Label::OPERATOR_CDAGGER;
		if (what_ == Label::OPERATOR_CDAGGER)
			l = Label::OPERATOR_C;
		if (what_ == Label::OPERATOR_SPLUS)
			l = Label::OPERATOR_SMINUS;
		if (what_ == Label::OPERATOR_SMINUS)
			l = Label::OPERATOR_SPLUS;
		return l;
	}

private:

	Label what_;
	SizeType x_;
};
}
#endif // LABELEDOPERATOR_H
