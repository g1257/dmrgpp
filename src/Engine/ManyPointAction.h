#ifndef MANYPOINTACTION_H
#define MANYPOINTACTION_H
#include "PredicateAwesome.h"

namespace Dmrg {

class ManyPointAction {

public:

	typedef PsimagLite::PredicateAwesome<> PredicateAwesomeType;

	ManyPointAction(bool hasNonTrivialAction, PsimagLite::String actionString)
	    : nonTrivial_(hasNonTrivialAction)
	    , actionString_(actionString)
	{ }

	bool operator()(SizeType s0, SizeType s1, SizeType s2, SizeType s3) const
	{
		if (!nonTrivial_)
			return true;

		PredicateAwesomeType pred(actionString_, "~");
		return pred.isTrue("%0", s0, "%1", s1, "%2", s2, "%3", s3);
	}

	bool operator()(SizeType s0, SizeType s1) const
	{
		if (!nonTrivial_)
			return true;

		PredicateAwesomeType pred(actionString_, "~");
		return pred.isTrue("%0", s0, "%1", s1);
	}

	bool nonTrivial() const { return nonTrivial_; }

private:

	bool nonTrivial_;
	PsimagLite::String actionString_;
};

}
#endif // MANYPOINTACTION_H
