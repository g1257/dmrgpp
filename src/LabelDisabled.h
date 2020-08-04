#ifndef LABEL_DISABLED_H
#define LABEL_DISABLED_H
#include "Vector.h"

namespace PsimagLite {

class LabelDisabled {

	typedef Vector<String>::Type VectorStringType;

public:

	void disable(String label)
	{
		data_.push_back(label);
	}

	bool operator()(String label) const
	{
		return (find(data_.begin(),data_.end(),label) != data_.end());
	}

private:

	VectorStringType data_;
};

}
#endif // LABEL_DISABLED_H
