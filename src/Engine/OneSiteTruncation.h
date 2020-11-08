#ifndef ONESITETRUNCATION_H
#define ONESITETRUNCATION_H

namespace Dmrg {

template<typename ModelType, typename VectorWithOffsetType>
class OneSiteTruncation {

public:

	OneSiteTruncation(const ModelType& model)
	    : model_(model) {}

	void update(const VectorWithOffsetType& psi)
	{}

private:

	const ModelType& model_;
};
}
#endif // ONESITETRUNCATION_H
