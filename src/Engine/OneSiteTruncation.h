#ifndef ONESITETRUNCATION_H
#define ONESITETRUNCATION_H

namespace Dmrg {

template<typename ModelType, typename VectorWithOffsetType>
class OneSiteTruncation {

public:

	typedef typename ModelType::MatrixType MatrixType;

	OneSiteTruncation(const ModelType& model)
	    : model_(model) {}

	void update(const VectorWithOffsetType& psi)
	{
		if (!model_.isOneSiteTruncationActive()) return;

		// compute U ...
		MatrixType U;
		computeU(U);

		// ... and send it to model
		model_.oneSiteTruncationUpdate(U);
	}

private:

	void computeU(MatrixType&)
	{
		std::cerr<<"OneSiteTruncation::computeU() unimplemented, sorry, todo, fixme\n";
	}

	const ModelType& model_;
};
}
#endif // ONESITETRUNCATION_H
