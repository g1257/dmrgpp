#ifndef MULTIPOINTINSITU_H
#define MULTIPOINTINSITU_H
#include "Braket.h"
#include "Observer.h"
#include "HelperForMultiPointInSitu.h"
#include "ManyPointAction.h"

namespace Dmrg {

template<typename VectorWithOffsetType, typename ModelType>
class MultiPointInSitu {

public:

	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
	typedef HelperForMultiPointInSitu<VectorWithOffsetType, LeftRightSuperType>
	HelperForMultiPointInSituType;
	typedef typename HelperForMultiPointInSituType::BogusInput BogusInputType;
	typedef typename HelperForMultiPointInSituType::MatrixType MatrixType;
	typedef Observer<HelperForMultiPointInSituType, ModelType> ObserverType;
	typedef Braket<ModelType> BraketType;

	MultiPointInSitu(const ModelType& model)
	    : model_(model),
	      bogusInput_(model.superGeometry().numberOfSites()),
	      observer_(bogusInput_, 0, 0, 0, model)
	{}


	void operator()(const BraketType& braket, SizeType centerOfOrtho) const
	{
		// TODO FIXME: Use ObservableLibrary instead of Observer
		if (braket.points() != 2)
			err("MultiPointInSitu: only two point for now\n");

		constexpr bool needsPrinting = true;
		ManyPointAction action(false, "");
		const SizeType n = model_.superGeometry().numberOfSites();
		MatrixType storage(n, n);

		for (SizeType site = 0; site < n; ++site) {
			if (site == centerOfOrtho) continue;
			BraketType braket2 = buildBraketWithSites(braket, site, centerOfOrtho);
			observer_.twoPoint(storage, braket2, needsPrinting, action);
		}
	}

private:

	BraketType buildBraketWithSites(const BraketType& braket, SizeType site0, SizeType site1) const
	{
		if (braket.points() != 2)
			err("MultiPointInSitu::buildBraketWithSites: expected two sites\n");

		if (site0 == site1)
			err("MultiPointInSitu::buildBraketWithSites: expected different sites\n");

		const PsimagLite::String op0 = braket.opName(0) + "[" + ttos(site0) + "]";
		const PsimagLite::String op1 = braket.opName(1) + "[" + ttos(site1) + "]";
		const PsimagLite::String bra = "<" + braket.bra().toString() + "|";
		const PsimagLite::String ket = "|" + braket.ket().toString() + ">";
		return BraketType(model_, bra + op0 + ";" + op1 + ket);
	}

	const ModelType& model_;
	BogusInputType bogusInput_;
	ObserverType observer_;
};
}
#endif // MULTIPOINTINSITU_H
