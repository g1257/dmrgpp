#ifndef MULTIPOINTINSITU_H
#define MULTIPOINTINSITU_H
#include "Braket.h"
#include "Observer.h"
#include "HelperForMultiPointInSitu.h"
#include "ManyPointAction.h"
#include "Wft/WaveFunctionTransfFactory.h"
#include "Checkpoint.h"

namespace Dmrg {

template<typename VectorWithOffsetType, typename ModelType>
class MultiPointInSitu {

public:

	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelType::ParametersType ParametersType;
	typedef typename ParametersType::OptionsType OptionsType;
	typedef WaveFunctionTransfFactory<LeftRightSuperType, VectorWithOffsetType, OptionsType>
	WaveFunctionTransfType;
	typedef Checkpoint<ModelType, WaveFunctionTransfType> CheckpointType;
	typedef HelperForMultiPointInSitu<CheckpointType> HelperForMultiPointInSituType;
	typedef typename HelperForMultiPointInSituType::BogusInput BogusInputType;
	typedef typename HelperForMultiPointInSituType::MatrixType MatrixType;
	typedef Observer<HelperForMultiPointInSituType, ModelType> ObserverType;
	typedef Braket<ModelType> BraketType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;

	MultiPointInSitu(const ModelType& model,
	                 const CheckpointType& checkpoint,
	                 const WaveFunctionTransfType& wft,
	                 ProgramGlobals::DirectionEnum dir)
	    : model_(model),
	      bogusInput_(model.superGeometry().numberOfSites(), checkpoint, wft, dir),
	      observer_(bogusInput_, 0, 0, 0, model)
	{
		if (seen_.size() == 0)
			seen_.resize(bogusInput_.numberOfSites(), false);
	}

	void operator()(const BraketType& braket, SizeType centerOfOrtho)
	{
		if (bogusInput_.direction() == ProgramGlobals::DirectionEnum::INFINITE) return;

		if (!everySiteSeen()) {
			seen_[centerOfOrtho] = true;
			return;
		}

		// TODO FIXME: Use ObservableLibrary instead of Observer
		if (braket.points() != 2)
			err("MultiPointInSitu: only two point for now\n");

		constexpr bool needsPrinting = true;
		ManyPointAction action(false, "");
		const SizeType n = model_.superGeometry().numberOfSites();
		MatrixType storage(n, n);

		SizeType start = 0;
		SizeType end = 0;
		if (bogusInput_.direction() == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			start = 0;
			end = centerOfOrtho;
		} else {
			start = centerOfOrtho + 1;;
			end = n - 1;
		}

		for (SizeType site = start; site < end; ++site) {
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

	static bool everySiteSeen()
	{
		for (auto it = seen_.begin(); it != seen_.end(); ++it)
			if (!*it) return false;
		return true;
	}

	const ModelType& model_;
	BogusInputType bogusInput_;
	ObserverType observer_;
	static VectorBoolType seen_;
};

template<typename T1, typename T2>
typename MultiPointInSitu<T1, T2>::VectorBoolType MultiPointInSitu<T1, T2>::seen_;
}
#endif // MULTIPOINTINSITU_H
