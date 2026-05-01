#ifndef LANCZOS_MODEL_SELECTOR_H
#define LANCZOS_MODEL_SELECTOR_H

#include "../Models/FeBasedSc/BasisFeAsBasedSc.h"
#include "../Models/FeBasedSc/BasisFeAsSpinOrbit.h"
#include "../Models/FeBasedSc/FeBasedSc.h"
#include "../Models/FermionSpinless/FermionSpinless.hh"
#include "../Models/Heisenberg/Heisenberg.h"
#include "../Models/HubbardOneOrbital/HubbardOneOrbital.h"
#include "../Models/HubbardOneOrbital/HubbardOneOrbitalRashbaSOC.h"
#include "../Models/Immm/Immm.h"
#include "../Models/Kitaev/Kitaev.h"
#include "../Models/TjMultiOrb/TjMultiOrb.h"
#include "LanczosModelBase.hpp"
#include "Vector.h"

namespace LanczosPlusPlus {

template <typename ComplexOrRealType, typename GeometryType, typename InputType>
class ModelSelector {

	typedef TjMultiOrb<ComplexOrRealType, GeometryType, InputType>        TjMultiOrbType;
	typedef Immm<ComplexOrRealType, GeometryType, InputType>              ImmmType;
	typedef HubbardOneOrbital<ComplexOrRealType, GeometryType, InputType> HubbardOneOrbitalType;
	typedef HubbardOneOrbitalRashbaSOC<ComplexOrRealType, GeometryType, InputType>
	                                         HubbardOneBandRashbaSOCType;
	typedef BasisFeAsBasedSc<GeometryType>   BasisFeAsBasedScType;
	typedef BasisFeAsSpinOrbit<GeometryType> BasisFeAsSpinOrbitType;
	typedef FeBasedSc<ComplexOrRealType, BasisFeAsBasedScType, InputType> FeBasedScType;
	typedef FeBasedSc<ComplexOrRealType, BasisFeAsSpinOrbitType, InputType>
	                                                                    FeBasedScSpinOrbitType;
	typedef Heisenberg<ComplexOrRealType, GeometryType, InputType>      HeisenbergType;
	typedef Kitaev<ComplexOrRealType, GeometryType, InputType>          KitaevType;
	typedef FermionSpinless<ComplexOrRealType, GeometryType, InputType> FermionSpinlessType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type          RealType;

public:

	typedef LanczosModelBase<ComplexOrRealType, GeometryType, InputType> ModelBaseType;

	/** @class hide_modelselector
	 - TargetElectronsUp=integer
	 - TargetElectronsDown=integer
	 - Model=string One of TjMultiOrb Heisenberg Immm
	                                           HubbardOneBand HubbardOneBandExtended
	 SuperHubbardExtended FeAsBasedSc FeAsBasedScExtended
	*/
	ModelSelector(InputType& io, const GeometryType& geometry)
	    : modelPtr_(0)
	{
		PsimagLite::String model("");
		io.readline(model, "Model=");

		SizeType nup         = 0;
		SizeType ndown       = 0;
		SizeType szPlusConst = 0;

		if (model != "Kitaev" && model != "HubbardOneBandRashbaSOC"
		    && model != "FermionSpinless") {
			try {
				io.readline(nup, "TargetElectronsUp=");
				io.readline(ndown, "TargetElectronsDown=");
			} catch (std::exception&) {
				io.readline(szPlusConst, "TargetSzPlusConst=");
			}
		}

		if (model == "HubbardOneBandRashbaSOC" || model == "FermionSpinless")
			io.readline(szPlusConst, "TargetElectronsTotal=");

		PsimagLite::Matrix<ComplexOrRealType> spinOrbit;
		try {
			io.read(spinOrbit, "SpinOrbit");
		} catch (std::exception&) { }

		if (model == "TjMultiOrb") {
			modelPtr_ = new TjMultiOrbType(nup, ndown, io, geometry);
		} else if (model == "Immm") {
			modelPtr_ = new ImmmType(nup, ndown, io, geometry);
		} else if (model == "HubbardOneBand" || model == "HubbardOneBandExtended"
		           || model == "SuperHubbardExtended" || model == "KaneMeleHubbard") {
			modelPtr_ = new HubbardOneOrbitalType(nup, ndown, io, geometry);
		} else if (model == "FeAsBasedSc" || model == "FeAsBasedScExtended") {
			if (spinOrbit.n_row() != 4)
				modelPtr_ = new FeBasedScType(nup, ndown, io, geometry);
			else
				modelPtr_ = new FeBasedScSpinOrbitType(nup, ndown, io, geometry);
		} else if (model == "Heisenberg") {
			modelPtr_ = new HeisenbergType(szPlusConst, io, geometry);
		} else if (model == "Kitaev") {
			modelPtr_ = new KitaevType(io, geometry);
		} else if (model == "HubbardOneBandRashbaSOC") {
			modelPtr_ = new HubbardOneBandRashbaSOCType(szPlusConst, io, geometry);
		} else if (model == "FermionSpinless") {
			modelPtr_ = new FermionSpinlessType(szPlusConst, io, geometry);
		} else {
			PsimagLite::String str("No known model " + model + "\n");
			throw PsimagLite::RuntimeError(str);
		}
	}

	~ModelSelector()
	{
		delete modelPtr_;
		modelPtr_ = 0;
	}

	const ModelBaseType& operator()() const { return *modelPtr_; }

private:

	ModelBaseType* modelPtr_;
};

} // namespace LanczosPlusPlus

#endif
