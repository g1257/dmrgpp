/*
Copyright (c) 2012-2015, UT-Battelle, LLC
All rights reserved

[MPS++, Version 0.1]
[by K. Al-Hassanieh, Oak Ridge National Laboratory]
[by J. Rincon, Oak Ridge National Laboratory]
[by G.A., Oak Ridge National Laboratory]

See full open source LICENSE under file LICENSE
in the root directory of this distribution.

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup Dmrg */
/*@{*/
/** \file ModelSelector.h
 */

#ifndef MODEL_SELECTOR_H
#define MODEL_SELECTOR_H

#include "ProgramGlobals.h"
#include "Utils.h"
#include <stdexcept>
// start headers: // DO NOT REMOVE MARK
#include "../Models/ExtendedHubbard1Orb/ExtendedHubbard1Orb.h"
#include "../Models/FeAsBasedScExtended/FeAsBasedScExtended.h"
#include "../Models/FeAsModel/ModelFeBasedSc.h"
#include "../Models/FermionSpinless/FermionSpinless.h"
#include "../Models/GaugeSpin/GaugeSpin.h"
#include "../Models/Graphene/Graphene.h"
#include "../Models/Heisenberg/HeisenbergMix.h"
#include "../Models/Heisenberg/ModelHeisenberg.h"
#include "../Models/HeisenbergAncillaC/HeisenbergAncillaC.h"
#include "../Models/HolsteinSpinlessThin/HolsteinSpinlessThin.h"
#include "../Models/HolsteinThin/HolsteinThin.h"
#include "../Models/HubbardAncilla/HubbardAncilla.h"
#include "../Models/HubbardAncillaExtended/HubbardAncillaExtended.h"
#include "../Models/HubbardHolstein/HubbardHolstein.h"
#include "../Models/HubbardHolsteinSpinless/HubbardHolsteinSpinless.h"
#include "../Models/HubbardMultiBand/ModelHubbardMultiBand.h"
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "../Models/Immm/Immm.h"
#include "../Models/IsingMultiOrb/ModelIsingMultiOrb.h"
#include "../Models/Kitaev/Kitaev.h"
#include "../Models/Kondo/Kondo.h"
#include "../Models/LiouvillianHeisenberg/LiouvillianHeisenberg.hh"
#include "../Models/SpinOrbital/SpinOrbitalModel.h"
#include "../Models/Su3/Su3Model.h"
#include "../Models/SuperExtendedHubbard1Orb/SuperExtendedHubbard1Orb.h"
#include "../Models/SuperHubbardExtended/SuperHubbardExtended.h"
#include "../Models/TjAncillaC/TjAncillaC.h"
#include "../Models/TjAncillaC2/TjAncillaC2.h"
#include "../Models/TjAncillaG/TjAncillaG.h"
#include "../Models/TjAnisotropic/TjAnisotropic.h"
#include "../Models/TjMultiOrb/TjMultiOrb.h"
#include "../Models/UlsOsu/UlsOsu.h"
// end models DO NOT REMOVE MARK

namespace Dmrg {

template <typename ModelBaseType> class ModelSelector {

	using ModelHelperType    = typename ModelBaseType::ModelHelperType;
	using SolverParamsType   = typename ModelBaseType::SolverParamsType;
	using SuperGeometryType  = typename ModelBaseType::SuperGeometryType;
	using InputValidatorType = typename ModelBaseType::InputValidatorType;
	using SparseMatrixType   = typename ModelHelperType::SparseMatrixType;

	// start models here:  DO NOT REMOVE MARK
	using ModelHeisenbergType         = ModelHeisenberg<ModelBaseType>;
	using ModelIsingMultiOrbType      = ModelIsingMultiOrb<ModelBaseType>;
	using ModelHubbardType            = ModelHubbard<ModelBaseType>;
	using HeisenbergAncillaCType      = HeisenbergAncillaC<ModelBaseType>;
	using ModelHubbardExtType         = ExtendedHubbard1Orb<ModelBaseType>;
	using ModelHubbardExtSuperType    = ExtendedSuperHubbard1Orb<ModelBaseType>;
	using FeBasedScType               = ModelFeBasedSc<ModelBaseType>;
	using FeBasedScExtType            = FeAsBasedScExtended<ModelBaseType>;
	using ImmmType                    = Immm<ModelBaseType>;
	using TjMultiOrbType              = TjMultiOrb<ModelBaseType>;
	using TjAnisotropicType           = TjAnisotropic<ModelBaseType>;
	using UlsOsuType                  = UlsOsu<ModelBaseType>;
	using TjAncillaC2Type             = TjAncillaC2<ModelBaseType>;
	using TjAncillaCType              = TjAncillaC<ModelBaseType>;
	using TjAncillaGType              = TjAncillaG<ModelBaseType>;
	using SuperHubbardExtendedType    = SuperHubbardExtended<ModelBaseType>;
	using HubbardAncillaType          = HubbardAncilla<ModelBaseType>;
	using HubbardAncillaExtendedType  = HubbardAncillaExtended<ModelBaseType>;
	using FermionSpinlessType         = FermionSpinless<ModelBaseType>;
	using KitaevType                  = Kitaev<ModelBaseType>;
	using ModelHubbardMultiBandType   = ModelHubbardMultiBand<ModelBaseType>;
	using HubbardHolsteinType         = HubbardHolstein<ModelBaseType>;
	using HubbardHolsteinSpinlessType = HubbardHolsteinSpinless<ModelBaseType>;
	using KondoType                   = Kondo<ModelBaseType>;
	using GrapheneType                = Graphene<ModelBaseType>;
	using HolsteinThinType            = HolsteinThin<ModelBaseType>;
	using HolsteinSpinlessThinType    = HolsteinSpinlessThin<ModelBaseType>;
	using GaugeSpinType               = GaugeSpin<ModelBaseType>;
	using HeisenbergMixType           = HeisenbergMix<ModelBaseType>;
	using SpinOrbitalModelType        = SpinOrbitalModel<ModelBaseType>;
	using Su3ModelType                = Su3Model<ModelBaseType>;
	using LiouvillianHeisenbergType   = LiouvillianHeisenberg<ModelBaseType>;
	// end models  DO NOT REMOVE MARK

public:

	ModelSelector(const PsimagLite::String& name)
	    : name_(name)
	    , model_(0)
	{ }

	~ModelSelector()
	{
		if (model_)
			delete model_;
	}

	ModelBaseType& operator()(const SolverParamsType&  solverParams,
	                          InputValidatorType&      io,
	                          const SuperGeometryType& geometry)
	{
		if (model_)
			return *model_;

		PsimagLite::String hdf5fileIfAny = findHdf5FileIfAny(solverParams);

		// named models start  DO NOT REMOVE MARK
		if (name_ == "Heisenberg") {
			model_ = new ModelHeisenbergType(solverParams, io, geometry, "");
		} else if (name_ == "HeisenbergAnisotropic") {
			model_ = new ModelHeisenbergType(solverParams, io, geometry, "Anisotropic");
		} else if (name_ == "Aklt") {
			model_ = new ModelHeisenbergType(solverParams, io, geometry, "Aklt");
		} else if (name_ == "Heisenberg2") {
			model_ = new ModelHeisenbergType(solverParams, io, geometry, "2");
		} else if (name_ == "IsingMultiOrb") {
			model_ = new ModelIsingMultiOrbType(solverParams, io, geometry, "");
		} else if (name_ == "HeisenbergAncillaC") {
			model_ = new HeisenbergAncillaCType(solverParams, io, geometry);
		} else if (name_ == "HubbardOneBandExtended") {
			model_ = new ModelHubbardExtType(solverParams, io, geometry, "");
		} else if (name_.substr(0, 27) == "HubbardOneBandExtendedSuper") {
			PsimagLite::String tmp = getExtension("HubbardOneBandExtendedSuper");
			model_ = new ModelHubbardExtSuperType(solverParams, io, geometry, tmp);
		} else if (name_.substr(0, 11) == "FeAsBasedSc" && name_ != "FeAsBasedScExtended") {
			PsimagLite::String tmp = getExtension("FeAsBasedSc");
			model_                 = new FeBasedScType(solverParams, io, geometry, tmp);
		} else if (name_.substr(0, 19) == "FeAsBasedScExtended") {
			model_ = new FeBasedScExtType(solverParams, io, geometry);
		} else if (name_ == "Immm") {
			model_ = new ImmmType(solverParams, io, geometry);
		} else if (name_ == "TjMultiOrb") {
			model_ = new TjMultiOrbType(solverParams, io, geometry);
		} else if (name_ == "TjAnisotropic") {
			model_ = new TjAnisotropicType(solverParams, io, geometry);
		} else if (name_ == "UlsOsu") {
			model_ = new UlsOsuType(solverParams, io, geometry);
		} else if (name_ == "TjAncillaC2") {
			model_ = new TjAncillaC2Type(solverParams, io, geometry);
		} else if (name_ == "TjAncillaC") {
			model_ = new TjAncillaCType(solverParams, io, geometry);
		} else if (name_ == "TjAncillaG") {
			model_ = new TjAncillaGType(solverParams, io, geometry);
		} else if (name_ == "SuperHubbardExtended") {
			model_ = new SuperHubbardExtendedType(solverParams, io, geometry);
		} else if (name_ == "KaneMeleHubbard") {
			model_ = new ModelHubbardType(solverParams, io, geometry, "");
		} else if (name_ == "HubbardAncilla") {
			model_ = new HubbardAncillaType(solverParams, io, geometry);
		} else if (name_ == "HubbardAncillaExtended") {
			model_ = new HubbardAncillaExtendedType(solverParams, io, geometry);
		} else if (name_ == "FermionSpinless") {
			model_ = new FermionSpinlessType(solverParams, io, geometry, "");
		} else if (name_ == "FermionSpinlessWithDelta") {
			model_ = new FermionSpinlessType(solverParams, io, geometry, "WithDelta");
		} else if (name_.substr(0, 6) == "Kitaev") {
			PsimagLite::String tmp = getExtension("Kitaev");
			model_                 = new KitaevType(solverParams, io, geometry, tmp);
		} else if (name_ == "ModelHubbardMultiBand") {
			model_ = new ModelHubbardMultiBandType(solverParams, io, geometry);
		} else if (name_ == "HubbardHolstein") {
			model_ = new HubbardHolsteinType(
			    solverParams, io, geometry, "", hdf5fileIfAny);
		} else if (name_ == "HubbardHolsteinSSH") {
			model_ = new HubbardHolsteinType(
			    solverParams, io, geometry, "SSH", hdf5fileIfAny);
		} else if (name_ == "HubbardHolsteinLRH") {
			model_ = new HubbardHolsteinType(
			    solverParams, io, geometry, "LRH", hdf5fileIfAny);
		} else if (name_ == "HolsteinThin") {
			model_ = new HolsteinThinType(solverParams, io, geometry, "");
		} else if (name_ == "HubbardHolsteinSpinless") {
			model_ = new HubbardHolsteinSpinlessType(
			    solverParams, io, geometry, "", hdf5fileIfAny);
		} else if (name_ == "HubbardHolsteinSpinlessSSH") {
			model_ = new HubbardHolsteinSpinlessType(
			    solverParams, io, geometry, "SSH", hdf5fileIfAny);
		} else if (name_ == "HubbardHolsteinSpinlessLRH") {
			model_ = new HubbardHolsteinSpinlessType(
			    solverParams, io, geometry, "LRH", hdf5fileIfAny);
		} else if (name_ == "HolsteinSpinlessThin") {
			model_ = new HolsteinSpinlessThinType(solverParams, io, geometry, "");
		} else if (name_.substr(0, 5) == "Kondo") {
			PsimagLite::String tmp = getExtension("Kondo");
			model_                 = new KondoType(solverParams, io, geometry, tmp);
		} else if (name_.substr(0, 14) == "HubbardOneBand") {
			PsimagLite::String tmp = getExtension("HubbardOneBand");
			model_ = new ModelHubbardType(solverParams, io, geometry, tmp);
		} else if (name_ == "Graphene") {
			model_ = new GrapheneType(solverParams, io, geometry);
		} else if (name_ == "GaugeSpin") {
			model_ = new GaugeSpinType(solverParams, io, geometry);
		} else if (name_ == "HeisenbergMix") {
			model_ = new HeisenbergMixType(solverParams, io, geometry);
		} else if (name_.substr(0, 11) == "SpinOrbital") {
			PsimagLite::String tmp = getExtension("SpinOrbital");
			model_ = new SpinOrbitalModelType(solverParams, io, geometry, tmp);
		} else if (name_ == "Su3Model") {
			model_ = new Su3ModelType(solverParams, io, geometry);
		} else if (name_ == "LiouvillianHeisenberg") {
			model_ = new LiouvillianHeisenbergType(solverParams, io, geometry);
			// end models  DO NOT REMOVE MARK
		} else {
			PsimagLite::String s(__FILE__);
			s += " Unknown model " + name_ + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		model_->postCtor();

		return *model_;
	}

private:

	std::string getExtension(const std::string& str) const
	{
		SizeType           l     = str.length();
		SizeType           namel = name_.length();
		PsimagLite::String tmp   = (namel == l) ? "" : name_.substr(l, namel - l);
		return tmp;
	}

	PsimagLite::String findHdf5FileIfAny(const SolverParamsType& solverParams) const
	{
		// check first for observe
		bool isObserve = solverParams.options.isSet("observe");
		if (isObserve)
			return solverParams.filename;

		// then for restart
		bool isRestart = solverParams.options.isSet("restart");
		if (isRestart)
			return solverParams.checkpoint.filename();

		return "";
	}

	PsimagLite::String name_;
	ModelBaseType*     model_;

}; // ModelSelector

} // namespace Dmrg

/*@}*/
#endif // MODEL_SELECTOR_H
