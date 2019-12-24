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

#define ALL_MODELS 1
#include <stdexcept>
#include "ProgramGlobals.h"
#include "Utils.h"
#include "../Models/Heisenberg/ModelHeisenberg.h"
#if ALL_MODELS
#include "../Models/HubbardOneBand/ModelHubbard.h"
#include "../Models/HeisenbergAncillaC/HeisenbergAncillaC.h"
#include "../Models/ExtendedHubbard1Orb/ExtendedHubbard1Orb.h"
#include "../Models/SuperExtendedHubbard1Orb/SuperExtendedHubbard1Orb.h"
#include "../Models/FeAsModel/ModelFeBasedSc.h"
#include "../Models/FeAsBasedScExtended/FeAsBasedScExtended.h"
#include "../Models/Immm/Immm.h"
#include "../Models/TjMultiOrb/TjMultiOrb.h"
#include "../Models/TjAnisotropic/TjAnisotropic.h"
#include "../Models/TjAncillaC2/TjAncillaC2.h"
#include "../Models/TjAncillaC/TjAncillaC.h"
#include "../Models/TjAncillaG/TjAncillaG.h"
#include "../Models/SuperHubbardExtended/SuperHubbardExtended.h"
#include "../Models/HubbardAncilla/HubbardAncilla.h"
#include "../Models/HubbardAncillaExtended/HubbardAncillaExtended.h"
#include "../Models/FermionSpinless/FermionSpinless.h"
#include "../Models/Kitaev/Kitaev.h"
#include "../Models/HubbardMultiBand/ModelHubbardMultiBand.h"
#include "../Models/HubbardHolstein/HubbardHolstein.h"
#include "../Models/HolsteinThin/HolsteinThin.h"
#include "../Models/Kondo/Kondo.h"
#include "../Models/UlsOsu/UlsOsu.h"
#include "../Models/Graphene/Graphene.h"

#endif

namespace Dmrg {

template<typename ModelBaseType>
class ModelSelector {

	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename ModelBaseType::SuperGeometryType SuperGeometryType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;

	// start models here:
	typedef ModelHeisenberg<ModelBaseType> ModelHeisenbergType;
#if ALL_MODELS
	typedef ModelHubbard<ModelBaseType> ModelHubbardType;
	typedef HeisenbergAncillaC<ModelBaseType> HeisenbergAncillaCType;
	typedef ExtendedHubbard1Orb<ModelBaseType> ModelHubbardExtType;
	typedef ExtendedSuperHubbard1Orb<ModelBaseType> ModelHubbardExtSuperType;
	typedef ModelFeBasedSc<ModelBaseType> FeBasedScType;
	typedef FeAsBasedScExtended<ModelBaseType> FeBasedScExtType;
	typedef Immm<ModelBaseType> ImmmType;
	typedef TjMultiOrb<ModelBaseType> TjMultiOrbType;
	typedef TjAnisotropic<ModelBaseType> TjAnisotropicType;
    typedef UlsOsu<ModelBaseType> UlsOsuType;
    typedef TjAncillaC2<ModelBaseType> TjAncillaC2Type;
	typedef TjAncillaC<ModelBaseType> TjAncillaCType;
	typedef TjAncillaG<ModelBaseType> TjAncillaGType;
	typedef SuperHubbardExtended<ModelBaseType> SuperHubbardExtendedType;
	typedef HubbardAncilla<ModelBaseType> HubbardAncillaType;
	typedef HubbardAncillaExtended<ModelBaseType> HubbardAncillaExtendedType;
	typedef FermionSpinless<ModelBaseType> FermionSpinlessType;
	typedef Kitaev<ModelBaseType> KitaevType;
	typedef ModelHubbardMultiBand<ModelBaseType> ModelHubbardMultiBandType;
	typedef HubbardHolstein<ModelBaseType> HubbardHolsteinType;
	typedef Kondo<ModelBaseType> KondoType;
	typedef Graphene<ModelBaseType> GrapheneType;
	typedef HolsteinThin<ModelBaseType> HolsteinThinType;

#endif
	// end models

public:

	ModelSelector(const PsimagLite::String& name)
	    : name_(name),model_(0)
	{}

	~ModelSelector()
	{
		if (model_) delete model_;
	}

	const ModelBaseType& operator()(const SolverParamsType& solverParams,
	                                InputValidatorType& io,
	                                const SuperGeometryType& geometry)
	{
		if (model_) return *model_;

		if (name_ == "Heisenberg") {
			model_ = new ModelHeisenbergType(solverParams,io,geometry,"");
		} else if (name_ == "HeisenbergAnisotropic") {
			model_ = new ModelHeisenbergType(solverParams,io,geometry,"Anisotropic");
		}
#if ALL_MODELS
		else if (name_ == "HubbardOneBand") {
			model_ = new ModelHubbardType(solverParams, io, geometry, "");
		} else if (name_ == "HeisenbergAncillaC") {
			model_ = new HeisenbergAncillaCType(solverParams,io,geometry);
		} else if (name_ == "HubbardOneBandExtended") {
			model_ = new ModelHubbardExtType(solverParams,io,geometry);
		} else if (name_ == "HubbardOneBandExtendedSuper") {
			model_ = new ModelHubbardExtSuperType(solverParams,io,geometry);
		} else if (name_ == "FeAsBasedSc") {
			model_ = new FeBasedScType(solverParams,io,geometry);
		} else if (name_ == "FeAsBasedScExtended") {
			model_ = new FeBasedScExtType(solverParams,io,geometry);
		} else if (name_ == "Immm") {
			model_ = new ImmmType(solverParams,io,geometry);
		} else if (name_ == "TjMultiOrb") {
			model_ = new TjMultiOrbType(solverParams,io,geometry);
		} else if (name_ == "TjAnisotropic") {
			model_ = new TjAnisotropicType(solverParams,io,geometry);
        } else if (name_ == "UlsOsu") {
            model_ = new UlsOsuType(solverParams, io, geometry);
		} else if (name_ == "TjAncillaC2") {
			model_ = new TjAncillaC2Type(solverParams,io,geometry);
		} else if (name_ == "TjAncillaC") {
			model_ = new TjAncillaCType(solverParams,io,geometry);
		} else if (name_ == "TjAncillaG") {
			model_ = new TjAncillaGType(solverParams,io,geometry);
		} else if (name_ == "SuperHubbardExtended") {
			model_ = new SuperHubbardExtendedType(solverParams,io,geometry);
		} else if (name_ == "KaneMeleHubbard") {
			model_ = new ModelHubbardType(solverParams,io,geometry,"");
		} else if (name_ == "HubbardAncilla") {
			model_ = new HubbardAncillaType(solverParams,io,geometry);
		} else if (name_ == "HubbardAncillaExtended") {
			model_ = new HubbardAncillaExtendedType(solverParams,io,geometry);
		} else if (name_ == "FermionSpinless") {
			model_ = new FermionSpinlessType(solverParams,io,geometry);
		} else if (name_ == "Kitaev") {
			model_ = new KitaevType(solverParams,io,geometry,"");
		} else if (name_ == "KitaevExtended") {
			model_ = new KitaevType(solverParams,io,geometry,"Extended");
		} else if (name_ == "KitaevWithGammas") {
			model_ = new KitaevType(solverParams,io,geometry,"WithGammas");
		} else if (name_ == "ModelHubbardMultiBand") {
			model_ = new ModelHubbardMultiBandType(solverParams,io,geometry);
		} else if (name_ == "HubbardHolstein") {
			model_ = new HubbardHolsteinType(solverParams, io, geometry, "");
		} else if (name_ == "HubbardHolsteinSSH") {
			model_ = new HubbardHolsteinType(solverParams, io, geometry, "SSH");
		} else if (name_ == "HolsteinThin") {
			model_ = new HolsteinThinType(solverParams, io, geometry, "");
		} else if (name_ == "Kondo") {
			model_ = new KondoType(solverParams, io, geometry, "");
		} else if (name_ == "KondoEx") {
			model_ = new KondoType(solverParams, io, geometry, "Ex");
		} else if (name_ == "HubbardOneBandRashbaSOC") {
			model_ = new ModelHubbardType(solverParams, io, geometry, "RashbaSOC");
		} else if (name_ == "Graphene") {
			model_ = new GrapheneType(solverParams, io, geometry);
		}
#endif
		else {
			PsimagLite::String s(__FILE__);
			s += " Unknown model " + name_ + "\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		model_->postCtor();

		return *model_;
	}

private:

	PsimagLite::String name_;
	ModelBaseType* model_;

}; // ModelSelector

} // namespace Dmrg

/*@}*/
#endif // MODEL_SELECTOR_H

