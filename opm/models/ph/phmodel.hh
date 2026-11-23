// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 NORCE Research AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \brief Generic compositional pressure-enthalpy model: N components in up to three phases, non-isothermal,
 *        with pressure, specific enthalpy and N-1 switched variables as primary variables.
 *
 * N + 1 equations per cell: mass of the components and energy. The temperature, the saturations and the
 * compositions are secondary variables found by the p-h flash of the property PhFlash, which is selected by the
 * property PhFlashKind (see phflash.hh): pure water (flow_steam), brine and CO2 (flow_co2store) and later the
 * cubic EOS flash. A binary sets the kind, for instance
 *   template<class TypeTag> struct PhFlashKind<TypeTag, TTag::MyProblem> { static constexpr PhFlashType value = PhFlashType::Steam; };
 */
#ifndef OPM_PH_MODEL_HH
#define OPM_PH_MODEL_HH

#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>

#include <opm/models/common/multiphasebasemodel.hh>

#include <opm/models/flash/flashextensivequantities.hh>

#include <opm/models/ph/co2storeflash.hh>
#include <opm/models/ph/phflash.hh>
#include <opm/models/ph/phindices.hh>
#include <opm/models/ph/phintensivequantities.hh>
#include <opm/models/ph/phlocalresidual.hh>
#include <opm/models/ph/phnewtonmethod.hh>
#include <opm/models/ph/phprimaryvariables.hh>
#include <opm/models/ph/phratevector.hh>
#include <opm/models/ph/steamflash.hh>

#include <string>

namespace Opm {

template <class TypeTag>
class PhModel;

}

namespace Opm::Properties {

namespace TTag {
struct PhModel { using InheritsFrom = std::tuple<MultiPhaseBaseModel>; };
} // namespace TTag

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PhModel>
{ using type = PhLocalResidual<TypeTag>; };

template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::PhModel>
{ using type = PhNewtonMethod<TypeTag>; };

template<class TypeTag>
struct Model<TypeTag, TTag::PhModel>
{ using type = PhModel<TypeTag>; };

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::PhModel>
{ using type = PhPrimaryVariables<TypeTag>; };

// the rate vectors are generic: one entry per equation
template<class TypeTag>
struct RateVector<TypeTag, TTag::PhModel>
{ using type = PhRateVector<TypeTag>; };

template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::PhModel>
{ using type = PhBoundaryRateVector<TypeTag>; };

// the flash class follows from the kind of flash, which a binary sets (PhFlashKind has no default)
template<class TypeTag>
struct PhFlash<TypeTag, TTag::PhModel>
{ using type = typename PhFlashSelector<getPropValue<TypeTag, Properties::PhFlashKind>(),
                                        GetPropType<TypeTag, Properties::Scalar>>::type; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PhModel>
{ using type = typename GetPropType<TypeTag, Properties::PhFlash>::FluidSystem; };

template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::PhModel>
{ using type = PhIntensiveQuantities<TypeTag>; };

// the advective flux is the multi-phase Darcy flux of the flux module; the
// energy conduction is done in the local residual
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::PhModel>
{ using type = FlashExtensiveQuantities<TypeTag>; };

template<class TypeTag>
struct Indices<TypeTag, TTag::PhModel>
{ using type = PhIndices<TypeTag, GetPropType<TypeTag, Properties::PhFlash>, /*PVOffset=*/0>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::PhModel>
{ static constexpr bool value = false; };

// the generic energy module is not used: the model has its own energy equation
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::PhModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableWater<TypeTag, TTag::PhModel>
{ static constexpr int value = true; };

} // namespace Opm::Properties

namespace Opm {

template <class TypeTag>
class PhModel : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Flash = GetPropType<TypeTag, Properties::PhFlash>;

public:
    explicit PhModel(Simulator& simulator) : ParentType(simulator) {}

    static void registerParameters()
    {
        ParentType::registerParameters();
        Parameters::SetDefault<Parameters::EnableIntensiveQuantityCache>(true);
    }

    std::string primaryVarName(unsigned pvIdx) const
    {
        if (pvIdx == static_cast<unsigned>(Indices::pressure0Idx)) {
            return Flash::numComponents == 1 ? "pressure_liquid" : "pressure_aqueous";
        }
        if (pvIdx == static_cast<unsigned>(Indices::enthalpy0Idx)) {
            return "enthalpy";
        }
        return "switch";
    }

    std::string eqName(unsigned eqIdx) const
    {
        if (eqIdx == static_cast<unsigned>(Indices::contiEnergyEqIdx)) {
            return "energy";
        }
        if (Flash::numComponents == 1) {
            return "continuity^H2O";
        }
        return "continuity^" + std::to_string(eqIdx - static_cast<unsigned>(Indices::conti0EqIdx));
    }
};

} // namespace Opm

#endif
