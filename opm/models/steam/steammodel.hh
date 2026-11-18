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
 * \brief Pure-water two-phase (liquid/steam) non-isothermal model with
 *        pressure and specific enthalpy as primary variables.
 *
 * Two equations per cell: mass of water and energy. The temperature and
 * the steam saturation are secondary variables found by WaterPHFlash from
 * (p, h), and boiling and condensation need no primary variable switching.
 */
#ifndef OPM_STEAM_MODEL_HH
#define OPM_STEAM_MODEL_HH

#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>

#include <opm/models/common/multiphasebasemodel.hh>

#include <opm/models/flash/flashextensivequantities.hh>

#include <opm/models/steam/steamindices.hh>
#include <opm/models/steam/steamintensivequantities.hh>
#include <opm/models/steam/steamlocalresidual.hh>
#include <opm/models/steam/steamnewtonmethod.hh>
#include <opm/models/steam/steamprimaryvariables.hh>
#include <opm/models/steam/steamratevector.hh>

#include <string>

namespace Opm {

template <class TypeTag>
class SteamModel;

}

namespace Opm::Properties {

namespace TTag {
struct SteamModel { using InheritsFrom = std::tuple<MultiPhaseBaseModel>; };
} // namespace TTag

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::SteamModel>
{ using type = SteamLocalResidual<TypeTag>; };

template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::SteamModel>
{ using type = SteamNewtonMethod<TypeTag>; };

template<class TypeTag>
struct Model<TypeTag, TTag::SteamModel>
{ using type = SteamModel<TypeTag>; };

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::SteamModel>
{ using type = SteamPrimaryVariables<TypeTag>; };

template<class TypeTag>
struct RateVector<TypeTag, TTag::SteamModel>
{ using type = SteamRateVector<TypeTag>; };

template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::SteamModel>
{ using type = SteamBoundaryRateVector<TypeTag>; };

template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::SteamModel>
{ using type = SteamIntensiveQuantities<TypeTag>; };

// the advective flux is the multi-phase Darcy flux of the flux module; the
// energy conduction is done in the local residual
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::SteamModel>
{ using type = FlashExtensiveQuantities<TypeTag>; };

template<class TypeTag>
struct Indices<TypeTag, TTag::SteamModel>
{ using type = SteamIndices<TypeTag, /*PVOffset=*/0>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::SteamModel>
{ static constexpr bool value = false; };

// the generic energy module is not used: the model has its own energy equation
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::SteamModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableWater<TypeTag, TTag::SteamModel>
{ static constexpr int value = true; };

} // namespace Opm::Properties

namespace Opm {

template <class TypeTag>
class SteamModel : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

public:
    explicit SteamModel(Simulator& simulator) : ParentType(simulator) {}

    static void registerParameters()
    {
        ParentType::registerParameters();
        Parameters::SetDefault<Parameters::EnableIntensiveQuantityCache>(true);
    }

    std::string primaryVarName(unsigned pvIdx) const
    {
        if (pvIdx == static_cast<unsigned>(Indices::pressure0Idx)) {
            return "pressure_liquid";
        }
        return "enthalpy";
    }

    std::string eqName(unsigned eqIdx) const
    {
        if (eqIdx == static_cast<unsigned>(Indices::conti0EqIdx)) {
            return "continuity^H2O";
        }
        return "energy";
    }
};

} // namespace Opm

#endif
