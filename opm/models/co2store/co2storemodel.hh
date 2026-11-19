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
 * \brief CO2 storage model: brine and CO2 in three phases (brine, liquid CO2, gas CO2), non-isothermal,
 *        with pressure, specific enthalpy and a switched third variable as primary variables.
 *
 * Three equations per cell: mass of brine, mass of CO2 and energy. The temperature, the saturations and
 * the compositions are secondary variables found by a pressure-enthalpy provider
 * (Opm::Co2BrinePressureEnthalpyProvider), see Opm::Co2StorePrimaryVariables for the switching.
 */
#ifndef OPM_CO2STORE_MODEL_HH
#define OPM_CO2STORE_MODEL_HH

#include <opm/material/constraintsolvers/Co2BrinePressureEnthalpyProvider.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>

#include <opm/models/common/enthalpyproperties.hh>
#include <opm/models/common/multiphasebasemodel.hh>

#include <opm/models/flash/flashextensivequantities.hh>

#include <opm/models/co2store/co2storeindices.hh>
#include <opm/models/co2store/co2storeintensivequantities.hh>
#include <opm/models/co2store/co2storelocalresidual.hh>
#include <opm/models/co2store/co2storenewtonmethod.hh>
#include <opm/models/co2store/co2storeprimaryvariables.hh>
#include <opm/models/steam/steamratevector.hh>

#include <string>

namespace Opm {

template <class TypeTag>
class Co2StoreModel;

}

namespace Opm::Properties {

namespace TTag {
struct Co2StoreModel { using InheritsFrom = std::tuple<MultiPhaseBaseModel>; };
} // namespace TTag

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Co2StoreModel>
{ using type = Co2StoreLocalResidual<TypeTag>; };

template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::Co2StoreModel>
{ using type = Co2StoreNewtonMethod<TypeTag>; };

template<class TypeTag>
struct Model<TypeTag, TTag::Co2StoreModel>
{ using type = Co2StoreModel<TypeTag>; };

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::Co2StoreModel>
{ using type = Co2StorePrimaryVariables<TypeTag>; };

// the rate vectors are generic: one entry per equation
template<class TypeTag>
struct RateVector<TypeTag, TTag::Co2StoreModel>
{ using type = SteamRateVector<TypeTag>; };

template<class TypeTag>
struct BoundaryRateVector<TypeTag, TTag::Co2StoreModel>
{ using type = SteamBoundaryRateVector<TypeTag>; };

// CO2 from the p-h tables, brine properties of the CO2 storage fluid system
template<class TypeTag>
struct PressureEnthalpyProvider<TypeTag, TTag::Co2StoreModel>
{ using type = Co2BrinePressureEnthalpyProvider<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct IntensiveQuantities<TypeTag, TTag::Co2StoreModel>
{ using type = Co2StoreIntensiveQuantities<TypeTag>; };

// the advective flux is the multi-phase Darcy flux of the flux module; the
// energy conduction is done in the local residual
template<class TypeTag>
struct ExtensiveQuantities<TypeTag, TTag::Co2StoreModel>
{ using type = FlashExtensiveQuantities<TypeTag>; };

template<class TypeTag>
struct Indices<TypeTag, TTag::Co2StoreModel>
{ using type = Co2StoreIndices<TypeTag, /*PVOffset=*/0>; };

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::Co2StoreModel>
{ static constexpr bool value = false; };

// the generic energy module is not used: the model has its own energy equation
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::Co2StoreModel>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableWater<TypeTag, TTag::Co2StoreModel>
{ static constexpr int value = true; };

} // namespace Opm::Properties

namespace Opm {

template <class TypeTag>
class Co2StoreModel : public MultiPhaseBaseModel<TypeTag>
{
    using ParentType = MultiPhaseBaseModel<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

public:
    explicit Co2StoreModel(Simulator& simulator) : ParentType(simulator) {}

    static void registerParameters()
    {
        ParentType::registerParameters();
        Parameters::SetDefault<Parameters::EnableIntensiveQuantityCache>(true);
    }

    std::string primaryVarName(unsigned pvIdx) const
    {
        if (pvIdx == static_cast<unsigned>(Indices::pressure0Idx)) {
            return "pressure_aqueous";
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
        return eqIdx == static_cast<unsigned>(Indices::conti0EqIdx) ? "continuity^brine" : "continuity^CO2";
    }
};

} // namespace Opm

#endif
