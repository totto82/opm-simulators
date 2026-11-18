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
 * \brief Newton method of the steam model with chopping of the (p, h) update.
 */
#ifndef OPM_STEAM_NEWTON_METHOD_HH
#define OPM_STEAM_NEWTON_METHOD_HH

#include <opm/material/constraintsolvers/WaterPHFlash.hpp>

#include <opm/models/common/enthalpynewtonlimits.hh>
#include <opm/models/nonlinear/newtonmethod.hh>

#include <algorithm>
#include <cmath>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

template <class TypeTag>
class SteamNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Flash = WaterPHFlash<Scalar>;

public:
    explicit SteamNewtonMethod(Simulator& simulator) : ParentType(simulator) {}

protected:
    friend ParentType;
    friend NewtonMethod<TypeTag>;

    void updatePrimaryVariables_(unsigned /* globalDofIdx */,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& /* currentResidual */)
    {
        const PrimaryVariables old = currentValue;
        nextValue = old;
        nextValue -= update;

        Scalar p = nextValue[Indices::pressure0Idx];
        Scalar h = nextValue[Indices::enthalpy0Idx];
        limitPressureEnthalpyStep(old[Indices::pressure0Idx], old[Indices::enthalpy0Idx], p, h, WaterLimits{});
        nextValue[Indices::pressure0Idx] = p;
        nextValue[Indices::enthalpy0Idx] = h;
    }

    //! Validity range and phase regions of pure water from the steam tables
    struct WaterLimits
    {
        Scalar minPressure() const { return Flash::minPressure(); }
        Scalar maxPressure() const { return Flash::maxPressure(); }
        Scalar minEnthalpy(Scalar p) const { return Flash::minEnthalpy(p); }
        Scalar maxEnthalpy(Scalar p) const { return Flash::maxEnthalpy(p); }

        // a step is limited to a fraction of the latent heat, so a cell cannot
        // jump across the whole two-phase region in one iteration
        Scalar maxEnthalpyStep(Scalar p) const
        {
            constexpr Scalar maxLatentFraction = 0.25;
            return maxLatentFraction * (Flash::saturatedVaporEnthalpy(p) - Flash::saturatedLiquidEnthalpy(p));
        }

        // 0: liquid, 1: two-phase, 2: vapour
        int region(Scalar p, Scalar h) const
        {
            if (h <= Flash::saturatedLiquidEnthalpy(p)) {
                return 0;
            }
            return h >= Flash::saturatedVaporEnthalpy(p) ? 2 : 1;
        }
    };
};

} // namespace Opm

#endif
