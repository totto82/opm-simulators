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

        // pressure: at most 20 % per iteration, inside the validity range of the tables
        constexpr Scalar maxRelPressureChange = 0.2;
        Scalar p = std::clamp(nextValue[Indices::pressure0Idx],
                              old[Indices::pressure0Idx] * (1.0 - maxRelPressureChange),
                              old[Indices::pressure0Idx] * (1.0 + maxRelPressureChange));
        const Scalar margin = 1.0e-3;
        p = std::clamp(p, Flash::minPressure() * (1.0 + margin), Flash::maxPressure() * (1.0 - margin));
        nextValue[Indices::pressure0Idx] = p;

        // enthalpy: a step is limited to a fraction of the latent heat, so a
        // cell cannot jump across the whole two-phase region in one iteration,
        // and kept inside the validity range at the new pressure
        constexpr Scalar maxLatentFraction = 0.25;
        const Scalar hl = Flash::saturatedLiquidEnthalpy(p);
        const Scalar hg = Flash::saturatedVaporEnthalpy(p);
        const Scalar maxStep = maxLatentFraction * (hg - hl);
        Scalar h = std::clamp(nextValue[Indices::enthalpy0Idx],
                              old[Indices::enthalpy0Idx] - maxStep,
                              old[Indices::enthalpy0Idx] + maxStep);
        const Scalar hmin = Flash::minEnthalpy(p);
        const Scalar hmax = Flash::maxEnthalpy(p);
        const Scalar hmargin = 1.0e-3 * (hmax - hmin);
        h = std::clamp(h, hmin + hmargin, hmax - hmargin);
        nextValue[Indices::enthalpy0Idx] = h;

        // a step that crosses a phase boundary (saturated liquid or vapour
        // enthalpy) is cut to land just beyond the boundary. A full step across
        // the kink of T(p,h) and S_g(p,h) otherwise makes the iteration jump back
        // and forth over the boundary.
        const Scalar p0 = old[Indices::pressure0Idx];
        const Scalar h0 = old[Indices::enthalpy0Idx];
        const int region0 = region_(p0, h0);
        if (region_(p, h) != region0) {
            const auto at = [&](Scalar alpha) { return region_(p0 + alpha * (p - p0), h0 + alpha * (h - h0)); };

            // first sample that has left the old region, then bisection
            constexpr int numSamples = 20;
            Scalar lo = 0.0;
            Scalar hi = 1.0;
            for (int i = 1; i <= numSamples; ++i) {
                const Scalar alpha = static_cast<Scalar>(i) / numSamples;
                if (at(alpha) != region0) {
                    hi = alpha;
                    break;
                }
                lo = alpha;
            }
            for (int it = 0; it < 40; ++it) {
                const Scalar mid = 0.5 * (lo + hi);
                if (at(mid) == region0) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            // land a little inside the new region
            constexpr Scalar overshoot = 0.1;
            const Scalar alpha = hi + overshoot * (1.0 - hi);
            nextValue[Indices::pressure0Idx] = p0 + alpha * (p - p0);
            nextValue[Indices::enthalpy0Idx] = h0 + alpha * (h - h0);
        }
    }

    // 0: liquid, 1: two-phase, 2: vapour
    static int region_(Scalar p, Scalar h)
    {
        if (h <= Flash::saturatedLiquidEnthalpy(p)) {
            return 0;
        }
        return h >= Flash::saturatedVaporEnthalpy(p) ? 2 : 1;
    }
};

} // namespace Opm

#endif
