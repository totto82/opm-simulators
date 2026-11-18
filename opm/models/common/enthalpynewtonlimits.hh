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
 * \brief Newton step limiting for models with pressure and specific enthalpy
 *        as primary variables.
 *
 * T(p,h) and the saturations are only continuous, not differentiable, across
 * phase boundaries, and the flash is only valid in a limited (p,h) range. The
 * step limiter keeps an update inside that range, limits its size, and cuts a
 * step that crosses a phase boundary so that it lands just beyond it. A full
 * step across the kink otherwise makes the iteration jump back and forth
 * over the boundary.
 *
 * What is specific to the fluid description is supplied by a policy
 * (see Opm::PressureEnthalpyLimitsPolicy): the validity range, the largest
 * enthalpy step and an identifier of the phase region of a (p,h) point.
 * For pure water the policy is built from the steam tables (saturated liquid
 * and vapour enthalpies); for a mixture the region would be the number of
 * phases of the flash.
 */
#ifndef OPM_ENTHALPY_NEWTON_LIMITS_HH
#define OPM_ENTHALPY_NEWTON_LIMITS_HH

#include <algorithm>
#include <concepts>

namespace Opm {

template <class Policy, class Scalar>
concept PressureEnthalpyLimitsPolicy = requires(const Policy& policy, Scalar p, Scalar h) {
    { policy.minPressure() } -> std::convertible_to<Scalar>;
    { policy.maxPressure() } -> std::convertible_to<Scalar>;
    //! lowest/highest enthalpy of the validity range at pressure p
    { policy.minEnthalpy(p) } -> std::convertible_to<Scalar>;
    { policy.maxEnthalpy(p) } -> std::convertible_to<Scalar>;
    //! largest allowed enthalpy change in one iteration at pressure p
    { policy.maxEnthalpyStep(p) } -> std::convertible_to<Scalar>;
    //! identifier of the phase region of (p,h); a change of the value means a phase boundary was crossed
    { policy.region(p, h) } -> std::convertible_to<int>;
};

/*!
 * \brief Limits the update (p0,h0) -> (p,h); the result is returned in p and h.
 *
 * - the relative pressure change is at most maxRelPressureChange, and p stays inside the validity range,
 * - the enthalpy change is at most policy.maxEnthalpyStep(p), and h stays inside the validity range,
 * - a step that leaves the phase region of (p0,h0) is cut to land just inside the next region
 *   (bisection on the segment, then an overshoot of a fraction of the remaining step).
 */
template <class Scalar, class Policy>
void limitPressureEnthalpyStep(const Scalar p0,
                               const Scalar h0,
                               Scalar& p,
                               Scalar& h,
                               const Policy& policy,
                               const Scalar maxRelPressureChange = 0.2,
                               const Scalar overshoot = 0.1)
{
    static_assert(PressureEnthalpyLimitsPolicy<Policy, Scalar>);

    // pressure: limited relative change, inside the validity range
    p = std::clamp(p, p0 * (1.0 - maxRelPressureChange), p0 * (1.0 + maxRelPressureChange));
    const Scalar margin = 1.0e-3;
    p = std::clamp(p, policy.minPressure() * (1.0 + margin), policy.maxPressure() * (1.0 - margin));

    // enthalpy: limited step, inside the validity range at the new pressure
    const Scalar maxStep = policy.maxEnthalpyStep(p);
    h = std::clamp(h, h0 - maxStep, h0 + maxStep);
    const Scalar hmin = policy.minEnthalpy(p);
    const Scalar hmax = policy.maxEnthalpy(p);
    const Scalar hmargin = 1.0e-3 * (hmax - hmin);
    h = std::clamp(h, hmin + hmargin, hmax - hmargin);

    // phase boundary
    const int region0 = policy.region(p0, h0);
    if (policy.region(p, h) == region0) {
        return;
    }

    const auto at = [&](Scalar alpha) { return policy.region(p0 + alpha * (p - p0), h0 + alpha * (h - h0)); };

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
    const Scalar alpha = hi + overshoot * (1.0 - hi);
    const Scalar pNew = p0 + alpha * (p - p0);
    const Scalar hNew = h0 + alpha * (h - h0);
    p = pNew;
    h = hNew;
}

} // namespace Opm

#endif
