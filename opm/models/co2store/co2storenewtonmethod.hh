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
 * \brief Newton method of the CO2 storage model: limits the updates of the primary variables,
 *        cuts a step that crosses a phase boundary of the CO2 and switches the meaning of the
 *        third variable when a phase appears or disappears.
 */
#ifndef OPM_CO2STORE_NEWTON_METHOD_HH
#define OPM_CO2STORE_NEWTON_METHOD_HH

#include <opm/common/Exceptions.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/models/common/enthalpyproperties.hh>
#include <opm/models/nonlinear/newtonmethod.hh>

#include <algorithm>
#include <cmath>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

template <class TypeTag>
class Co2StoreNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Provider = GetPropType<TypeTag, Properties::PressureEnthalpyProvider>;
    using Eval1 = DenseAd::Evaluation<Scalar, 1>;

public:
    explicit Co2StoreNewtonMethod(Simulator& simulator) : ParentType(simulator) {}

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

        int meaning = old.meaning();
        const Scalar p0 = old[Indices::pressure0Idx];
        const Scalar h0 = old[Indices::enthalpy0Idx];
        const Scalar s0 = old[Indices::switchIdx];

        // --- size and range of the update ---

        // pressure: at most 20 % per iteration, inside the validity range
        constexpr Scalar maxRelPressureChange = 0.2;
        constexpr Scalar margin = 1.0e-3;
        Scalar p = std::clamp(nextValue[Indices::pressure0Idx], p0 * (1.0 - maxRelPressureChange), p0 * (1.0 + maxRelPressureChange));
        p = std::clamp(p, Provider::minPressure() * (1.0 + margin), Provider::maxPressure() * (1.0 - margin));

        // enthalpy: limited step inside the range of the temperature limits at the new pressure
        const Scalar maxStep = Provider::maxEnthalpyStep(p);
        Scalar h = std::clamp(nextValue[Indices::enthalpy0Idx], h0 - maxStep, h0 + maxStep);
        const Scalar hmin = Provider::minEnthalpy(p), hmax = Provider::maxEnthalpy(p);
        h = std::clamp(h, hmin + margin * (hmax - hmin), hmax - margin * (hmax - hmin));

        // third variable: the aqueous saturation changes by at most 0.2 per iteration and stays in [0,1];
        // a mass fraction is not negative
        constexpr Scalar maxSaturationChange = 0.2;
        Scalar s = nextValue[Indices::switchIdx];
        if (meaning == Provider::bothPhases) {
            s = std::clamp(s, std::max(Scalar(0.0), s0 - maxSaturationChange), std::min(Scalar(1.0), s0 + maxSaturationChange));
        }
        else {
            s = std::max(s, Scalar(0.0));
        }

        // --- the state of the new variables must exist: otherwise take a part of the step ---
        const auto stateExists = [&](Scalar a) {
            try {
                Provider::solve(Eval1(p0 + a * (p - p0)), Eval1(h0 + a * (h - h0)), Eval1(s0 + a * (s - s0)), meaning);
                return true;
            }
            catch (const NumericalProblem&) {
                return false;
            }
        };
        Scalar scale = 1.0;
        for (int i = 0; i < 12 && !stateExists(scale); ++i) {
            scale *= 0.5;
        }
        if (scale < 1.0) {
            p = p0 + scale * (p - p0);
            h = h0 + scale * (h - h0);
            s = s0 + scale * (s - s0);
        }

        // --- a step that crosses a phase boundary of the CO2 is cut to land just beyond it ---
        // (T and the phase split are only continuous across the edges of the two-phase region of the CO2)
        const auto regionAt = [&](Scalar a) {
            try {
                return Provider::solve(Eval1(p0 + a * (p - p0)), Eval1(h0 + a * (h - h0)), Eval1(s0 + a * (s - s0)),
                                       meaning).regionId();
            }
            catch (const NumericalProblem&) {
                return -1;
            }
        };
        const int region0 = regionAt(0.0);
        if (region0 >= 0 && regionAt(1.0) != region0) {
            constexpr int numSamples = 20;
            Scalar lo = 0.0, hi = 1.0;
            for (int i = 1; i <= numSamples; ++i) {
                const Scalar a = static_cast<Scalar>(i) / numSamples;
                if (regionAt(a) != region0) { hi = a; break; }
                lo = a;
            }
            for (int it = 0; it < 40; ++it) {
                const Scalar mid = 0.5 * (lo + hi);
                (regionAt(mid) == region0 ? lo : hi) = mid;
            }
            constexpr Scalar overshoot = 0.1;
            const Scalar a = hi + overshoot * (1.0 - hi);
            if (regionAt(a) >= 0) {
                const Scalar pn = p0 + a * (p - p0), hn = h0 + a * (h - h0), sn = s0 + a * (s - s0);
                p = pn; h = hn; s = sn;
            }
        }

        // --- a phase appears or disappears: switch the meaning of the third variable ---
        try {
            Provider::adaptVariables(p, h, s, meaning);
        }
        catch (const NumericalProblem&) {
            // the state with the other meaning does not exist: keep the meaning, the next iteration decides
        }

        nextValue[Indices::pressure0Idx] = p;
        nextValue[Indices::enthalpy0Idx] = h;
        nextValue[Indices::switchIdx] = s;
        nextValue.setMeaning(meaning);
    }
};

} // namespace Opm

#endif
