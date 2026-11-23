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
 * \brief Flash of brine and CO2 for Opm::PhModel: the p-h tables of CO2 and the brine properties of the CO2 storage
 *        fluid system (Opm::Co2BrinePressureEnthalpyProvider of opm-common). Two components (brine, CO2); the primary
 *        variables are (p, h, s) where the meaning of s (aqueous saturation, dissolved CO2, or brine in the CO2)
 *        is switched when a phase appears or disappears.
 */
#ifndef OPM_PH_CO2STORE_FLASH_HH
#define OPM_PH_CO2STORE_FLASH_HH

#include <opm/common/Exceptions.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/material/constraintsolvers/Co2BrinePressureEnthalpyProvider.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/Co2StoreFluidSystem.hpp>

#include <opm/models/ph/phflash.hh>
#include <opm/models/wellbore/Co2WellFluid.hpp>

#include <algorithm>
#include <stdexcept>

namespace Opm {

template <class Scalar>
struct Co2StoreFlash
{
    using Provider = Co2BrinePressureEnthalpyProvider<Scalar>;
    using FluidSystem = Co2StoreFluidSystem<Scalar>;

    static constexpr int numComponents = 2;
    static constexpr int numSwitchVariables = 1;
    static constexpr bool oilPhaseEnabled = true;
    static constexpr bool hasCo2Output = true;

    // meaning of the switched variable
    static constexpr int bothPhases = Provider::bothPhases;

    static constexpr const char* saturationKeyword = "SWAT";
    static constexpr bool saturationRequired = true;

    static Scalar minPressure() { return Provider::minPressure(); }
    static Scalar maxPressure() { return Provider::maxPressure(); }
    static Scalar minEnthalpy(Scalar p) { return Provider::minEnthalpy(p); }
    static Scalar maxEnthalpy(Scalar p) { return Provider::maxEnthalpy(p); }
    static Scalar maxEnthalpyStep(Scalar p) { return Provider::maxEnthalpyStep(p); }

    template <class Evaluation>
    static auto solve(const Evaluation& p, const Evaluation& h, const Evaluation& s, int meaning)
    { return Provider::template solve<Evaluation>(p, h, s, meaning); }

    /*!
     * \brief Limits the Newton update (p0, h0, s0) -> (p, h, s): size and range of the update, a state that exists,
     *        a cut at the phase boundaries of the CO2 and the switch of the meaning of s.
     */
    static void limitUpdate(Scalar p0, Scalar h0, Scalar s0, int meaning0,
                            Scalar& p, Scalar& h, Scalar& s, int& meaning)
    {
        using Eval1 = DenseAd::Evaluation<Scalar, 1>;
        meaning = meaning0;

        // --- size and range of the update ---

        // pressure: at most 20 % per iteration, inside the validity range
        constexpr Scalar maxRelPressureChange = 0.2;
        constexpr Scalar margin = 1.0e-3;
        p = std::clamp(p, p0 * (1.0 - maxRelPressureChange), p0 * (1.0 + maxRelPressureChange));
        p = std::clamp(p, minPressure() * (1.0 + margin), maxPressure() * (1.0 - margin));

        // enthalpy: limited step inside the range of the temperature limits at the new pressure
        const Scalar maxStep = maxEnthalpyStep(p);
        h = std::clamp(h, h0 - maxStep, h0 + maxStep);
        const Scalar hmin = minEnthalpy(p), hmax = maxEnthalpy(p);
        h = std::clamp(h, hmin + margin * (hmax - hmin), hmax - margin * (hmax - hmin));

        // third variable: the aqueous saturation changes by at most 0.2 per iteration and stays in [0,1];
        // a mass fraction is not negative
        constexpr Scalar maxSaturationChange = 0.2;
        if (meaning == bothPhases) {
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
    }

    //! primary variables of the initial state from PRESSURE, TEMPI and SWAT
    static void initialVariables(Scalar p, Scalar T, Scalar Sw, Scalar& h, Scalar& s, int& meaning)
    { Provider::initialVariables(p, T, Sw, /*X=*/0.0, /*Y=*/0.0, h, s, meaning); }

    //! enthalpy of the injected stream (CO2 if gas, otherwise brine) at pressure p and temperature T
    template <class Evaluation>
    static Evaluation injectionEnthalpy(bool gas, const Evaluation& p, const Evaluation& T)
    { return Provider::injectionEnthalpy(gas, p, T); }

    //! the component of a source of gas (CO2) or water (brine)
    static unsigned sourceComponentIndex(bool gas)
    { return gas ? static_cast<unsigned>(FluidSystem::CO2Idx) : static_cast<unsigned>(FluidSystem::BrineIdx); }

    // --- wells ---

    template <class S>
    using WellFluid = Co2WellFluid<S>;

    //! the stream a well injects: the composition {CO2, brine}
    template <class Fluid>
    static typename Fluid::Stream injectionStream(bool gas)
    { return gas ? typename Fluid::Stream{1.0, 0.0} : typename Fluid::Stream{0.0, 1.0}; }

    //! equation (component) index of the component k of the wellbore fluid (0: CO2, 1: brine)
    static unsigned wellComponentIndex(std::size_t k)
    { return k == 0 ? static_cast<unsigned>(FluidSystem::CO2Idx) : static_cast<unsigned>(FluidSystem::BrineIdx); }
    static const char* wellComponentName(std::size_t k) { return k == 0 ? "CO2" : "brine"; }
    static constexpr const char* gasName = "CO2";
    static constexpr const char* waterName = "brine";

    static void checkDeck(const EclipseState& eclState)
    {
        const auto& runspec = eclState.runspec();
        if (!eclState.getSimulationConfig().isEnthalpy() || !runspec.co2Storage()) {
            throw std::runtime_error("flow_ph with the CO2 storage flash requires the ENTHALPY and CO2STORE keywords in the "
                                     "RUNSPEC section (and no THERMAL or TEMP).");
        }
        if (!runspec.phases().active(Phase::OIL) || !runspec.phases().active(Phase::GAS) ||
            !runspec.phases().active(Phase::WATER)) {
            throw std::runtime_error("flow_ph with the CO2 storage flash requires the phases OIL (liquid CO2), GAS and WATER.");
        }
    }
};

template <class Scalar>
struct PhFlashSelector<PhFlashType::Co2Store, Scalar>
{ using type = Co2StoreFlash<Scalar>; };

} // namespace Opm

#endif
