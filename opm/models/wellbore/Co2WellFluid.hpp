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
 *
 * \brief The fluid of Opm::WellBore for CO2 and brine: the CO2 storage provider
 *        (Opm::Co2BrinePressureEnthalpyProvider) behind a flash from the pressure, the mixture
 *        enthalpy and the overall CO2 mass fraction Z.
 *
 * The provider takes a third variable s with a meaning (aqueous saturation, dissolved CO2
 * fraction, or brine fraction of the CO2). The flow in a well carries the overall composition
 * instead, which is the same as the dissolved fraction X in an aqueous-only state, one minus the
 * brine fraction Y in a CO2-only state, and fixes the aqueous saturation in between. The flash finds the
 * state that holds the composition:
 *  1. aqueous only, s = Z, if Z is at most the solubility at the temperature of that state,
 *  2. CO2 only, s = 1 - Z, if the brine content 1 - Z is at most the saturated value,
 *  3. otherwise both phases: for a given temperature (the enthalpy u of the CO2 component) the compositions, densities
 *     and enthalpies of the equilibrium phases are fixed and the aqueous saturation that gives Z is explicit, so that the
 *     state is a secant iteration on one scalar (Co2BrinePressureEnthalpyProvider::bothPhasesFromComposition(), about
 *     30 us). The iteration on the aqueous saturation of the provider is the fallback.
 * The regimes are tried in the order of the last result, since a well is marched through states that are close to each other.
 *
 * The wellbore model sees two flowing "phases": the aqueous phase (its "liquid") and the CO2-rich phase
 * (its "vapour"), which is liquid CO2, gas CO2 or both. The CO2-rich phase has the density of a homogeneous
 * mixture of its liquid and gas parts. The interfacial tension is a constant, it only enters the drift velocity.
 * The components are 0: CO2 and 1: brine (mass fractions).
 */
#ifndef OPM_CO2_WELL_FLUID_HPP
#define OPM_CO2_WELL_FLUID_HPP

#include <opm/common/Exceptions.hpp>

#include <opm/material/constraintsolvers/Co2BrinePressureEnthalpyProvider.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidsystems/Co2StoreFluidSystem.hpp>
#include <opm/models/wellbore/WellBore.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace Opm {

template <class Scalar>
class Co2WellFluid
{
    using Provider = Co2BrinePressureEnthalpyProvider<Scalar>;
    using FS = Co2StoreFluidSystem<Scalar>;
    using Eval1 = DenseAd::Evaluation<Scalar, 1>;
    template <class E>
    using ProviderState = typename Provider::template State<E>;

    static constexpr unsigned oilIdx = FS::oilPhaseIdx;      // liquid CO2
    static constexpr unsigned gasIdx = FS::gasPhaseIdx;      // gas CO2
    static constexpr unsigned waterIdx = FS::waterPhaseIdx;  // brine
    static constexpr unsigned co2Comp = FS::CO2Idx;

public:
    static constexpr std::size_t numComponents = 2;   //!< 0: CO2, 1: brine
    using Composition = std::array<Scalar, 2>;
    //! the composition (mass fractions) of what is injected
    using Stream = std::array<Scalar, 2>;

    //! interfacial tension of brine and CO2 [N/m]
    static constexpr Scalar surfaceTension() { return 0.04; }

    struct State
    {
        Scalar temperature;
        Scalar density;       //!< homogeneous (no slip) mixture density
        Scalar viscosity;     //!< mass-fraction harmonic mean of the phase viscosities
        Scalar quality;       //!< mass fraction of the CO2-rich phase
        Scalar voidFraction;  //!< volume fraction of the CO2-rich phase without slip
        Scalar densityLiquid; //!< aqueous phase
        Scalar densityVapor;  //!< CO2-rich phase (homogeneous mixture of its liquid and gas)
        Scalar surfaceTension;
        int region;           //!< meaning * 3 + CO2 region, see Co2BrinePressureEnthalpyProvider
        int meaning;
        Scalar aqueousSaturation;
        Scalar dissolvedCo2;  //!< CO2 mass fraction of the aqueous phase
    };

    static Scalar minPressure() { return Provider::minPressure(); }
    static Scalar maxPressure() { return Provider::maxPressure(); }

    //! the flash result: the primary variables of the provider and its state
    struct Flash
    {
        int meaning;
        Scalar s;
        ProviderState<Scalar> state;
    };

    /*!
     * \brief The state at pressure p, mixture enthalpy h and overall CO2 mass fraction z[0].
     *
     * Throws Opm::NumericalProblem if the state is outside the validity range of the provider.
     */
    static Flash flash(Scalar p, Scalar h, const Composition& z)
    {
        const Scalar Z = std::clamp(z[0], Scalar(0.0), Scalar(1.0));
        if (!(p >= minPressure() && p <= maxPressure())) {
            throw NumericalProblem("Co2WellFluid: pressure " + std::to_string(p) + " Pa outside the range");
        }
        const auto equilibrium = [&](Scalar T, Scalar& X, Scalar& Y) {
            Eval1 Xe, Ye;
            FS::equilibriumMassFractions(Eval1(T), Eval1(p), Xe, Ye);
            X = Xe.value();
            Y = Ye.value();
        };

        // The three regimes are tried in the order of the last result (a well is marched through states that are
        // close to each other): at most one of them holds the state.
        static thread_local int last = Provider::bothPhases;
        std::array<int, 3> order = {Provider::bothPhases, Provider::aqueousOnly, Provider::co2Only};
        std::swap(order[0], order[last == Provider::bothPhases ? 0 : (last == Provider::aqueousOnly ? 1 : 2)]);
        for (const int meaning : order) {
            if (meaning == Provider::aqueousOnly) {
                try {
                    auto st = Provider::template solve<Scalar>(p, h, Z, Provider::aqueousOnly);
                    Scalar X, Y;
                    equilibrium(st.temperature(), X, Y);
                    if (Z <= X * (1.0 + 1e-12)) {
                        last = meaning;
                        return {Provider::aqueousOnly, Z, std::move(st)};
                    }
                }
                catch (const NumericalProblem&) {
                }
            }
            else if (meaning == Provider::co2Only) {
                try {
                    auto st = Provider::template solve<Scalar>(p, h, 1.0 - Z, Provider::co2Only);
                    Scalar X, Y;
                    equilibrium(st.temperature(), X, Y);
                    if (1.0 - Z <= Y * (1.0 + 1e-12)) {
                        last = meaning;
                        return {Provider::co2Only, 1.0 - Z, std::move(st)};
                    }
                }
                catch (const NumericalProblem&) {
                }
            }
            else {
                Scalar Sw;
                ProviderState<Scalar> st;
                if (Provider::bothPhasesFromComposition(p, h, Z, Sw, st)) {
                    last = meaning;
                    return {Provider::bothPhases, Sw, std::move(st)};
                }
            }
        }
        // nothing found by the regimes above: the iteration on the aqueous saturation of the provider
        auto result = bothPhasesFlash_(p, h, Z);
        last = result.meaning;
        return result;
    }

    static State state(Scalar p, Scalar h, const Composition& z)
    {
        const auto f = flash(p, h, z);
        const auto& st = f.state;
        State s;
        s.temperature = st.temperature();
        s.meaning = f.meaning;
        s.region = st.regionId();

        std::array<Scalar, 3> mass;
        Scalar total = 0.0;
        for (const unsigned a : {oilIdx, gasIdx, waterIdx}) {
            mass[a] = st.saturation(a) * st.density(a);
            total += mass[a];
        }
        const Scalar mCo2Rich = mass[oilIdx] + mass[gasIdx];
        s.quality = mCo2Rich / total;
        s.densityLiquid = st.density(waterIdx);
        s.densityVapor = mCo2Rich > 0.0
            ? mCo2Rich / (mass[oilIdx] / st.density(oilIdx) + mass[gasIdx] / st.density(gasIdx))
            : st.density(gasIdx);
        s.aqueousSaturation = st.saturation(waterIdx);
        s.voidFraction = 1.0 - st.saturation(waterIdx);
        s.density = 1.0 / (s.quality / s.densityVapor + (1.0 - s.quality) / s.densityLiquid);
        Scalar invMu = 0.0;
        for (const unsigned a : {oilIdx, gasIdx, waterIdx}) {
            invMu += mass[a] / total / st.viscosity(a);
        }
        s.viscosity = 1.0 / invMu;
        s.surfaceTension = surfaceTension();
        s.dissolvedCo2 = st.massFraction(waterIdx, co2Comp);
        return s;
    }

    /*!
     * \brief Enthalpy [J/kg] of the fluid with composition z at rest at temperature T (found from the flash).
     */
    static Scalar equilibriumEnthalpy(Scalar p, Scalar T, const Composition& z)
    {
        Scalar lo = Provider::minEnthalpy(p);
        Scalar hi = Provider::maxEnthalpy(p);
        const Scalar mid = 0.5 * (lo + hi);
        const auto f = [&](Scalar h) {
            try {
                return flash(p, h, z).state.temperature() - T;
            }
            catch (const NumericalProblem&) {
                return h < mid ? -1.0e3 : 1.0e3;   // outside the temperature range of the provider
            }
        };
        Scalar flo = f(lo);
        Scalar fhi = f(hi);
        if (!(flo <= 0.0 && fhi >= 0.0)) {
            throw NumericalProblem("Co2WellFluid: the temperature " + std::to_string(T) + " K is outside the range");
        }
        for (int it = 0; it < 100 && hi - lo > 1e-6 * (1.0 + std::abs(hi)); ++it) {
            const Scalar m = 0.5 * (lo + hi);
            (f(m) < 0.0 ? lo : hi) = m;
        }
        return 0.5 * (lo + hi);
    }

    //! Enthalpy [J/kg] of the injected stream at temperature T and pressure p
    static Scalar injectionEnthalpy(Scalar p, Scalar T, const Stream& z)
    {
        if (z[0] >= 1.0) {
            return Provider::template injectionEnthalpy<Scalar>(true, p, T);
        }
        if (z[0] <= 0.0) {
            return Provider::template injectionEnthalpy<Scalar>(false, p, T);
        }
        return equilibriumEnthalpy(p, T, z);
    }

    static Composition injectionComposition(const Stream& z) { return z; }

    //! Corey relative permeability exponent of the reservoir cell, kr = S^exponent of every phase
    struct RelativePermeability
    {
        Scalar exponent{2.0};
    };

    /*!
     * \brief The inflow of a producer perforation from a reservoir cell.
     *
     * The cell (pressure, temperature, aqueous saturation Sw, and the dissolved CO2 mass fraction X
     * for Sw >= 1 or the brine fraction Y for Sw <= 0, see the provider) delivers
     *   productivity = WI sum_a rho_a kr_a / mu_a,  enthalpy and composition weighted by the same mass rates,
     * and injectivity = WI sum_a kr_a / mu_a. WI is the well index [m3].
     */
    static WellBorePerforation<Scalar, 2> perforation(Scalar md, Scalar pRes, Scalar TRes, Scalar Sw, Scalar X,
                                                      Scalar Y, Scalar wellIndex, const RelativePermeability& kr = {})
    {
        Scalar h, s;
        int meaning;
        Provider::initialVariables(pRes, TRes, Sw, X, Y, h, s, meaning);
        const auto st = Provider::template solve<Scalar>(pRes, h, s, meaning);
        Scalar productivity = 0.0, hSum = 0.0, co2Sum = 0.0, mobility = 0.0;
        for (const unsigned a : {oilIdx, gasIdx, waterIdx}) {
            const Scalar lambda = std::pow(st.saturation(a), kr.exponent) / st.viscosity(a);
            const Scalar c = wellIndex * st.density(a) * lambda;
            productivity += c;
            hSum += c * st.enthalpy(a);
            co2Sum += c * st.massFraction(a, co2Comp);
            mobility += wellIndex * lambda;
        }
        WellBorePerforation<Scalar, 2> perf;
        perf.md = md;
        perf.pressure = pRes;
        perf.productivity = productivity;
        perf.inflowEnthalpy = hSum / productivity;
        perf.inflowComposition = {co2Sum / productivity, 1.0 - co2Sum / productivity};
        perf.injectivity = mobility;
        return perf;
    }

private:
    // overall CO2 mass fraction of a state
    template <class E>
    static E overallCo2_(const ProviderState<E>& st)
    {
        E co2 = 0.0, total = 0.0;
        for (const unsigned a : {oilIdx, gasIdx, waterIdx}) {
            const E m = st.saturation(a) * st.density(a);
            total += m;
            co2 += m * st.massFraction(a, co2Comp);
        }
        return co2 / total;
    }

    static Flash bothPhasesFlash_(Scalar p, Scalar h, Scalar Z)
    {
        // g(Sw) = overall CO2 mass fraction - Z decreases with the aqueous saturation
        const auto g = [&](Scalar Sw) -> std::optional<std::pair<Scalar, Scalar>> {
            try {
                const auto st = Provider::template solve<Eval1>(Eval1(p), Eval1(h), Eval1::createVariable(Sw, 0),
                                                                Provider::bothPhases);
                const Eval1 z = overallCo2_<Eval1>(st);
                return std::pair{z.value() - Z, z.derivative(0)};
            }
            catch (const NumericalProblem&) {
                return std::nullopt;
            }
        };
        Scalar lo = 0.0, hi = 1.0;
        Scalar Sw = 0.5;
        std::optional<Scalar> lastValid;
        for (int it = 0; it < 80; ++it) {
            const auto r = g(Sw);
            Scalar next;
            if (!r) {
                // outside the range of the provider: move back towards the last valid saturation
                if (!lastValid) {
                    // nothing valid yet: probe the other saturations
                    next = it % 2 == 0 ? 0.5 * Sw : Sw + 0.5 * (1.0 - Sw);
                }
                else {
                    (Sw > *lastValid ? hi : lo) = Sw;
                    next = 0.5 * (Sw + *lastValid);
                }
                Sw = next;
                continue;
            }
            lastValid = Sw;
            if (std::abs(r->first) < 1e-13) {
                break;
            }
            // g decreases with Sw
            (r->first > 0.0 ? lo : hi) = Sw;
            next = r->second < 0.0 ? Sw - r->first / r->second : 0.5 * (lo + hi);
            if (!(next > lo && next < hi)) {
                next = 0.5 * (lo + hi);
            }
            if (std::abs(next - Sw) < 1e-15) {
                break;
            }
            Sw = next;
        }
        if (!lastValid) {
            throw NumericalProblem("Co2WellFluid: no state for p = " + std::to_string(p) + " Pa, h = " +
                                   std::to_string(h) + " J/kg, CO2 mass fraction " + std::to_string(Z));
        }
        auto st = Provider::template solve<Scalar>(p, h, Sw, Provider::bothPhases);
        if (!(std::abs(overallCo2_<Scalar>(st) - Z) < 1e-8)) {
            throw NumericalProblem("Co2WellFluid: the aqueous saturation did not converge for p = " +
                                   std::to_string(p) + " Pa, h = " + std::to_string(h) + " J/kg");
        }
        return {Provider::bothPhases, Sw, std::move(st)};
    }
};

} // namespace Opm

#endif
