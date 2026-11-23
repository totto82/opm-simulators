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
 * \brief The fluid of Opm::WellBore for water and steam: a steam pressure-enthalpy provider
 *        (raw IAPWS or tabulated) with the homogeneous mixture properties of the wellbore model,
 *        and the inflow description of a producer perforation from a reservoir cell.
 */
#ifndef OPM_STEAM_WELL_FLUID_HPP
#define OPM_STEAM_WELL_FLUID_HPP

#include <opm/material/constraintsolvers/WaterPHFlash.hpp>
#include <opm/material/fluidsystems/SteamFluidSystem.hpp>
#include <opm/models/wellbore/WellBore.hpp>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace Opm {

template <class Scalar, class Provider>
class SteamWellFluid
{
    using FluidSystem = SteamFluidSystem<Scalar>;
    using Flash = WaterPHFlash<Scalar>;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

public:
    //! one component
    static constexpr std::size_t numComponents = 1;
    using Composition = std::array<Scalar, 1>;
    //! what is injected: steam (true) or liquid water (false)
    using Stream = bool;

    struct State
    {
        Scalar temperature;
        Scalar density;       //!< homogeneous mixture density [kg/m3]
        Scalar viscosity;     //!< mixture viscosity, harmonic in the mass fractions [Pa s]
        Scalar quality;       //!< steam mass fraction
        Scalar voidFraction;  //!< steam volume fraction without slip
        Scalar densityLiquid; //!< [kg/m3]
        Scalar densityVapor;  //!< [kg/m3]
        Scalar surfaceTension;//!< [N/m]
        int region;           //!< 0: liquid, 1: two-phase, 2: steam
    };

    static Scalar minPressure() { return Provider::minPressure(); }
    static Scalar maxPressure() { return Provider::maxPressure(); }

    static State state(Scalar p, Scalar h, const Composition&) { return state(p, h); }

    //! Properties of the fluid at (p, h). Throws std::invalid_argument outside the validity range.
    static State state(Scalar p, Scalar h)
    {
        if (!(p >= Provider::minPressure() && p <= Provider::maxPressure())) {
            throw std::invalid_argument("SteamWellFluid: pressure " + std::to_string(p) + " Pa outside the range");
        }
        if (!(h >= Provider::minEnthalpy(p) && h <= Provider::maxEnthalpy(p))) {
            throw std::invalid_argument("SteamWellFluid: enthalpy " + std::to_string(h) + " J/kg outside the range");
        }
        const auto s = Provider::solve(p, h);
        const Scalar x = s.steamQuality();
        const Scalar rhoL = s.density(waterPhaseIdx);
        const Scalar rhoG = s.density(gasPhaseIdx);
        State st;
        st.temperature = s.temperature();
        st.quality = x;
        st.voidFraction = s.saturation(gasPhaseIdx);
        st.densityLiquid = rhoL;
        st.densityVapor = rhoG;
        st.surfaceTension = surfaceTension(st.temperature);
        st.density = 1.0 / (x / rhoG + (1.0 - x) / rhoL);
        st.viscosity = 1.0 / (x / s.viscosity(gasPhaseIdx) + (1.0 - x) / s.viscosity(waterPhaseIdx));
        st.region = s.regionId();
        return st;
    }

    //! Surface tension of water against its vapour [N/m] (IAPWS 1994), zero at and above the critical temperature
    static Scalar surfaceTension(Scalar T)
    {
        constexpr Scalar Tc = 647.096;
        const Scalar tau = 1.0 - T / Tc;
        return tau > 0.0 ? 0.2358 * std::pow(tau, 1.256) * (1.0 - 0.625 * tau) : 0.0;
    }

    //! h minus the saturated enthalpy of the boundary: 0 is liquid / two-phase (h - hl), 1 is two-phase / steam (h - hv)
    static Scalar phaseChangeMargin(Scalar p, Scalar h, int boundary)
    { return h - (boundary == 0 ? Flash::saturatedLiquidEnthalpy(p) : Flash::saturatedVaporEnthalpy(p)); }

    /*!
     * Enthalpy of the injected stream: liquid water or steam at temperature T and pressure p. Liquid water
     * above the saturation temperature and steam below it do not exist at p: Throws std::invalid_argument
     * (the stream is not clipped to the saturation temperature as for the SOURCE keyword, because a wellbore
     * solution with a different temperature than specified would be wrong).
     */
    static Scalar injectionEnthalpy(Scalar p, Scalar T, bool vapour)
    {
        const Scalar Tsat = Flash::saturationTemperature(p);
        constexpr Scalar tolerance = 1.0e-6;   // K
        if ((!vapour && T > Tsat + tolerance) || (vapour && T < Tsat - tolerance)) {
            throw std::invalid_argument(std::string(vapour ? "steam" : "liquid water") + " at " + std::to_string(T) +
                                        " K does not exist at " + std::to_string(p) + " Pa (saturation temperature " +
                                        std::to_string(Tsat) + " K)");
        }
        return Provider::injectionEnthalpy(vapour ? gasPhaseIdx : waterPhaseIdx, p, T);
    }

    static Composition injectionComposition(bool) { return {1.0}; }

    //! water at rest at temperature T: liquid (at the saturation temperature at most)
    static Scalar equilibriumEnthalpy(Scalar p, Scalar T, const Composition&)
    { return Provider::injectionEnthalpy(waterPhaseIdx, p, T); }

    static Scalar phaseChangeMargin(Scalar p, Scalar h, const Composition&, int boundary)
    { return phaseChangeMargin(p, h, boundary); }

    //! Corey relative permeabilities of the reservoir cell
    struct RelativePermeability
    {
        Scalar liquidExponent{2.0};
        Scalar vaporExponent{2.0};
    };

    /*!
     * \brief The inflow of a producer perforation from a reservoir cell at (pRes, hRes).
     *
     * The cell fluid flows with the mobilities of its liquid and steam,
     * productivity = WI sum_a rho_a kr_a / mu_a, at the mobility weighted enthalpy.
     * The injectivity is WI sum_a kr_a / mu_a. WI is the well index [m3] (WI * mobility is m3/(s Pa)).
     */
    static WellBorePerforation<Scalar> perforation(Scalar md, Scalar pRes, Scalar hRes, Scalar wellIndex,
                                                   const RelativePermeability& kr = {})
    {
        const auto s = Provider::solve(pRes, hRes);
        const Scalar Sg = s.saturation(gasPhaseIdx);
        const Scalar lamL = std::pow(1.0 - Sg, kr.liquidExponent) / s.viscosity(waterPhaseIdx);
        const Scalar lamG = std::pow(Sg, kr.vaporExponent) / s.viscosity(gasPhaseIdx);
        const Scalar cL = wellIndex * s.density(waterPhaseIdx) * lamL;
        const Scalar cG = wellIndex * s.density(gasPhaseIdx) * lamG;
        WellBorePerforation<Scalar> perf;
        perf.md = md;
        perf.pressure = pRes;
        perf.productivity = cL + cG;
        perf.inflowEnthalpy = (cL * s.enthalpy(waterPhaseIdx) + cG * s.enthalpy(gasPhaseIdx)) / (cL + cG);
        perf.injectivity = wellIndex * (lamL + lamG);
        return perf;
    }
};

} // namespace Opm

#endif
