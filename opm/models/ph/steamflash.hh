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
 * \brief Flash of pure water (liquid and steam) for Opm::PhModel: the IAPWS-IF97 provider of opm-common,
 *        raw or tabulated. One component, no switched variable: the primary variables are (p, h).
 */
#ifndef OPM_PH_STEAM_FLASH_HH
#define OPM_PH_STEAM_FLASH_HH

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/material/constraintsolvers/SteamPressureEnthalpyProvider.hpp>
#include <opm/material/constraintsolvers/TabulatedSteamPressureEnthalpyProvider.hpp>
#include <opm/material/fluidsystems/SteamFluidSystem.hpp>

#include <opm/models/common/enthalpynewtonlimits.hh>
#include <opm/models/ph/phflash.hh>
#include <opm/models/wellbore/SteamWellFluid.hpp>

#include <stdexcept>

namespace Opm {

template <class Scalar, class ProviderT>
struct SteamFlash
{
    using Provider = ProviderT;
    using FluidSystem = SteamFluidSystem<Scalar>;

    static constexpr int numComponents = 1;
    static constexpr int numSwitchVariables = 0;
    static constexpr bool oilPhaseEnabled = false;
    static constexpr bool hasCo2Output = false;

    //! the deck keyword of the initial saturation and whether it is required
    static constexpr const char* saturationKeyword = "SGAS";
    static constexpr bool saturationRequired = false;

    static Scalar minPressure() { return Provider::minPressure(); }
    static Scalar maxPressure() { return Provider::maxPressure(); }
    static Scalar minEnthalpy(Scalar p) { return Provider::minEnthalpy(p); }
    static Scalar maxEnthalpy(Scalar p) { return Provider::maxEnthalpy(p); }
    static Scalar maxEnthalpyStep(Scalar p) { return Provider::maxEnthalpyStep(p); }

    // the switched variable and its meaning do not exist
    template <class Evaluation>
    static auto solve(const Evaluation& p, const Evaluation& h, const Evaluation& /*s*/, int /*meaning*/)
    { return Provider::template solve<Evaluation>(p, h); }

    //! limits the Newton update (p0, h0) -> (p, h)
    static void limitUpdate(Scalar p0, Scalar h0, Scalar /*s0*/, int /*meaning0*/,
                            Scalar& p, Scalar& h, Scalar& /*s*/, int& /*meaning*/)
    { limitPressureEnthalpyStep(p0, h0, p, h, Provider{}); }

    //! primary variables of the initial state from PRESSURE, TEMPI and SGAS
    static void initialVariables(Scalar p, Scalar T, Scalar Sg, Scalar& h, Scalar& s, int& meaning)
    {
        h = Provider::initialEnthalpy(p, T, Sg);
        s = 0.0;
        meaning = 0;
    }

    //! enthalpy of the injected stream (steam if gas) at pressure p and temperature T
    template <class Evaluation>
    static Evaluation injectionEnthalpy(bool gas, const Evaluation& p, const Evaluation& T)
    {
        return Provider::injectionEnthalpy(gas ? FluidSystem::gasPhaseIdx : FluidSystem::waterPhaseIdx, p, T);
    }

    //! the component of a source of water or steam
    static unsigned sourceComponentIndex(bool /*gas*/) { return 0; }

    // --- wells ---

    //! the fluid of the wellbore: always the tabulated provider (the raw IAPWS flash is too slow for the root searches)
    template <class S>
    using WellFluid = SteamWellFluid<S, TabulatedSteamPressureEnthalpyProvider<S>>;

    //! the stream a well injects: for the steam fluid whether it is steam
    template <class Fluid>
    static typename Fluid::Stream injectionStream(bool gas) { return gas; }

    //! equation (component) index of the component k of the wellbore fluid
    static unsigned wellComponentIndex(std::size_t /*k*/) { return 0; }
    static const char* wellComponentName(std::size_t /*k*/) { return "water"; }
    static constexpr const char* gasName = "steam";
    static constexpr const char* waterName = "water";

    static void checkDeck(const EclipseState& eclState)
    {
        const auto& simConfig = eclState.getSimulationConfig();
        if (!simConfig.isEnthalpy() || !simConfig.hasSTEAM()) {
            throw std::runtime_error("flow_ph with the steam flash requires the ENTHALPY and STEAM keywords in the "
                                     "RUNSPEC section (and no THERMAL or TEMP).");
        }
        if (eclState.runspec().phases().active(Phase::OIL)) {
            throw std::runtime_error("flow_ph with the steam flash does not support the OIL phase.");
        }
    }
};

template <class Scalar>
struct PhFlashSelector<PhFlashType::Steam, Scalar>
{ using type = SteamFlash<Scalar, SteamPressureEnthalpyProvider<Scalar>>; };

template <class Scalar>
struct PhFlashSelector<PhFlashType::SteamTabulated, Scalar>
{ using type = SteamFlash<Scalar, TabulatedSteamPressureEnthalpyProvider<Scalar>>; };

} // namespace Opm

#endif
