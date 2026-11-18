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
 * \brief Primary variable and equation indices of the steam model: pure water
 *        in two phases with pressure and specific enthalpy as primary variables.
 */
#ifndef OPM_STEAM_INDICES_HH
#define OPM_STEAM_INDICES_HH

namespace Opm {

template <class TypeTag, int PVOffset>
struct SteamIndices
{
    static constexpr bool waterEnabled = true;   // liquid water
    static constexpr bool gasEnabled = true;     // steam
    static constexpr bool oilEnabled = false;
    static constexpr int numPhases = 3;          // oil slot inactive

    static constexpr int numEq = 2;

    //! pressure of the liquid phase [Pa]
    static constexpr int pressure0Idx = PVOffset;
    //! specific enthalpy of the fluid mixture [J/kg]
    static constexpr int enthalpy0Idx = PVOffset + 1;
    //! under the name the flow linear solver uses
    static constexpr int pressureSwitchIdx = pressure0Idx;

    //! mass conservation of water
    static constexpr int conti0EqIdx = PVOffset;
    //! energy conservation
    static constexpr int contiEnergyEqIdx = PVOffset + 1;

    static constexpr bool enableMICP = false;
};

} // namespace Opm

#endif
