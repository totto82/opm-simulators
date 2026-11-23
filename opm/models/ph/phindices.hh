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
 * rief Indices of the generic compositional pressure-enthalpy model: N components in three phase slots
 *        (aqueous, liquid, gas) with pressure, specific enthalpy and N-1 switched variables as primary variables.
 */
#ifndef OPM_PH_INDICES_HH
#define OPM_PH_INDICES_HH

namespace Opm {

//! Flash: the flash class (see phflash.hh)
template <class TypeTag, class Flash, int PVOffset>
struct PhIndices
{
    static constexpr bool waterEnabled = true;
    static constexpr bool gasEnabled = true;
    //! the oil slot is the liquid CO2 of the CO2 storage flash; the steam flash has no such phase
    static constexpr bool oilEnabled = Flash::oilPhaseEnabled;
    static constexpr int numPhases = 3;
    static constexpr int numComponents = Flash::numComponents;

    static constexpr int numEq = numComponents + 1;

    //! pressure of the aqueous (liquid) phase [Pa]
    static constexpr int pressure0Idx = PVOffset;
    //! specific enthalpy of the fluid mixture [J/kg]
    static constexpr int enthalpy0Idx = PVOffset + 1;
    //! first switched variable: the meaning is stored with the primary variables, see PhPrimaryVariables
    static constexpr int switchIdx = PVOffset + 2;
    //! under the name the flow linear solver uses
    static constexpr int pressureSwitchIdx = pressure0Idx;

    //! mass conservation of the components: conti0EqIdx + compIdx with the component indices of the fluid system
    static constexpr int conti0EqIdx = PVOffset;
    //! energy conservation
    static constexpr int contiEnergyEqIdx = PVOffset + numComponents;

    static constexpr bool enableMICP = false;
};

} // namespace Opm

#endif
