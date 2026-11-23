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
 * \brief The flash of the generic compositional pressure-enthalpy model (Opm::PhModel).
 *
 * The model has the primary variables (p, h, s_1, ..., s_{N-1}) and the equations of the N components and of
 * the energy. Everything that depends on the fluid is behind the "flash" class selected by the property
 * PhFlashKind:
 *
 *   PhFlashType::Steam, SteamTabulated  pure water (N = 1), IAPWS-IF97 raw or tabulated
 *   PhFlashType::Co2Store               brine and CO2 (N = 2), the third variable is the switched variable
 *   PhFlashType::Cubic                  (later) the general cubic EOS flash
 *
 * A flash class F provides (see SteamFlash and Co2StoreFlash):
 *   F::FluidSystem, F::numComponents, F::numSwitchVariables (= N - 1)
 *   F::minPressure(), maxPressure(), minEnthalpy(p), maxEnthalpy(p), maxEnthalpyStep(p)
 *   F::solve(p, h, s, meaning)    state with saturation(phase), density(phase), enthalpy(phase), viscosity(phase),
 *                                 temperature(), regionId() and, if N > 1, massFraction(phase, comp)
 *   F::limitUpdate(...)           the Newton update limiter (range, size, phase boundaries, variable switching)
 *   F::initialVariables(...)      primary variables of the initial state of the deck
 *   F::injectionEnthalpy(...), F::sourceComponentIndex(...)   sources
 *   F::checkDeck(...)             what the deck must have for this flash
 */
#ifndef OPM_PH_FLASH_HH
#define OPM_PH_FLASH_HH

#include <opm/models/utils/propertysystem.hh>

namespace Opm {

enum class PhFlashType { Steam, SteamTabulated, Co2Store, Cubic };

//! The flash class of a flash type; the specialisations are in steamflash.hh and co2storeflash.hh
template <PhFlashType kind, class Scalar>
struct PhFlashSelector;

} // namespace Opm

namespace Opm::Properties {

//! The kind of flash (a value of Opm::PhFlashType)
template <class TypeTag, class MyTypeTag>
struct PhFlashKind { using type = UndefinedProperty; };

//! The flash class (the default follows from PhFlashKind)
template <class TypeTag, class MyTypeTag>
struct PhFlash { using type = UndefinedProperty; };

} // namespace Opm::Properties

#endif
