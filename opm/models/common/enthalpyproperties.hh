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
 * \brief Property of models that have specific enthalpy as a primary variable.
 */
#ifndef OPM_ENTHALPY_PROPERTIES_HH
#define OPM_ENTHALPY_PROPERTIES_HH

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Properties {

//! The pressure-enthalpy provider (see Opm::PressureEnthalpyProvider in opm-common):
//! the flash, the validity range and the initialisation of the fluid.
template <class TypeTag, class MyTypeTag>
struct PressureEnthalpyProvider { using type = UndefinedProperty; };

} // namespace Opm::Properties

#endif
