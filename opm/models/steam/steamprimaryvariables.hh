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
 * \brief Primary variables of the steam model: (p, h).
 */
#ifndef OPM_STEAM_PRIMARY_VARIABLES_HH
#define OPM_STEAM_PRIMARY_VARIABLES_HH

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>

#include <ostream>

namespace Opm {

template <class TypeTag>
class SteamPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

public:
    SteamPrimaryVariables() : ParentType()
    { Opm::Valgrind::SetDefined(*this); }

    SteamPrimaryVariables(const SteamPrimaryVariables&) = default;
    SteamPrimaryVariables& operator=(const SteamPrimaryVariables&) = default;
    using ParentType::operator=;

    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState, const MaterialLawParams&, bool = false)
    { assignNaive(fluidState); }

    /*!
     * \brief Assign (p, h) from a fluid state carrying the liquid pressure,
     *        the saturations, densities and phase enthalpies.
     */
    template <class FluidState>
    void assignNaive(const FluidState& fs)
    {
        (*this) = 0.0;
        constexpr unsigned w = FluidSystem::waterPhaseIdx;
        constexpr unsigned g = FluidSystem::gasPhaseIdx;
        const Scalar mw = getValue(fs.saturation(w)) * getValue(fs.density(w));
        const Scalar mg = getValue(fs.saturation(g)) * getValue(fs.density(g));
        (*this)[Indices::pressure0Idx] = getValue(fs.pressure(w));
        (*this)[Indices::enthalpy0Idx] =
            (mw * getValue(fs.enthalpy(w)) + mg * getValue(fs.enthalpy(g))) / (mw + mg);
    }

    void print(std::ostream& os) const
    {
        os << "(p = " << (*this)[Indices::pressure0Idx]
           << ", h = " << (*this)[Indices::enthalpy0Idx] << ")" << std::flush;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        using FV = Dune::FieldVector<Scalar, getPropValue<TypeTag, Properties::NumEq>()>;
        serializer(static_cast<FV&>(*this));
    }
};

} // namespace Opm

#endif
