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
 * \brief Primary variables of the generic pressure-enthalpy model: pressure, specific enthalpy and, if the flash has
 *        components beyond the first, a switched variable whose meaning is stored here (see Opm::Co2StoreFlash).
 */
#ifndef OPM_PH_PRIMARY_VARIABLES_HH
#define OPM_PH_PRIMARY_VARIABLES_HH

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/ph/phflash.hh>

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>

#include <ostream>
#include <stdexcept>

namespace Opm {

template <class TypeTag>
class PhPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Flash = GetPropType<TypeTag, Properties::PhFlash>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

public:
    PhPrimaryVariables() : ParentType()
    { Opm::Valgrind::SetDefined(*this); }

    PhPrimaryVariables(const PhPrimaryVariables&) = default;
    PhPrimaryVariables& operator=(const PhPrimaryVariables&) = default;
    using ParentType::operator=;

    //! meaning of the switched variable (the values are those of the flash; 0 if there is none)
    int meaning() const { return meaning_; }
    void setMeaning(int meaning) { meaning_ = meaning; }

    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState, const MaterialLawParams&, bool = false)
    { assignNaive(fluidState); }

    /*!
     * \brief Assign (p, h) from a fluid state carrying the liquid pressure, the saturations, densities and phase
     *        enthalpies. Only for one component: the other flashes are initialised explicitly (PRESSURE, TEMPI, ...).
     */
    template <class FluidState>
    void assignNaive(const FluidState& fs)
    {
        if constexpr (Flash::numComponents == 1) {
            (*this) = 0.0;
            constexpr unsigned w = FluidSystem::waterPhaseIdx;
            constexpr unsigned g = FluidSystem::gasPhaseIdx;
            const Scalar mw = getValue(fs.saturation(w)) * getValue(fs.density(w));
            const Scalar mg = getValue(fs.saturation(g)) * getValue(fs.density(g));
            (*this)[Indices::pressure0Idx] = getValue(fs.pressure(w));
            (*this)[Indices::enthalpy0Idx] =
                (mw * getValue(fs.enthalpy(w)) + mg * getValue(fs.enthalpy(g))) / (mw + mg);
        }
        else {
            throw std::logic_error("This pressure-enthalpy flash is initialised explicitly (PRESSURE, TEMPI, SWAT, ...)");
        }
    }

    void print(std::ostream& os) const
    {
        os << "(p = " << (*this)[Indices::pressure0Idx]
           << ", h = " << (*this)[Indices::enthalpy0Idx];
        if constexpr (Flash::numSwitchVariables > 0) {
            os << ", s = " << (*this)[Indices::switchIdx]
               << ", meaning " << meaning_;
        }
        os << ")" << std::flush;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        using FV = Dune::FieldVector<Scalar, getPropValue<TypeTag, Properties::NumEq>()>;
        serializer(static_cast<FV&>(*this));
        serializer(meaning_);
    }

private:
    int meaning_{0};
};

} // namespace Opm

#endif
