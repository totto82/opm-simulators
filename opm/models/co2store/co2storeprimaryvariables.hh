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
 * \brief Primary variables of the CO2 storage model: pressure, specific enthalpy and a third
 *        variable whose meaning is switched, see Opm::Co2BrinePressureEnthalpyProvider.
 */
#ifndef OPM_CO2STORE_PRIMARY_VARIABLES_HH
#define OPM_CO2STORE_PRIMARY_VARIABLES_HH

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>

#include <ostream>
#include <stdexcept>

namespace Opm {

template <class TypeTag>
class Co2StorePrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

public:
    Co2StorePrimaryVariables() : ParentType()
    { Opm::Valgrind::SetDefined(*this); }

    Co2StorePrimaryVariables(const Co2StorePrimaryVariables&) = default;
    Co2StorePrimaryVariables& operator=(const Co2StorePrimaryVariables&) = default;
    using ParentType::operator=;

    //! meaning of the third variable, one of Provider::bothPhases, aqueousOnly, co2Only
    int meaning() const { return meaning_; }
    void setMeaning(int meaning) { meaning_ = meaning; }

    template <class FluidState>
    void assignMassConservative(const FluidState&, const MaterialLawParams&, bool = false)
    { throw std::logic_error("The CO2 storage model is initialised explicitly (PRESSURE, TEMPI, SWAT, ...)"); }

    template <class FluidState>
    void assignNaive(const FluidState&)
    { throw std::logic_error("The CO2 storage model is initialised explicitly (PRESSURE, TEMPI, SWAT, ...)"); }

    void print(std::ostream& os) const
    {
        os << "(p = " << (*this)[Indices::pressure0Idx]
           << ", h = " << (*this)[Indices::enthalpy0Idx]
           << ", s = " << (*this)[Indices::switchIdx]
           << ", meaning " << meaning_ << ")" << std::flush;
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
