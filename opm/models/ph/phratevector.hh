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
 * \brief Rate vectors of the ph model. Units: mass [kg/(s m3)] and
 *        energy [J/(s m3)].
 */
#ifndef OPM_PH_RATE_VECTOR_HH
#define OPM_PH_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <opm/material/common/Valgrind.hpp>

namespace Opm {

template <class TypeTag>
class PhRateVector
    : public Dune::FieldVector<GetPropType<TypeTag, Properties::Evaluation>,
                               getPropValue<TypeTag, Properties::NumEq>()>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using ParentType = Dune::FieldVector<Evaluation, numEq>;

public:
    PhRateVector() : ParentType() { Valgrind::SetUndefined(*this); }
    explicit PhRateVector(const Evaluation& value) : ParentType(value) {}
    PhRateVector(const PhRateVector&) = default;
    PhRateVector& operator=(const PhRateVector&) = default;
    using ParentType::operator=;
};

template <class TypeTag>
class PhBoundaryRateVector : public PhRateVector<TypeTag>
{
    using ParentType = PhRateVector<TypeTag>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    PhBoundaryRateVector() : ParentType() {}
    explicit PhBoundaryRateVector(const Evaluation& value) : ParentType(value) {}
    PhBoundaryRateVector(const PhBoundaryRateVector&) = default;
    PhBoundaryRateVector& operator=(const PhBoundaryRateVector&) = default;

    void setNoFlow()
    {
        for (unsigned i = 0; i < (*this).size(); ++i) {
            (*this)[i] = Evaluation(0.0);
        }
    }
};

} // namespace Opm

#endif
