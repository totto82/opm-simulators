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
 * \brief Newton method of the generic pressure-enthalpy model: the update of the primary variables is limited by
 *        the flash (range, size, cut at phase boundaries and, for flashes with a switched variable, the switch of its
 *        meaning), see Opm::SteamFlash::limitUpdate and Opm::Co2StoreFlash::limitUpdate.
 */
#ifndef OPM_PH_NEWTON_METHOD_HH
#define OPM_PH_NEWTON_METHOD_HH

#include <opm/models/nonlinear/newtonmethod.hh>
#include <opm/models/ph/phflash.hh>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

template <class TypeTag>
class PhNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Flash = GetPropType<TypeTag, Properties::PhFlash>;

public:
    explicit PhNewtonMethod(Simulator& simulator) : ParentType(simulator) {}

protected:
    friend ParentType;
    friend NewtonMethod<TypeTag>;

    void updatePrimaryVariables_(unsigned /* globalDofIdx */,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& /* currentResidual */)
    {
        const PrimaryVariables old = currentValue;
        nextValue = old;
        nextValue -= update;

        Scalar p = nextValue[Indices::pressure0Idx];
        Scalar h = nextValue[Indices::enthalpy0Idx];
        if constexpr (Flash::numSwitchVariables == 0) {
            Scalar s = 0.0;
            int meaning = old.meaning();
            Flash::limitUpdate(old[Indices::pressure0Idx], old[Indices::enthalpy0Idx], 0.0, meaning, p, h, s, meaning);
        }
        else {
            Scalar s = nextValue[Indices::switchIdx];
            int meaning = old.meaning();
            Flash::limitUpdate(old[Indices::pressure0Idx], old[Indices::enthalpy0Idx], old[Indices::switchIdx], old.meaning(),
                               p, h, s, meaning);
            nextValue[Indices::switchIdx] = s;
            nextValue.setMeaning(meaning);
        }
        nextValue[Indices::pressure0Idx] = p;
        nextValue[Indices::enthalpy0Idx] = h;
    }
};

} // namespace Opm

#endif
