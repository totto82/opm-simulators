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
/*!
 * \file
 * \brief flow_steam with the tabulated steam provider instead of the raw IAPWS flash.
 */
#include "config.h"

#include <opm/material/constraintsolvers/TabulatedSteamPressureEnthalpyProvider.hpp>

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include "flow_steam.hpp"

namespace Opm::Properties {

namespace TTag {
struct FlowSteamTabProblem {
    using InheritsFrom = std::tuple<FlowSteamProblem>;
};
}

template<class TypeTag>
struct PressureEnthalpyProvider<TypeTag, TTag::FlowSteamTabProblem>
{ using type = TabulatedSteamPressureEnthalpyProvider<GetPropType<TypeTag, Properties::Scalar>>; };

} // namespace Opm::Properties

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowSteamTabProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    const auto ret = mainObject->runStatic<TypeTag>();
    mainObject.reset();
    return ret;
}
