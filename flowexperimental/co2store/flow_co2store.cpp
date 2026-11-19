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
#include "config.h"

#include <opm/simulators/flow/FlowGenericProblem_impl.hpp>

#include "flow_co2store.hpp"

int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::FlowCo2StoreProblem;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    const auto ret = mainObject->runStatic<TypeTag>();
    mainObject.reset();
    return ret;
}
