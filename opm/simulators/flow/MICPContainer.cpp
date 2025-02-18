// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/simulators/flow/MICPContainer.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/data/Solution.hpp>

#include <array>
#include <string>
#include <tuple>


namespace Opm {

template<class Scalar>
void MICPContainer<Scalar>::
allocate(const unsigned bufferSize)
{
    cMicrobes_.resize(bufferSize, 0.0);
    cOxygen_.resize(bufferSize, 0.0);
    cUrea_.resize(bufferSize, 0.0);
    cBiofilm_.resize(bufferSize, 0.0);
    cCalcite_.resize(bufferSize, 0.0);

    allocated_ = true;
}

template<class Scalar>
void MICPContainer<Scalar>::
outputRestart(data::Solution& sol)
{
    if (!this->allocated_) {
        return;
    }

    using DataEntry =
        std::tuple<std::string, UnitSystem::measure, std::vector<Scalar>&>;

    auto solutionVectors = std::array {
        DataEntry{"BIOFILM",  UnitSystem::measure::identity, cBiofilm_},
        DataEntry{"CALCITE",  UnitSystem::measure::identity, cCalcite_},
        DataEntry{"MICROBES", UnitSystem::measure::density,  cMicrobes_},
        DataEntry{"OXYGEN",   UnitSystem::measure::density,  cOxygen_},
        DataEntry{"UREA",     UnitSystem::measure::density,  cUrea_},
    };

    std::for_each(solutionVectors.begin(), solutionVectors.end(),
                  [&sol](auto& entry)
                  {
                      if (!std::get<2>(entry).empty()) {
                          sol.insert(std::get<std::string>(entry),
                          std::get<UnitSystem::measure>(entry),
                          std::move(std::get<2>(entry)),
                          data::TargetType::RESTART_OPM_EXTENDED);
                      }
                   });

    allocated_ = false;
}

template class MICPContainer<double>;

#if FLOW_INSTANTIATE_FLOAT
template class MICPContainer<float>;
#endif

} // namespace Opm
