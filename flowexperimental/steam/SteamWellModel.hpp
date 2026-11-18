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
 * \brief Placeholder well model of the steam model. There are no wells yet:
 *        injection and production are specified with the SOURCE keyword. All
 *        well interfaces needed by the Flow stack are provided as no-ops.
 */
#ifndef OPM_STEAM_WELL_MODEL_HPP
#define OPM_STEAM_WELL_MODEL_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/output/data/Wells.hpp>

#include <opm/input/eclipse/Schedule/Events.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm {

class RegionVariableCollection;
struct NewtonIterationContext;

//! stands in for a well in the container the linear solver loops over; never instantiated
struct SteamWellPlaceholder
{
    std::vector<int> cells() const { return {}; }
    std::string name() const { return {}; }
    template <class X, class Y>
    void apply(const X&, Y&) const {}
};

template <class TypeTag>
class SteamWellModel : public BaseAuxiliaryModule<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    static const int numEq = Indices::numEq;
    using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
    using BVector = Dune::BlockVector<VectorBlockType>;
    using NeighborSet = std::set<unsigned>;

    SteamWellModel(Simulator& simulator, const NewtonIterationContext&)
        : simulator_(simulator)
    {}

    static void registerParameters() {}

    unsigned numDofs() const override { return 0; }
    void addNeighbors(std::vector<NeighborSet>&) const override {}
    void linearize(SparseMatrixAdapter&, GlobalEqVector&) override {}
    void applyInitial() override {}

    template <class Restarter> void serialize(Restarter&) {}
    template <class Restarter> void deserialize(Restarter&) {}
    template <class Serializer> void serializeOp(Serializer&) {}

    void beginEpisode() {}
    void beginReportStep(unsigned) {}
    void beginTimeStep() {}
    void beginIteration() {}
    void restoreLastValidState() {}
    void updateFailed() {}
    void advanceTimeLevel() {}
    void init() {}
    void endIteration() const {}
    void endTimeStep() {}
    void endEpisode() {}

    void computeTotalRatesForDof(RateVector&, unsigned) const {}

    [[nodiscard]] data::Wells wellData() const { return data::Wells{}; }
    [[nodiscard]] data::WellBlockAveragePressures wellBlockAveragePressures() const
    { return data::WellBlockAveragePressures{}; }
    [[nodiscard]] data::GroupAndNetworkValues groupAndNetworkData(const int&) const
    { return data::GroupAndNetworkValues{}; }
    [[nodiscard]] WellTestState wellTestState() const { return WellTestState{}; }

    void reportIntervalConnectionOilProduction(const double, const std::size_t,
                                               RegionVariableCollection&) const {}

    void recoverWellSolutionAndUpdateWellState(const BVector&) {}

    bool addMatrixContributions() const { return false; }
    const Schedule& schedule() const { return simulator_.vanguard().schedule(); }
    auto begin() const { return wells_.begin(); }
    auto end() const { return wells_.end(); }
    bool empty() const { return true; }

    int numLocalWellsEnd() const { return 0; }

    template <class PressureMatrix, class Vector>
    void addWellPressureEquations(PressureMatrix&, const Vector&, const bool, const bool) const {}
    template <class PressureMatrix>
    void addWellPressureEquationsStruct(PressureMatrix&) const {}
    template <class PressureMatrix, class Vector>
    void addWellPressureEquationsDomain(PressureMatrix&, const Vector&, const bool, const bool, const int) const
    { OPM_THROW(std::logic_error, "NLDD is not supported by the steam model"); }

    const std::vector<std::vector<int>>& well_local_cells() const
    { OPM_THROW(std::logic_error, "NLDD is not supported by the steam model"); }
    const std::map<std::string, int>& well_domain() const
    { OPM_THROW(std::logic_error, "NLDD is not supported by the steam model"); }

    const SimulatorReportSingle& lastReport() const { return last_report_; }
    void prepareDeserialize(const int) {}
    const std::map<std::string, double>& wellOpenTimes() const { return well_open_times_; }
    const std::map<std::string, double>& wellCloseTimes() const { return well_close_times_; }
    const WellGroupEvents& reportStepStartEvents() const { return report_step_start_events_; }
    bool forceShutWellByName(const std::string&, double, bool) { return false; }

    template <class ReservoirCouplingSlave> void setReservoirCouplingSlave(ReservoirCouplingSlave*) {}
    template <class ReservoirCouplingMaster> void setReservoirCouplingMaster(ReservoirCouplingMaster*) {}

    bool getWellConvergence() const { return true; }
    bool updateWellStatusAtConvergence() { return false; }

private:
    Simulator& simulator_;
    std::vector<std::shared_ptr<SteamWellPlaceholder>> wells_;
    SimulatorReportSingle last_report_{};
    std::map<std::string, double> well_open_times_;
    std::map<std::string, double> well_close_times_;
    WellGroupEvents report_step_start_events_{};
};

} // namespace Opm

#endif
