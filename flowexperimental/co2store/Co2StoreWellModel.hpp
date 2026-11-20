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
 * \brief Well model of the CO2 storage model: wells with a wellbore (opm-common Opm::WellBore with
 *        Opm::Co2WellFluid) in which CO2 and brine flow and change phase between the perforations and the
 *        reference depth (CO2 coming out of solution, CO2 boiling, ...).
 *
 * It is the model of Opm::Co2StoreWellModel for brine and CO2, see there for the principle. Wells come from the Schedule
 * (WELSPECS, COMPDAT, WCONINJE, WCONPROD, WINJTEMP):
 *  - the reference depth of the well (WELSPECS) is the depth of the control. Put it at the surface (0): the
 *    wellbore then runs from the deepest connection to the surface and the pressure of a producer (the BHP item of
 *    WCONPROD, control mode BHP) is the pressure at the surface. The wellbore is vertical (measured depth = depth -
 *    reference depth), its diameter is twice the largest connection radius, and it is adiabatic with homogeneous flow.
 *  - an injector has control mode RESV: the reservoir volume rate is the volume of the injected stream at the
 *    injection temperature (WINJTEMP) and the pressure at the reference depth. Type GAS injects pure CO2, type WATER
 *    injects brine without CO2.
 *  - a producer produces the fluid of its cells with the mobilities of the three phases; the composition at the
 *    surface follows from the mixing along the wellbore.
 *
 * At the start of every Newton iteration the wellbore of each open well is solved from the current cell states
 * (warm-started from the previous solution). The connections enter the residuals of the components (brine, CO2) and of the
 * energy as sources with derivatives wrt the pressure and the enthalpy of the cell,
 *   producer: sum_a WI rho_a lambda_a (p_cell - p_well)^+  of mass, carrying the component fractions and h_a of the phases,
 *   injector: WI sum_a lambda_a (p_well - p_cell)^+ rho_well  of mass, carrying the composition and h_well of the stream,
 * with the well pressure following the cell pressure with the sensitivity of the well solution (see Co2StoreWellModel).
 *
 * The temperature and pressure ranges of the CO2 storage provider (280-400 K, 0.5-100 MPa) apply to the whole wellbore.
 * Not supported yet: parallel runs, group and rate controls, BHP limits of injectors, injected brine with dissolved
 * CO2, wells in the restart and summary files.
 */
#ifndef OPM_CO2STORE_WELL_MODEL_HPP
#define OPM_CO2STORE_WELL_MODEL_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/output/data/Wells.hpp>

#include <opm/material/wellbore/Co2WellFluid.hpp>
#include <opm/material/wellbore/WellBore.hpp>

#include <opm/input/eclipse/Schedule/Events.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/models/common/enthalpyenergyterms.hh>
#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

class RegionVariableCollection;
struct NewtonIterationContext;

//! stands in for a well in the container the linear solver loops over; never instantiated
struct Co2StoreWellPlaceholder
{
    std::vector<int> cells() const { return {}; }
    std::string name() const { return {}; }
    template <class X, class Y>
    void apply(const X&, Y&) const {}
};

template <class TypeTag>
class Co2StoreWellModel : public BaseAuxiliaryModule<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using WellFluid = Co2WellFluid<Scalar>;
    using WellBoreType = WellBore<Scalar, WellFluid>;
    using WellPerforation = WellBorePerforation<Scalar, 2>;
    using Stream = typename WellFluid::Stream;

    static constexpr unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;     // liquid CO2
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;     // gas CO2
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx; // brine
    static constexpr unsigned numPhases = FluidSystem::numPhases;
    static constexpr unsigned co2CompIdx = FluidSystem::CO2Idx;
    static constexpr unsigned brineCompIdx = FluidSystem::BrineIdx;

public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    static const int numEq = Indices::numEq;
    using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
    using BVector = Dune::BlockVector<VectorBlockType>;
    using NeighborSet = std::set<unsigned>;

    Co2StoreWellModel(Simulator& simulator, const NewtonIterationContext&)
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

    void beginEpisode() { buildWells_(); }
    void beginReportStep(unsigned) {}
    void beginTimeStep() {}
    void beginIteration()
    {
        for (auto& well : wellData_) {
            solveWell_(well);
        }
    }
    void restoreLastValidState() {}
    void updateFailed() {}
    void advanceTimeLevel() {}
    void init() {}
    void endIteration() const {}
    void endTimeStep() {}
    void endEpisode() { report_(); }

    void computeTotalRatesForDof(RateVector&, unsigned) const {}

    /*!
     * \brief Adds the connections of the wells in a cell to the source (per volume, in the units of
     *        the residual) with derivatives wrt the pressure and the enthalpy of the cell.
     */
    void addConnectionSources(RateVector& rate, unsigned cell, Scalar volume) const
    {
        const auto it = cellConnections_.find(cell);
        if (it == cellConnections_.end()) {
            return;
        }
        const auto& iq = simulator_.model().intensiveQuantities(cell, /*timeIdx=*/0);
        const auto& fs = iq.fluidState();
        const Evaluation& p = fs.pressure(waterPhaseIdx);
        constexpr Scalar energyScale = Parameters::EnthalpyEnergyScaling::value;

        for (const auto& [w, c] : it->second) {
            const auto& well = wellData_[w];
            const auto& conn = well.connections[c];
            // the pressure of the well fluid follows the cell pressure with the sensitivity of the well solution
            Evaluation pWell = conn.wellPressure + conn.sensitivity * (p - conn.cellPressure);
            // what the cell delivers (with its derivatives) decides the wellbore solution
            Evaluation coefficient = 0.0;
            if (well.injector) {
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    coefficient += conn.wellIndex * iq.mobility(phaseIdx);
                }
            }
            else {
                Evaluation hSum = 0.0, co2Sum = 0.0;
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    const Evaluation c_a = conn.wellIndex * fs.density(phaseIdx) * iq.mobility(phaseIdx);
                    coefficient += c_a;
                    hSum += c_a * fs.enthalpy(phaseIdx);
                    co2Sum += c_a * iq.massFraction(phaseIdx, co2CompIdx);
                }
                if (coefficient > 0.0) {
                    if (conn.enthalpySensitivity != 0.0) {
                        pWell += conn.enthalpySensitivity * (hSum / coefficient - conn.inflowEnthalpy);
                    }
                    if (conn.co2Sensitivity != 0.0) {
                        pWell += conn.co2Sensitivity * (co2Sum / coefficient - conn.inflowCo2);
                    }
                }
            }
            if (conn.coefficientSensitivity != 0.0) {
                pWell += conn.coefficientSensitivity * (coefficient - conn.coefficient);
            }
            if (well.injector) {
                const Evaluation dp = pWell - p;
                if (!(dp > 0.0)) {
                    continue;
                }
                Evaluation lambda = 0.0;
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    lambda += iq.mobility(phaseIdx);
                }
                const Evaluation mass = conn.wellIndex * lambda * dp * conn.wellDensity;
                // the stream is pure: its composition is that of the well (CO2 first, brine second)
                rate[Indices::conti0EqIdx + co2CompIdx] += mass * well.stream[0] / volume;
                rate[Indices::conti0EqIdx + brineCompIdx] += mass * well.stream[1] / volume;
                rate[Indices::contiEnergyEqIdx] += mass * conn.wellEnthalpy / volume * energyScale;
            }
            else {
                const Evaluation dp = p - pWell;
                if (!(dp > 0.0)) {
                    continue;
                }
                Evaluation co2 = 0.0;
                Evaluation brine = 0.0;
                Evaluation energy = 0.0;
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    const Evaluation c_a = conn.wellIndex * fs.density(phaseIdx) * iq.mobility(phaseIdx) * dp;
                    co2 += c_a * iq.massFraction(phaseIdx, co2CompIdx);
                    brine += c_a * iq.massFraction(phaseIdx, brineCompIdx);
                    energy += c_a * fs.enthalpy(phaseIdx);
                }
                rate[Indices::conti0EqIdx + co2CompIdx] -= co2 / volume;
                rate[Indices::conti0EqIdx + brineCompIdx] -= brine / volume;
                rate[Indices::contiEnergyEqIdx] -= energy / volume * energyScale;
            }
        }
    }

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
    { OPM_THROW(std::logic_error, "NLDD is not supported by the CO2 storage model"); }

    const std::vector<std::vector<int>>& well_local_cells() const
    { OPM_THROW(std::logic_error, "NLDD is not supported by the CO2 storage model"); }
    const std::map<std::string, int>& well_domain() const
    { OPM_THROW(std::logic_error, "NLDD is not supported by the CO2 storage model"); }

    const SimulatorReportSingle& lastReport() const { return last_report_; }
    void prepareDeserialize(const int) {}
    const std::map<std::string, double>& wellOpenTimes() const { return well_open_times_; }
    const std::map<std::string, double>& wellCloseTimes() const { return well_close_times_; }
    const WellGroupEvents& reportStepStartEvents() const { return report_step_start_events_; }
    bool forceShutWellByName(const std::string&, double, bool) { return false; }

    template <class ReservoirCouplingSlave> void setReservoirCouplingSlave(ReservoirCouplingSlave*) {}
    template <class ReservoirCouplingMaster> void setReservoirCouplingMaster(ReservoirCouplingMaster*) {}

    bool getWellConvergence() const
    {
        return std::none_of(wellData_.begin(), wellData_.end(), [](const auto& well) { return well.failed; });
    }
    bool updateWellStatusAtConvergence() { return false; }

private:
    struct ConnectionData
    {
        unsigned cell{0};
        Scalar md{0};                 //!< measured depth below the reference depth [m]
        Scalar wellIndex{0};          //!< [m3]
        Scalar radius{0};
        Scalar wellPressure{std::numeric_limits<Scalar>::max()};   //!< pressure of the well fluid at the connection [Pa]
        Scalar cellPressure{0};       //!< pressure of the cell when the well was solved [Pa]
        //! d(wellPressure)/d(cell pressure) of this connection: 1 for a well with a rate control (its pressure
        //! follows the cell), about 0 for a producer with a pressure control
        Scalar sensitivity{0};
        //! the inflow (producer) or outflow (injector) coefficient, the enthalpy and the CO2 mass fraction of what a producer draws
        //! when the well was solved, and the sensitivity of the pressure of the well fluid to them: the wellbore solution
        //! depends on what the cells deliver (the enthalpy and the composition decide how the fluid expands and how heavy the column is,
        //! the mobilities decide the rate)
        Scalar coefficient{0};
        Scalar coefficientSensitivity{0};
        Scalar inflowEnthalpy{0};
        Scalar enthalpySensitivity{0};
        Scalar inflowCo2{0};
        Scalar co2Sensitivity{0};
        Scalar wellDensity{0};        //!< injector: density of the well fluid at the connection
        Scalar wellEnthalpy{0};       //!< injector: enthalpy of the well fluid at the connection
        Scalar massRate{0};           //!< last rate through the connection [kg/s], positive in the flow direction
    };

    struct WellData
    {
        std::string name;
        bool injector{false};
        Stream stream{0.0, 1.0};      //!< injector: mass fractions {CO2, brine} of the stream: CO2 (type GAS) or brine (type WATER)
        Scalar referenceDepth{0};
        Scalar pressureTarget{0};     //!< producer: pressure at the reference depth [Pa]
        Scalar volumeRate{0};         //!< injector: volume rate at the reference depth [m3/s]
        Scalar temperature{0};        //!< injector: [K]
        Scalar diameter{0.2};
        std::vector<ConnectionData> connections;
        WellBoreResult<Scalar> result;
        bool flowing{false};
        bool rateReached{true};       //!< injector: false if the weight of the column alone injects more than the rate
        bool failed{false};           //!< the wellbore could not be solved in the last iteration
    };

    //! the wells of the current report step
    void buildWells_()
    {
        wellData_.clear();
        cellConnections_.clear();
        const auto& vanguard = simulator_.vanguard();
        const auto& schedule = vanguard.schedule();
        const auto& summaryState = vanguard.summaryState();
        const std::size_t reportStep = static_cast<std::size_t>(std::max(0, simulator_.episodeIndex()));

        for (const auto& well : schedule.getWells(reportStep)) {
            if (well.getStatus() != WellStatus::OPEN) {
                continue;
            }
            WellData wd;
            wd.name = well.name();
            wd.injector = well.isInjector();
            wd.referenceDepth = well.getRefDepth();
            Scalar maxRadius = 0.0;
            for (const auto& conn : well.getConnections()) {
                if (conn.state() != Connection::State::OPEN) {
                    continue;
                }
                const int cell = vanguard.compressedIndexForInterior(static_cast<int>(conn.global_index()));
                if (cell < 0) {
                    continue;   // the cell is not on this process
                }
                ConnectionData cd;
                cd.cell = static_cast<unsigned>(cell);
                const Scalar below = conn.depth() - wd.referenceDepth;
                if (below < -1.0) {
                    OPM_THROW(std::runtime_error, "Well " + wd.name + ": the connection at depth " +
                              std::to_string(conn.depth()) + " m is above the reference depth " +
                              std::to_string(wd.referenceDepth) + " m");
                }
                // a connection at the reference depth gets a wellbore of a millimetre
                cd.md = std::max(below, Scalar(1.0e-3));
                cd.wellIndex = conn.CF();
                cd.radius = conn.rw();
                maxRadius = std::max(maxRadius, cd.radius);
                wd.connections.push_back(cd);
            }
            if (wd.connections.empty()) {
                continue;
            }
            wd.diameter = maxRadius > 0.0 ? 2.0 * maxRadius : 0.2;

            if (wd.injector) {
                const auto controls = well.injectionControls(summaryState);
                if (controls.cmode != WellInjectorCMode::RESV) {
                    OPM_THROW(std::runtime_error, "Well " + wd.name + ": the CO2 storage model only supports RESV control "
                              "for injectors (the volume at the reference depth)");
                }
                if (controls.injector_type != InjectorType::WATER && controls.injector_type != InjectorType::GAS) {
                    OPM_THROW(std::runtime_error, "Well " + wd.name + ": the injector must be of type WATER (brine) or GAS (CO2)");
                }
                wd.stream = controls.injector_type == InjectorType::GAS ? Stream{1.0, 0.0} : Stream{0.0, 1.0};
                wd.volumeRate = controls.reservoir_rate;
                wd.temperature = well.inj_temperature();
            }
            else {
                const auto controls = well.productionControls(summaryState);
                if (controls.cmode != WellProducerCMode::BHP) {
                    OPM_THROW(std::runtime_error, "Well " + wd.name + ": the CO2 storage model only supports BHP control "
                              "for producers (the pressure at the reference depth)");
                }
                wd.pressureTarget = controls.bhp_limit;
            }
            wellData_.push_back(std::move(wd));
        }

        for (std::size_t w = 0; w < wellData_.size(); ++w) {
            for (std::size_t c = 0; c < wellData_[w].connections.size(); ++c) {
                cellConnections_[wellData_[w].connections[c].cell].emplace_back(w, c);
            }
            OpmLog::info("Well " + wellData_[w].name + (wellData_[w].injector ? " (injector of " +
                         std::string(wellData_[w].stream[0] > 0.5 ? "CO2" : "brine") + ")" : " (producer)") + ": " +
                         std::to_string(wellData_[w].connections.size()) + " connections, reference depth " +
                         std::to_string(wellData_[w].referenceDepth) + " m");
        }
    }

    //! solves the wellbore of a well from the current states of its cells
    void solveWell_(WellData& well)
    {
        std::vector<WellPerforation> perforations;
        Scalar mdMax = 0.0;
        for (const auto& conn : well.connections) {
            const auto& iq = simulator_.model().intensiveQuantities(conn.cell, /*timeIdx=*/0);
            const auto& fs = iq.fluidState();
            WellPerforation perf;
            perf.md = conn.md;
            perf.pressure = getValue(fs.pressure(waterPhaseIdx));
            Scalar productivity = 0.0, hSum = 0.0, co2Sum = 0.0, mobility = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                const Scalar lambda = getValue(iq.mobility(phaseIdx));
                const Scalar c = conn.wellIndex * getValue(fs.density(phaseIdx)) * lambda;
                productivity += c;
                hSum += c * getValue(fs.enthalpy(phaseIdx));
                co2Sum += c * getValue(iq.massFraction(phaseIdx, co2CompIdx));
                mobility += conn.wellIndex * lambda;
            }
            perf.productivity = productivity;
            perf.injectivity = mobility;
            perf.inflowEnthalpy = productivity > 0.0 ? hSum / productivity : getValue(fs.enthalpy(waterPhaseIdx));
            const Scalar zCo2 = productivity > 0.0 ? co2Sum / productivity : 0.0;
            perf.inflowComposition = {zCo2, 1.0 - zCo2};   // CO2 first, brine second
            perforations.push_back(perf);
            mdMax = std::max(mdMax, conn.md);
        }

        auto geometry = WellBoreGeometry<Scalar>::vertical(mdMax, well.diameter);
        geometry.maxStepLength = 20.0;
        // the solution of the previous iteration is the starting point of the root search of the wellbore
        const auto hintIt = pressureHints_.find(well.name);
        const std::optional<Scalar> hint = hintIt != pressureHints_.end() ? std::optional<Scalar>(hintIt->second) : std::nullopt;
        const auto solve = [&](const std::vector<WellPerforation>& perfs, std::optional<Scalar> pressureHint) {
            const WellBoreType wellbore(geometry, perfs);
            return well.injector ? wellbore.solveInjector(well.volumeRate, well.temperature, well.stream, pressureHint)
                                 : wellbore.solveProducer(well.pressureTarget, pressureHint);
        };
        const auto nodePressure = [](const WellBoreResult<Scalar>& r, Scalar md) {
            const auto it = std::find_if(r.profile.begin(), r.profile.end(), [&](const auto& pt) {
                return std::abs(pt.md - md) < 1e-6 * (1.0 + md);
            });
            assert(it != r.profile.end());
            return it->pressure;
        };
        well.failed = false;
        try {
            well.result = solve(perforations, hint);
        }
        catch (const std::exception& e) {
            // The states of a Newton iteration can be outside what the wellbore can handle (for example a cell
            // pressure beyond the range of the fluid). The well is switched off for this iteration and
            // the iteration cannot converge while it is failed (getWellConvergence()).
            well.failed = true;
            pressureHints_.erase(well.name);
            std::string message = "Well " + well.name + ": the wellbore could not be solved (" + e.what() + "), cell pressures";
            for (const auto& perf : perforations) {
                message += " " + std::to_string(perf.pressure / 1e5) + " bar";
            }
            OpmLog::debug(message);
            well.flowing = false;
            well.rateReached = true;
            for (auto& conn : well.connections) {
                conn.wellPressure = well.injector ? 0.0 : std::numeric_limits<Scalar>::max();
                conn.wellDensity = 0.0;
                conn.wellEnthalpy = 0.0;
                conn.sensitivity = 0.0;
                conn.coefficientSensitivity = 0.0;
                conn.enthalpySensitivity = 0.0;
                conn.co2Sensitivity = 0.0;
                conn.cellPressure = 0.0;
                conn.massRate = 0.0;
            }
            return;
        }
        well.flowing = well.result.flowing;
        well.rateReached = well.result.rateReached;
        if (well.flowing && well.rateReached) {
            pressureHints_[well.name] = well.injector ? well.result.surfacePressure : well.result.bottomHolePressure;
        }
        else {
            pressureHints_.erase(well.name);
        }

        for (std::size_t c = 0; c < well.connections.size(); ++c) {
            auto& conn = well.connections[c];
            const auto it = std::find_if(well.result.profile.begin(), well.result.profile.end(), [&](const auto& pt) {
                return std::abs(pt.md - conn.md) < 1e-6 * (1.0 + conn.md);
            });
            assert(it != well.result.profile.end());
            conn.cellPressure = perforations[c].pressure;
            conn.sensitivity = 0.0;
            conn.coefficient = well.injector ? perforations[c].injectivity : perforations[c].productivity;
            conn.coefficientSensitivity = 0.0;
            conn.inflowEnthalpy = perforations[c].inflowEnthalpy;
            conn.enthalpySensitivity = 0.0;
            conn.inflowCo2 = perforations[c].inflowComposition[0];
            conn.co2Sensitivity = 0.0;
            if (!well.flowing) {
                conn.wellPressure = well.injector ? 0.0 : std::numeric_limits<Scalar>::max();
                conn.wellDensity = 0.0;
                conn.wellEnthalpy = 0.0;
                continue;
            }
            conn.wellPressure = it->pressure;
            conn.wellDensity = it->density;
            conn.wellEnthalpy = it->enthalpy;

            // sensitivity of the pressure of the well fluid at this connection to its own cell pressure
            // (the cross terms between the cells of a well are neglected)
            const Scalar dp = 1.0e-5 * perforations[c].pressure;
            auto perturbed = perforations;
            perturbed[c].pressure += dp;
            try {
                const auto r = solve(perturbed, well.injector ? well.result.surfacePressure : well.result.bottomHolePressure);
                if (r.flowing) {
                    conn.sensitivity = (nodePressure(r, conn.md) - conn.wellPressure) / dp;
                }
            }
            catch (const std::exception&) {
                conn.sensitivity = 0.0;
            }
            // sensitivities to the other descriptors of the connection, each from the wellbore solution with a perturbed perforation
            const Scalar hint2 = well.injector ? well.result.surfacePressure : well.result.bottomHolePressure;
            const auto sensitivityTo = [&](auto&& perturb, Scalar change, Scalar& result) {
                auto perturbed = perforations;
                perturb(perturbed[c], change);
                try {
                    const auto r = solve(perturbed, hint2);
                    if (r.flowing) {
                        result = (nodePressure(r, conn.md) - conn.wellPressure) / change;
                    }
                }
                catch (const std::exception&) {
                    result = 0.0;
                }
            };
            {
                // a relative change of 0.1% of the coefficient
                const Scalar change = 1.0e-3 * conn.coefficient;
                if (change > 0.0) {
                    sensitivityTo([&](WellPerforation& perf, Scalar d) { (well.injector ? perf.injectivity : perf.productivity) += d; },
                                  change, conn.coefficientSensitivity);
                }
            }
            if (!well.injector) {
                sensitivityTo([](WellPerforation& perf, Scalar d) { perf.inflowEnthalpy += d; }, 1.0e3, conn.enthalpySensitivity);
                sensitivityTo([](WellPerforation& perf, Scalar d) {
                    perf.inflowComposition[0] += d;
                    perf.inflowComposition[1] -= d;
                }, perforations[c].inflowComposition[0] < 0.99 ? 0.01 : -0.01, conn.co2Sensitivity);
            }
        }
        for (std::size_t c = 0; c < well.connections.size(); ++c) {
            // perforation rates are in the order of the perforations, which is that of the connections
            well.connections[c].massRate = well.result.perforationRate[c];
        }
    }

    //! one line per well at the end of a report step
    void report_() const
    {
        for (const auto& well : wellData_) {
            const auto& r = well.result;
            std::string line = "Well " + well.name + ": " + (well.injector ? "injecting " : "producing ") +
                std::to_string(r.massFlow) + " kg/s, pressure at the reference depth " +
                std::to_string(r.surfacePressure / 1e5) + " bar, bottom-hole pressure " +
                std::to_string(r.bottomHolePressure / 1e5) + " bar";
            if (r.flowing && !r.profile.empty()) {
                // composition at the reference depth: CO2 fraction of the flow and the rates of the components
                const Scalar zCo2 = well.injector ? well.stream[0] : r.profile.front().composition[0];
                line += ", CO2 " + std::to_string(r.massFlow * zCo2) + " kg/s, brine " + std::to_string(r.massFlow * (1.0 - zCo2)) + " kg/s";
                if (!well.injector) {
                    line += ", CO2-rich phase " + std::to_string(r.profile.front().steamQuality) + " of the mass at the reference depth";
                }
            }
            if (!well.rateReached) {
                line += ", the weight of the fluid column alone injects more than the specified rate";
            }
            if (r.phaseChangeMd) {
                line += ", phase change at " + std::to_string(*r.phaseChangeMd + well.referenceDepth) + " m depth";
            }
            OpmLog::info(line);
        }
    }

    Simulator& simulator_;
    std::vector<WellData> wellData_;
    std::map<unsigned, std::vector<std::pair<std::size_t, std::size_t>>> cellConnections_;
    //! last solution of every well: BHP of a producer, surface pressure of an injector
    std::map<std::string, Scalar> pressureHints_;
    std::vector<std::shared_ptr<Co2StoreWellPlaceholder>> wells_;
    SimulatorReportSingle last_report_{};
    std::map<std::string, double> well_open_times_;
    std::map<std::string, double> well_close_times_;
    WellGroupEvents report_step_start_events_{};
};

} // namespace Opm

#endif
