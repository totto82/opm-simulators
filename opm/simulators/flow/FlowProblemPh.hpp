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
 * \copydoc Opm::FlowProblemPh
 */
#ifndef OPM_FLOW_PROBLEM_PH_HPP
#define OPM_FLOW_PROBLEM_PH_HPP

#include <opm/models/common/enthalpyenergyterms.hh>
#include <opm/models/ph/phflash.hh>
#include <opm/models/ph/phlocalresidual.hh>
#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/OutputPhModule.hpp>
#include <opm/simulators/flow/NonlinearSystemPh.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>
#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/input/eclipse/Schedule/Source.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Opm {

template <class TypeTag>
class FlowProblemPh;

}

namespace Opm::Properties {

namespace TTag {

struct FlowBaseProblemPh {
    using InheritsFrom = std::tuple<FlowBaseProblem>;
};

}

template<class TypeTag>
struct NonlinearSystem<TypeTag, TTag::FlowBaseProblemPh>
{ using type = NonlinearSystemPh<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FlowBaseProblemPh>
{ using type = FlowProblemPh<TypeTag>; };

template<class TypeTag>
struct TracerModel<TypeTag, TTag::FlowBaseProblemPh>
{ using type = ::Opm::TracerModel<TypeTag>; };

// the energy equation is always on: the thermal laws of the problem are used
template<class TypeTag>
struct EnergyModuleType<TypeTag, TTag::FlowBaseProblemPh>
{ static constexpr EnergyModules value = EnergyModules::FullyImplicitThermal; };

template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::FlowBaseProblemPh>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using Traits = ThreePhaseMaterialTraits<Scalar,
                                            /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                            /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                            /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx,
                                            /* hysteresis */ true,
                                            /* endpointscaling */ true>;

public:
    using EclMaterialLawManager = ::Opm::EclMaterialLaw::Manager<Traits>;
    using type = typename EclMaterialLawManager::MaterialLaw;
};

template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::FlowBaseProblemPh>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;
    using type = typename EclThermalLawManager::SolidEnergyLaw;
};

template<class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::FlowBaseProblemPh>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;
    using type = typename EclThermalLawManager::ThermalConductionLaw;
};

template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::FlowBaseProblemPh>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableConvectiveMixing<TypeTag, TTag::FlowBaseProblemPh>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::FlowBaseProblemPh>
{ static constexpr bool value = false; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup PhSimulator
 *
 * \brief Flow problem of the generic pressure-enthalpy model (Opm::PhModel): the components of the flash
 *        (pure water, or brine and CO2, ...) with pressure, specific enthalpy and the switched variables of the
 *        flash as primary variables.
 *
 * The model is initialised explicitly (PRESSURE, TEMPI and the saturation of the flash: SGAS for water, SWAT for
 * CO2 storage) and driven by the SOURCE keyword and by wells. What depends on the fluid (deck requirements, initial
 * variables, source components and enthalpies) comes from the flash class, see phflash.hh.
 */
template <class TypeTag>
class FlowProblemPh : public FlowProblem<TypeTag>
{
    using FlowProblemType = FlowProblem<TypeTag>;

    using typename FlowProblemType::Scalar;
    using typename FlowProblemType::Simulator;
    using typename FlowProblemType::GridView;
    using typename FlowProblemType::FluidSystem;
    using typename FlowProblemType::Vanguard;

    using FlowProblemType::dim;
    using FlowProblemType::dimWorld;
    using FlowProblemType::numPhases;

    using FlowProblemType::gasPhaseIdx;
    using FlowProblemType::oilPhaseIdx;
    using FlowProblemType::waterPhaseIdx;

    using typename FlowProblemType::Indices;
    using typename FlowProblemType::PrimaryVariables;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using typename FlowProblemType::Evaluation;
    using typename FlowProblemType::IntensiveQuantities;
    using typename FlowProblemType::MaterialLaw;
    using typename FlowProblemType::RateVector;

    using Flash = GetPropType<TypeTag, Properties::PhFlash>;
    using EclWriterType = EclWriter<TypeTag, OutputPhModule<TypeTag>>;

    // pressure, temperature and saturations of the initial state of a cell
    using InitialFluidState = SimpleModularFluidState<Scalar, numPhases, /*numComponents=*/Flash::numComponents, FluidSystem,
                                                      /*storePressure=*/true, /*storeTemperature=*/true,
                                                      /*storeComposition=*/false, /*storeFugacity=*/false,
                                                      /*storeSaturation=*/true, /*storeDensity=*/false,
                                                      /*storeViscosity=*/false, /*storeEnthalpy=*/false>;

public:
    using FlowProblemType::porosity;
    using FlowProblemType::pvtRegionIndex;

    static void registerParameters()
    {
        FlowProblemType::registerParameters();
        Parameters::SetDefault<Parameters::InputSkipMode>("300");
        // the mass and energy residuals are in kg/(m3 s) and MJ/(m3 s): the flow default
        // of 1e-2 would accept errors of the order of the storage term itself
        Parameters::SetDefault<Parameters::NewtonTolerance<Scalar>>(1e-7);

        EclWriterType::registerParameters();
        GetPropType<TypeTag, Properties::WellModel>::registerParameters();
    }

    explicit FlowProblemPh(Simulator& simulator)
        : FlowProblemType(simulator)
    {
        eclWriter_ = std::make_unique<EclWriterType>(simulator);
        enableEclOutput_ = Parameters::Get<Parameters::EnableEclOutput>();
    }

    void finishInit()
    {
        Flash::checkDeck(this->simulator().vanguard().eclState());
        FlowProblemType::finishInit();

        auto& simulator = this->simulator();
        const bool transmissibilitiesFinished =
            this->prepareTransmissibilityOutput_(*eclWriter_, enableEclOutput_);
        if (!transmissibilitiesFinished) {
            this->finishTransmissibilities_();
        }

        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();
        this->initializeSimulatorTime_();

        this->initFluidSystem_();
        if (eclState.getSimulationConfig().rock_config().active()) {
            throw std::runtime_error("Rock compaction (ROCKCOMP) is not supported by the pressure-enthalpy model");
        }
        this->initializeModelProperties_();

        if (enableEclOutput_) {
            eclWriter_->writeInit();
        }

        const auto& initconfig = eclState.getInitConfig();
        if (initconfig.restartRequested()) {
            throw std::logic_error("Restarting is not supported by the pressure-enthalpy model yet");
        }
        this->readInitialCondition_();

        FlowProblemType::updatePffDofData_();

        if (this->enableVtkOutput_() && eclState.getIOConfig().initOnly()) {
            simulator.setTimeStepSize(0.0);
            FlowProblemType::writeOutput(true);
        }

        simulator.startNextEpisode(schedule.seconds(1));
        simulator.setEpisodeIndex(0);
        simulator.setTimeStepIndex(0);
    }

    void initialSolutionApplied() override
    {
        FlowProblemType::initialSolutionApplied();
        this->thresholdPressures_.finishInit();
    }

    void endTimeStep() override
    {
        FlowProblemType::endTimeStep();

        this->eclWriter_->mutableOutputModule().invalidateLocalData();

        const auto& grid = this->simulator().vanguard().gridView().grid();
        using GridType = std::remove_cv_t<std::remove_reference_t<decltype(grid)>>;
        constexpr bool isCpGrid = std::is_same_v<GridType, Dune::CpGrid>;
        if (!isCpGrid || (grid.maxLevel() == 0)) {
            this->eclWriter_->evalSummaryState(! this->episodeWillBeOver());
        }
    }

    void writeReports(const SimulatorTimer& timer)
    {
        if (enableEclOutput_) {
            eclWriter_->writeReports(timer);
        }
    }

    void writeOutput(bool verbose) override
    {
        FlowProblemType::writeOutput(verbose);

        if (! this->enableEclOutput_) {
            return;
        }

        const auto isSubStep = !this->episodeWillBeOver();

        if (!isSubStep || Parameters::Get<Parameters::EnableWriteAllSolutions>()) {
            auto localCellData = data::Solution {};

            this->eclWriter_->writeOutput(std::move(localCellData), isSubStep,
                                          this->simulator().vanguard().schedule()
                                         .exitStatus().has_value());
        }
    }

    //! no-flow on all boundaries
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned /* timeIdx */) const
    {
        if (!context.intersection(spaceIdx).boundary())
            return;

        values.setNoFlow();

        if (this->nonTrivialBoundaryConditions()) {
            throw std::logic_error("boundary conditions are not supported by the pressure-enthalpy model yet");
        }
    }

    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        values[Indices::pressure0Idx] = initialPressure_[globalDofIdx];
        values[Indices::enthalpy0Idx] = initialEnthalpy_[globalDofIdx];
        if constexpr (Flash::numSwitchVariables > 0) {
            values[Indices::switchIdx] = initialSwitch_[globalDofIdx];
            values.setMeaning(initialMeaning_[globalDofIdx]);
        }
    }

    /*!
     * \brief Mass and energy sources from the SOURCE keyword.
     *
     * A positive GAS rate injects the gas component (steam or CO2) and a positive WATER rate injects the water
     * component (liquid water or brine), at the given temperature (default: the cell temperature; steam is not
     * injected below, and liquid not above, the saturation temperature). A negative rate removes the component at its
     * mean enthalpy in the cell. An explicit HRATE overrides the energy rate.
     */
    void addToSourceDense(RateVector& rate, unsigned globalDofIdx, unsigned /*timeIdx*/) const override
    {
        // the connections of the wells
        this->wellModel().addConnectionSources(rate, globalDofIdx, this->model().dofTotalVolume(globalDofIdx));

        const auto& source = this->simulator().vanguard().schedule()[this->episodeIndex()].source();
        std::array<int, 3> ijk;
        this->simulator().vanguard().cartesianCoordinate(globalDofIdx, ijk);
        if (!source.hasSource(ijk)) {
            return;
        }

        const Scalar volume = this->model().dofTotalVolume(globalDofIdx);
        const auto& iq = this->simulator().model().intensiveQuantities(globalDofIdx, /*timeIdx=*/0);
        const auto& fs = iq.fluidState();
        const Evaluation p = fs.pressure(waterPhaseIdx);

        constexpr Scalar energyScale = Parameters::EnthalpyEnergyScaling::value;

        for (const auto& [comp, isGas] : {std::pair{SourceComponent::WATER, false},
                                          std::pair{SourceComponent::GAS, true}}) {
            const unsigned compIdx = Flash::sourceComponentIndex(isGas);
            const Scalar massRate = source.rate(ijk, comp) / volume;
            if (massRate == 0.0) {
                continue;
            }
            rate[Indices::conti0EqIdx + compIdx] += massRate;

            Evaluation energyRate;
            if (const auto hrate = source.hrate(ijk, comp)) {
                energyRate = hrate.value() / volume;
            }
            else if (massRate < 0.0) {
                // mean enthalpy of the component in the cell: mass of the component in each phase times the phase enthalpy
                Evaluation mass = 0.0, energy = 0.0;
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }
                    const Evaluation m = fs.saturation(phaseIdx) * fs.density(phaseIdx) * iq.massFraction(phaseIdx, compIdx);
                    mass += m;
                    energy += m * fs.enthalpy(phaseIdx);
                }
                energyRate = massRate * energy / mass;
            }
            else {
                const Scalar T = source.temperature(ijk, comp).value_or(getValue(fs.temperature(0)));
                energyRate = massRate * Flash::injectionEnthalpy(isGas, p, Evaluation(T));
            }
            rate[Indices::contiEnergyEqIdx] += energyRate * energyScale;
        }
        if (source.rate(ijk, SourceComponent::OIL) != 0.0) {
            throw std::logic_error("The OIL source component is not supported by the pressure-enthalpy model: use GAS and WATER");
        }
    }

    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }
    std::vector<InitialFluidState>& initialFluidStates()
    { return initialFluidStates_; }
    const std::vector<InitialFluidState>& initialFluidStates() const
    { return initialFluidStates_; }

    const EclWriterType& eclWriter() const { return *eclWriter_; }
    EclWriterType& eclWriter() { return *eclWriter_; }

    void setSubStepReport(const SimulatorReportSingle& report)
    { return eclWriter_->setSubStepReport(report); }

    void setSimulationReport(const SimulatorReport& report)
    { return eclWriter_->setSimulationReport(report); }

    void finalizeOutput()
    {
        OPM_TIMEBLOCK(finalizeOutput);
        eclWriter_.reset();
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<FlowProblemType&>(*this));
        serializer(*eclWriter_);
    }

protected:
    void updateExplicitQuantities_(int /* episodeIdx*/, int /* timeStepSize */, bool /* first_step_after_restart */) override
    {}

    void readEquilInitialCondition_() override
    {
        throw std::runtime_error(std::string("The pressure-enthalpy model does not support EQUIL: initialise with "
                                             "PRESSURE, TEMPI and ") + Flash::saturationKeyword);
    }

    /*!
     * \brief Initialise from PRESSURE (pressure of the aqueous phase), TEMPI and the saturation of the flash
     *        (Flash::saturationKeyword; see the flash for its meaning, for CO2 storage SWAT = 1 is brine without CO2,
     *        0 < SWAT < 1 brine and CO2 in equilibrium and SWAT = 0 dry CO2).
     */
    void readExplicitInitialCondition_() override
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& fp = vanguard.eclState().fieldProps();
        if (!fp.has_double("PRESSURE")) {
            throw std::runtime_error("The pressure-enthalpy model requires PRESSURE in the deck");
        }
        const bool hasSaturation = fp.has_double(Flash::saturationKeyword);
        if (Flash::saturationRequired && !hasSaturation) {
            throw std::runtime_error(std::string("The pressure-enthalpy model requires ") + Flash::saturationKeyword +
                                     " in the deck");
        }
        const auto& pressure = fp.get_double("PRESSURE");
        const auto& temperature = fp.get_double("TEMPI");
        const std::vector<double>* saturation = hasSaturation ? &fp.get_double(Flash::saturationKeyword) : nullptr;

        const std::size_t numDof = this->model().numGridDof();
        initialPressure_.resize(numDof);
        initialEnthalpy_.resize(numDof);
        initialSwitch_.resize(numDof);
        initialMeaning_.resize(numDof);
        initialFluidStates_.resize(numDof);

        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            std::array<int, 3> ijk;
            vanguard.cartesianCoordinate(dofIdx, ijk);
            const auto cell = fmt::format("cell ({}, {}, {})", ijk[0] + 1, ijk[1] + 1, ijk[2] + 1);

            const Scalar p = pressure[dofIdx];
            const Scalar T = temperature[dofIdx];
            const Scalar sat = saturation ? (*saturation)[dofIdx] : 0.0;
            Scalar h = 0.0, s = 0.0;
            int meaning = 0;
            try {
                Flash::initialVariables(p, T, sat, h, s, meaning);
            }
            catch (const std::exception& e) {
                throw std::invalid_argument(fmt::format("{} in {}", e.what(), cell));
            }
            initialPressure_[dofIdx] = p;
            initialEnthalpy_[dofIdx] = h;
            initialSwitch_[dofIdx] = s;
            initialMeaning_[dofIdx] = meaning;

            // secondary variables of the initial state
            const auto st = Flash::solve(Evaluation(p), Evaluation(h), Evaluation(s), meaning);
            auto& ifs = initialFluidStates_[dofIdx];
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                ifs.setPressure(phaseIdx, p);
                ifs.setTemperature(phaseIdx, getValue(st.temperature()));
                ifs.setSaturation(phaseIdx, getValue(st.saturation(phaseIdx)));
            }
        }
    }

    void handleSolventBC(const BCState::BCFace&, RateVector&) const override {}
    void handlePolymerBC(const BCState::BCFace&, RateVector&) const override {}
    void handleMicrBC(const BCState::BCFace&, RateVector&) const override {}
    void handleOxygBC(const BCState::BCFace&, RateVector&) const override {}
    void handleUreaBC(const BCState::BCFace&, RateVector&) const override {}

private:
    std::unique_ptr<EclWriterType> eclWriter_;
    bool enableEclOutput_{true};
    std::vector<Scalar> initialPressure_;
    std::vector<Scalar> initialEnthalpy_;
    std::vector<Scalar> initialSwitch_;
    std::vector<int> initialMeaning_;
    std::vector<InitialFluidState> initialFluidStates_;
};

} // namespace Opm

#endif
