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
 * \brief Intensive quantities of the CO2 storage model.
 *
 * (p, h, s) are the primary variables (s only if the flash has more than one component). The flash turns them into the
 * temperature, the saturations of the phases, the phase properties and the compositions,
 * with correct derivatives wrt (p, h, s).
 */
#ifndef OPM_PH_INTENSIVE_QUANTITIES_HH
#define OPM_PH_INTENSIVE_QUANTITIES_HH

#include <dune/common/fmatrix.hh>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#include <opm/models/ph/phflash.hh>
#include <opm/models/ph/phindices.hh>

#include <array>

namespace Opm {

template <class TypeTag>
class PhIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;

    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidEnergyLaw = GetPropType<TypeTag, Properties::SolidEnergyLaw>;
    using ThermalConductionLaw = GetPropType<TypeTag, Properties::ThermalConductionLaw>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    using Flash = GetPropType<TypeTag, Properties::PhFlash>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = FluidSystem::numComponents };
    enum { dimWorld = GridView::dimensionworld };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    // pressure, temperature, saturation, density, viscosity and enthalpy per phase. The mass
    // fractions of the components are kept separately. The thermal laws ask the state for the fluid system.
    class FluidState
        : public SimpleModularFluidState<Evaluation, numPhases, numComponents, FluidSystem,
                                         /*storePressure=*/true,
                                         /*storeTemperature=*/true,
                                         /*storeComposition=*/false,
                                         /*storeFugacity=*/false,
                                         /*storeSaturation=*/true,
                                         /*storeDensity=*/true,
                                         /*storeViscosity=*/true,
                                         /*storeEnthalpy=*/true>
    {
    public:
        const FluidSystem& fluidSystem() const
        {
            static const FluidSystem fluidSystem;
            return fluidSystem;
        }
    };

    PhIntensiveQuantities() = default;
    PhIntensiveQuantities(const PhIntensiveQuantities&) = default;
    PhIntensiveQuantities& operator=(const PhIntensiveQuantities&) = default;

    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();

        const Evaluation p = priVars.makeEvaluation(Indices::pressure0Idx, timeIdx);
        const Evaluation h = priVars.makeEvaluation(Indices::enthalpy0Idx, timeIdx);
        Evaluation s{};
        if constexpr (Flash::numSwitchVariables > 0) {
            s = priVars.makeEvaluation(Indices::switchIdx, timeIdx);
        }

        // p-h flash: temperature, saturations, phase properties and compositions
        const auto st = Flash::solve(p, h, s, priVars.meaning());
        regionId_ = st.regionId();

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setSaturation(phaseIdx, st.saturation(phaseIdx));
            fluidState_.setPressure(phaseIdx, p);
            fluidState_.setTemperature(phaseIdx, st.temperature());
            if constexpr (numComponents > 1) {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    massFraction_[phaseIdx][compIdx] = st.massFraction(phaseIdx, compIdx);
                }
            }
        }

        // capillary pressure: the primary pressure is the pressure of the aqueous phase
        const MaterialLawParams& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        std::array<Evaluation, numPhases> pC{};
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                fluidState_.setPressure(phaseIdx, p + (pC[phaseIdx] - pC[waterPhaseIdx]));
            }
        }

        // phase properties
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setDensity(phaseIdx, st.density(phaseIdx));
            fluidState_.setEnthalpy(phaseIdx, st.enthalpy(phaseIdx));
            fluidState_.setViscosity(phaseIdx, st.viscosity(phaseIdx));
        }

        // relative permeabilities and mobilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            mobility_[phaseIdx] = FluidSystem::phaseIsActive(phaseIdx)
                ? relativePermeability_[phaseIdx] / st.viscosity(phaseIdx)
                : Evaluation(0.0);
        }

        // porosity, scaled with pressure by the rock compressibility
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        const Scalar rockCompressibility = problem.rockCompressibility(elemCtx, dofIdx, timeIdx);
        if (rockCompressibility > 0.0) {
            const Scalar rockRefPressure = problem.rockReferencePressure(elemCtx, dofIdx, timeIdx);
            const Evaluation x = rockCompressibility * (p - rockRefPressure);
            porosity_ *= 1.0 + x + 0.5 * x * x;
        }

        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // rock energy and thermal conductivity
        const auto& solidEnergyParams = problem.solidEnergyLawParams(elemCtx, dofIdx, timeIdx);
        rockInternalEnergy_ = SolidEnergyLaw::solidInternalEnergy(solidEnergyParams, fluidState_);
        const auto& thermalCondParams = problem.thermalConductionLawParams(elemCtx, dofIdx, timeIdx);
        totalThermalConductivity_ = ThermalConductionLaw::thermalConductivity(thermalCondParams, fluidState_);
        rockFraction_ = problem.rockFraction(elemCtx.globalSpaceIndex(dofIdx, timeIdx), timeIdx);
    }

    const FluidState& fluidState() const { return fluidState_; }

    //! mass fraction of a component in a phase (always one for a single component)
    const Evaluation& massFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        if constexpr (numComponents > 1) {
            return massFraction_[phaseIdx][compIdx];
        }
        else {
            static const Evaluation one(1.0);
            return one;
        }
    }

    //! id of the phase region of the flash
    int regionId() const { return regionId_; }

    const DimMatrix& intrinsicPermeability() const { return intrinsicPerm_; }
    const Evaluation& relativePermeability(unsigned phaseIdx) const { return relativePermeability_[phaseIdx]; }
    const Evaluation& mobility(unsigned phaseIdx) const { return mobility_[phaseIdx]; }
    const Evaluation& mobility(unsigned phaseIdx, FaceDir::DirEnum) const { return mobility_[phaseIdx]; }
    Scalar rockCompTransMultiplier() const { return 1.0; }
    const Evaluation& porosity() const { return porosity_; }

    const Evaluation& rockInternalEnergy() const { return rockInternalEnergy_; }
    const Evaluation& totalThermalConductivity() const { return totalThermalConductivity_; }
    Scalar rockFraction() const { return rockFraction_; }

    bool phaseIsPresent(unsigned phaseIdx) const
    { return saturationForOutput(phaseIdx) > Scalar{0}; }

    //! saturation for reporting
    Scalar saturationForOutput(unsigned phaseIdx) const
    { return FluidSystem::phaseIsActive(phaseIdx) ? getValue(fluidState_.saturation(phaseIdx)) : Scalar{0}; }

private:
    int regionId_{0};
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation rockInternalEnergy_;
    Evaluation totalThermalConductivity_;
    Scalar rockFraction_{0};
    std::array<Evaluation, numPhases> relativePermeability_;
    std::array<Evaluation, numPhases> mobility_;
    std::array<std::array<Evaluation, (numComponents > 1 ? numComponents : 0)>, numPhases> massFraction_;
};

} // namespace Opm

#endif
