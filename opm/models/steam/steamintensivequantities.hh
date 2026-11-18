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
 * \brief Intensive quantities of the steam model.
 *
 * (p, h) are the primary variables. WaterPHFlash turns them into the
 * temperature, the steam saturation and the phase properties, so T and S_g
 * are secondary variables with correct derivatives wrt (p, h).
 */
#ifndef OPM_STEAM_INTENSIVE_QUANTITIES_HH
#define OPM_STEAM_INTENSIVE_QUANTITIES_HH

#include <dune/common/fmatrix.hh>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/constraintsolvers/WaterPHFlash.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#include <opm/models/steam/steamindices.hh>

#include <array>

namespace Opm {

template <class TypeTag>
class SteamIntensiveQuantities
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
    using Flash = WaterPHFlash<Scalar>;
    using Water = H2O<Scalar>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { dimWorld = GridView::dimensionworld };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    // one component: pressure, temperature, saturation, density, viscosity and
    // enthalpy per phase (no composition). The thermal laws ask the state for
    // the fluid system.
    class FluidState
        : public SimpleModularFluidState<Evaluation, numPhases, /*numComponents=*/1, FluidSystem,
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

    SteamIntensiveQuantities() = default;
    SteamIntensiveQuantities(const SteamIntensiveQuantities&) = default;
    SteamIntensiveQuantities& operator=(const SteamIntensiveQuantities&) = default;

    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();

        const Evaluation p = priVars.makeEvaluation(Indices::pressure0Idx, timeIdx);
        const Evaluation h = priVars.makeEvaluation(Indices::enthalpy0Idx, timeIdx);

        // p-h flash: temperature, steam quality/saturation and phase properties
        const auto st = Flash::solve(p, h);
        region_ = st.region;
        steamQuality_ = st.steamQuality;

        // phases that are absent are evaluated on the saturation line
        const Evaluation& Tsat = st.saturationTemperature;
        const Evaluation Tl = st.region == WaterRegion::Vapor ? Tsat : st.temperature;
        const Evaluation Tg = st.region == WaterRegion::Liquid ? Tsat : st.temperature;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setSaturation(phaseIdx, Evaluation(0.0));
            fluidState_.setPressure(phaseIdx, p);
            fluidState_.setTemperature(phaseIdx, st.temperature);
        }
        fluidState_.setSaturation(gasPhaseIdx, st.steamSaturation);
        fluidState_.setSaturation(waterPhaseIdx, 1.0 - st.steamSaturation);

        // capillary pressure: the primary pressure is the liquid pressure
        const MaterialLawParams& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        std::array<Evaluation, numPhases> pC{};
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                fluidState_.setPressure(phaseIdx, p + (pC[phaseIdx] - pC[waterPhaseIdx]));
            }
        }

        // phase properties
        fluidState_.setDensity(waterPhaseIdx, st.liquidDensity);
        fluidState_.setDensity(gasPhaseIdx, st.vaporDensity);
        fluidState_.setEnthalpy(waterPhaseIdx, st.liquidEnthalpy);
        fluidState_.setEnthalpy(gasPhaseIdx, st.vaporEnthalpy);
        // inactive oil slot: keep the state defined
        fluidState_.setDensity(oilPhaseIdx, st.liquidDensity);
        fluidState_.setEnthalpy(oilPhaseIdx, st.liquidEnthalpy);

        const Evaluation muL = Water::liquidViscosity(Tl, p);
        const Evaluation muG = Water::gasViscosity(Tg, p);
        fluidState_.setViscosity(waterPhaseIdx, muL);
        fluidState_.setViscosity(gasPhaseIdx, muG);
        fluidState_.setViscosity(oilPhaseIdx, muL);

        // relative permeabilities and mobilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        mobility_[oilPhaseIdx] = 0.0;
        mobility_[waterPhaseIdx] = relativePermeability_[waterPhaseIdx] / muL;
        mobility_[gasPhaseIdx] = relativePermeability_[gasPhaseIdx] / muG;

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
    WaterRegion region() const { return region_; }
    const Evaluation& steamQuality() const { return steamQuality_; }

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
    WaterRegion region_{WaterRegion::Liquid};
    Evaluation steamQuality_;
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation rockInternalEnergy_;
    Evaluation totalThermalConductivity_;
    Scalar rockFraction_{0};
    std::array<Evaluation, numPhases> relativePermeability_;
    std::array<Evaluation, numPhases> mobility_;
};

} // namespace Opm

#endif
