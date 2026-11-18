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
 * \brief Terms of the energy equation for models with specific enthalpy as a
 *        primary variable (the ENTHALPY formulation).
 *
 *   d/dt [ phi sum_a S_a rho_a u_a + (1-phi) u_R ] + div( sum_a rho_a h_a v_a - lambda grad T ) = q_e
 *
 * with the internal energy u_a = h_a - p_a / rho_a.
 *
 * The terms do not depend on the primary variables of the model or on how the
 * fluid state was obtained. Whether h or T is a primary variable only changes
 * which derivatives the AD evaluations carry; these come from the flash that
 * fills the fluid state (for example a p-h flash for pure water, a p-h flash
 * for a compositional mixture, or a temperature solve for a black-oil fluid).
 * A model uses these functions from its local residual, next to its own mass
 * balance equations, and supplies the energy equation index.
 *
 * Requirements on the model (checked at compile time where possible, see
 * Opm::EnthalpyFluidState):
 *  - IntensiveQuantities: fluidState(), porosity(), rockFraction(),
 *    rockInternalEnergy(), totalThermalConductivity().
 *  - Fluid state: saturation(phase), density(phase), enthalpy(phase),
 *    pressure(phase) and temperature(phase). Density and enthalpy must be on
 *    the same basis (both per mass, or both per mole) so that rho*h and
 *    p/rho*rho are energy per bulk volume. The energy equation is then in
 *    J/m3 (per mass) or in the molar equivalent.
 *  - ExtensiveQuantities: upstreamIndex(phase), volumeFlux(phase).
 *  - Problem: thermalHalfTransmissibilityIn/Out(elemCtx, scvfIdx, timeIdx).
 *  - FluidSystem: numPhases, phaseIsActive(phase).
 *
 * All phases that are active in the fluid system take part in the sums; a
 * phase that is absent in a cell must report zero saturation.
 */
#ifndef OPM_ENTHALPY_ENERGY_TERMS_HH
#define OPM_ENTHALPY_ENERGY_TERMS_HH

#include <opm/material/common/MathToolbox.hpp>

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <concepts>
#include <type_traits>

namespace Opm {

namespace Parameters {
//! Scaling of the energy equation so that its residual is of the order of the mass residuals
struct EnthalpyEnergyScaling { static constexpr double value = 1.0e-6; };
}

//! The fluid state accessors that the enthalpy energy terms need.
template <class FluidState>
concept EnthalpyFluidState = requires(const FluidState& fs, unsigned phaseIdx) {
    fs.saturation(phaseIdx);
    fs.density(phaseIdx);
    fs.enthalpy(phaseIdx);
    fs.pressure(phaseIdx);
    fs.temperature(phaseIdx);
};

template <class TypeTag>
class EnthalpyEnergyTerms
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using Toolbox = MathToolbox<Evaluation>;

    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    //! Factor applied to every term of the energy equation
    static constexpr Scalar scaling()
    { return Parameters::EnthalpyEnergyScaling::value; }

    /*!
     * \brief Adds the (scaled) energy storage of the fluids and the rock to storage[energyEqIdx].
     */
    template <class LhsEval, class Storage>
    static void addStorage(Storage& storage,
                           const IntensiveQuantities& iq,
                           unsigned energyEqIdx)
    {
        static_assert(EnthalpyFluidState<std::remove_cvref_t<decltype(iq.fluidState())>>,
                      "the fluid state does not provide what the enthalpy energy terms need");

        const auto& fs = iq.fluidState();
        const LhsEval poro = Toolbox::template decay<LhsEval>(iq.porosity());

        LhsEval energy = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const LhsEval rho = Toolbox::template decay<LhsEval>(fs.density(phaseIdx));
            const LhsEval S = Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx));
            const LhsEval h = Toolbox::template decay<LhsEval>(fs.enthalpy(phaseIdx));
            const LhsEval p = Toolbox::template decay<LhsEval>(fs.pressure(phaseIdx));
            // internal energy u = h - p/rho
            energy += poro * S * rho * (h - p / rho);
        }
        energy += iq.rockFraction() * Toolbox::template decay<LhsEval>(iq.rockInternalEnergy());

        storage[energyEqIdx] += energy * scaling();
    }

    /*!
     * \brief Adds the (scaled) advective energy flux sum_a rho_a h_a v_a to flux[energyEqIdx].
     *
     * Only the energy equation is touched; the model adds the mass fluxes itself.
     */
    static void addAdvectiveFlux(RateVector& flux,
                                 const ElementContext& elemCtx,
                                 unsigned scvfIdx,
                                 unsigned timeIdx,
                                 unsigned energyEqIdx)
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const unsigned focusDofIdx = elemCtx.focusDofIndex();

        Evaluation energyFlux = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const auto upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const auto& upFs = elemCtx.intensiveQuantities(upIdx, timeIdx).fluidState();
            // only the focus cell carries derivatives
            if (upIdx == focusDofIdx) {
                energyFlux += upFs.enthalpy(phaseIdx) * upFs.density(phaseIdx) * extQuants.volumeFlux(phaseIdx);
            }
            else {
                energyFlux += Toolbox::value(upFs.enthalpy(phaseIdx)) * Toolbox::value(upFs.density(phaseIdx))
                    * extQuants.volumeFlux(phaseIdx);
            }
        }
        flux[energyEqIdx] += energyFlux * scaling();
    }

    /*!
     * \brief Adds the (scaled) heat conduction -lambda grad T, with the thermal half
     *        transmissibilities of the problem, to flux[energyEqIdx].
     */
    static void addConductiveFlux(RateVector& flux,
                                  const ElementContext& elemCtx,
                                  unsigned scvfIdx,
                                  unsigned timeIdx,
                                  unsigned energyEqIdx)
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);
        const unsigned inIdx = scvf.interiorIndex();
        const unsigned exIdx = scvf.exteriorIndex();
        const unsigned focusIdx = elemCtx.focusDofIndex();
        const auto& inIq = elemCtx.intensiveQuantities(inIdx, timeIdx);
        const auto& exIq = elemCtx.intensiveQuantities(exIdx, timeIdx);

        const Scalar inAlpha = elemCtx.problem().thermalHalfTransmissibilityIn(elemCtx, scvfIdx, timeIdx);
        const Scalar outAlpha = elemCtx.problem().thermalHalfTransmissibilityOut(elemCtx, scvfIdx, timeIdx);

        // all phases of a cell are in thermal equilibrium, take the first active one
        unsigned refPhase = 0;
        while (refPhase < numPhases - 1 && !FluidSystem::phaseIsActive(refPhase)) {
            ++refPhase;
        }

        const auto pick = [&](unsigned idx, const Evaluation& e) {
            return idx == focusIdx ? e : Evaluation(Toolbox::value(e));
        };
        const Evaluation Tin = pick(inIdx, inIq.fluidState().temperature(refPhase));
        const Evaluation Tex = pick(exIdx, exIq.fluidState().temperature(refPhase));
        const Evaluation inH = pick(inIdx, inIq.totalThermalConductivity()) * inAlpha;
        const Evaluation exH = pick(exIdx, exIq.totalThermalConductivity()) * outAlpha;

        if (inH > 0 && exH > 0) {
            const Evaluation H = 1.0 / (1.0 / inH + 1.0 / exH);
            flux[energyEqIdx] += -(Tex - Tin) * H / scvf.area() * scaling();
        }
    }
};

} // namespace Opm

#endif
