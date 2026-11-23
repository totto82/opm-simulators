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
 * \brief Local residual of the CO2 storage model: mass of the components (brine and CO2) and energy.
 *
 *   d/dt [ phi sum_a S_a rho_a X_a^k ]  + div(sum_a rho_a X_a^k v_a) = q_k       (k = 1..N, X = 1 for N = 1)
 *   d/dt [ phi sum_a S_a rho_a u_a + (1-phi) u_R ] + div(sum_a rho_a h_a v_a - lambda grad T) = q_e
 *
 * The energy terms are the shared ones of Opm::EnthalpyEnergyTerms.
 */
#ifndef OPM_PH_LOCAL_RESIDUAL_HH
#define OPM_PH_LOCAL_RESIDUAL_HH

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/common/enthalpyenergyterms.hh>
#include <opm/models/ph/phindices.hh>

namespace Opm {

template <class TypeTag>
class PhLocalResidual : public GetPropType<TypeTag, Properties::DiscLocalResidual>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using EnergyTerms = EnthalpyEnergyTerms<TypeTag>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = Indices::numComponents };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiEnergyEqIdx = Indices::contiEnergyEqIdx };

    using Toolbox = MathToolbox<Evaluation>;

public:
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        const IntensiveQuantities& iq = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = iq.fluidState();
        const LhsEval poro = Toolbox::template decay<LhsEval>(iq.porosity());

        storage = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const LhsEval rhoS = Toolbox::template decay<LhsEval>(fs.density(phaseIdx))
                * Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx));
            if constexpr (numComponents == 1) {
                storage[conti0EqIdx] += poro * rhoS;
            }
            else {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    storage[conti0EqIdx + compIdx] += poro * rhoS * Toolbox::template decay<LhsEval>(iq.massFraction(phaseIdx, compIdx));
                }
            }
        }
        EnergyTerms::template addStorage<LhsEval>(storage, iq, contiEnergyEqIdx);
    }

    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        flux = 0.0;
        addAdvectiveMassFlux(flux, elemCtx, scvfIdx, timeIdx);
        EnergyTerms::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx, contiEnergyEqIdx);
        EnergyTerms::addConductiveFlux(flux, elemCtx, scvfIdx, timeIdx, contiEnergyEqIdx);
        Valgrind::CheckDefined(flux);
    }

    //! mass fluxes of the components: the upstream phase mass flux times the mass fractions
    void addAdvectiveMassFlux(RateVector& flux,
                              const ElementContext& elemCtx,
                              unsigned scvfIdx,
                              unsigned timeIdx) const
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const unsigned focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const auto upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
            const auto& fs = up.fluidState();
            if (upIdx == focusDofIdx) {
                const Evaluation massFlux = fs.density(phaseIdx) * extQuants.volumeFlux(phaseIdx);
                if constexpr (numComponents == 1) {
                    flux[conti0EqIdx] += massFlux;
                }
                else {
                    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                        flux[conti0EqIdx + compIdx] += massFlux * up.massFraction(phaseIdx, compIdx);
                    }
                }
            }
            else {
                const Evaluation massFlux = Toolbox::value(fs.density(phaseIdx)) * extQuants.volumeFlux(phaseIdx);
                if constexpr (numComponents == 1) {
                    flux[conti0EqIdx] += massFlux;
                }
                else {
                    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                        flux[conti0EqIdx + compIdx] += massFlux * Toolbox::value(up.massFraction(phaseIdx, compIdx));
                    }
                }
            }
        }
    }

    void computeSource(RateVector& source,
                       const ElementContext& elemCtx,
                       unsigned dofIdx,
                       unsigned timeIdx) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
        Valgrind::CheckDefined(source);
    }
};

} // namespace Opm

#endif
