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
 *
 * \brief Steady-state wellbore model with phase behaviour: the fluid in the well flows
 *        between the perforations and the surface, and may change phase on the way (hot
 *        water that flashes to steam while it rises, steam that condenses while it is
 *        injected).
 *
 * The well is a 1D path from the surface (measured depth 0) to the deepest perforation, with
 * a true vertical depth, an inner diameter and a roughness for every interval. The control is
 * applied at the surface:
 *  - a producer has a tubing head pressure; the bottom-hole pressure is a result,
 *  - an injector has a volume rate at the surface pressure and the injection temperature; the
 *    surface pressure (and the bottom-hole pressure) are results.
 *
 * Along the path, for the flow direction s (up for a producer, down for an injector), with the
 * mass flow mdot, the mass flux G = mdot/A and the elevation z (up positive):
 *   mass:      mdot is the sum of the perforation inflows (outflows) that were passed,
 *   energy:    h + kinetic energy + g z is constant in an interval apart from the heat loss (mixing at
 *              perforations uses the static enthalpies),
 *   momentum:  dp/ds = - rho_m g dz/ds - f/D G^2/(2 rho) - d(momentum flux)/ds.
 * The phase state, the temperature and the phase densities follow from the (p, h) flash of the
 * fluid, as for the reservoir models. The phases flow with the same velocity (homogeneous flow,
 * the default) or with slip according to the Zuber-Findlay drift-flux model,
 *   u_g = C0 j + u_d,   alpha = j_g / (C0 j + u_d),
 * with j the total superficial velocity, the bubble rise velocity u_d = 1.53 (g sigma (rho_l - rho_g) / rho_l^2)^(1/4)
 * times the sine of the inclination, and the distribution parameter C0. u_d changes sign with the
 * direction of the flow (the gas rises against a downward flow), and C0 -> 1 and u_d -> 0 for a pure
 * vapour so that alpha = 1 at a quality of 1 (a regularisation: the correlation is for intermediate
 * void fractions). The in-situ mixture density rho_m = alpha rho_g + (1 - alpha) rho_l goes into the
 * hydrostatic term and the kinetic energy and momentum flux follow from the phase velocities; the
 * friction uses the no-slip density and viscosity and the Churchill friction factor.
 *
 * The well can lose heat to the formation, q' (T - T_ambient(z)) per length with an overall heat
 * transfer coefficient q' [W/(m K)] and an ambient temperature that is linear in the true vertical
 * depth. The loss is implicit in the temperature at the end of a step, so that the fluid cannot
 * cool below the ambient temperature however low the flow. A well without flow is in thermal
 * equilibrium with the formation (liquid at the ambient temperature).
 *
 * The inflow performance of a perforation does not depend on how the reservoir cell was obtained:
 *  producer: mdot_i = productivity_i (p_res,i - p_well,i)^+ at the enthalpy inflowEnthalpy_i
 *  injector: mdot_i = injectivity_i rho_well,i (p_well,i - p_res,i)^+ at the enthalpy of the well fluid
 * where productivity [kg/(s Pa)] and injectivity [m3/(s Pa)] include the well index and the
 * mobilities of the cell. A standalone run builds them from a reservoir state, a model that
 * couples the well to a simulator takes them from the cells.
 *
 * Marching along the well is explicit for a given bottom-hole (producer) or surface (injector)
 * pressure, so that the well is solved by a scalar root search on that pressure. A producer
 * takes the first solution from the high pressure side, which is the stable one if the
 * inflow performance curve and the outflow curve intersect twice.
 *
 * The fluid is a class with
 *   static constexpr std::size_t numComponents; using Composition = std::array<Scalar, numComponents>;
 *   using Stream = ...;   // what is injected (the phase of the injected steam, or its composition)
 *   struct State { temperature, density (no slip), viscosity, quality (steam mass fraction), voidFraction (no slip),
 *                  densityLiquid, densityVapor, surfaceTension; int region; };
 *   static State state(Scalar p, Scalar h, const Composition& z);   // throws if (p,h,z) is not valid
 *   static Scalar phaseChangeMargin(Scalar p, Scalar h, const Composition& z, int boundary);   // optional
 *                           // changes sign at the boundary between the regions boundary and boundary + 1
 *   static Scalar injectionEnthalpy(Scalar p, Scalar T, const Stream&);
 *   static Composition injectionComposition(const Stream&);
 *   static Scalar equilibriumEnthalpy(Scalar p, Scalar T, const Composition&);   // fluid at rest at temperature T
 *   static Scalar minPressure(); static Scalar maxPressure();
 * see Opm::SteamWellFluid.
 */
#ifndef OPM_WELL_BORE_HPP
#define OPM_WELL_BORE_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm {

//! The path of the well from the surface down
template <class Scalar>
struct WellBoreGeometry
{
    struct Point
    {
        Scalar md;    //!< measured depth [m]
        Scalar tvd;   //!< true vertical depth [m]
    };

    std::vector<Point> points;       //!< from the surface (md = 0) down, md increasing
    std::vector<Scalar> diameter;    //!< inner diameter [m] of the interval that ends at points[k + 1]
    std::vector<Scalar> roughness;   //!< absolute roughness [m] of the same interval
    Scalar maxStepLength{10.0};      //!< the intervals are subdivided so that no step is longer [m]

    //! a vertical well of the given length with a constant diameter
    static WellBoreGeometry vertical(Scalar length, Scalar diameter, Scalar roughness = 4.5e-5)
    {
        WellBoreGeometry g;
        g.points = {{0.0, 0.0}, {length, length}};
        g.diameter = {diameter};
        g.roughness = {roughness};
        return g;
    }
};

//! Heat loss to the formation and slip between the phases
template <class Scalar>
struct WellBoreOptions
{
    struct HeatLoss
    {
        Scalar coefficient{0.0};            //!< overall coefficient per length q' [W/(m K)], 0: adiabatic
        Scalar surfaceTemperature{288.15};  //!< ambient temperature at the surface [K]
        Scalar gradient{0.03};              //!< increase of the ambient temperature with the depth [K/m]
    } heatLoss;

    struct Slip
    {
        bool driftFlux{false};              //!< false: homogeneous flow
        Scalar distributionParameter{1.2};  //!< C0
        Scalar driftCoefficient{1.53};      //!< u_d = coefficient (g sigma (rho_l - rho_g) / rho_l^2)^(1/4)
    } slip;
};

template <class Scalar, std::size_t numComponents = 1>
struct WellBorePerforation
{
    Scalar md{0};                //!< measured depth [m]
    Scalar pressure{0};          //!< reservoir pressure [Pa]
    Scalar productivity{0};      //!< producer: mass rate per pressure difference [kg/(s Pa)]
    Scalar inflowEnthalpy{0};    //!< producer: enthalpy of the inflow [J/kg]
    Scalar injectivity{0};       //!< injector: volume rate per pressure difference [m3/(s Pa)]
    //! producer: mass fractions of the components of the inflow
    std::array<Scalar, numComponents> inflowComposition{};
};

template <class Scalar>
struct WellBoreProfilePoint
{
    Scalar md, tvd;
    Scalar pressure, temperature, enthalpy;
    Scalar massFlow;       //!< [kg/s]
    Scalar steamQuality;   //!< mass fraction of the vapour
    Scalar voidFraction;   //!< volume fraction of the vapour in the pipe (with slip)
    Scalar density;        //!< in-situ mixture density [kg/m3]
    Scalar velocity;       //!< total superficial velocity [m/s]
    Scalar kineticEnergy;  //!< kinetic energy of the flow per mass [J/kg]
    Scalar heatLossRate{0};//!< heat lost to the formation per length at this point [W/m]
    std::vector<Scalar> composition;   //!< mass fractions of the components
    int region;
};

template <class Scalar>
struct WellBoreResult
{
    bool flowing{false};
    Scalar surfacePressure{0};        //!< [Pa]
    Scalar bottomHolePressure{0};     //!< pressure at the deepest perforation [Pa]
    Scalar massFlow{0};               //!< at the surface [kg/s]
    Scalar surfaceEnthalpy{0};        //!< [J/kg]
    Scalar surfaceTemperature{0};     //!< [K]
    //! per perforation: mass rate [kg/s] (producer: into the well, injector: into the reservoir)
    std::vector<Scalar> perforationRate;
    //! per perforation: enthalpy of the fluid that crosses it [J/kg]
    std::vector<Scalar> perforationEnthalpy;
    //! per perforation: pressure of the fluid in the well at the perforation [Pa]
    std::vector<Scalar> perforationPressure;
    //! the profile from the surface down (substeps included)
    std::vector<WellBoreProfilePoint<Scalar>> profile;
    //! where the fluid first changes phase region in the direction of the flow (flashing in a
    //! producer, condensation or the end of the superheat in an injector)
    std::optional<Scalar> phaseChangeMd;
    std::optional<Scalar> phaseChangePressure;
    //! injector: mass flow that was not taken up by the perforations (zero at the solution)
    Scalar remainingMassFlow{0};
    //! injector: false if the perforations take more than the specified rate even at the lowest surface pressure
    //! (the weight of the fluid column alone is enough; there is no choke in the model). The result is then the
    //! state at the lowest surface pressure.
    bool rateReached{true};
};

template <class Scalar, class Fluid>
class WellBore
{
public:
    using Geometry = WellBoreGeometry<Scalar>;
    static constexpr std::size_t numComponents = Fluid::numComponents;
    using Composition = std::array<Scalar, numComponents>;
    using Stream = typename Fluid::Stream;
    using Perforation = WellBorePerforation<Scalar, numComponents>;
    using Result = WellBoreResult<Scalar>;
    using ProfilePoint = WellBoreProfilePoint<Scalar>;
    using Options = WellBoreOptions<Scalar>;

    static constexpr Scalar gravity = 9.80665;

    /*!
     * The perforations must lie on the path. The path is cut after the deepest perforation.
     */
    WellBore(const Geometry& geometry, std::vector<Perforation> perforations, const Options& options = {})
        : perforations_(std::move(perforations))
        , options_(options)
    {
        if (geometry.points.size() < 2 || geometry.diameter.size() + 1 != geometry.points.size()
            || geometry.roughness.size() + 1 != geometry.points.size()) {
            throw std::invalid_argument("WellBore: inconsistent geometry");
        }
        if (perforations_.empty()) {
            throw std::invalid_argument("WellBore: no perforations");
        }
        geometry_ = geometry;
        maxStep_ = geometry.maxStepLength;
        const Scalar mdLast = geometry.points.back().md;
        Scalar mdDeepest = 0.0;
        for (const auto& perf : perforations_) {
            if (!(perf.md > 0.0) || perf.md > mdLast * (1.0 + 1e-12)) {
                throw std::invalid_argument("WellBore: a perforation is not on the well path");
            }
            mdDeepest = std::max(mdDeepest, perf.md);
        }

        // the nodes: the points of the path and the perforations, cut at the deepest perforation
        std::vector<Scalar> mds;
        for (const auto& pt : geometry.points) {
            if (pt.md < mdDeepest * (1.0 - 1e-12)) {
                mds.push_back(pt.md);
            }
        }
        mds.push_back(mdDeepest);
        for (const auto& perf : perforations_) {
            mds.push_back(perf.md);
        }
        std::sort(mds.begin(), mds.end());
        mds.erase(std::unique(mds.begin(), mds.end(),
                              [](Scalar a, Scalar b) { return std::abs(a - b) < 1e-9 * (1.0 + std::abs(a)); }),
                  mds.end());

        for (const Scalar md : mds) {
            Node node;
            node.md = md;
            // the interval of the path that contains md
            std::size_t k = 0;
            while (k + 2 < geometry.points.size() && geometry.points[k + 1].md < md * (1.0 - 1e-12)) {
                ++k;
            }
            const auto& a = geometry.points[k];
            const auto& b = geometry.points[k + 1];
            node.tvd = a.tvd + (b.tvd - a.tvd) * (md - a.md) / (b.md - a.md);
            node.diameter = geometry.diameter[k];
            node.roughness = geometry.roughness[k];
            nodes_.push_back(node);
        }
        for (std::size_t i = 0; i < perforations_.size(); ++i) {
            for (std::size_t k = 0; k < nodes_.size(); ++k) {
                if (std::abs(nodes_[k].md - perforations_[i].md) < 1e-9 * (1.0 + nodes_[k].md)) {
                    nodes_[k].perforation = static_cast<int>(i);
                }
            }
        }
    }

    const std::vector<Perforation>& perforations() const { return perforations_; }

    //! number of marches along the well (producerAt and injectorAt calls) so far, to measure the cost of a solve
    std::size_t marchCount() const { return *marchCount_; }

    /*!
     * \brief Producer with the given tubing head pressure.
     *
     * If the well cannot flow (the pressure at the surface of a well without flow is below
     * the tubing head pressure) the result has flowing = false and no flow. bottomPressureHint is an
     * optional guess of the bottom-hole pressure that saves marches (the result does not depend on it).
     */
    Result solveProducer(Scalar tubingHeadPressure, std::optional<Scalar> bottomPressureHint = std::nullopt) const
    {
        Scalar pHigh = 0.0;
        for (const auto& perf : perforations_) {
            pHigh = std::max(pHigh, perf.pressure);
        }
        const auto f = [&](Scalar pb) {
            try {
                return producerAt(pb, false).surfacePressure - tubingHeadPressure;
            }
            catch (const std::exception&) {
                return -huge_;   // too low a pressure for a valid state
            }
        };

        // A hint (the solution of a similar problem, for example the previous Newton iteration of a simulator) is
        // used to bracket the root with a few marches. Anything unexpected falls back to the scan below.
        if (bottomPressureHint && *bottomPressureHint > Fluid::minPressure() && *bottomPressureHint < pHigh) {
            if (const auto pb = warmProducer_(f, *bottomPressureHint, pHigh, tubingHeadPressure)) {
                return producerAt(*pb, true);
            }
        }

        // Scan down from the pressure at which nothing flows. f is the surface pressure that the
        // well delivers minus the tubing head pressure. With heat loss f is negative at low rates (a cold
        // column is heavy) and positive when the flow keeps the fluid hot: the well can only flow
        // once f has been positive. The stable operating point is where f falls to zero again
        // (the outflow curve rises with the rate), the first such point from the high pressure side.
        const Scalar step = 0.005 * pHigh;
        const Scalar pMin = Fluid::minPressure();
        Scalar upper = pHigh;
        Scalar lower = pHigh;
        Scalar fUpper = f(pHigh);
        Scalar fLower = fUpper;
        bool seenPositive = fUpper > 0.0;
        bool found = false;
        while (lower - step > pMin) {
            upper = lower;
            fUpper = fLower;
            lower -= step;
            fLower = f(lower);
            if (fLower > 0.0) {
                seenPositive = true;
            }
            else if (seenPositive) {
                found = true;
                break;
            }
        }
        if (!seenPositive) {
            // dead well: the pressure cannot be delivered, whatever a shallower perforation would still let in
            Result r = producerAt(pHigh, true);
            r.flowing = false;
            r.massFlow = 0.0;
            std::fill(r.perforationRate.begin(), r.perforationRate.end(), Scalar(0));
            return r;
        }
        if (!found) {
            throw std::runtime_error("WellBore: no producer solution above the minimum pressure");
        }
        const Scalar pb = brent_(f, lower, upper, fLower, fUpper, 1e-10 * pHigh);
        if (!(std::abs(f(pb)) < 1e-6 * tubingHeadPressure)) {
            throw std::runtime_error("WellBore: the producer solution is outside the validity range of the fluid");
        }
        return producerAt(pb, true);
    }

    /*!
     * \brief Injector with the given volume rate [m3/s] of the stream (steam or liquid for water,
     *        a composition for CO2 and brine) at the injection temperature and the surface pressure.
     *
     * The mass rate is the volume rate times the density at the surface pressure, so that for
     * steam the problem can have several solutions: one at a very low surface pressure, where the
     * little steam condenses and fills the well with water, and the solutions where the steam reaches
     * the reservoir. The solution with the highest surface pressure is taken. If the perforations take more than
     * the rate even at the lowest surface pressure, the result is that state with rateReached = false.
     */
    Result solveInjector(Scalar surfaceVolumeRate, Scalar temperature, const Stream& stream,
                         std::optional<Scalar> surfacePressureHint = std::nullopt) const
    {
        // the remaining flow at the bottom; not valid if the march leaves the range of the fluid
        const auto remaining = [&](Scalar ps) -> std::optional<Scalar> {
            try {
                return injectorAt(ps, surfaceVolumeRate, temperature, stream, false).remainingMassFlow;
            }
            catch (const std::exception&) {
                return std::nullopt;
            }
        };
        const auto f = [&](Scalar ps) {
            const auto v = remaining(ps);
            return v ? *v : huge_;
        };

        // warm start from the hint, see solveProducer()
        if (surfacePressureHint && *surfacePressureHint > Fluid::minPressure() && *surfacePressureHint < Fluid::maxPressure()) {
            if (const auto ps = warmInjector_(remaining, *surfacePressureHint)) {
                const Result r = injectorAt(*ps, surfaceVolumeRate, temperature, stream, true);
                if (std::abs(r.remainingMassFlow) < 1e-6 * r.massFlow) {
                    return r;
                }
            }
        }

        // Scan the surface pressure upwards. The remaining flow falls through zero where the surface pressure is
        // right; only two neighbouring valid evaluations count as a bracket, and states where the stream does not
        // exist (for example liquid water above its saturation temperature) are invalid.
        const Scalar pMax = Fluid::maxPressure();
        Scalar prev = Fluid::minPressure() * 1.05;
        auto fPrev = remaining(prev);
        std::optional<std::pair<Scalar, Scalar>> firstValid;
        if (fPrev) {
            firstValid = std::pair{prev, *fPrev};
        }
        bool found = false;
        Scalar lower = prev, upper = prev, fLower = 0.0, fUpper = 0.0;
        // Without a root so far, a remaining flow that has fallen steadily over this many steps (a factor of 1.02^25
        // = 1.6 in the pressure) is the gravity-driven injector: the scan stops instead of marching to the maximum pressure.
        constexpr int stepsOfDecline = 25;
        int falling = 0;
        while (prev < pMax) {
            const Scalar ps = std::min(pMax, prev * 1.02);
            const auto fv = remaining(ps);
            if (fv && !firstValid) {
                firstValid = std::pair{ps, *fv};
            }
            falling = (!found && fPrev && fv && *fv <= 0.0 && *fv < *fPrev) ? falling + 1 : 0;
            if (falling >= stepsOfDecline) {
                break;
            }
            if (fPrev && fv && *fPrev > 0.0 && *fv <= 0.0) {
                found = true;
                lower = prev;
                upper = ps;
                fLower = *fPrev;
                fUpper = *fv;
            }
            prev = ps;
            fPrev = fv;
        }
        if (!found) {
            if (firstValid && firstValid->second <= 0.0) {
                // the weight of the fluid column alone injects more than the rate, there is no choke in the model:
                // the state at the lowest surface pressure at which the stream exists
                Result r = injectorAt(firstValid->first, surfaceVolumeRate, temperature, stream, true);
                r.rateReached = false;
                return r;
            }
            throw std::runtime_error("WellBore: the injection rate cannot be taken up below the maximum pressure");
        }
        // f decreases through zero: use the negated function so that brent_ sees an increasing one
        const auto g = [&](Scalar ps) { return -f(ps); };
        const Scalar ps = brent_(g, lower, upper, -fLower, -fUpper, 1e-10 * upper);
        const Result r = injectorAt(ps, surfaceVolumeRate, temperature, stream, true);
        if (!(std::abs(r.remainingMassFlow) < 1e-6 * r.massFlow)) {
            throw std::runtime_error("WellBore: the injector solution is outside the validity range of the fluid");
        }
        return r;
    }

    //! A change of one descriptor of a perforation, for the sensitivities of the solution
    struct Perturbation
    {
        std::size_t perforation{0};
        std::function<void(Perforation&)> apply;   //!< changes the perforation by `size` in some variable
        Scalar size{1.0};                          //!< the size of that change
    };

    //! The same well with other perforations (the path and the options are kept)
    WellBore withPerforations(std::vector<Perforation> perforations) const
    {
        WellBore other(geometry_, std::move(perforations), options_);
        other.marchCount_ = marchCount_;
        return other;
    }

    /*!
     * \brief Sensitivities of the pressure of the fluid at the perforations of a producer to descriptors of the perforations.
     *
     * \param base the result of solveProducer(thp) (with the profile)
     * The solution is the pressure pb at the bottom with p_top(pb; d) = thp. With the implicit function theorem
     *   dpb/dd = - (dp_top/dd) / (dp_top/dpb),
     * so that the sensitivities cost one march for dp_top/dpb and one per perturbation, instead of a solution of the
     * well for every perturbation. The pressure at a perforation i is a function of pb and d, and
     *   dp_i/dd = (partial p_i / partial d) + (partial p_i / partial pb) dpb/dd.
     * \return result[k][i]: derivative of the pressure at perforation i wrt perturbation k; nothing if the well does not flow
     *         or the derivative of the surface pressure wrt pb is too small to be trusted (the caller solves the perturbed wells).
     */
    std::optional<std::vector<std::vector<Scalar>>> producerSensitivities(const Result& base,
                                                                         const std::vector<Perturbation>& perturbations) const
    {
        if (!base.flowing || base.perforationPressure.size() != perforations_.size()) {
            return std::nullopt;
        }
        const Scalar pb = base.bottomHolePressure;
        const Scalar dpb = 1.0e-6 * pb;
        Result up;
        try {
            up = producerAt(pb + dpb, false);
        }
        catch (const std::exception&) {
            return std::nullopt;
        }
        const Scalar kTop = (up.surfacePressure - base.surfacePressure) / dpb;   // dp_top / dpb
        if (!(kTop > 1.0e-3)) {
            return std::nullopt;
        }
        std::vector<Scalar> dNodeDpb(perforations_.size());
        for (std::size_t i = 0; i < perforations_.size(); ++i) {
            dNodeDpb[i] = (up.perforationPressure[i] - base.perforationPressure[i]) / dpb;
        }
        std::vector<std::vector<Scalar>> result;
        for (const auto& pert : perturbations) {
            auto perforations = perforations_;
            pert.apply(perforations[pert.perforation]);
            Result r;
            try {
                r = withPerforations(std::move(perforations)).producerAt(pb, false);
            }
            catch (const std::exception&) {
                return std::nullopt;
            }
            const Scalar dTopDd = (r.surfacePressure - base.surfacePressure) / pert.size;
            const Scalar dpbDd = -dTopDd / kTop;
            std::vector<Scalar> row(perforations_.size());
            for (std::size_t i = 0; i < perforations_.size(); ++i) {
                row[i] = (r.perforationPressure[i] - base.perforationPressure[i]) / pert.size + dNodeDpb[i] * dpbDd;
            }
            result.push_back(std::move(row));
        }
        return result;
    }

    /*!
     * \brief Sensitivities of the pressure of the fluid at the perforations of an injector, as producerSensitivities().
     *
     * The solution is the surface pressure ps at which the perforations take up the rate, remaining(ps; d) = 0.
     * \param base the result of solveInjector(...) with rateReached
     */
    std::optional<std::vector<std::vector<Scalar>>> injectorSensitivities(const Result& base, Scalar volumeRate,
                                                                         Scalar temperature, const Stream& stream,
                                                                         const std::vector<Perturbation>& perturbations) const
    {
        if (!base.flowing || !base.rateReached || base.perforationPressure.size() != perforations_.size()) {
            return std::nullopt;
        }
        const Scalar ps = base.surfacePressure;
        const Scalar dps = 1.0e-6 * ps;
        Result up;
        try {
            up = injectorAt(ps + dps, volumeRate, temperature, stream, false);
        }
        catch (const std::exception&) {
            return std::nullopt;
        }
        const Scalar kRemaining = (up.remainingMassFlow - base.remainingMassFlow) / dps;   // d remaining / d ps
        if (!(kRemaining < -1.0e-9)) {
            return std::nullopt;
        }
        std::vector<Scalar> dNodeDps(perforations_.size());
        for (std::size_t i = 0; i < perforations_.size(); ++i) {
            dNodeDps[i] = (up.perforationPressure[i] - base.perforationPressure[i]) / dps;
        }
        std::vector<std::vector<Scalar>> result;
        for (const auto& pert : perturbations) {
            auto perforations = perforations_;
            pert.apply(perforations[pert.perforation]);
            Result r;
            try {
                r = withPerforations(std::move(perforations)).injectorAt(ps, volumeRate, temperature, stream, false);
            }
            catch (const std::exception&) {
                return std::nullopt;
            }
            const Scalar dRemainingDd = (r.remainingMassFlow - base.remainingMassFlow) / pert.size;
            const Scalar dpsDd = -dRemainingDd / kRemaining;
            std::vector<Scalar> row(perforations_.size());
            for (std::size_t i = 0; i < perforations_.size(); ++i) {
                row[i] = (r.perforationPressure[i] - base.perforationPressure[i]) / pert.size + dNodeDps[i] * dpsDd;
            }
            result.push_back(std::move(row));
        }
        return result;
    }

    /*!
     * \brief March up from the bottom with the given pressure at the deepest perforation.
     */
    Result producerAt(Scalar bottomPressure, bool record) const
    {
        ++*marchCount_;
        Result r;
        r.perforationRate.assign(perforations_.size(), 0.0);
        r.perforationEnthalpy.assign(perforations_.size(), 0.0);
        r.perforationPressure.assign(perforations_.size(), 0.0);
        std::vector<ProfilePoint> prof;

        Scalar p = bottomPressure;
        Scalar h = 0.0;
        Composition z{};
        Scalar mdot = 0.0;
        bool hDefined = false;
        const int last = static_cast<int>(nodes_.size()) - 1;
        for (int k = last; k >= 0; --k) {
            const Node& node = nodes_[k];
            if (k < last) {
                const Node& below = nodes_[k + 1];
                march_(p, h, z, mdot, below, node, below.diameter, below.roughness, record ? &prof : nullptr);
            }
            if (node.perforation >= 0) {
                const auto& perf = perforations_[node.perforation];
                const Scalar m = std::max(0.0, perf.productivity * (perf.pressure - p));
                if (!hDefined) {
                    h = perf.inflowEnthalpy;   // nominal until there is flow
                    z = perf.inflowComposition;
                    hDefined = true;
                }
                if (m > 0.0) {
                    if (mdot > 0.0) {
                        h = (mdot * h + m * perf.inflowEnthalpy) / (mdot + m);
                        for (std::size_t c = 0; c < numComponents; ++c) {
                            z[c] = (mdot * z[c] + m * perf.inflowComposition[c]) / (mdot + m);
                        }
                    }
                    else {
                        h = perf.inflowEnthalpy;
                        z = perf.inflowComposition;
                    }
                    mdot += m;
                }
                r.perforationRate[node.perforation] = m;
                r.perforationEnthalpy[node.perforation] = perf.inflowEnthalpy;
                r.perforationPressure[node.perforation] = p;
            }
            if (record) {
                Scalar sinTheta, direction;
                if (k > 0) { orientation_(nodes_[k], nodes_[k - 1], sinTheta, direction); }
                else { orientation_(nodes_[1], nodes_[0], sinTheta, direction); }
                prof.push_back(point_(node, p, h, z, mdot, sinTheta, direction));
            }
        }
        r.flowing = mdot > 0.0;
        r.bottomHolePressure = bottomPressure;
        r.surfacePressure = p;
        r.massFlow = mdot;
        r.surfaceEnthalpy = h;
        const auto s = Fluid::state(p, h, z);
        r.surfaceTemperature = s.temperature;
        if (record) {
            r.profile.assign(prof.rbegin(), prof.rend());   // from the surface down
            findPhaseChange_(r, /*producerDirection=*/true);
        }
        return r;
    }

    /*!
     * \brief March down from the surface with the given surface pressure.
     */
    Result injectorAt(Scalar surfacePressure, Scalar volumeRate, Scalar temperature, const Stream& stream,
                      bool record) const
    {
        ++*marchCount_;
        Result r;
        r.perforationRate.assign(perforations_.size(), 0.0);
        r.perforationEnthalpy.assign(perforations_.size(), 0.0);
        r.perforationPressure.assign(perforations_.size(), 0.0);
        std::vector<ProfilePoint> prof;

        Scalar p = surfacePressure;
        Scalar h = Fluid::injectionEnthalpy(p, temperature, stream);
        const Composition z = Fluid::injectionComposition(stream);
        const auto s0 = Fluid::state(p, h, z);
        Scalar mdot = volumeRate * s0.density;
        r.massFlow = mdot;
        r.surfaceEnthalpy = h;
        r.surfaceTemperature = s0.temperature;
        r.surfacePressure = surfacePressure;

        for (std::size_t k = 0; k < nodes_.size(); ++k) {
            const Node& node = nodes_[k];
            if (k > 0) {
                const Node& above = nodes_[k - 1];
                march_(p, h, z, mdot, above, node, node.diameter, node.roughness, record ? &prof : nullptr);
            }
            if (record) {
                // the flow that arrives at the node
                Scalar sinTheta, direction;
                if (k > 0) { orientation_(nodes_[k - 1], nodes_[k], sinTheta, direction); }
                else { orientation_(nodes_[0], nodes_[1], sinTheta, direction); }
                prof.push_back(point_(node, p, h, z, mdot, sinTheta, direction));
            }
            if (node.perforation >= 0) {
                const auto& perf = perforations_[node.perforation];
                const auto s = Fluid::state(p, h, z);
                const Scalar out = perf.injectivity * s.density * std::max(0.0, p - perf.pressure);
                r.perforationRate[node.perforation] = out;
                r.perforationEnthalpy[node.perforation] = h;
                r.perforationPressure[node.perforation] = p;
                mdot -= out;
                if (mdot <= 0.0 && k + 1 < nodes_.size()) {
                    break;   // everything is gone before the bottom: the pressure is too high
                }
            }
        }
        r.flowing = true;
        r.remainingMassFlow = mdot;
        r.bottomHolePressure = p;
        if (record) {
            r.profile = std::move(prof);
            findPhaseChange_(r, /*producerDirection=*/false);
        }
        return r;
    }

private:
    struct Node
    {
        Scalar md{0}, tvd{0};
        Scalar diameter{0}, roughness{0};   //!< of the interval above the node
        int perforation{-1};
    };

    static constexpr Scalar huge_ = 1.0e30;

    //! The mixture properties of a flowing state
    struct Mix
    {
        Scalar alpha;      //!< vapour volume fraction in the pipe
        Scalar rho;        //!< in-situ mixture density
        Scalar kinetic;    //!< kinetic energy per mass
        Scalar momentum;   //!< momentum flux G u [Pa]
        Scalar velocity;   //!< total superficial velocity
    };

    /*!
     * \param sinTheta  sine of the inclination from the horizontal
     * \param direction +1 if the flow goes up, -1 if it goes down (irrelevant for horizontal flow)
     */
    Mix mix_(const typename Fluid::State& s, Scalar G, Scalar sinTheta, Scalar direction) const
    {
        Mix m;
        const Scalar x = s.quality;
        const Scalar rhoL = s.densityLiquid;
        const Scalar rhoG = s.densityVapor;
        const Scalar jG = x * G / rhoG;
        const Scalar jL = (1.0 - x) * G / rhoL;
        const Scalar j = jG + jL;
        m.velocity = j;
        if (x <= 0.0) {
            m.alpha = 0.0;
        }
        else if (x >= 1.0) {
            m.alpha = 1.0;
        }
        else if (!(j > 0.0)) {
            // no flow: the equilibrium state without slip; with drift flux the gas rises out, which is the limit alpha -> 0
            m.alpha = options_.slip.driftFlux ? 0.0 : s.voidFraction;
        }
        else if (!options_.slip.driftFlux) {
            m.alpha = jG / j;
        }
        else {
            const Scalar beta = jG / j;
            const Scalar c0 = 1.0 + (options_.slip.distributionParameter - 1.0) * (1.0 - beta);
            const Scalar sigma = std::max(s.surfaceTension, 0.0);
            const Scalar ud = options_.slip.driftCoefficient
                * std::pow(gravity * sigma * std::max(rhoL - rhoG, 0.0) / (rhoL * rhoL), 0.25)
                * sinTheta * direction * std::sqrt(1.0 - beta);
            m.alpha = std::clamp(jG / std::max(c0 * j + ud, jG), 0.0, 1.0);
        }
        m.rho = m.alpha * rhoG + (1.0 - m.alpha) * rhoL;
        const Scalar uG = m.alpha > 0.0 ? jG / m.alpha : 0.0;
        const Scalar uL = m.alpha < 1.0 ? jL / (1.0 - m.alpha) : 0.0;
        m.kinetic = 0.5 * (x * uG * uG + (1.0 - x) * uL * uL);
        m.momentum = G * (x * uG + (1.0 - x) * uL);
        return m;
    }

    Scalar ambientTemperature_(Scalar tvd) const
    { return options_.heatLoss.surfaceTemperature + options_.heatLoss.gradient * tvd; }

    //! sine of the inclination and the flow direction of an interval, from a to b
    static void orientation_(const Node& a, const Node& b, Scalar& sinTheta, Scalar& direction)
    {
        const Scalar length = std::abs(b.md - a.md);
        const Scalar rise = -(b.tvd - a.tvd);
        sinTheta = length > 0.0 ? std::min(1.0, std::abs(rise) / length) : 1.0;
        direction = rise >= 0.0 ? 1.0 : -1.0;
    }

    ProfilePoint point_(const Node& node, Scalar p, Scalar h, const Composition& z, Scalar mdot, Scalar sinTheta,
                        Scalar direction) const
    {
        return pointAt_(node.md, node.tvd, p, h, z, mdot, node.diameter, sinTheta, direction);
    }

    ProfilePoint pointAt_(Scalar md, Scalar tvd, Scalar p, Scalar h, const Composition& z, Scalar mdot,
                          Scalar diameter, Scalar sinTheta, Scalar direction) const
    {
        const auto s = Fluid::state(p, h, z);
        const Scalar area = 0.25 * M_PI * diameter * diameter;
        const Mix m = mix_(s, mdot / area, sinTheta, direction);
        ProfilePoint pt;
        pt.md = md;
        pt.tvd = tvd;
        pt.pressure = p;
        pt.temperature = s.temperature;
        pt.enthalpy = h;
        pt.massFlow = mdot;
        pt.steamQuality = s.quality;
        pt.voidFraction = m.alpha;
        pt.density = m.rho;
        pt.velocity = m.velocity;
        pt.kineticEnergy = m.kinetic;
        pt.heatLossRate = mdot > 0.0 ? options_.heatLoss.coefficient * (s.temperature - ambientTemperature_(tvd)) : 0.0;
        pt.region = s.region;
        pt.composition.assign(z.begin(), z.end());
        return pt;
    }

    //! Churchill (1977) Darcy friction factor
    static Scalar frictionFactor_(Scalar Re, Scalar relativeRoughness)
    {
        Re = std::max(Re, 1.0);
        const Scalar a = std::pow(2.457 * std::log(1.0 / (std::pow(7.0 / Re, 0.9) + 0.27 * relativeRoughness)), 16);
        const Scalar b = std::pow(37530.0 / Re, 16);
        return 8.0 * std::pow(std::pow(8.0 / Re, 12) + 1.0 / std::pow(a + b, 1.5), 1.0 / 12.0);
    }

    /*!
     * Advances the state from node a to node b in the direction of the flow. The elevation
     * increases if b is above a. Steps of at most maxStep_ are taken, each solved by iteration
     * for the (p, h) at its end.
     */
    void march_(Scalar& p, Scalar& h, const Composition& z, Scalar mdot, const Node& a, const Node& b, Scalar diameter,
                Scalar roughness, std::vector<ProfilePoint>* prof) const
    {
        const Scalar length = std::abs(b.md - a.md);
        const int n = std::max(1, static_cast<int>(std::ceil(length / maxStep_ - 1e-9)));
        const Scalar dl = length / n;
        const Scalar rise = -(b.tvd - a.tvd) / n;   // elevation gain per step
        const Scalar area = 0.25 * M_PI * diameter * diameter;
        const Scalar G = mdot / area;
        Scalar sinTheta, direction;
        orientation_(a, b, sinTheta, direction);
        const Scalar ua = options_.heatLoss.coefficient;

        for (int step = 1; step <= n; ++step) {
            const Scalar tvdEnd = a.tvd + (b.tvd - a.tvd) * step / n;
            const auto sa = Fluid::state(p, h, z);
            const Mix ma = mix_(sa, G, sinTheta, direction);
            Scalar pb = p;
            Scalar hb = h;
            if (mdot <= 0.0) {
                // no flow: a static column, in thermal equilibrium with the formation if the well loses heat.
                // Without heat loss it is the limit of a vanishing flow: the enthalpy follows the potential energy
                // (a jump at zero flow would give the well a spurious solution there).
                if (ua <= 0.0) {
                    hb = h - gravity * rise;
                }
                for (int it = 0; it < 10; ++it) {
                    if (ua > 0.0) {
                        hb = Fluid::equilibriumEnthalpy(pb, ambientTemperature_(tvdEnd), z);
                    }
                    const auto sb = Fluid::state(pb, hb, z);
                    pb = p - 0.5 * (ma.rho + sb.density) * gravity * rise;
                }
            }
            else {
                bool converged = false;
                for (int it = 0; it < 200 && !converged; ++it) {
                    const auto sb = Fluid::state(pb, hb, z);
                    const Mix mb = mix_(sb, G, sinTheta, direction);
                    const Scalar rhoNoSlip = 0.5 * (sa.density + sb.density);
                    const Scalar mu = 0.5 * (sa.viscosity + sb.viscosity);
                    const Scalar f = frictionFactor_(G * diameter / mu, roughness / diameter);
                    const Scalar dpFriction = f * dl / diameter * G * G / (2.0 * rhoNoSlip);
                    const Scalar dpGravity = 0.5 * (ma.rho + mb.rho) * gravity * rise;
                    const Scalar dpAcceleration = mb.momentum - ma.momentum;
                    const Scalar pNew = p - dpGravity - dpFriction - dpAcceleration;
                    // heat loss, implicit in the temperature at the end of the step
                    const Scalar heatLoss = ua * dl * (sb.temperature - ambientTemperature_(tvdEnd));
                    const Scalar hNew = h - gravity * rise + (ma.kinetic - mb.kinetic) - heatLoss / mdot;
                    if (!(pNew > 0.0) || !std::isfinite(hNew)) {
                        throw std::runtime_error("WellBore: the pressure falls to zero (choked flow?)");
                    }
                    converged = std::abs(pNew - pb) < 1e-11 * p && std::abs(hNew - hb) < 1e-9 * std::abs(h) + 1e-7;
                    Scalar omega = it < 30 ? 1.0 : 0.5;
                    if (ua > 0.0) {
                        // Newton on the stiff temperature term: d(heat loss/mdot)/dh = ua dl / (mdot cp)
                        const Scalar dh = 1.0e3;
                        const Scalar dTdh = (Fluid::state(pb, hb + dh, z).temperature - sb.temperature) / dh;
                        omega /= 1.0 + ua * dl / mdot * std::max(dTdh, 0.0);
                    }
                    pb += (pNew - pb) * (it < 30 ? 1.0 : 0.5);
                    hb += omega * (hNew - hb);
                }
                if (!converged) {
                    throw std::runtime_error("WellBore: the step did not converge");
                }
            }
            p = pb;
            h = hb;
            if (prof && step < n) {
                const Scalar f = static_cast<Scalar>(step) / n;
                prof->push_back(pointAt_(a.md + (b.md - a.md) * f, tvdEnd, p, h, z, mdot, diameter, sinTheta, direction));
            }
        }
    }

    //! The first change of the phase region in the direction of the flow, located by the margin to the boundary
    void findPhaseChange_(Result& r, bool producerDirection) const
    {
        // the profile is from the surface down: a producer flows up (look from the bottom), an injector down
        const auto& prof = r.profile;
        if (prof.size() < 2) {
            return;
        }
        const auto locate = [&](std::size_t i, std::size_t j) {
            // between the points i (upstream) and j with different regions
            Scalar w = 0.5;   // the middle of the step if the fluid has no margin to the boundary
            if constexpr (requires(const Composition& zc) { Fluid::phaseChangeMargin(Scalar(0), Scalar(0), zc, 0); }) {
                const int boundary = std::min(prof[i].region, prof[j].region);
                Composition zi, zj;
                std::copy(prof[i].composition.begin(), prof[i].composition.end(), zi.begin());
                std::copy(prof[j].composition.begin(), prof[j].composition.end(), zj.begin());
                const Scalar mi = Fluid::phaseChangeMargin(prof[i].pressure, prof[i].enthalpy, zi, boundary);
                const Scalar mj = Fluid::phaseChangeMargin(prof[j].pressure, prof[j].enthalpy, zj, boundary);
                w = std::abs(mi - mj) > 0.0 ? std::clamp(mi / (mi - mj), Scalar(0), Scalar(1)) : 0.5;
            }
            r.phaseChangeMd = prof[i].md + w * (prof[j].md - prof[i].md);
            r.phaseChangePressure = prof[i].pressure + w * (prof[j].pressure - prof[i].pressure);
        };
        if (producerDirection) {
            for (std::size_t i = prof.size() - 1; i > 0; --i) {
                if (prof[i].region != prof[i - 1].region) {
                    locate(i, i - 1);
                    return;
                }
            }
        }
        else {
            for (std::size_t i = 0; i + 1 < prof.size(); ++i) {
                if (prof[i].region != prof[i + 1].region) {
                    locate(i, i + 1);
                    return;
                }
            }
        }
    }

    /*!
     * Secant iteration for a root of f from a guess, where f increases through zero at the root. The root is only
     * accepted if the slope is positive (the stable root of a producer, the root of the negated remaining flow of an
     * injector) and the pressure stays within (lo, hi); the caller checks the residual. Values that are not
     * finite or are sentinels of invalid states give up.
     */
    template <class F>
    static std::optional<Scalar> secant_(const F& f, Scalar x0, Scalar lo, Scalar hi, Scalar xtol)
    {
        constexpr Scalar invalid = 1.0e20;
        Scalar f0 = f(x0);
        if (!(std::abs(f0) < invalid)) {
            return std::nullopt;
        }
        // the second point: towards the root of an increasing function
        Scalar x1 = f0 > 0.0 ? x0 * (1.0 - 1.0e-4) : x0 * (1.0 + 1.0e-4);
        if (x1 <= lo || x1 >= hi) {
            return std::nullopt;
        }
        Scalar f1 = f(x1);
        for (int it = 0; it < 12; ++it) {
            if (!(std::abs(f1) < invalid)) {
                return std::nullopt;
            }
            const Scalar slope = (f1 - f0) / (x1 - x0);
            if (!(slope > 0.0)) {
                return std::nullopt;
            }
            const Scalar x2 = x1 - f1 / slope;
            if (!(x2 > lo && x2 < hi)) {
                return std::nullopt;
            }
            if (std::abs(x2 - x1) < xtol) {
                return x2;
            }
            x0 = x1;
            f0 = f1;
            x1 = x2;
            f1 = f(x1);
        }
        return std::nullopt;
    }

    //! The stable root of the producer from a guess: f(pb) = pressure at the surface - head pressure rises with pb.
    template <class F>
    std::optional<Scalar> warmProducer_(const F& f, Scalar guess, Scalar pHigh, Scalar head) const
    {
        const Scalar pMin = Fluid::minPressure();
        // secant iteration: a few marches if the guess is good; its result is checked by the delivered pressure
        if (const auto pb = secant_(f, guess, pMin, pHigh, 1e-10 * pHigh)) {
            if (std::abs(f(*pb)) < 1e-6 * head) {
                return pb;
            }
        }
        Scalar step = 1.0e-5 * guess;
        Scalar a = guess;       // the end of the bracket that was evaluated last
        Scalar fa = f(a);
        const bool above = fa > 0.0;   // the root is below the guess if f is positive there
        for (int k = 0; k < 14; ++k) {
            const Scalar b = above ? a - step : std::min(a + step, pHigh);
            if (b <= pMin) {
                return std::nullopt;
            }
            const Scalar fb = f(b);
            if (above ? fb <= 0.0 : fb > 0.0) {
                // the sign changes between a and b: f > 0 at the higher pressure
                const Scalar lo = above ? b : a, hi = above ? a : b;
                const Scalar flo = above ? fb : fa, fhi = above ? fa : fb;
                const Scalar pb = brent_(f, lo, hi, flo, fhi, 1e-10 * pHigh);
                return std::abs(f(pb)) < 1e-6 * head ? std::optional<Scalar>(pb) : std::nullopt;
            }
            if (!above && b >= pHigh) {
                return std::nullopt;   // f stays non-positive up to the pressure of no flow: nothing flows
            }
            a = b;
            fa = fb;
            step *= 4.0;
        }
        return std::nullopt;
    }

    //! The root of the remaining flow of an injector near a guess; the remaining flow falls through zero.
    //! A pressure at which the march is not valid is read as beyond the root in the direction of the search.
    template <class R>
    std::optional<Scalar> warmInjector_(const R& remaining, Scalar guess) const
    {
        const Scalar pMin = Fluid::minPressure();
        const Scalar pMax = Fluid::maxPressure();
        const auto fGuess = remaining(guess);
        if (!fGuess) {
            return std::nullopt;
        }
        {
            // secant iteration on the negated remaining flow (increasing); the caller checks the remaining flow
            const auto g = [&](Scalar ps) { const auto w = remaining(ps); return w ? -*w : (ps > guess ? 1.0e30 : -1.0e30); };
            if (const auto ps = secant_(g, guess, pMin, pMax, 1e-10 * guess)) {
                return ps;
            }
        }
        Scalar a = guess;
        Scalar fa = *fGuess;
        Scalar step = 1.0e-5 * guess;
        const bool upwards = fa > 0.0;   // too little taken up: raise the pressure
        for (int k = 0; k < 14; ++k) {
            const Scalar b = upwards ? a + step : a - step;
            if (b <= pMin || b >= pMax) {
                return std::nullopt;
            }
            const auto v = remaining(b);
            // invalid: too high a pressure when going up (negative remaining flow), too low when going down
            const Scalar fb = v ? *v : (upwards ? -huge_ : huge_);
            if (upwards ? fb <= 0.0 : fb > 0.0) {
                const Scalar lo = upwards ? a : b, hi = upwards ? b : a;
                const Scalar flo = upwards ? fa : fb, fhi = upwards ? fb : fa;
                // the remaining flow decreases with the pressure: use its negative so that it increases
                const auto g = [&](Scalar ps) { const auto w = remaining(ps); return w ? -*w : (ps > guess ? huge_ : -huge_); };
                return brent_(g, lo, hi, -flo, -fhi, 1e-10 * hi);
            }
            a = b;
            fa = fb;
            step *= 4.0;
        }
        return std::nullopt;
    }

    //! root of f in [a, b] (f(a), f(b) of opposite signs), Brent's method
    template <class F>
    static Scalar brent_(const F& f, Scalar a, Scalar b, Scalar fa, Scalar fb, Scalar tol)
    {
        if (fa == 0.0) {
            return a;
        }
        if (fb == 0.0) {
            return b;
        }
        Scalar c = a, fc = fa, d = b - a, e = d;
        for (int it = 0; it < 200; ++it) {
            if ((fb > 0.0) == (fc > 0.0)) {
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }
            if (std::abs(fc) < std::abs(fb)) {
                a = b; b = c; c = a;
                fa = fb; fb = fc; fc = fa;
            }
            const Scalar tol1 = 2.0 * std::numeric_limits<Scalar>::epsilon() * std::abs(b) + 0.5 * tol;
            const Scalar xm = 0.5 * (c - b);
            if (std::abs(xm) <= tol1 || fb == 0.0) {
                return b;
            }
            if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
                Scalar p, q;
                const Scalar s = fb / fa;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                }
                else {
                    const Scalar qq = fa / fc;
                    const Scalar r = fb / fc;
                    p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
                    q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0) {
                    q = -q;
                }
                p = std::abs(p);
                if (2.0 * p < std::min(3.0 * xm * q - std::abs(tol1 * q), std::abs(e * q))) {
                    e = d;
                    d = p / q;
                }
                else {
                    d = xm;
                    e = d;
                }
            }
            else {
                d = xm;
                e = d;
            }
            a = b;
            fa = fb;
            b += std::abs(d) > tol1 ? d : (xm > 0 ? tol1 : -tol1);
            fb = f(b);
        }
        return b;
    }

    std::vector<Node> nodes_;
    std::vector<Perforation> perforations_;
    Geometry geometry_;
    Options options_;
    Scalar maxStep_{10.0};
    //! shared with the wells made by withPerforations(), so that the cost of the sensitivities is counted
    std::shared_ptr<std::size_t> marchCount_{std::make_shared<std::size_t>(0)};
};

} // namespace Opm

#endif
