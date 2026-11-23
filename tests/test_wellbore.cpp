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

#define BOOST_TEST_MODULE WellBore
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/SteamPressureEnthalpyProvider.hpp>
#include <opm/material/constraintsolvers/TabulatedSteamPressureEnthalpyProvider.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/constraintsolvers/WaterPHFlash.hpp>
#include <opm/models/wellbore/SteamWellFluid.hpp>
#include <opm/models/wellbore/WellBore.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using Scalar = double;
using Tab = Opm::TabulatedSteamPressureEnthalpyProvider<Scalar>;
using Raw = Opm::SteamPressureEnthalpyProvider<Scalar>;
using Fluid = Opm::SteamWellFluid<Scalar, Tab>;
using Well = Opm::WellBore<Scalar, Fluid>;
using Flash = Opm::WaterPHFlash<Scalar>;
using Geometry = Opm::WellBoreGeometry<Scalar>;
using Perf = Opm::WellBorePerforation<Scalar>;

namespace {

constexpr Scalar bar = 1.0e5;
constexpr Scalar g = Well::gravity;

Scalar celsius(Scalar c) { return c + 273.15; }

void printProfile(const Opm::WellBoreResult<Scalar>& r, int every)
{
    std::cout << "   md[m]   p[bar]   T[C]    h[kJ/kg]  quality  void    rho[kg/m3]  v[m/s]  region\n";
    for (std::size_t i = 0; i < r.profile.size(); i += every) {
        const auto& pt = r.profile[i];
        std::cout << std::setw(8) << pt.md << std::setw(9) << std::setprecision(5) << pt.pressure / bar
                  << std::setw(8) << std::setprecision(4) << pt.temperature - 273.15 << std::setw(11)
                  << std::setprecision(6) << pt.enthalpy / 1e3 << std::setw(9) << std::setprecision(3)
                  << pt.steamQuality << std::setw(8) << pt.voidFraction << std::setw(11) << std::setprecision(5)
                  << pt.density << std::setw(9) << std::setprecision(3) << pt.velocity << std::setw(4) << pt.region
                  << "\n";
    }
}

// total energy per mass, h + v^2/2 + g z (z up)
Scalar totalEnergy(const Opm::WellBoreProfilePoint<Scalar>& pt)
{ return pt.enthalpy + pt.kineticEnergy - g * pt.tvd; }

// hydrostatic integration from the surface down (RK4) at the enthalpy hTop + g*tvd (adiabatic, no flow)
Scalar hydrostaticPressure(Scalar pTop, Scalar hTop, Scalar tvd)
{
    const int n = 400;
    const Scalar dz = tvd / n;
    Scalar p = pTop;
    const auto rho = [&](Scalar pp, Scalar z) { return Fluid::state(pp, hTop + g * z).density; };
    for (int i = 0; i < n; ++i) {
        const Scalar z = i * dz;
        const Scalar k1 = rho(p, z) * g;
        const Scalar k2 = rho(p + 0.5 * dz * k1, z + 0.5 * dz) * g;
        const Scalar k3 = rho(p + 0.5 * dz * k2, z + 0.5 * dz) * g;
        const Scalar k4 = rho(p + dz * k3, z + dz) * g;
        p += dz / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    return p;
}

} // namespace

BOOST_AUTO_TEST_CASE(FluidAgainstRawProvider)
{
    using RawFluid = Opm::SteamWellFluid<Scalar, Raw>;
    for (const auto& [p, h] : {std::pair{50 * bar, 1.2e6}, std::pair{20 * bar, 1.8e6}, std::pair{5 * bar, 2.9e6},
                               std::pair{60 * bar, 4.0e5}}) {
        const auto a = Fluid::state(p, h);
        const auto b = RawFluid::state(p, h);
        BOOST_CHECK_CLOSE(a.density, b.density, 0.05);
        BOOST_CHECK_SMALL(a.temperature - b.temperature, 0.1);
        BOOST_CHECK_SMALL(a.quality - b.quality, 1e-4);
        BOOST_CHECK_EQUAL(a.region, b.region);
    }
}

BOOST_AUTO_TEST_CASE(SingleLiquidPhaseProducerIsHydrostatic)
{
    const Scalar H = 1000.0;
    const Scalar hRes = Raw::initialEnthalpy(150 * bar, celsius(80.0), 0.0);
    const Well well(Geometry::vertical(H, 0.3), {Perf{H, 150 * bar, 1e-5, hRes, 0.0}});
    const Scalar thp = 40 * bar;
    auto r = well.solveProducer(thp);
    BOOST_REQUIRE(r.flowing);
    BOOST_CHECK_CLOSE(r.surfacePressure, thp, 1e-6);
    BOOST_CHECK_CLOSE(r.perforationRate[0], 1e-5 * (150 * bar - r.bottomHolePressure), 1e-6);
    BOOST_CHECK_CLOSE(r.massFlow, r.perforationRate[0], 1e-9);
    // an independent integration of the hydrostatic column from the top: the difference is the friction
    const Scalar hTop = hRes - g * H;
    const Scalar pBottom = hydrostaticPressure(thp, hTop, H);
    std::cout << "single phase producer: rate " << r.massFlow << " kg/s, BHP " << r.bottomHolePressure / bar
              << " bar, hydrostatic " << pBottom / bar << " bar, friction etc "
              << (r.bottomHolePressure - pBottom) / bar << " bar\n";
    BOOST_CHECK_SMALL(r.bottomHolePressure - pBottom, 0.1 * bar);
    BOOST_CHECK_GT(r.bottomHolePressure, pBottom);   // friction is positive
    BOOST_CHECK(!r.phaseChangeMd);
}

BOOST_AUTO_TEST_CASE(HotWaterFlashesInTheWell)
{
    const Scalar H = 1500.0;
    const Scalar pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(300.0), 0.0);
    Geometry geo = Geometry::vertical(H, 0.2);
    geo.maxStepLength = 10.0;
    const Well well(geo, {Perf{H, pRes, 2e-6, hRes, 0.0}});
    const Scalar thp = 55 * bar;
    const auto r = well.solveProducer(thp);
    BOOST_REQUIRE(r.flowing);
    std::cout << "hot water producer: rate " << r.massFlow << " kg/s, BHP " << r.bottomHolePressure / bar
              << " bar, THP " << r.surfacePressure / bar << " bar, surface T " << r.surfaceTemperature - 273.15
              << " C, surface quality " << r.profile.front().steamQuality << "\n";
    BOOST_REQUIRE(r.phaseChangeMd);
    std::cout << "flashing at md " << *r.phaseChangeMd << " m, p = " << *r.phaseChangePressure / bar << " bar\n";
    printProfile(r, 15);

    BOOST_CHECK_CLOSE(r.surfacePressure, thp, 1e-6);
    BOOST_CHECK_GT(*r.phaseChangeMd, 0.0);
    BOOST_CHECK_LT(*r.phaseChangeMd, H);
    // the fluid is hot water at the bottom and steam and water at the top
    BOOST_CHECK_EQUAL(r.profile.back().region, 0);
    BOOST_CHECK_EQUAL(r.profile.front().region, 1);
    BOOST_CHECK_GT(r.profile.front().steamQuality, 0.05);

    // from the bottom up: the pressure falls, the steam quality and the velocity rise, the region changes once
    int changes = 0;
    for (std::size_t i = r.profile.size() - 1; i > 0; --i) {
        const auto& lo = r.profile[i];       // deeper
        const auto& up = r.profile[i - 1];
        BOOST_CHECK_LT(up.pressure, lo.pressure);
        BOOST_CHECK_GE(up.steamQuality, lo.steamQuality - 1e-12);
        BOOST_CHECK_GE(up.velocity, lo.velocity * (1.0 - 1e-9));
        changes += up.region != lo.region;
        // two-phase mixtures are at the saturation temperature
        if (up.region == 1) {
            BOOST_CHECK_SMALL(up.temperature - Flash::saturationTemperature(up.pressure), 0.05);
        }
    }
    BOOST_CHECK_EQUAL(changes, 1);
    // the phase change is where the profile switches region
    for (const auto& pt : r.profile) {
        if (pt.md > *r.phaseChangeMd + 10.0) {
            BOOST_CHECK_EQUAL(pt.region, 0);
        }
        if (pt.md < *r.phaseChangeMd - 10.0) {
            BOOST_CHECK_EQUAL(pt.region, 1);
        }
    }

    // energy: the total energy h + v^2/2 + g z is constant along the profile (one perforation)
    const Scalar E0 = totalEnergy(r.profile.back());
    Scalar maxDev = 0.0;
    for (const auto& pt : r.profile) {
        maxDev = std::max(maxDev, std::abs(totalEnergy(pt) - E0));
    }
    std::cout << "max deviation of the total energy " << maxDev << " J/kg of " << E0 << "\n";
    BOOST_CHECK_SMALL(maxDev, 2.0);

    // flashing makes the column lighter: less pressure drop than for the same fluid kept liquid
    const Scalar pLiquid = r.surfacePressure + 0.9 * 700.0 * g * H;
    BOOST_CHECK_LT(r.bottomHolePressure, pLiquid + 20 * bar);
}

BOOST_AUTO_TEST_CASE(PhaseChangeDepthFromTheHydrostaticColumn)
{
    // the depth where p = psat(T) of the rising water, from an independent integration from the bottom up
    const Scalar H = 1500.0;
    const Scalar pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(300.0), 0.0);
    const Well well(Geometry::vertical(H, 0.2), {Perf{H, pRes, 2e-6, hRes, 0.0}});
    const auto r = well.solveProducer(55 * bar);
    BOOST_REQUIRE(r.phaseChangeMd);
    // at the phase change h = hl(p): check with the flash that the interpolated point is on the saturation line
    // within the enthalpy changes over a step (g*10 m = 98 J/kg)
    const auto& prof = r.profile;
    std::size_t k = 0;
    while (prof[k].md < *r.phaseChangeMd) {
        ++k;
    }
    // prof[k] is below the change, prof[k-1] above
    const Scalar mLow = Fluid::phaseChangeMargin(prof[k].pressure, prof[k].enthalpy, 0);
    const Scalar mUp = Fluid::phaseChangeMargin(prof[k - 1].pressure, prof[k - 1].enthalpy, 0);
    BOOST_CHECK_LT(mLow, 0.0);
    BOOST_CHECK_GE(mUp, 0.0);
    BOOST_CHECK_SMALL(mUp - mLow, 1.0e4);   // a step crosses the line by less than 10 kJ/kg
}

BOOST_AUTO_TEST_CASE(ThreePerforationsMassAndEnergyBalance)
{
    const Scalar H = 1200.0;
    const Scalar hot = Raw::initialEnthalpy(130 * bar, celsius(280.0), 0.0);
    const Scalar hotter = Raw::initialEnthalpy(140 * bar, celsius(300.0), 0.0);
    const Scalar warm = Raw::initialEnthalpy(120 * bar, celsius(260.0), 0.0);
    std::vector<Perf> perfs = {Perf{500.0, 120 * bar, 2e-6, warm, 0.0},
                               Perf{900.0, 130 * bar, 3e-6, hot, 0.0},
                               Perf{H, 140 * bar, 4e-6, hotter, 0.0}};
    Geometry geo = Geometry::vertical(H, 0.2);
    const Well well(geo, perfs);
    const auto r = well.solveProducer(8 * bar);
    BOOST_REQUIRE(r.flowing);
    Scalar sum = 0.0, energy = 0.0;
    for (std::size_t i = 0; i < perfs.size(); ++i) {
        std::cout << "perforation at " << perfs[i].md << " m: " << r.perforationRate[i] << " kg/s\n";
        BOOST_CHECK_GT(r.perforationRate[i], 0.0);
        sum += r.perforationRate[i];
        energy += r.perforationRate[i] * (r.perforationEnthalpy[i] - g * perfs[i].md);
    }
    BOOST_CHECK_CLOSE(sum, r.massFlow, 1e-8);
    const auto& top = r.profile.front();
    const Scalar etop = r.massFlow * (top.enthalpy + top.kineticEnergy);
    std::cout << "energy flux at the surface " << etop / 1e6 << " MW, from the perforations " << energy / 1e6
              << " MW\n";
    BOOST_CHECK_SMALL(etop - energy, 1e-3 * std::abs(energy));
    // the perforation nodes appear in the profile
    for (const auto& perf : perfs) {
        const auto it = std::find_if(r.profile.begin(), r.profile.end(),
                                     [&](const auto& pt) { return std::abs(pt.md - perf.md) < 1e-6; });
        BOOST_CHECK(it != r.profile.end());
    }
}

BOOST_AUTO_TEST_CASE(DeadWell)
{
    const Scalar H = 1000.0;
    const Scalar hRes = Raw::initialEnthalpy(150 * bar, celsius(80.0), 0.0);
    const Well well(Geometry::vertical(H, 0.3), {Perf{H, 150 * bar, 1e-5, hRes, 0.0}});
    // the reservoir cannot lift the water (98 bar of column) to a tubing head pressure of 100 bar
    const auto r = well.solveProducer(100 * bar);
    BOOST_CHECK(!r.flowing);
    BOOST_CHECK_EQUAL(r.massFlow, 0.0);
}

BOOST_AUTO_TEST_CASE(LiquidInjector)
{
    const Scalar H = 1000.0;
    const Scalar pRes = 150 * bar;
    const Well well(Geometry::vertical(H, 0.15), {Perf{600.0, pRes, 0, 0, 1.0e-8}, Perf{H, pRes + 2 * bar, 0, 0, 1.0e-8}});
    const Scalar Q = 0.01;   // m3/s at the surface
    const auto r = well.solveInjector(Q, celsius(20.0), /*vapour=*/false);
    BOOST_REQUIRE(r.flowing);
    const Scalar rho0 = Fluid::state(r.surfacePressure, r.surfaceEnthalpy).density;
    std::cout << "liquid injector: surface pressure " << r.surfacePressure / bar << " bar, BHP "
              << r.bottomHolePressure / bar << " bar, mass flow " << r.massFlow << " kg/s\n";
    BOOST_CHECK_CLOSE(r.massFlow, Q * rho0, 1e-9);
    BOOST_CHECK_CLOSE(r.perforationRate[0] + r.perforationRate[1], r.massFlow, 1e-6);
    BOOST_CHECK_SMALL(r.remainingMassFlow, 1e-6 * r.massFlow);
    // the hydrostatic pressure at the bottom from the surface pressure, plus friction
    const Scalar pBottom = hydrostaticPressure(r.surfacePressure, r.surfaceEnthalpy, H);
    BOOST_CHECK_SMALL(r.bottomHolePressure - pBottom, 0.5 * bar);
    BOOST_CHECK_GT(r.bottomHolePressure, r.surfacePressure);
    // the surface pressure was found so that the perforations take up everything
    BOOST_CHECK(!r.phaseChangeMd);
}

BOOST_AUTO_TEST_CASE(InjectorWithTooMuchHeadSaysSo)
{
    // 98 bar of water column against a 80 bar reservoir: gravity alone injects more than the rate (no choke)
    const Scalar H = 1000.0;
    const Well well(Geometry::vertical(H, 0.15), {Perf{H, 80 * bar, 0, 0, 1.0e-8}});
    const auto before = well.marchCount();
    const auto r = well.solveInjector(0.01, celsius(20.0), false);
    // the scan stops once the remaining flow has been falling for a long stretch: far fewer marches than up to the maximum pressure
    BOOST_CHECK_LT(well.marchCount() - before, 120u);
    BOOST_CHECK(!r.rateReached);
    BOOST_CHECK_LT(r.surfacePressure, 5e3);                      // the lowest pressure at which liquid water of 20 C exists (2.3 kPa)
    BOOST_CHECK_GT(r.perforationRate[0], r.massFlow);              // it takes more than the rate
    BOOST_CHECK_LT(r.remainingMassFlow, 0.0);
    // a well that can hold the rate reaches it
    const Well other(Geometry::vertical(H, 0.15), {Perf{H, 150 * bar, 0, 0, 1.0e-8}});
    BOOST_CHECK(other.solveInjector(0.01, celsius(20.0), false).rateReached);
}

BOOST_AUTO_TEST_CASE(SteamInjector)
{
    const Scalar H = 800.0;
    const Scalar pRes = 40 * bar;
    const Well well(Geometry::vertical(H, 0.2), {Perf{H, pRes, 0, 0, 1.0e-6}});
    const Scalar Q = 0.2;   // m3/s of steam at the surface
    const auto r = well.solveInjector(Q, celsius(300.0), /*vapour=*/true);
    BOOST_REQUIRE(r.flowing);
    std::cout << "steam injector: surface pressure " << r.surfacePressure / bar << " bar, BHP "
              << r.bottomHolePressure / bar << " bar, mass flow " << r.massFlow << " kg/s\n";
    printProfile(r, 20);
    BOOST_CHECK_CLOSE(r.perforationRate[0], r.massFlow, 1e-6);
    // adiabatic: the steam stays dry, and the total energy is conserved
    const Scalar E0 = totalEnergy(r.profile.front());
    for (const auto& pt : r.profile) {
        BOOST_CHECK_EQUAL(pt.region, 2);
        BOOST_CHECK_SMALL(totalEnergy(pt) - E0, 2.0);
    }
    BOOST_CHECK_GT(r.bottomHolePressure, r.surfacePressure);
}

namespace {

using Options = Opm::WellBoreOptions<Scalar>;

// integral of the heat loss rate along the profile (trapezoid in the measured depth) [W]
Scalar integratedHeatLoss(const Opm::WellBoreResult<Scalar>& r)
{
    Scalar sum = 0.0;
    for (std::size_t i = 0; i + 1 < r.profile.size(); ++i) {
        sum += 0.5 * (r.profile[i].heatLossRate + r.profile[i + 1].heatLossRate)
            * (r.profile[i + 1].md - r.profile[i].md);
    }
    return sum;
}

// the drift-flux void fraction of a profile point, computed from its state
Scalar driftFluxVoid(const Opm::WellBoreProfilePoint<Scalar>& pt, Scalar diameter, Scalar c0Param, Scalar direction)
{
    const auto s = Fluid::state(pt.pressure, pt.enthalpy);
    const Scalar G = pt.massFlow / (0.25 * M_PI * diameter * diameter);
    const Scalar jG = s.quality * G / s.densityVapor;
    const Scalar jL = (1.0 - s.quality) * G / s.densityLiquid;
    const Scalar j = jG + jL;
    const Scalar beta = jG / j;
    const Scalar c0 = 1.0 + (c0Param - 1.0) * (1.0 - beta);
    const Scalar ud = 1.53 * std::pow(g * s.surfaceTension * (s.densityLiquid - s.densityVapor)
                                      / (s.densityLiquid * s.densityLiquid), 0.25)
        * direction * std::sqrt(1.0 - beta);
    return std::min(1.0, jG / (c0 * j + ud));   // the gas cannot occupy more than the pipe
}

} // namespace

BOOST_AUTO_TEST_CASE(DeviatedAndHorizontalPaths)
{
    const Scalar hRes = Raw::initialEnthalpy(150 * bar, celsius(80.0), 0.0);
    const Scalar thp = 40 * bar;
    // the same true vertical depth with a longer path, and a horizontal section at the end
    Geometry dev;
    dev.points = {{0.0, 0.0}, {300.0, 300.0}, {1500.0, 1000.0}};
    dev.diameter = {0.3, 0.3};
    dev.roughness = {4.5e-5, 4.5e-5};
    Geometry hor;
    hor.points = {{0.0, 0.0}, {1000.0, 1000.0}, {2000.0, 1000.0}};
    hor.diameter = {0.3, 0.3};
    hor.roughness = {4.5e-5, 4.5e-5};
    const auto vert = Well(Geometry::vertical(1000.0, 0.3), {Perf{1000.0, 150 * bar, 1e-5, hRes, 0.0}}).solveProducer(thp);
    const auto rd = Well(dev, {Perf{1500.0, 150 * bar, 1e-5, hRes, 0.0}}).solveProducer(thp);
    const auto rh = Well(hor, {Perf{2000.0, 150 * bar, 1e-5, hRes, 0.0}}).solveProducer(thp);
    const Scalar hydro = hydrostaticPressure(thp, hRes - g * 1000.0, 1000.0);
    std::cout << "BHP vertical " << vert.bottomHolePressure / bar << " bar, deviated " << rd.bottomHolePressure / bar
              << " bar, with a horizontal section " << rh.bottomHolePressure / bar << " bar, hydrostatic (TVD 1000 m) "
              << hydro / bar << " bar\n";
    // the same hydrostatic head, and only the friction of the longer path on top
    BOOST_CHECK_GT(rd.bottomHolePressure, vert.bottomHolePressure);
    BOOST_CHECK_GT(rh.bottomHolePressure, rd.bottomHolePressure);
    BOOST_CHECK_SMALL(rd.bottomHolePressure - hydro, 0.1 * bar);
    BOOST_CHECK_SMALL(rh.bottomHolePressure - hydro, 0.15 * bar);
    // the pressure changes along the horizontal section by friction only
    for (const auto& pt : rh.profile) {
        if (pt.md > 1000.0) {
            BOOST_CHECK_SMALL(pt.pressure - rh.bottomHolePressure, 0.1 * bar);
            BOOST_CHECK_LE(pt.pressure, rh.bottomHolePressure);
        }
    }
}

BOOST_AUTO_TEST_CASE(HeatLossEnergyBalance)
{
    const Scalar H = 1500.0;
    const Scalar pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(300.0), 0.0);
    Options opt;
    opt.heatLoss.coefficient = 0.5;   // W/(m K)
    opt.heatLoss.surfaceTemperature = celsius(15.0);
    opt.heatLoss.gradient = 0.03;
    const std::vector<Perf> perf = {Perf{H, pRes, 2e-6, hRes, 0.0}};
    const Well adiabatic(Geometry::vertical(H, 0.2), perf);
    const Well losing(Geometry::vertical(H, 0.2), perf, opt);
    const auto ra = adiabatic.solveProducer(45 * bar);
    const auto rl = losing.solveProducer(45 * bar);
    BOOST_REQUIRE(ra.flowing && rl.flowing);
    BOOST_REQUIRE(ra.phaseChangeMd && rl.phaseChangeMd);
    std::cout << "adiabatic: " << ra.massFlow << " kg/s, surface " << ra.surfaceTemperature - 273.15 << " C, quality "
              << ra.profile.front().steamQuality << ", flash at " << *ra.phaseChangeMd << " m\n";
    std::cout << "heat loss: " << rl.massFlow << " kg/s, surface " << rl.surfaceTemperature - 273.15 << " C, quality "
              << rl.profile.front().steamQuality << ", flash at " << *rl.phaseChangeMd << " m, loss "
              << integratedHeatLoss(rl) / 1e3 << " kW\n";
    BOOST_CHECK_CLOSE(rl.surfacePressure, 45 * bar, 1e-6);
    BOOST_CHECK_LT(rl.profile.front().steamQuality, ra.profile.front().steamQuality);
    // both are two-phase at the surface: the temperature is the saturation temperature at the same pressure
    BOOST_CHECK_SMALL(rl.surfaceTemperature - ra.surfaceTemperature, 1e-6);
    BOOST_CHECK_GT(integratedHeatLoss(rl), 0.0);
    for (const auto* res : {&ra, &rl}) {
        int last = res->profile.back().region;
        for (std::size_t i = res->profile.size() - 1; i > 0; --i) {
            if (res->profile[i - 1].region != last) {
                std::cout << "  region " << last << " -> " << res->profile[i - 1].region << " between md "
                          << res->profile[i].md << " and " << res->profile[i - 1].md << "\n";
                last = res->profile[i - 1].region;
            }
        }
    }
    // the energy that leaves the fluid is the heat that was lost
    const auto& top = rl.profile.front();
    const auto& bottom = rl.profile.back();
    const Scalar dE = rl.massFlow * (totalEnergy(bottom) - totalEnergy(top));
    BOOST_CHECK_SMALL(dE - integratedHeatLoss(rl), 0.02 * integratedHeatLoss(rl));
    // zero coefficient is the adiabatic well
    Options none;
    const auto r0 = Well(Geometry::vertical(H, 0.2), perf, none).solveProducer(45 * bar);
    BOOST_CHECK_CLOSE(r0.bottomHolePressure, ra.bottomHolePressure, 1e-9);
}

BOOST_AUTO_TEST_CASE(SlowFlowCoolsToTheAmbientTemperature)
{
    const Scalar H = 1500.0;
    const Scalar pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(120.0), 0.0);
    Options opt;
    opt.heatLoss.coefficient = 30.0;
    opt.heatLoss.surfaceTemperature = celsius(15.0);
    opt.heatLoss.gradient = 0.03;
    const Well well(Geometry::vertical(H, 0.2), {Perf{H, pRes, 1e-7, hRes, 0.0}}, opt);
    // the cold column (147 bar) leaves little for the wellhead: the flow is too slow to keep the fluid warm
    const auto r = well.solveProducer(2 * bar);
    BOOST_REQUIRE(r.flowing);
    std::cout << "slow producer: " << r.massFlow << " kg/s, surface T " << r.surfaceTemperature - 273.15 << " C\n";
    // the heat loss is implicit: the fluid never gets colder than the formation, and reaches it
    for (const auto& pt : r.profile) {
        BOOST_CHECK_GE(pt.temperature, opt.heatLoss.surfaceTemperature + opt.heatLoss.gradient * pt.tvd - 1e-3);
    }
    BOOST_CHECK_SMALL(r.surfaceTemperature - opt.heatLoss.surfaceTemperature, 0.1);
}

BOOST_AUTO_TEST_CASE(SteamCondensesInTheInjector)
{
    const Scalar H = 800.0;
    const Scalar pRes = 40 * bar;
    Options opt;
    opt.heatLoss.coefficient = 15.0;
    opt.heatLoss.surfaceTemperature = celsius(15.0);
    opt.heatLoss.gradient = 0.03;
    const std::vector<Perf> perf = {Perf{H, pRes, 0, 0, 1.0e-6}};
    const Well well(Geometry::vertical(H, 0.2), perf, opt);
    const auto r = well.solveInjector(0.2, celsius(300.0), /*vapour=*/true);
    BOOST_REQUIRE(r.flowing);
    BOOST_REQUIRE(r.phaseChangeMd);
    std::cout << "steam injector with heat loss: surface " << r.surfacePressure / bar << " bar, BHP "
              << r.bottomHolePressure / bar << " bar, condensation starts at md " << *r.phaseChangeMd
              << " m, bottom quality " << r.profile.back().steamQuality << ", loss " << integratedHeatLoss(r) / 1e3
              << " kW\n";
    printProfile(r, 20);
    BOOST_CHECK_GT(*r.phaseChangeMd, 0.0);
    BOOST_CHECK_LT(*r.phaseChangeMd, H);
    BOOST_CHECK_LT(r.profile.back().steamQuality, 0.9);
    BOOST_CHECK_EQUAL(r.profile.front().region, 2);
    BOOST_CHECK_EQUAL(r.profile.back().region, 1);
    BOOST_CHECK_CLOSE(r.perforationRate[0], r.massFlow, 1e-6);
    // downward: the quality never increases, the temperature of the wet steam is the saturation temperature
    for (std::size_t i = 0; i + 1 < r.profile.size(); ++i) {
        BOOST_CHECK_LE(r.profile[i + 1].steamQuality, r.profile[i].steamQuality + 1e-12);
        if (r.profile[i + 1].region == 1) {
            BOOST_CHECK_SMALL(r.profile[i + 1].temperature - Flash::saturationTemperature(r.profile[i + 1].pressure),
                              0.05);
        }
    }
    // the energy lost by the fluid is the heat loss
    const Scalar dE = r.massFlow * (totalEnergy(r.profile.front()) - totalEnergy(r.profile.back()));
    BOOST_CHECK_SMALL(dE - integratedHeatLoss(r), 0.02 * integratedHeatLoss(r));
}

BOOST_AUTO_TEST_CASE(DriftFluxSlip)
{
    const Scalar H = 1500.0;
    const Scalar pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(300.0), 0.0);
    const std::vector<Perf> perf = {Perf{H, pRes, 2e-6, hRes, 0.0}};
    Options homogeneous;
    Options noSlip;   // drift flux with C0 = 1 and no drift velocity is the homogeneous model
    noSlip.slip.driftFlux = true;
    noSlip.slip.distributionParameter = 1.0;
    noSlip.slip.driftCoefficient = 0.0;
    Options slip;
    slip.slip.driftFlux = true;
    const auto rh = Well(Geometry::vertical(H, 0.2), perf, homogeneous).solveProducer(40 * bar);
    const auto rn = Well(Geometry::vertical(H, 0.2), perf, noSlip).solveProducer(40 * bar);
    const auto rs = Well(Geometry::vertical(H, 0.2), perf, slip).solveProducer(40 * bar);
    BOOST_REQUIRE(rh.flowing && rn.flowing && rs.flowing);
    BOOST_CHECK_CLOSE(rn.bottomHolePressure, rh.bottomHolePressure, 1e-8);
    BOOST_CHECK_CLOSE(rn.massFlow, rh.massFlow, 1e-8);
    std::cout << "producer, homogeneous: BHP " << rh.bottomHolePressure / bar << " bar, " << rh.massFlow
              << " kg/s, void at the surface " << rh.profile.front().voidFraction << "\n";
    std::cout << "producer, drift flux:  BHP " << rs.bottomHolePressure / bar << " bar, " << rs.massFlow
              << " kg/s, void at the surface " << rs.profile.front().voidFraction << "\n";
    // the gas rises faster than the liquid: less gas in the pipe and a heavier column
    BOOST_CHECK_GT(rs.bottomHolePressure, rh.bottomHolePressure);
    BOOST_CHECK_LT(rs.profile.front().voidFraction,
                   Fluid::state(rs.profile.front().pressure, rs.profile.front().enthalpy).voidFraction);
    // every two-phase point follows the drift-flux relation
    int twoPhase = 0;
    for (const auto& pt : rs.profile) {
        if (pt.region == 1) {
            ++twoPhase;
            BOOST_CHECK_CLOSE(pt.voidFraction, driftFluxVoid(pt, 0.2, 1.2, +1.0), 1e-8);
            const Scalar rhoM = pt.voidFraction * Fluid::state(pt.pressure, pt.enthalpy).densityVapor
                + (1.0 - pt.voidFraction) * Fluid::state(pt.pressure, pt.enthalpy).densityLiquid;
            BOOST_CHECK_CLOSE(pt.density, rhoM, 1e-8);
        }
        else {
            BOOST_CHECK_EQUAL(pt.voidFraction, pt.region == 0 ? 0.0 : 1.0);
        }
    }
    BOOST_CHECK_GT(twoPhase, 10);
    // energy is still conserved with the kinetic energy of the two velocities
    const Scalar E0 = totalEnergy(rs.profile.back());
    for (const auto& pt : rs.profile) {
        BOOST_CHECK_SMALL(totalEnergy(pt) - E0, 2.0);
    }

    // the same relation in a downward flow: the drift velocity changes sign
    Options inj;
    inj.slip.driftFlux = true;
    inj.heatLoss.coefficient = 15.0;
    const Well injector(Geometry::vertical(800.0, 0.2), {Perf{800.0, 40 * bar, 0, 0, 1.0e-6}}, inj);
    const auto ri = injector.solveInjector(0.2, celsius(300.0), true);
    int wet = 0;
    for (const auto& pt : ri.profile) {
        if (pt.region == 1) {
            ++wet;
            BOOST_CHECK_CLOSE(pt.voidFraction, driftFluxVoid(pt, 0.2, 1.2, -1.0), 1e-8);
        }
    }
    BOOST_CHECK_GT(wet, 5);
}

BOOST_AUTO_TEST_CASE(ProducerNearTheFlowThreshold)
{
    // a reservoir pressure slightly below what is needed to lift the fluid gives a dead well, and the rate
    // grows continuously from zero as the pressure rises above the threshold (no jump at zero flow)
    const Scalar H = 1500.0;
    Scalar previous = 0.0;
    int flowing = 0;
    for (Scalar pr = 134.0; pr <= 140.0; pr += 0.25) {
        const Scalar hRes = Raw::initialEnthalpy(pr * bar, celsius(300.0), 0.0);
        const Well well(Geometry::vertical(H, 0.2), {Perf{H, pr * bar, 2e-4, hRes, 0.0}});
        const auto r = well.solveProducer(55 * bar);
        BOOST_CHECK_GE(r.massFlow, previous - 1e-9);
        if (r.flowing) {
            ++flowing;
            BOOST_CHECK_CLOSE(r.surfacePressure, 55 * bar, 1e-6);
        }
        else {
            BOOST_CHECK_EQUAL(r.massFlow, 0.0);
        }
        previous = r.massFlow;
    }
    BOOST_CHECK_GT(flowing, 5);
    BOOST_CHECK_LT(flowing, 24);
}

BOOST_AUTO_TEST_CASE(HotWaterInjectorNeedsASurfacePressureAboveSaturation)
{
    // A column of cold water (about 147 bar over 1500 m) is heavier than the 140 bar of the reservoir, so cold water
    // runs in by gravity. Water at 150 C is lighter and can only exist as a liquid above 4.76 bar: the solution is
    // a surface pressure above the saturation pressure, not the (invalid) vacuum state.
    const Scalar H = 1500.0, pRes = 140 * bar;
    const Well well(Geometry::vertical(H, 0.2), {Perf{H, pRes, 0, 0, 7.0e-8}});
    const Scalar Q = 0.0035;   // m3/s
    const auto hot = well.solveInjector(Q, celsius(150.0), false);
    BOOST_REQUIRE(hot.rateReached);
    std::cout << "150 C injector: surface pressure " << hot.surfacePressure / bar << " bar, BHP "
              << hot.bottomHolePressure / bar << " bar\n";
    BOOST_CHECK_GT(hot.surfacePressure, Opm::H2O<Scalar>::vaporPressure(celsius(150.0)));
    BOOST_CHECK_CLOSE(hot.surfaceTemperature, celsius(150.0), 1e-6);
    BOOST_CHECK_CLOSE(hot.perforationRate[0], hot.massFlow, 1e-6);
    BOOST_CHECK_GT(hot.bottomHolePressure, pRes);
    const auto cold = well.solveInjector(Q, celsius(20.0), false);
    BOOST_CHECK(!cold.rateReached);
}

BOOST_AUTO_TEST_CASE(WarmStartedRootSearches)
{
    // a hint saves marches and never changes the result
    const Scalar H = 1500.0, pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(300.0), 0.0);
    const std::vector<Perf> perf = {Perf{H, pRes, 2e-6, hRes, 0.0}};
    const Well well(Geometry::vertical(H, 0.2), perf);
    const Scalar thp = 55 * bar;

    const auto marches = [&](const Well& w, auto&& solve) {
        const auto before = w.marchCount();
        auto r = solve();
        return std::pair{std::move(r), w.marchCount() - before};
    };
    const auto [cold, nCold] = marches(well, [&] { return well.solveProducer(thp); });
    std::cout << "producer, cold search: " << nCold << " marches\n";
    for (const Scalar factor : {1.0, 1.0 + 1e-4, 1.0 - 1e-4, 1.0 + 2e-3, 1.0 - 2e-3, 1.02, 0.98, 0.9, 1.05}) {
        const auto [warm, n] = marches(well, [&] { return well.solveProducer(thp, cold.bottomHolePressure * factor); });
        std::cout << "  hint x" << factor << ": " << n << " marches\n";
        BOOST_CHECK_CLOSE(warm.bottomHolePressure, cold.bottomHolePressure, 1e-7);
        BOOST_CHECK_CLOSE(warm.massFlow, cold.massFlow, 1e-6);
        if (std::abs(factor - 1.0) <= 2e-3) {
            BOOST_CHECK_LE(n, nCold / 2);
        }
    }
    // a hint outside the pressures of the well, or a nonsense one, is ignored
    for (const Scalar hint : {1.0, 1e3, 1e9, 0.0}) {
        const auto warm = well.solveProducer(thp, hint);
        BOOST_CHECK_CLOSE(warm.bottomHolePressure, cold.bottomHolePressure, 1e-7);
    }

    // an injector
    const Well injector(Geometry::vertical(1000.0, 0.15), {Perf{600.0, 150 * bar, 0, 0, 1.0e-8}, Perf{1000.0, 152 * bar, 0, 0, 1.0e-8}});
    const auto [coldI, nColdI] = marches(injector, [&] { return injector.solveInjector(0.01, celsius(20.0), false); });
    std::cout << "injector, cold search: " << nColdI << " marches\n";
    for (const Scalar factor : {1.0, 1.0 + 1e-4, 1.0 - 1e-4, 1.0 + 2e-3, 1.0 - 2e-3, 1.1, 0.9, 3.0, 0.01}) {
        const auto [warm, n] = marches(injector, [&] { return injector.solveInjector(0.01, celsius(20.0), false, coldI.surfacePressure * factor); });
        std::cout << "  hint x" << factor << ": " << n << " marches\n";
        BOOST_CHECK_CLOSE(warm.surfacePressure, coldI.surfacePressure, 1e-6);
        BOOST_CHECK_CLOSE(warm.massFlow, coldI.massFlow, 1e-8);
        // 3 and 0.01 times the pressure are outside what the fluid can do: the full search takes over
        if (factor == 3.0 || factor == 0.01 || factor == 1.1) {
            continue;   // 1.1: the bottom hole would be above the pressure range of the tabulated steam
        }
        BOOST_CHECK_LE(n, nColdI / 10);
    }
}

BOOST_AUTO_TEST_CASE(WarmStartWithTwoRoots)
{
    // with heat loss a producer has two operating points; a hint in the region of no flow or near the unstable
    // one must still give the stable solution of the cold search
    const Scalar H = 1500.0, pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(300.0), 0.0);
    Options opt;
    opt.heatLoss.coefficient = 0.5;
    opt.heatLoss.surfaceTemperature = celsius(15.0);
    opt.heatLoss.gradient = 0.03;
    const Well well(Geometry::vertical(H, 0.2), {Perf{H, pRes, 2e-6, hRes, 0.0}}, opt);
    const auto cold = well.solveProducer(45 * bar);
    BOOST_REQUIRE(cold.flowing);
    for (const Scalar hint : {0.9999 * pRes, 0.999 * pRes, 0.99 * pRes, 0.95 * pRes, cold.bottomHolePressure, 0.9 * cold.bottomHolePressure}) {
        const auto warm = well.solveProducer(45 * bar, hint);
        BOOST_CHECK_CLOSE(warm.bottomHolePressure, cold.bottomHolePressure, 1e-6);
        BOOST_CHECK_CLOSE(warm.massFlow, cold.massFlow, 1e-5);
    }
    // a dead well stays dead
    BOOST_CHECK(!well.solveProducer(120 * bar, 0.97 * pRes).flowing);
}

BOOST_AUTO_TEST_CASE(SensitivitiesFromTheImplicitFunctionTheorem)
{
    // the derivative of the pressure of the fluid at the perforations wrt the descriptors of the perforations, from one
    // march per descriptor, against finite differences of the fully solved well
    const Scalar H = 1500.0;
    const Scalar h1 = Raw::initialEnthalpy(150 * bar, celsius(300.0), 0.0);
    const Scalar h2 = Raw::initialEnthalpy(145 * bar, celsius(290.0), 0.0);
    const std::vector<Perf> perfs = {Perf{900.0, 145 * bar, 3e-6, h2, 0.0}, Perf{H, 150 * bar, 2e-6, h1, 0.0}};
    const Well well(Geometry::vertical(H, 0.2), perfs);
    const Scalar thp = 40 * bar;
    const auto base = well.solveProducer(thp);
    std::cout << "two-perforation producer: BHP " << base.bottomHolePressure / bar << " bar, surface " << base.surfacePressure / bar << " bar, rates " << base.perforationRate[0] << " " << base.perforationRate[1] << "\n";
    BOOST_REQUIRE_CLOSE(base.surfacePressure, thp, 1e-6);
    BOOST_REQUIRE(base.flowing);
    BOOST_REQUIRE_EQUAL(base.perforationPressure.size(), 2u);

    struct Descriptor { const char* name; Scalar size; std::function<void(Perf&)> change; };
    const Scalar dp = 1.0e-5 * 150 * bar;
    const std::vector<Descriptor> descriptors = {
        {"pressure", dp, [&](Perf& p) { p.pressure += dp; }},
        {"productivity", 2e-9, [](Perf& p) { p.productivity += 2e-9; }},
        {"enthalpy", 1.0e3, [](Perf& p) { p.inflowEnthalpy += 1.0e3; }},
    };
    std::vector<Well::Perturbation> perturbations;
    for (std::size_t c = 0; c < 2; ++c) {
        for (const auto& d : descriptors) {
            perturbations.push_back({c, [&d](Perf& p) { d.change(p); }, d.size});
        }
    }
    const auto before = well.marchCount();
    const auto ift = well.producerSensitivities(base, perturbations);
    const auto marches = well.marchCount() - before;
    BOOST_REQUIRE(ift);
    std::cout << "producer: " << perturbations.size() << " perturbations of 2 perforations in " << marches << " marches\n";
    BOOST_CHECK_EQUAL(marches, perturbations.size() + 1);
    // finite differences of the fully solved wells; a sensitivity is compared with 2% of the largest one of its kind
    // (a perforation that takes no flow has none: the finite difference there is noise of the root finder)
    const auto compare = [&](const std::vector<Descriptor>& descs, const std::vector<Perf>& base_perfs, const Geometry& geo,
                             const auto& solveFor, const auto& pressuresOf, const std::vector<std::vector<Scalar>>& ift,
                             const std::vector<Scalar>& basePressures) {
        std::vector<std::vector<Scalar>> fd;
        std::vector<Scalar> scale(descs.size(), 0.0);
        for (std::size_t c = 0; c < 2; ++c) {
            for (std::size_t j = 0; j < descs.size(); ++j) {
                auto changed = base_perfs;
                descs[j].change(changed[c]);
                const auto p = pressuresOf(solveFor(Well(geo, changed)));
                std::vector<Scalar> row(2);
                for (std::size_t i = 0; i < 2; ++i) {
                    row[i] = (p[i] - basePressures[i]) / descs[j].size;
                    scale[j] = std::max(scale[j], std::abs(row[i]));
                }
                fd.push_back(row);
            }
        }
        for (std::size_t c = 0, k = 0; c < 2; ++c) {
            for (std::size_t j = 0; j < descs.size(); ++j, ++k) {
                for (std::size_t i = 0; i < 2; ++i) {
                    std::cout << "  perforation " << c << " " << descs[j].name << " -> pressure at " << i << ": IFT " << ift[k][i]
                              << ", finite difference " << fd[k][i] << "\n";
                    BOOST_CHECK_SMALL(ift[k][i] - fd[k][i], 0.02 * scale[j] + 1e-9);
                }
            }
        }
    };
    compare(descriptors, perfs, Geometry::vertical(H, 0.2),
            [&](const Well& w) { return w.solveProducer(thp, base.bottomHolePressure); },
            [](const auto& r) { return r.perforationPressure; }, *ift, base.perforationPressure);

    // an injector with two perforations
    const std::vector<Perf> iperfs = {Perf{600.0, 150 * bar, 0, 0, 1.0e-8}, Perf{1000.0, 152 * bar, 0, 0, 2.0e-8}};
    const Well injector(Geometry::vertical(1000.0, 0.15), iperfs);
    const auto ibase = injector.solveInjector(0.01, celsius(20.0), false);
    BOOST_REQUIRE(ibase.rateReached);
    const Scalar dpi = 1.0e-5 * 150 * bar;
    const std::vector<Descriptor> idescriptors = {
        {"pressure", dpi, [&](Perf& p) { p.pressure += dpi; }},
        {"injectivity", 1e-11, [](Perf& p) { p.injectivity += 1e-11; }},
    };
    std::vector<Well::Perturbation> iperturbations;
    for (std::size_t c = 0; c < 2; ++c) {
        for (const auto& d : idescriptors) {
            iperturbations.push_back({c, [&d](Perf& p) { d.change(p); }, d.size});
        }
    }
    const auto ibefore = injector.marchCount();
    const auto iift = injector.injectorSensitivities(ibase, 0.01, celsius(20.0), false, iperturbations);
    BOOST_REQUIRE(iift);
    BOOST_CHECK_EQUAL(injector.marchCount() - ibefore, iperturbations.size() + 1);
    std::cout << "injector:\n";
    compare(idescriptors, iperfs, Geometry::vertical(1000.0, 0.15),
            [&](const Well& w) { return w.solveInjector(0.01, celsius(20.0), false, ibase.surfacePressure); },
            [](const auto& r) { return r.perforationPressure; }, *iift, ibase.perforationPressure);
    // no sensitivities for a well that does not flow
    BOOST_CHECK(!well.producerSensitivities(Well(Geometry::vertical(H, 0.2), perfs).solveProducer(300 * bar), perturbations));
}

BOOST_AUTO_TEST_CASE(DeadWellWithSeveralPerforations)
{
    // the wellhead pressure cannot be delivered: nothing flows, even if the pressure of a shallower perforation is
    // above the pressure of the well there
    const Scalar H = 1500.0;
    const Scalar h1 = Raw::initialEnthalpy(150 * bar, celsius(300.0), 0.0);
    const Scalar h2 = Raw::initialEnthalpy(145 * bar, celsius(290.0), 0.0);
    const Well well(Geometry::vertical(H, 0.2), {Perf{900.0, 145 * bar, 3e-6, h2, 0.0}, Perf{H, 150 * bar, 2e-6, h1, 0.0}});
    const auto r = well.solveProducer(55 * bar);
    BOOST_CHECK(!r.flowing);
    BOOST_CHECK_EQUAL(r.massFlow, 0.0);
    BOOST_CHECK_EQUAL(r.perforationRate[0], 0.0);
    BOOST_CHECK_EQUAL(r.perforationRate[1], 0.0);
}

BOOST_AUTO_TEST_CASE(DriftFluxProducerNearTheZeroFlowJump)
{
    // With drift flux the surface pressure of a vanishing flow (about 42 bar) is far below that of the static column
    // (52 bar): f jumps at zero flow. A tubing head pressure between the two used to end in an exception, because the
    // bracket of the first crossing holds the jump and no root.
    const Scalar H = 1500.0;
    const Scalar pRes = 150 * bar;
    const Scalar hRes = Raw::initialEnthalpy(pRes, celsius(300.0), 0.0);
    Options drift;
    drift.slip.driftFlux = true;
    const Well well(Geometry::vertical(H, 0.2), {Perf{H, pRes, 2e-6, hRes, 0.0}}, drift);
    const auto r = well.solveProducer(45 * bar);
    BOOST_REQUIRE(r.flowing);
    BOOST_CHECK_CLOSE(r.surfacePressure, 45 * bar, 1e-6);
    BOOST_CHECK_GT(r.massFlow, 1.0);   // the operating point of the outflow curve, not the jump at zero flow
    // above the highest surface pressure that a flow can deliver (49.06 bar) the well is dead
    const auto dead = well.solveProducer(50 * bar);
    BOOST_CHECK(!dead.flowing);
    BOOST_CHECK_EQUAL(dead.massFlow, 0.0);
}
