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

#define BOOST_TEST_MODULE Co2WellBore
#include <boost/test/unit_test.hpp>

#include <opm/material/constraintsolvers/Co2BrinePressureEnthalpyProvider.hpp>
#include <opm/material/components/CO2PHTables.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/models/wellbore/Co2WellFluid.hpp>
#include <opm/models/wellbore/WellBore.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

using Scalar = double;
using Fluid = Opm::Co2WellFluid<Scalar>;
using Provider = Opm::Co2BrinePressureEnthalpyProvider<Scalar>;
using FS = Opm::Co2StoreFluidSystem<Scalar>;
using Well = Opm::WellBore<Scalar, Fluid>;
using Geometry = Opm::WellBoreGeometry<Scalar>;
using Perf = Opm::WellBorePerforation<Scalar, 2>;
using Eval1 = Opm::DenseAd::Evaluation<Scalar, 1>;

namespace {

constexpr Scalar bar = 1.0e5;
constexpr Scalar g = Well::gravity;

// solubility X and the saturated brine content Y at (T, p)
void equilibrium(Scalar T, Scalar p, Scalar& X, Scalar& Y)
{
    Eval1 Xe, Ye;
    FS::equilibriumMassFractions(Eval1(T), Eval1(p), Xe, Ye);
    X = Xe.value();
    Y = Ye.value();
}

Scalar totalEnergy(const Opm::WellBoreProfilePoint<Scalar>& pt)
{ return pt.enthalpy + pt.kineticEnergy - g * pt.tvd; }

void printProfile(const Opm::WellBoreResult<Scalar>& r, int every)
{
    std::cout << "   md[m]   p[bar]   T[C]    h[kJ/kg]  CO2-rich  void   rho[kg/m3]  CO2 frac  region\n";
    for (std::size_t i = 0; i < r.profile.size(); i += every) {
        const auto& pt = r.profile[i];
        std::cout << std::setw(8) << pt.md << std::setw(9) << std::setprecision(5) << pt.pressure / bar
                  << std::setw(8) << std::setprecision(4) << pt.temperature - 273.15 << std::setw(11)
                  << std::setprecision(6) << pt.enthalpy / 1e3 << std::setw(10) << std::setprecision(3)
                  << pt.steamQuality << std::setw(8) << pt.voidFraction << std::setw(11) << std::setprecision(5)
                  << pt.density << std::setw(10) << std::setprecision(4) << pt.composition[0] << std::setw(4)
                  << pt.region << "\n";
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(FlashRoundTrip)
{
    struct Case { Scalar p, T, Sw, Xfraction, Yfraction; };
    const std::vector<Case> cases = {
        {200 * bar, 330.0, 1.0, 0.5, 0.0},     // undersaturated brine
        {200 * bar, 330.0, 1.0, 0.999, 0.0},   // almost saturated brine
        {150 * bar, 340.0, 0.3, 0.0, 0.0},     // both phases, supercritical CO2
        {80 * bar, 310.0, 0.8, 0.0, 0.0},      // both phases
        {60 * bar, 285.0, 0.5, 0.0, 0.0},      // liquid CO2 and brine
        {100 * bar, 350.0, 0.0, 0.0, 0.5},     // CO2 with some water
        {50 * bar, 300.0, 0.0, 0.0, 0.0},      // dry CO2
    };
    for (const auto& c : cases) {
        Scalar X, Y;
        equilibrium(c.T, c.p, X, Y);
        Scalar h, s;
        int meaning;
        Provider::initialVariables(c.p, c.T, c.Sw, c.Xfraction * X, c.Yfraction * Y, h, s, meaning);
        const auto st = Provider::solve<Scalar>(c.p, h, s, meaning);
        // the overall CO2 mass fraction of the state
        Scalar co2 = 0.0, total = 0.0;
        for (unsigned a = 0; a < 3; ++a) {
            const Scalar m = st.saturation(a) * st.density(a);
            total += m;
            co2 += m * st.massFraction(a, FS::CO2Idx);
        }
        const Scalar Z = co2 / total;
        const auto f = Fluid::flash(c.p, h, {Z, 1.0 - Z});
        std::cout << "p " << c.p / bar << " bar, T " << c.T << " K: Z = " << Z << ", meaning " << meaning << " -> "
                  << f.meaning << ", s " << s << " -> " << f.s << "\n";
        BOOST_CHECK_EQUAL(f.meaning, meaning);
        BOOST_CHECK_SMALL(f.s - s, 1e-9);
        BOOST_CHECK_SMALL(f.state.temperature() - st.temperature(), 1e-6);
        for (unsigned a = 0; a < 3; ++a) {
            BOOST_CHECK_CLOSE(f.state.density(a), st.density(a), 1e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(WellFluidState)
{
    const Scalar p = 150 * bar, T = 340.0;
    Scalar X, Y;
    equilibrium(T, p, X, Y);
    Scalar h, s;
    int meaning;
    Provider::initialVariables(p, T, 0.3, 0.0, 0.0, h, s, meaning);
    const auto st = Provider::solve<Scalar>(p, h, s, meaning);
    Scalar co2 = 0.0, total = 0.0;
    for (unsigned a = 0; a < 3; ++a) {
        co2 += st.saturation(a) * st.density(a) * st.massFraction(a, FS::CO2Idx);
        total += st.saturation(a) * st.density(a);
    }
    const auto ws = Fluid::state(p, h, {co2 / total, 1.0 - co2 / total});
    // the CO2-rich phase is the non-aqueous part: its mass fraction, volume fraction and density
    BOOST_CHECK_CLOSE(ws.voidFraction, 1.0 - s, 1e-8);
    BOOST_CHECK_CLOSE(ws.aqueousSaturation, s, 1e-8);
    BOOST_CHECK_CLOSE(ws.densityLiquid, st.density(2), 1e-8);
    // mixture density = volume weighted
    const Scalar rho = s * st.density(2) + (1.0 - s) * ws.densityVapor;
    BOOST_CHECK_CLOSE(ws.density, rho, 1e-6);
    // gas and liquid CO2 are one pseudo phase: its density is between those of the two
    BOOST_CHECK_GT(ws.quality, 0.0);
    BOOST_CHECK_LT(ws.quality, 1.0);
}

BOOST_AUTO_TEST_CASE(DissolvedCo2ComesOutOfSolution)
{
    // brine with 90% of the CO2 that it can dissolve at the reservoir pressure: the CO2 exsolves
    // where the pressure of the rising brine falls to the bubble point pressure
    const Scalar H = 1500.0, pRes = 200 * bar, TRes = 340.0;
    Scalar Xeq, Yeq;
    equilibrium(TRes, pRes, Xeq, Yeq);
    const Scalar X = 0.9 * Xeq;
    Geometry geo = Geometry::vertical(H, 0.2);
    geo.maxStepLength = 25.0;
    const Well well(geo, {Fluid::perforation(H, pRes, TRes, 1.0, X, 0.0, 1.0e-12)});
    std::cout << "perforation: productivity " << well.perforations()[0].productivity << " kg/(s Pa), CO2 fraction "
              << well.perforations()[0].inflowComposition[0] << " (X = " << X << ")\n";
    // choose the productivity of the reservoir so that the rate is a few kg/s
    auto perf = well.perforations()[0];
    perf.productivity = 1e-5;
    const Well w(geo, {perf});
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = w.solveProducer(40 * bar);
    const auto t1 = std::chrono::steady_clock::now();
    BOOST_REQUIRE(r.flowing);
    std::cout << "brine producer: rate " << r.massFlow << " kg/s, BHP " << r.bottomHolePressure / bar << " bar, THP "
              << r.surfacePressure / bar << " bar, surface T " << r.surfaceTemperature - 273.15 << " C, solve "
              << std::chrono::duration<double>(t1 - t0).count() << " s\n";
    BOOST_REQUIRE(r.phaseChangeMd);
    std::cout << "CO2 exsolves at md " << *r.phaseChangeMd << " m, p = " << *r.phaseChangePressure / bar << " bar\n";
    printProfile(r, 6);

    BOOST_CHECK_CLOSE(r.surfacePressure, 40 * bar, 1e-6);
    // the composition does not change with one perforation
    for (const auto& pt : r.profile) {
        BOOST_CHECK_CLOSE(pt.composition[0], X, 1e-9);
    }
    // aqueous only below, two phases above; the region changes once
    int changes = 0;
    for (std::size_t i = r.profile.size() - 1; i > 0; --i) {
        changes += r.profile[i].region != r.profile[i - 1].region;
        BOOST_CHECK_GE(r.profile[i - 1].steamQuality, r.profile[i].steamQuality - 1e-12);
    }
    BOOST_CHECK_EQUAL(changes, 1);
    BOOST_CHECK_EQUAL(r.profile.back().steamQuality, 0.0);
    BOOST_CHECK_GT(r.profile.front().steamQuality, 0.0);
    // the bubble point: at the point below the change X is at most the solubility, above it is more
    for (std::size_t i = 0; i < r.profile.size(); ++i) {
        Scalar Xe, Ye;
        equilibrium(r.profile[i].temperature, r.profile[i].pressure, Xe, Ye);
        if (r.profile[i].steamQuality == 0.0) {
            BOOST_CHECK_LE(X, Xe * (1.0 + 1e-9));
        }
        else {
            BOOST_CHECK_GE(X, Xe * (1.0 - 1e-9));
        }
    }
    // the phase change is within a step (25 m) of where the solubility equals X
    const auto bubble = [&](const Opm::WellBoreProfilePoint<Scalar>& pt) {
        Scalar Xe, Ye;
        equilibrium(pt.temperature, pt.pressure, Xe, Ye);
        return Xe - X;
    };
    std::size_t k = 0;
    while (r.profile[k].md < *r.phaseChangeMd) {
        ++k;
    }
    // the solubility is above X in the point below the change and below X in the point above it
    BOOST_CHECK_GE(bubble(r.profile[k]), 0.0);
    BOOST_CHECK_LE(bubble(r.profile[k - 1]), 0.0);
    BOOST_CHECK_LE(r.profile[k].md - r.profile[k - 1].md, 25.0 + 1e-9);
    // energy
    const Scalar E0 = totalEnergy(r.profile.back());
    Scalar maxDev = 0.0;
    for (const auto& pt : r.profile) {
        maxDev = std::max(maxDev, std::abs(totalEnergy(pt) - E0));
    }
    std::cout << "max deviation of the total energy " << maxDev << " J/kg of " << E0 << "\n";
    BOOST_CHECK_SMALL(maxDev, 5.0);
}

BOOST_AUTO_TEST_CASE(PureCo2Injector)
{
    const Scalar H = 1000.0, pRes = 125 * bar;
    Geometry geo = Geometry::vertical(H, 0.15);
    geo.maxStepLength = 25.0;
    Perf perf;
    perf.md = H;
    perf.pressure = pRes;
    perf.injectivity = 2.0e-8;
    const Well well(geo, {perf});
    const Scalar T = 283.15;   // liquid CO2 at 10 C
    const Scalar Q = 0.02;     // m3/s at the surface pressure
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = well.solveInjector(Q, T, Fluid::Stream{1.0, 0.0});
    const auto t1 = std::chrono::steady_clock::now();
    BOOST_REQUIRE(r.flowing);
    std::cout << "CO2 injector: surface pressure " << r.surfacePressure / bar << " bar, BHP " << r.bottomHolePressure / bar
              << " bar, mass flow " << r.massFlow << " kg/s, bottom T " << r.profile.back().temperature - 273.15
              << " C, solve " << std::chrono::duration<double>(t1 - t0).count() << " s\n";
    printProfile(r, 8);
    const auto s0 = Fluid::state(r.surfacePressure, r.surfaceEnthalpy, {1.0, 0.0});
    BOOST_CHECK_CLOSE(r.massFlow, Q * s0.density, 1e-9);
    BOOST_CHECK_CLOSE(r.perforationRate[0], r.massFlow, 1e-6);
    BOOST_CHECK_CLOSE(r.surfaceTemperature, T, 1e-6);
    BOOST_CHECK_GT(r.bottomHolePressure, r.surfacePressure);
    // pure CO2 all the way: no brine, CO2 only, no phase change
    for (const auto& pt : r.profile) {
        BOOST_CHECK_CLOSE(pt.composition[0], 1.0, 1e-12);
        BOOST_CHECK_EQUAL(pt.steamQuality, 1.0);
    }
    // the liquid column weighs about rho g H
    const Scalar rhoMean = 0.5 * (r.profile.front().density + r.profile.back().density);
    BOOST_CHECK_SMALL(r.bottomHolePressure - r.surfacePressure - rhoMean * g * H, 0.5 * bar);
    // the CO2 warms by compression (the potential energy goes into the enthalpy): energy is conserved
    const Scalar E0 = totalEnergy(r.profile.front());
    for (const auto& pt : r.profile) {
        BOOST_CHECK_SMALL(totalEnergy(pt) - E0, 5.0);
    }
    BOOST_CHECK_GT(r.profile.back().temperature, r.profile.front().temperature);
}

BOOST_AUTO_TEST_CASE(BrineAndCo2PerforationsMix)
{
    const Scalar H = 1000.0;
    Scalar Xeq, Yeq;
    equilibrium(330.0, 130 * bar, Xeq, Yeq);
    auto brine = Fluid::perforation(H, 130 * bar, 330.0, 1.0, 0.8 * Xeq, 0.0, 1.0e-12);
    brine.productivity = 5.0e-6;
    auto co2 = Fluid::perforation(600.0, 95 * bar, 315.0, 0.0, 0.0, 0.0, 1.0e-12);
    co2.productivity = 2.0e-6;
    Geometry geo = Geometry::vertical(H, 0.2);
    geo.maxStepLength = 25.0;
    const Well well(geo, {co2, brine});
    const auto r = well.solveProducer(25 * bar);
    BOOST_REQUIRE(r.flowing);
    std::cout << "mixed producer: rate " << r.massFlow << " kg/s (CO2 perforation " << r.perforationRate[0]
              << ", brine perforation " << r.perforationRate[1] << "), BHP " << r.bottomHolePressure / bar
              << " bar, surface T " << r.surfaceTemperature - 273.15 << " C\n";
    printProfile(r, 6);
    BOOST_CHECK_GT(r.perforationRate[0], 0.0);
    BOOST_CHECK_GT(r.perforationRate[1], 0.0);
    BOOST_CHECK_CLOSE(r.perforationRate[0] + r.perforationRate[1], r.massFlow, 1e-8);

    // the composition of the mix and its energy flux
    const Scalar zTop = (r.perforationRate[0] * co2.inflowComposition[0] + r.perforationRate[1] * brine.inflowComposition[0])
        / r.massFlow;
    BOOST_CHECK_CLOSE(r.profile.front().composition[0], zTop, 1e-8);
    const auto& top = r.profile.front();
    const Scalar etop = r.massFlow * (top.enthalpy + top.kineticEnergy);
    const Scalar efeed = r.perforationRate[0] * (co2.inflowEnthalpy - g * 600.0)
        + r.perforationRate[1] * (brine.inflowEnthalpy - g * H);
    std::cout << "energy flux at the surface " << etop / 1e6 << " MW, from the perforations " << efeed / 1e6 << " MW\n";
    BOOST_CHECK_SMALL(etop - efeed, 2e-3 * std::abs(efeed));
    // above the CO2 perforation there is more CO2 than below it
    for (const auto& pt : r.profile) {
        if (pt.md < 600.0 - 1e-6) {
            BOOST_CHECK_CLOSE(pt.composition[0], zTop, 1e-8);
        }
        else if (pt.md > 600.0 + 1e-6) {
            BOOST_CHECK_CLOSE(pt.composition[0], brine.inflowComposition[0], 1e-8);
        }
    }
}

BOOST_AUTO_TEST_CASE(Co2BoilsInTheWell)
{
    // dense CO2 expands on its way up, the temperature falls and it enters the two-phase region of CO2
    const Scalar H = 1000.0;
    auto co2 = Fluid::perforation(H, 140 * bar, 300.5, 0.0, 0.0, 0.0, 1.0e-12);
    co2.productivity = 1.0e-5;
    Geometry geo = Geometry::vertical(H, 0.2);
    geo.maxStepLength = 25.0;
    const Well well(geo, {co2});
    const Scalar thp = 45 * bar;
    const auto r = well.solveProducer(thp);
    BOOST_REQUIRE(r.flowing);
    BOOST_REQUIRE(r.phaseChangeMd);
    std::cout << "CO2 producer: rate " << r.massFlow << " kg/s, BHP " << r.bottomHolePressure / bar << " bar, the CO2 starts to boil at md "
              << *r.phaseChangeMd << " m (" << *r.phaseChangePressure / bar << " bar), surface T "
              << r.surfaceTemperature - 273.15 << " C\n";
    BOOST_CHECK_CLOSE(r.surfacePressure, thp, 1e-6);
    // liquid-like CO2 below, two-phase CO2 (the temperature is the saturation temperature) at the surface
    BOOST_CHECK_EQUAL(r.profile.back().region, 6);
    BOOST_CHECK_EQUAL(r.profile.front().region, 7);
    const auto& tables = Opm::CO2PHTables::instance();
    BOOST_CHECK_SMALL(r.surfaceTemperature - tables.saturationTemperature(Opm::DenseAd::Evaluation<Scalar, 1>(thp)).value(), 0.05);
    int changes = 0;
    for (std::size_t i = r.profile.size() - 1; i > 0; --i) {
        changes += r.profile[i].region != r.profile[i - 1].region;
        BOOST_CHECK_LT(r.profile[i - 1].pressure, r.profile[i].pressure);
    }
    BOOST_CHECK_EQUAL(changes, 1);
    // the density falls as the CO2 boils
    BOOST_CHECK_LT(r.profile.front().density, 0.8 * r.profile.back().density);
    const Scalar E0 = totalEnergy(r.profile.back());
    for (const auto& pt : r.profile) {
        BOOST_CHECK_SMALL(totalEnergy(pt) - E0, 5.0);
    }
    // a higher wellhead pressure keeps the CO2 dense, and at 65 bar the reservoir cannot lift the column
    const auto r2 = well.solveProducer(55 * bar);
    BOOST_REQUIRE(r2.flowing);
    BOOST_CHECK(!r2.phaseChangeMd);
    BOOST_CHECK(!well.solveProducer(65 * bar).flowing);
}

BOOST_AUTO_TEST_CASE(Co2InjectorWithHeatLoss)
{
    // the formation (15 C + 0.03 K/m) is warmer than the injected CO2 at depth and heats it
    const Scalar H = 1000.0;
    Opm::WellBoreOptions<Scalar> opt;
    opt.heatLoss.coefficient = 2.0;
    opt.heatLoss.surfaceTemperature = 288.15;
    opt.heatLoss.gradient = 0.03;
    Perf perf;
    perf.md = H;
    perf.pressure = 125 * bar;
    perf.injectivity = 2.0e-8;
    Geometry geo = Geometry::vertical(H, 0.15);
    geo.maxStepLength = 25.0;
    const Well adiabatic(geo, {perf});
    const Well heated(geo, {perf}, opt);
    const auto ra = adiabatic.solveInjector(0.02, 283.15, Fluid::Stream{1.0, 0.0});
    const auto rh = heated.solveInjector(0.02, 283.15, Fluid::Stream{1.0, 0.0});
    std::cout << "CO2 injector bottom T: adiabatic " << ra.profile.back().temperature - 273.15 << " C, with heat from the formation "
              << rh.profile.back().temperature - 273.15 << " C\n";
    BOOST_CHECK_GT(rh.profile.back().temperature, ra.profile.back().temperature + 0.3);
    for (const auto& pt : rh.profile) {
        BOOST_CHECK_LT(pt.heatLossRate, 0.0);   // heat flows into the well
        BOOST_CHECK_LT(pt.temperature, 288.15 + 0.03 * pt.tvd);
    }
    BOOST_CHECK_CLOSE(rh.perforationRate[0], rh.massFlow, 1e-6);
}

BOOST_AUTO_TEST_CASE(WarmStartedSolvesOfCo2Wells)
{
    using Clock = std::chrono::steady_clock;
    const Scalar H = 1500.0, pRes = 200 * bar, TRes = 340.0;
    Scalar Xeq, Yeq;
    equilibrium(TRes, pRes, Xeq, Yeq);
    Geometry geo = Geometry::vertical(H, 0.2);
    geo.maxStepLength = 25.0;
    auto perf = Fluid::perforation(H, pRes, TRes, 1.0, 0.9 * Xeq, 0.0, 1.0e-12);
    perf.productivity = 1e-5;
    const Well producer(geo, {perf});
    const auto seconds = [](auto&& f) {
        const auto t0 = Clock::now();
        f();
        return std::chrono::duration<double>(Clock::now() - t0).count();
    };
    Opm::WellBoreResult<Scalar> cold, warm;
    const auto n0 = producer.marchCount();
    const double tCold = seconds([&] { cold = producer.solveProducer(40 * bar); });
    const auto nCold = producer.marchCount() - n0;
    const auto n1 = producer.marchCount();
    const double tWarm = seconds([&] { warm = producer.solveProducer(40 * bar, cold.bottomHolePressure * (1.0 + 1e-3)); });
    const auto nWarm = producer.marchCount() - n1;
    std::cout << "exsolving brine producer: cold " << tCold << " s (" << nCold << " marches), hint 0.1% off " << tWarm << " s ("
              << nWarm << " marches)\n";
    BOOST_CHECK_CLOSE(warm.bottomHolePressure, cold.bottomHolePressure, 1e-7);
    BOOST_CHECK_CLOSE(warm.massFlow, cold.massFlow, 1e-6);
    BOOST_CHECK_LE(nWarm, nCold / 2);
    BOOST_CHECK_LT(tWarm, 0.6 * tCold);

    // a CO2 injector
    Perf iperf;
    iperf.md = 1000.0;
    iperf.pressure = 125 * bar;
    iperf.injectivity = 2.0e-8;
    Geometry igeo = Geometry::vertical(1000.0, 0.15);
    igeo.maxStepLength = 25.0;
    const Well injector(igeo, {iperf});
    const auto m0 = injector.marchCount();
    const auto coldI = injector.solveInjector(0.02, 283.15, Fluid::Stream{1.0, 0.0});
    const auto mCold = injector.marchCount() - m0;
    const auto m1 = injector.marchCount();
    const auto warmI = injector.solveInjector(0.02, 283.15, Fluid::Stream{1.0, 0.0}, coldI.surfacePressure * (1.0 - 1e-3));
    const auto mWarm = injector.marchCount() - m1;
    std::cout << "CO2 injector: cold " << mCold << " marches, hint 0.1% off " << mWarm << " marches\n";
    BOOST_CHECK_CLOSE(warmI.surfacePressure, coldI.surfacePressure, 1e-6);
    BOOST_CHECK_CLOSE(warmI.massFlow, coldI.massFlow, 1e-8);
    BOOST_CHECK_LE(mWarm, mCold / 10);
}

BOOST_AUTO_TEST_CASE(FlashFromTheOverallCompositionOverTheRange)
{
    // states with both phases over the pressures, temperatures and saturations of the model (and CO2 in its
    // two-phase region), turned into their overall composition and back
    const auto& tables = Opm::CO2PHTables::instance();
    int count = 0, dome = 0;
    Scalar maxS = 0.0, maxT = 0.0, maxRho = 0.0;
    const auto check = [&](Scalar p, Scalar h, Scalar sw, int meaning) {
        const auto st = Provider::solve<Scalar>(p, h, sw, meaning);
        Scalar co2 = 0.0, total = 0.0;
        for (unsigned a = 0; a < 3; ++a) {
            co2 += st.saturation(a) * st.density(a) * st.massFraction(a, FS::CO2Idx);
            total += st.saturation(a) * st.density(a);
        }
        const Scalar Z = co2 / total;
        const auto f = Fluid::flash(p, h, {Z, 1.0 - Z});
        BOOST_REQUIRE_EQUAL(f.meaning, meaning);
        maxS = std::max(maxS, std::abs(f.s - sw));
        maxT = std::max(maxT, std::abs(f.state.temperature() - st.temperature()));
        for (unsigned a = 0; a < 3; ++a) {
            maxRho = std::max(maxRho, std::abs(f.state.density(a) / st.density(a) - 1.0));
        }
        ++count;
    };
    for (const Scalar p : {20.0, 45.0, 60.0, 80.0, 120.0, 200.0, 500.0}) {
        for (const Scalar T : {285.0, 300.0, 320.0, 350.0, 390.0}) {
            for (const Scalar Sw : {0.02, 0.2, 0.5, 0.8, 0.98}) {
                Scalar h, s;
                int meaning;
                try {
                    Provider::initialVariables(p * bar, T, Sw, 0.0, 0.0, h, s, meaning);
                }
                catch (const std::exception&) {
                    continue;   // CO2 at the saturation temperature or outside the tables
                }
                check(p * bar, h, s, meaning);
            }
        }
    }
    // CO2 in its two-phase region: the enthalpy of the CO2 component between those of the saturated phases
    for (const Scalar p : {42.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0}) {
        const Scalar pp = p * bar;
        const Scalar hl = tables.saturatedLiquidEnthalpy(Eval1(pp)).value();
        const Scalar hv = tables.saturatedVaporEnthalpy(Eval1(pp)).value();
        for (const Scalar q : {0.1, 0.5, 0.9}) {
            for (const Scalar Sw : {0.1, 0.5, 0.9}) {
                const Scalar u = hl + q * (hv - hl);
                const Scalar h = Provider::mixtureEnthalpy<Scalar>(pp, u, Sw, Provider::bothPhases);
                try {
                    check(pp, h, Sw, Provider::bothPhases);
                    ++dome;
                }
                catch (const std::exception&) {
                }
            }
        }
    }
    std::cout << count << " states (" << dome << " in the CO2 dome): max |dSw| " << maxS << ", max |dT| " << maxT
              << " K, max relative density difference " << maxRho << "\n";
    BOOST_CHECK_GT(count, 90);
    BOOST_CHECK_GT(dome, 20);
    BOOST_CHECK_SMALL(maxS, 1e-8);
    BOOST_CHECK_SMALL(maxT, 1e-6);
    BOOST_CHECK_SMALL(maxRho, 1e-7);
}

BOOST_AUTO_TEST_CASE(Co2SensitivitiesFromTheImplicitFunctionTheorem)
{
    // pressure, productivity, enthalpy and CO2 fraction of an exsolving brine producer; pressure and injectivity of a CO2 injector
    const Scalar H = 1500.0, pRes = 200 * bar, TRes = 340.0;
    Scalar Xeq, Yeq;
    equilibrium(TRes, pRes, Xeq, Yeq);
    Geometry geo = Geometry::vertical(H, 0.2);
    geo.maxStepLength = 25.0;
    auto perf = Fluid::perforation(H, pRes, TRes, 1.0, 0.9 * Xeq, 0.0, 1.0e-12);
    perf.productivity = 1e-5;
    const Well producer(geo, {perf});
    const Scalar thp = 40 * bar;
    const auto base = producer.solveProducer(thp);
    BOOST_REQUIRE(base.flowing);

    const Scalar dpr = 1.0e-5 * pRes;
    struct Descriptor { const char* name; Scalar size; std::function<void(Perf&)> change; };
    const std::vector<Descriptor> descriptors = {
        {"pressure", dpr, [&](Perf& p) { p.pressure += dpr; }},
        {"productivity", 1e-8, [](Perf& p) { p.productivity += 1e-8; }},
        {"enthalpy", 1.0e3, [](Perf& p) { p.inflowEnthalpy += 1.0e3; }},
        {"CO2 fraction", 0.001, [](Perf& p) { p.inflowComposition[0] += 0.001; p.inflowComposition[1] -= 0.001; }},
    };
    std::vector<Well::Perturbation> perturbations;
    for (const auto& d : descriptors) {
        perturbations.push_back({0, [&d](Perf& p) { d.change(p); }, d.size});
    }
    const auto before = producer.marchCount();
    const auto ift = producer.producerSensitivities(base, perturbations);
    const auto marches = producer.marchCount() - before;
    BOOST_REQUIRE(ift);
    std::cout << "CO2 producer: " << descriptors.size() << " sensitivities in " << marches << " marches\n";
    BOOST_CHECK_EQUAL(marches, descriptors.size() + 1);
    for (std::size_t k = 0; k < descriptors.size(); ++k) {
        auto changed = std::vector<Perf>{perf};
        descriptors[k].change(changed[0]);
        const auto r = Well(geo, changed).solveProducer(thp, base.bottomHolePressure);
        const Scalar fd = (r.perforationPressure[0] - base.perforationPressure[0]) / descriptors[k].size;
        std::cout << "  " << descriptors[k].name << ": IFT " << (*ift)[k][0] << ", finite difference " << fd << "\n";
        BOOST_CHECK_SMALL((*ift)[k][0] - fd, 0.03 * std::abs(fd) + 1e-6);
    }
    // the composition changes the weight of the column: more CO2, lighter column, lower pressure at the bottom
    BOOST_CHECK_LT((*ift)[3][0], 0.0);

    Perf iperf;
    iperf.md = 1000.0;
    iperf.pressure = 125 * bar;
    iperf.injectivity = 2.0e-8;
    Geometry igeo = Geometry::vertical(1000.0, 0.15);
    igeo.maxStepLength = 25.0;
    const Well injector(igeo, {iperf});
    const auto ibase = injector.solveInjector(0.02, 283.15, Fluid::Stream{1.0, 0.0});
    BOOST_REQUIRE(ibase.rateReached);
    const Scalar dpi = 1.0e-5 * 125 * bar;
    const std::vector<Descriptor> idescriptors = {
        {"pressure", dpi, [&](Perf& p) { p.pressure += dpi; }},
        {"injectivity", 2e-11, [](Perf& p) { p.injectivity += 2e-11; }},
    };
    std::vector<Well::Perturbation> iperturbations;
    for (const auto& d : idescriptors) {
        iperturbations.push_back({0, [&d](Perf& p) { d.change(p); }, d.size});
    }
    const auto iift = injector.injectorSensitivities(ibase, 0.02, 283.15, Fluid::Stream{1.0, 0.0}, iperturbations);
    BOOST_REQUIRE(iift);
    for (std::size_t k = 0; k < idescriptors.size(); ++k) {
        auto changed = std::vector<Perf>{iperf};
        idescriptors[k].change(changed[0]);
        const auto r = Well(igeo, changed).solveInjector(0.02, 283.15, Fluid::Stream{1.0, 0.0}, ibase.surfacePressure);
        const Scalar fd = (r.perforationPressure[0] - ibase.perforationPressure[0]) / idescriptors[k].size;
        std::cout << "  injector " << idescriptors[k].name << ": IFT " << (*iift)[k][0] << ", finite difference " << fd << "\n";
        BOOST_CHECK_SMALL((*iift)[k][0] - fd, 0.03 * std::abs(fd) + 1e-6);
    }
}
