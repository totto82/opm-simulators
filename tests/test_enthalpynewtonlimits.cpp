/*
  Copyright 2026 NORCE Research AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#define BOOST_TEST_MODULE EnthalpyNewtonLimits
#include <boost/test/unit_test.hpp>

#include <opm/models/common/enthalpynewtonlimits.hh>

namespace {

// Toy fluid: validity range p in [1e5, 1e7], h in [0, 3e6]; phase boundaries at
// hl(p) = 1e6 + 0.02 p and hg(p) = 2e6 + 0.02 p, as for a liquid/two-phase/vapour system.
struct ToyLimits
{
    double minPressure() const { return 1.0e5; }
    double maxPressure() const { return 1.0e7; }
    double minEnthalpy(double) const { return 0.0; }
    double maxEnthalpy(double) const { return 3.0e6; }
    double maxEnthalpyStep(double) const { return 0.25 * (hg(0) - hl(0)); }
    int region(double p, double h) const
    {
        if (h <= hl(p)) { return 0; }
        return h >= hg(p) ? 2 : 1;
    }
    static double hl(double p) { return 1.0e6 + 0.02 * p; }
    static double hg(double p) { return 2.0e6 + 0.02 * p; }
};

static_assert(Opm::PressureEnthalpyLimitsPolicy<ToyLimits, double>);

} // namespace

BOOST_AUTO_TEST_CASE(StepInsideRegionIsUnchanged)
{
    double p = 2.0e6;
    double h = 0.5e6 + 0.1e6;
    Opm::limitPressureEnthalpyStep(2.0e6, 0.5e6, p, h, ToyLimits{});
    BOOST_CHECK_EQUAL(p, 2.0e6);
    BOOST_CHECK_EQUAL(h, 0.6e6);
}

BOOST_AUTO_TEST_CASE(PressureChangeIsLimitedAndInRange)
{
    double p = 1.0e7;
    double h = 0.5e6;
    Opm::limitPressureEnthalpyStep(2.0e6, 0.5e6, p, h, ToyLimits{});
    BOOST_CHECK_CLOSE(p, 2.4e6, 1e-10); // +20 %

    p = -5.0e6;
    h = 0.5e6;
    Opm::limitPressureEnthalpyStep(2.0e5, 0.5e6, p, h, ToyLimits{});
    BOOST_CHECK_GE(p, 1.0e5); // not below the validity range
}

BOOST_AUTO_TEST_CASE(EnthalpyStepIsLimitedAndInRange)
{
    ToyLimits lim;
    double p = 2.0e6;
    double h = 2.9e6;
    Opm::limitPressureEnthalpyStep(2.0e6, 2.5e6, p, h, lim);
    BOOST_CHECK_LE(h, 2.5e6 + lim.maxEnthalpyStep(p) + 1e-6);
    BOOST_CHECK_LE(h, 3.0e6);
}

BOOST_AUTO_TEST_CASE(StepAcrossBoundaryLandsJustBeyondIt)
{
    ToyLimits lim;
    const double p0 = 2.0e6;
    const double h0 = ToyLimits::hl(p0) - 0.2e6; // liquid
    double p = p0;
    double h = ToyLimits::hl(p0) + 0.2e6;        // two-phase after the step
    Opm::limitPressureEnthalpyStep(p0, h0, p, h, lim);

    BOOST_CHECK_EQUAL(lim.region(p, h), 1);          // inside the new region
    BOOST_CHECK_LT(h, ToyLimits::hl(p0) + 0.2e6);    // but the step was cut
    BOOST_CHECK_GT(h, ToyLimits::hl(p0));            // and it did cross the boundary
    // only a small distance beyond the boundary: the overshoot is 10 % of the remaining step
    BOOST_CHECK_LT(h - ToyLimits::hl(p), 0.1 * 0.2e6 + 1.0);
}

BOOST_AUTO_TEST_CASE(StepAcrossWholeTwoPhaseRegionStopsInside)
{
    ToyLimits lim;
    const double p0 = 2.0e6;
    double p = p0;
    double h = ToyLimits::hg(p0) + 0.2e6;           // from liquid to beyond the vapour boundary
    Opm::limitPressureEnthalpyStep(p0, ToyLimits::hl(p0) - 0.01e6, p, h, lim);
    // limited to 25 % of the width, so it cannot reach the vapour boundary in one step
    BOOST_CHECK_LE(lim.region(p, h), 1);
}
