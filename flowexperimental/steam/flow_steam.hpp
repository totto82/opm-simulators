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
#ifndef FLOW_STEAM_HPP
#define FLOW_STEAM_HPP

#include <opm/material/fluidsystems/SteamFluidSystem.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/nonlinear/newtonmethod.hh>
#include <opm/models/steam/steammodel.hh>

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowProblemSteam.hpp>
#include <opm/simulators/flow/NewTranFluxModule.hpp>

#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <flowexperimental/comp/CompAquiferModel.hpp>
#include <flowexperimental/comp/EmptyModel.hpp>
#include <flowexperimental/steam/SteamWellModel.hpp>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

template <class TypeTag, class MyTypeTag>
struct WellModel;

template <class TypeTag, class MyTypeTag>
struct TracerModel;

namespace TTag {

struct FlowSteamProblem {
    using InheritsFrom = std::tuple<FlowBaseProblemSteam, SteamModel>;
};

}

template<class TypeTag>
struct SparseMatrixAdapter<TypeTag, TTag::FlowSteamProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using Block = MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Linear::IstlSparseMatrixAdapter<Block>;
};

template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::FlowSteamProblem>
{ using type = TTag::FlowIstlSolver; };

template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::FlowSteamProblem>
{ using type = ISTLSolver<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FlowSteamProblem>
{ using type = FlowProblemSteam<TypeTag>; };

template<class TypeTag>
struct AquiferModel<TypeTag, TTag::FlowSteamProblem> {
    using type = CompAquiferModel<TypeTag>;
};

template<class TypeTag>
struct WellModel<TypeTag, TTag::FlowSteamProblem> {
    using type = SteamWellModel<TypeTag>;
};

template<class TypeTag>
struct TracerModel<TypeTag, TTag::FlowSteamProblem> {
    using type = EmptyModel<TypeTag>;
};

//! Fluxes from the transmissibilities, so multipliers and edits apply
template<class TypeTag>
struct FluxModule<TypeTag, TTag::FlowSteamProblem>
{ using type = NewTranFluxModule<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FlowSteamProblem>
{ using type = Opm::SteamFluidSystem<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct Stencil<TypeTag, TTag::FlowSteamProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar, GridView>;
};

#define OPM_STEAM_DISABLE(PROP) \
    template<class TypeTag> \
    struct PROP<TypeTag, TTag::FlowSteamProblem> { static constexpr bool value = false; };

OPM_STEAM_DISABLE(EnableMech)
OPM_STEAM_DISABLE(EnableDisgasInWater)
OPM_STEAM_DISABLE(EnableApiTracking)
OPM_STEAM_DISABLE(EnableSaltPrecipitation)
OPM_STEAM_DISABLE(EnablePolymerMW)
OPM_STEAM_DISABLE(EnablePolymer)
OPM_STEAM_DISABLE(EnableDispersion)
OPM_STEAM_DISABLE(EnableBrine)
OPM_STEAM_DISABLE(EnableVapwat)
OPM_STEAM_DISABLE(EnableSolvent)
OPM_STEAM_DISABLE(EnableFoam)
OPM_STEAM_DISABLE(EnableExtbo)
OPM_STEAM_DISABLE(EnableBioeffects)
OPM_STEAM_DISABLE(EnableDebuggingChecks)

#undef OPM_STEAM_DISABLE

} // namespace Opm::Properties

#endif
