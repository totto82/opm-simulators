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
#ifndef FLOW_PH_HPP
#define FLOW_PH_HPP


#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/nonlinear/newtonmethod.hh>
#include <opm/models/ph/phmodel.hh>
#include <opm/models/ph/phflash.hh>

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowProblemPh.hpp>
#include <opm/simulators/flow/NewTranFluxModule.hpp>

#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <flowexperimental/comp/CompAquiferModel.hpp>
#include <flowexperimental/comp/EmptyModel.hpp>
#include <flowexperimental/ph/PhWellModel.hpp>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

template <class TypeTag, class MyTypeTag>
struct WellModel;

template <class TypeTag, class MyTypeTag>
struct TracerModel;

namespace TTag {

struct FlowPhProblem {
    using InheritsFrom = std::tuple<FlowBaseProblemPh, PhModel>;
};

}

template<class TypeTag>
struct SparseMatrixAdapter<TypeTag, TTag::FlowPhProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using Block = MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Linear::IstlSparseMatrixAdapter<Block>;
};

template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::FlowPhProblem>
{ using type = TTag::FlowIstlSolver; };

template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::FlowPhProblem>
{ using type = ISTLSolver<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FlowPhProblem>
{ using type = FlowProblemPh<TypeTag>; };

template<class TypeTag>
struct AquiferModel<TypeTag, TTag::FlowPhProblem> {
    using type = CompAquiferModel<TypeTag>;
};

template<class TypeTag>
struct WellModel<TypeTag, TTag::FlowPhProblem> {
    using type = PhWellModel<TypeTag>;
};

template<class TypeTag>
struct TracerModel<TypeTag, TTag::FlowPhProblem> {
    using type = EmptyModel<TypeTag>;
};

//! Fluxes from the transmissibilities, so multipliers and edits apply
template<class TypeTag>
struct FluxModule<TypeTag, TTag::FlowPhProblem>
{ using type = NewTranFluxModule<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FlowPhProblem>
{ using type = typename GetPropType<TypeTag, Properties::PhFlash>::FluidSystem; };

template<class TypeTag>
struct Stencil<TypeTag, TTag::FlowPhProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar, GridView>;
};

#define OPM_PH_DISABLE(PROP) \
    template<class TypeTag> \
    struct PROP<TypeTag, TTag::FlowPhProblem> { static constexpr bool value = false; };

OPM_PH_DISABLE(EnableMech)
OPM_PH_DISABLE(EnableDisgasInWater)
OPM_PH_DISABLE(EnableApiTracking)
OPM_PH_DISABLE(EnableSaltPrecipitation)
OPM_PH_DISABLE(EnablePolymerMW)
OPM_PH_DISABLE(EnablePolymer)
OPM_PH_DISABLE(EnableDispersion)
OPM_PH_DISABLE(EnableBrine)
OPM_PH_DISABLE(EnableVapwat)
OPM_PH_DISABLE(EnableSolvent)
OPM_PH_DISABLE(EnableFoam)
OPM_PH_DISABLE(EnableExtbo)
OPM_PH_DISABLE(EnableBioeffects)
OPM_PH_DISABLE(EnableDebuggingChecks)

#undef OPM_PH_DISABLE

} // namespace Opm::Properties

#endif
