# Explicit detached export example. Export does not run a solver or open a UI.
using LineCableModels, Gmsh

destination = isempty(ARGS) ? joinpath(@__DIR__, "onelab-toy") : abspath(only(ARGS))
copper = Material(kind=:conductor, rho=1.72e-8)
wire = build(CableDesign, "two-wire-export",
    terminal(:core, core(copper; r=0.005)))
system = build(LineCableSystem, [wire, wire], [(0.0,0.1),(0.2,-0.1)];
    connections=[Dict(:core=>1),Dict(:core=>2)])
problem = LineParametersProblem(system; frequencies=[50.0,10000.0],
    earth_props=homogeneous(rho=100.0,eps_r=10.0))
formulation = Formulation(:LineCableModelsFEM; options=(
    reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
# Eight PML layers make this a small execution/preservation toy, not a
# converged engineering reference. Export's production default is unchanged.
entry = export_data(:onelab,problem,formulation;
    file_name=joinpath(destination,"study.pro"),mesh_options=(pml_layers=8,))
println(entry)
