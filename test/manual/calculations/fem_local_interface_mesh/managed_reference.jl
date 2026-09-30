# Isolate construction-route differences from localization: use the captured
# pre-edit managed field builder in a separate module, with current equations.
# This does not replace a production method or change any frozen reference.
isdefined(@__MODULE__,:ROOT) || include("qualify.jl")
module ManagedReferenceControl
using LineCableModels, Gmsh
const gmsh = Gmsh.gmsh
const FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
using .FEM: FEMResolvedModel, FEMGeometry, FEMMeshPlan,
    _configure_conductor_mesh!, _entity_boundary, _interface_footprint
include(joinpath(Main.ROOT,"managed-reference-configure.jl"))
end

let
    problem, form = fixture(:mixed,1e6,.99)
    target = joinpath(ROOT,"mixed-f1.0e6-gamma0.99/managed-reference-terminal-shift")
    prepare_bundle(target,problem,form)
    if !isfile(joinpath(target,"study.msh"))
        model = FEM._resolved_fem_model(problem,form,
            computation_options(LineCableModelsFEM,ComputationOptions(;MESH...)))
        session = FEM._start_gmsh(2)
        try
            geometry = FEM._build_geometry!(model,"managed-reference")
            ManagedReferenceControl._configure_mesh!(model,geometry)
            gmsh.model.mesh.generate(2)
            gmsh.write(joinpath(target,"study.msh"))
            say("MANAGED REFERENCE nodes=",length(first(gmsh.model.mesh.get_nodes())))
        finally
            FEM._finish_gmsh(session)
        end
    end
    solve!(target)
    candidate = joinpath(ROOT,"feature-both-media-execution/f1.0e6-gamma-0.99/detached")
    native_reference = joinpath(ROOT,"mixed-f1.0e6-gamma0.99/baseline-terminal-shift")
    worst = 0.; route_change = 0.; flips = 0
    for quantity in ("Z","Y"), part in (real,imag)
        a = part.(native_matrix(target,quantity))
        b = part.(native_matrix(candidate,quantity))
        c = part.(native_matrix(native_reference,quantity))
        worst = max(worst,maximum(abs.((a-b)./a)))
        route_change = max(route_change,maximum(abs.((a-c)./c)))
        flips += count(sign.(a) .!= sign.(b))
    end
    result = Dict("managed_localization_relative_change"=>worst,
        "original_managed_vs_original_native_relative_change"=>route_change,
        "new_sign_changes"=>flips,"relative_tolerance"=>.02,
        "passed"=>worst<=.02 && flips==0)
    record(joinpath(target,"comparison.toml"),result)
    say("MANAGED CONSTRUCTION CONTROL ",result)
end
