# Serial qualification, using the already qualified native solve/recovery
# harness. No production source writes and no changed conductor targets.
using LineCableModels, Gmsh, JSON3, TOML, LinearAlgebra, Printf, Dates, SHA
include("prototype.jl")
install_prototype!()

function load_native_harness()
    old=joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-conductance-cost/diagnostics.jl")
    source=read(old,String)
    root=joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-physical-mesh")
    source=replace(source,"const ROOT=abspath(@__DIR__)"=>"const ROOT="*repr(root),
        "pml_layers=(96,96,96)"=>"pml_layers=(72,72,72)",
        "saved_mesh=false,quadrature=12"=>"saved_mesh=false,quadrature=12,mesh_file=nothing",
        "mesh=joinpath(SOURCE,\"mesh\",@sprintf(\"frequency_%04d.msh\",k))"=>
            "mesh=mesh_file===nothing ? joinpath(SOURCE,\"mesh\",@sprintf(\"frequency_%04d.msh\",k)) : mesh_file",
        "if !saved_mesh\n"=>"if !saved_mesh && mesh_file===nothing\n")
    # Imports/definitions only. Its old executable entry point is not run.
    source=first(split(source,"if abspath(PROGRAM_FILE)==@__FILE__"))
    Base.include_string(Main,source,"qualified_native_diagnostics.jl")
end
load_native_harness()

if isempty(ARGS) || only(ARGS)=="a"
    say("STEP A: six-strip interpolation distribution, (72,72,72); domain/thickness/stretch retained")
    say("Predicted corner triangles: 41472; reference 138240. No conductor or physical-domain size changes.")
    for (label,k) in (("step-a-low",1),("step-a-mid",4))
        execute(label,k)
    end
    say("STEP A COMPLETE; inspect mandatory G signs and full components before domain reduction")
elseif only(ARGS)=="diagnose"
    say("Step A rejected: all four G signs wrong at both probes. Domain reduction and full batch withheld.")
    say("DIAGNOSE 3/8: identical candidate mesh, 13-point native triangle quadrature")
    execute("diagnose-quadrature13",1;quadrature=13,
        mesh_file=joinpath(ROOT,"step-a-low/detached/study.msh"))
    say("DIAGNOSE 4/8: new side only; original top144 and bottom96")
    PML_DIAGNOSTIC_DESIGN[]=:side_only
    execute("diagnose-side-only",1;mesh_options=(pml_layers=(72,144,96),))
    say("DIAGNOSE 5/8: new top only; original side144 and bottom96")
    PML_DIAGNOSTIC_DESIGN[]=:top_only
    execute("diagnose-top-only",1;mesh_options=(pml_layers=(144,72,96),))
    say("DIAGNOSIS COMPLETE: five of eight exploratory coupled solves used")
elseif only(ARGS)=="guarded"
    PML_DIAGNOSTIC_DESIGN[]=:guarded
    say("DIAGNOSE 6/8: interpolation plus coefficient-variation density; native split at b*u^3=1")
    say("Side/top 109 cells, bottom76 at 0.1 Hz; maximum realized side log-stretch change 0.201 (reference 0.155)")
    result=execute("diagnose-coefficient-low",1)
    if result["all_G_signs_match"]
        say("DIAGNOSE 7/8: coefficient distribution at second probe")
        execute("diagnose-coefficient-mid",4)
    else
        say("COEFFICIENT CANDIDATE REJECTED: mandatory low-frequency G signs failed")
    end
elseif only(ARGS)=="strict"
    PML_DIAGNOSTIC_DESIGN[]=:strict_guard
    say("LAST TWO EXPLORATORY CONTROLS: coefficient-variation prescription 0.12 instead of 0.15")
    say("The previous guard reduced G12 absolute error by a factor 63 but retained its wrong sign.")
    say("This step tests the predicted second-order reduction in that error; physical domain stays fixed.")
    execute("diagnose-strict-low",1)
    execute("diagnose-strict-mid",4)
    say("EIGHT-SOLVE EXPLORATORY LIMIT REACHED; no further native solves authorized by this stage")
else
    error("unknown prescribed stage")
end
