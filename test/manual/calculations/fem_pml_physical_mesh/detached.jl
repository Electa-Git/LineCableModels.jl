# Independent native entry point, compared separately with remeshed managed
# results and the saved managed-worker execution on the exact same mesh.
root=joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh")
source=read(joinpath(dirname(root),"pml-conductance-cost/detached.jl"),String)
source=replace(source,"root=abspath(@__DIR__)"=>"root="*repr(root),
    "phase2-balanced144-mid"=>"diagnose-strict-mid")
Base.include_string(Main,source,"retained_detached_assessment.jl")
