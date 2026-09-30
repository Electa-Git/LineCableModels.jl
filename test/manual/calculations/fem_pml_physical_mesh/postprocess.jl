# Reuse the qualified public plotting path; only the evidence root and label
# differ. No new solver calls and no data clipping.
root=joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh")
source=read(joinpath(dirname(root),"pml-conductance-cost/postprocess.jl"),String)
source=replace(source,"root=abspath(@__DIR__)"=>"root="*repr(root),
    "FEM 144/144/96"=>"FEM physical PML mesh")
Base.include_string(Main,source,"retained_public_postprocessing.jl")
