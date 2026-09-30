# Feature verification after scientific qualification. Both public APIs use
# the implemented mesh owner; no in-memory geometry override is installed.
using LineCableModels, Gmsh, TOML, JSON3, LinearAlgebra, Dates, Test
const ROOT=joinpath(pwd(),".linecablemodels/fem/pml-physical-mesh")
const OUT=joinpath(ROOT,"feature-execution")
mkpath(OUT)
original=LineCableModels.ImportExport.deserialize_value(JSON3.read(read(
    joinpath(dirname(ROOT),"runs/run-WFpqS5/input/problem.json"),String),Dict{String,Any}))
problem=LineParametersProblem(original.system;frequencies=[1e6],
    temperature=original.temperature,earth_props=original.earth_props)
form=Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,
    reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
mesh=(domain_skin_depths=24.,pml_resolution=(interpolation_cells=72,coefficient_change=.12),
    mesh_size_factor=3.,exterior_mesh_size_factor=8.,volume_quadrature=12,
    conductor_geometry_tolerance=1e-3,conductor_skin_depth_elements=6.,
    conductor_mesh_growth=sqrt(1.25),conductor_skin_depths=5.,conductor_thickness_elements=4)
println(Dates.now()," FEATURE EXECUTION: public managed compute, one 1 MHz case");flush(stdout)
measurement=@timed compute(problem,form;options=(;mesh...,mesh_policy=:reuse,
    keep_run_directory=true,trace=true,timing=true,output_basis=:pul,
    resume_run_directory=:latest,frequency_workers=1,solver_threads=1,
    getdp_executable=joinpath(ROOT,"getdp-live"),getdp_verbosity=4))
result=measurement.value
runpath=details(result).data.fem.run.run_directory
entry=joinpath(OUT,"detached",basename(runpath),"study.pro")
isfile(entry) || export_data(:onelab,problem,form;file_name=entry,mesh_options=mesh)
output=joinpath(dirname(entry),"results/f0001-quasi-fw-b0000")
if !isfile(joinpath(output,"completed.txt"))
    println(Dates.now()," FEATURE EXECUTION: detached pro, exact managed mesh");flush(stdout)
    executable=joinpath(ROOT,"getdp-live")
    native_mesh=joinpath(runpath,"mesh/model.msh")
    open(joinpath(dirname(entry),"solver.log"),"w") do io
        run(pipeline(`$executable $entry -msh $native_mesh -solve LineCableModelsFEM -setnumber FrequencyIndex 1 -setnumber Physics 1 -setnumber PlotFieldMaps 0 -nt 1 -v 4`;stdout=io,stderr=io))
    end
end
function native_matrix(path)
    matrix=zeros(ComplexF64,2,2)
    for line in readlines(path)[3:end]
        row=split(line,'\t')
        matrix[parse(Int,row[1]),parse(Int,row[2])]=complex(parse(Float64,row[5]),parse(Float64,row[6]))
    end
    matrix
end
@testset "Implemented managed and detached execution" begin
    @test all(real.(Y(result)).<0)
    for (name,managed) in (("Z",Z(result)[:,:,1]),("Y",Y(result)[:,:,1]))
        detached=native_matrix(joinpath(output,"matrices/$name.tsv"))
        for component in (real,imag)
            scale=maximum(abs,component.(managed))
            @test all(abs.(component.(detached-managed)).<=2e-9.*abs.(component.(managed)).+100eps(Float64)*scale)
        end
    end
    # Public resume must recognize the complete physical prescription.
    reused=compute(problem,form;options=(;mesh...,mesh_policy=:reuse,
        keep_run_directory=true,trace=true,timing=true,output_basis=:pul,
        resume_run_directory=runpath,frequency_workers=1,solver_threads=1,
        getdp_executable=joinpath(ROOT,"getdp-live"),getdp_verbosity=4))
    @test details(reused).data.fem.run.reused
    @test Z(reused)==Z(result)
    @test Y(reused)==Y(result)
end
open(joinpath(OUT,"complete.toml"),"w") do io
    TOML.print(io,Dict("run_directory"=>runpath,"detached_entry"=>entry,
        "wall_seconds"=>measurement.time,"compile_seconds"=>measurement.compile_time,
        "same_mesh_parity"=>true,"managed_resume"=>true))
end
println(Dates.now()," FEATURE EXECUTION COMPLETE: managed compute, detached pro and explicit resume passed");flush(stdout)
