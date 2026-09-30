# Manual scientific regression, excluded from automatic test discovery.
# Thirty-six frequency solves. No mesh selection, retries, clipping or correction.
using LineCableModels, Gmsh, Test, Printf
include(joinpath(@__DIR__, "../../support/scenarios.jl"))

function run_fem_pml_signs(output; pml_layers=(144,144,96))
    mkpath(output)
    reductions=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    fem=Formulation(:LineCableModelsFEM;options=(;reductions...,physics=:quasi_fw))
    analytical=Formulation(earth_impedance=formula(:unified;options=(Γ=0.,)),
        earth_admittance=formula(:unified;options=(Γ=0.,));options=reductions)
    options=(;pml_layers,pml_grading=(192/191)*log(1536),domain_skin_depths=24.,
        mesh_size_factor=3.,exterior_mesh_size_factor=8.,volume_quadrature=12,
        conductor_geometry_tolerance=1e-3,conductor_skin_depth_elements=6.,
        conductor_mesh_growth=sqrt(1.25),conductor_skin_depths=5.,conductor_thickness_elements=4,
        frequency_workers=4,solver_threads=1,mesh_policy=:reuse,resume_run_directory=:latest,
        keep_run_directory=true,trace=true,timing=true,gmsh_verbosity=2,getdp_verbosity=4)
    problems=Pair{String,LineParametersProblem}[]
    copper=Material(MaterialsLibrary(add_defaults=true),:copper)
    for radius in (.001,.01,.085)
        design=build(CableDesign,"two_bare_wires",Stack(Group(:core,
            Region(:core_metal,Disk(radius),copper))))
        system=build(LineCableSystem,[design,design],[Pose2(0.,1.),Pose2(1.,1.)];
            connections=[Dict(:core=>1),Dict(:core=>2)],system_id="two_bare_wires",line_length=1.)
        # Preserve the ordinary study's entire frequency definition and mesh plan.
        problem=LineParametersProblem(system;temperature=20.,frequencies=collect(10.0.^range(-1,6;length=10)),
            earth_props=homogeneous(rho=.1,eps_r=1.,mu_r=1.))
        push!(problems,"radius-$radius"=>problem)
    end
    for layout in (:all_earth,:air_1)
        problem=CurrentScenarios.three_bare_wires_problem(;
            heights=getproperty(CurrentScenarios.three_bare_wires_layouts,layout),
            frequencies=[.1,21.544346900318832,1e6],name="three_bare_wires_$layout")
        push!(problems,"three-$layout"=>problem)
    end
    for (label,problem) in problems
        println("START ",label);flush(stdout)
        reference=compute(problem,analytical)
        result=compute(problem,fem;options)
        open(joinpath(output,label*".csv"),"w") do io
            println(io,"quantity,frequency_hz,receiver,source,real,imaginary,reference_real,reference_imaginary")
            for (q,a,b) in (("Z",Z(result),Z(reference)),("Y",Y(result),Y(reference))),
                    k in eachindex(problem.frequencies),j in axes(a,2),i in axes(a,1)
                @printf(io,"%s,%.17g,%d,%d,%.17g,%.17g,%.17g,%.17g\n",q,problem.frequencies[k],i,j,
                    real(a[i,j,k]),imag(a[i,j,k]),real(b[i,j,k]),imag(b[i,j,k]))
            end
        end
        @testset "$label conductance signs" begin
            # This requirement belongs to these scientific fixtures, not compute.
            @test sign.(real.(Y(result))) == sign.(real.(Y(reference)))
        end
        println("DONE ",label," ",details(result).data.fem.run.run_directory);flush(stdout)
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_fem_pml_signs(isempty(ARGS) ? joinpath(pwd(),".linecablemodels/fem/pml-sign-regression") : abspath(only(ARGS)))
end
