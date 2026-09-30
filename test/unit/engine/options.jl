@testitem "Engine / option grammar / owner dispatch contract" tags=[:unit] setup=[
    UseEngineSupport
] begin
    const Grammar=LineCableModels.Grammar
    const Engine=LineCableModels.Engine

    struct UnregisteredFormulation<:Grammar.AbstractFormulation end

    @test LineCableModels.FormulationOptions === Grammar.FormulationOptions !== NamedTuple
    @test LineCableModels.ComputationOptions === Grammar.ComputationOptions !== NamedTuple
    @test LineCableModels.formulation_options === Grammar.formulation_options
    @test LineCableModels.computation_options === Grammar.computation_options
    @test parentmodule(Grammar.formulation_options) === Grammar
    @test parentmodule(Grammar.computation_options) === Grammar

    formulation_owner=LineParametersFormulation
    computation_type=LineCableModelsCoaxial
    @test hasmethod(
        Grammar.formulation_options,
        Tuple{Type{LineParametersFormulation}, FormulationOptions}
    )
    @test hasmethod(
        Grammar.computation_options,
        Tuple{Type{LineCableModelsCoaxial}, ComputationOptions}
    )
    # Passive selections cannot invoke live normalization by accident.
    @test_throws MethodError Grammar.formulation_options((;))
    retained_options = (reduce_bundle=false, kron_reduction=true,
        ideal_transposition=false)
    @test_throws MethodError Grammar.formulation_options((leaf=formulation_owner =>
        (options=retained_options,),))
    @test_throws MethodError Grammar.computation_options((;))
    @test_throws MethodError Grammar.formulation_options(:analytical, (;))
    @test_throws MethodError Grammar.computation_options(:analytical, ComputationOptions((;)))
    @test_throws MethodError Grammar.formulation_options(
        UnregisteredFormulation,
        (;)
    )
    @test_throws MethodError Grammar.computation_options(
        UnregisteredFormulation, ComputationOptions((;)))
    @test_throws MethodError Grammar.formulation_options(
        formulation_owner, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(
        computation_type, Dict{Symbol, Any}())
    @test_throws MethodError Grammar.computation_options(computation_type, nothing)

    formulation=@inferred Grammar.formulation_options(formulation_owner, FormulationOptions())
    @test formulation.data == (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = true
    )
    @test_throws ArgumentError Grammar.formulation_options(
        formulation_owner, FormulationOptions(unknown = true))

    default_execution=@inferred Grammar.computation_options(computation_type, ComputationOptions((;)))
    @test default_execution.data.output_basis == Val(:pul)
    @test default_execution.data.trace == Val(false)
    @test default_execution.data.on_result === nothing
    execution=Grammar.computation_options(
        computation_type, ComputationOptions((
            verbosity = (default = 1, NLsolve = 0),
            output_basis = :total,
            trace = true
        )))
    @test execution.data == (
        verbosity = (default = 1, NLsolve = 0),
        output_basis = Val(:total),
        trace = Val(true),
        on_result = nothing,
        timing = false
    )
    @test Engine.verbosity(execution, :NLsolve) == 0
    @test Engine.verbosity(execution, :unlisted) == 1
    @test_throws ArgumentError Grammar.computation_options(
        computation_type, ComputationOptions((unknown = true,)))
    @test_throws ArgumentError Grammar.computation_options(
        computation_type, ComputationOptions((output_basis = :unknown,)))
    for retired_basis in (:per_length, :per_lenght, :per_unit_length)
        @test_throws ArgumentError Grammar.computation_options(
            computation_type, ComputationOptions((output_basis = retired_basis,)))
    end
end

@testitem "Engine / FEM computation option ownership" tags=[:unit] begin
    using LineCableModels

    formulation = LineCableModelsFEM(options=(physics=:quasi_fw,))
    @test formulation.options isa FormulationOptions
    @test formulation.options.data.physics === Symbol("quasi-fw")
    @test LineCableModelsFEM().options.data.physics === Symbol("quasi-fw")
    for value in (:quasi_fw, "quasi-fw", Symbol("quasi-fw"))
        @test LineCableModelsFEM(options=(physics=value,)).options.data.physics === Symbol("quasi-fw")
    end
    for value in (:quasi_tem, "quasi-tem", Symbol("quasi-tem"), :fullwave, 0)
        @test_throws ArgumentError LineCableModelsFEM(options=(physics=value,))
    end
    @test fieldnames(typeof(formulation)) == (:methods, :options, :definitions)
    @test !haskey(NamedTuple(formulation), :execution)
    # Execution settings cannot enter through the formulation, or vice versa.
    @test_throws ArgumentError LineCableModelsFEM(options=(frequency_workers=4,))
    @test_throws ArgumentError LineCableModelsFEM(options=(domain_skin_depths=1.5,))
    @test_throws MethodError LineCableModelsFEM(fem_options=(ui=true,))
    @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((physics=:quasi_fw,)))
    @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((unknown=true,)))
    defaults = @inferred computation_options(LineCableModelsFEM, ComputationOptions((;)))
    @test defaults isa ComputationOptions
    @test defaults.data.domain_skin_depths === 2.0
    @test defaults.data.pml_thickness === nothing
    @test defaults.data.pml_thickness_factor == 1.0
    @test defaults.data.mesh_size_factor == 1.0
    @test defaults.data.exterior_mesh_size_factor == 1.0
    @test defaults.data.interface_refinement_factor == 1.0
    @test computation_options(LineCableModelsFEM,
        ComputationOptions(interface_refinement_factor=2)).data.interface_refinement_factor === 2.0
    for value in (0, .5, Inf, NaN, true, "2")
        @test_throws ArgumentError computation_options(LineCableModelsFEM,
            ComputationOptions(interface_refinement_factor=value))
    end
    @test defaults.data.pml_layers == (128,128,128)
    @test defaults.data.pml_grading == ntuple(_ -> (192/191)*log(1536), 3)
    @test defaults.data.pml_resolution === nothing
    physical = computation_options(LineCableModelsFEM,
        ComputationOptions(pml_resolution=(interpolation_cells=72,coefficient_change=.12)))
    @test physical.data.pml_resolution == (interpolation_cells=72,coefficient_change=.12)
    @test physical.data.pml_layers === physical.data.pml_grading === nothing
    mesh_record = ComputationOptions(pml_resolution=physical.data.pml_resolution,
        pml_layers=physical.data.pml_layers,pml_grading=physical.data.pml_grading)
    @test computation_options(LineCableModelsFEM,mesh_record).data.pml_resolution == physical.data.pml_resolution
    @test_throws ArgumentError computation_options(LineCableModelsFEM,
        ComputationOptions(pml_resolution=(;),pml_layers=96))
    @test_throws ArgumentError computation_options(LineCableModelsFEM,
        ComputationOptions(pml_resolution=(;),pml_grading=2.))
    for invalid in (true, (), (unknown=1,), (interpolation_cells=0,),
            (interpolation_cells=true,), (interpolation_cells=typemax(Cint),),
            (coefficient_change=0.,), (coefficient_change=Inf,), (coefficient_change=true,))
        @test_throws ArgumentError computation_options(LineCableModelsFEM,
            ComputationOptions(pml_resolution=invalid))
    end
    for (layers, grading) in ((1,0), ((Int16(3),4,5),(0,1f0,2.)))
        request = ComputationOptions(; pml_layers=layers, pml_grading=grading)
        resolved = computation_options(LineCableModelsFEM, request)
        @test resolved.data.pml_layers isa NTuple{3,Int}
        @test resolved.data.pml_grading isa NTuple{3,Float64}
        @test resolved.data.pml_layers == (layers isa Tuple ? layers : (layers,layers,layers))
        @test resolved.data.pml_grading == (grading isa Tuple ? grading : (grading,grading,grading))
    end
    @test defaults.data.pml_reflection == 1e-10
    @test defaults.data.volume_quadrature == 12
    @test defaults.data.physical_volume_quadrature === nothing
    @test defaults.data.pml_element_family === :triangle
    @test defaults.data.pml_quadrature == 9
    discrete = computation_options(LineCableModelsFEM, ComputationOptions(
        physical_volume_quadrature=Int32(3), pml_element_family=:quadrangle,
        pml_quadrature=Int32(16)))
    @test discrete.data.physical_volume_quadrature === 3
    @test discrete.data.pml_element_family === :quadrangle
    @test discrete.data.pml_quadrature === 16
    for invalid in ((physical_volume_quadrature=2,), (physical_volume_quadrature=3.,),
            (pml_element_family=:quad,), (pml_element_family="quadrangle",),
            (pml_quadrature=7,), (pml_quadrature=9.,))
        @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions(invalid))
    end
    @test defaults.data.conductor_geometry_tolerance == 1e-3
    @test defaults.data.conductor_skin_depth_elements == 6.0
    @test defaults.data.conductor_mesh_growth == sqrt(1.25)
    @test defaults.data.conductor_skin_depths == 5.0
    @test defaults.data.conductor_thickness_elements == 4
    for name in (:conductor_geometry_tolerance, :conductor_skin_depth_elements,
            :conductor_mesh_growth, :conductor_skin_depths)
        for value in (0, -1, Inf, NaN, true, "6")
            @test_throws ArgumentError computation_options(LineCableModelsFEM,
                ComputationOptions(NamedTuple{(name,)}((value,))))
        end
    end
    for controls in ((conductor_mesh_growth=.9,), (conductor_thickness_elements=0,),
            (conductor_thickness_elements=true,), (conductor_thickness_elements=4.5,))
        @test_throws ArgumentError computation_options(LineCableModelsFEM,
            ComputationOptions(controls))
    end
    @test computation_options(LineCableModelsFEM,
        ComputationOptions(pml_thickness=3)).data.pml_thickness == (3.,3.,3.)
    @test computation_options(LineCableModelsFEM,
        ComputationOptions(pml_thickness=(2,3,4))).data.pml_thickness == (2.,3.,4.)
    @test_throws ArgumentError LineCableModelsFEM(options=(pml_layers=32,))
    @test computation_options(LineCableModelsFEM, ComputationOptions((;domain_skin_depths=1.5))).data.domain_skin_depths === 1.5
    callback = (problem, index, result) -> nothing
    raw = (frequency_workers=Int32(4), solver_threads=Int16(2),
        on_result=callback, trace=true, mesh_policy=:remesh,
        getdp_executable=SubString("/tmp/getdp", 1), verbosity=(default=1,))
    configured = computation_options(LineCableModelsFEM, ComputationOptions(raw))
    @test keys(configured.data) == keys(defaults.data)
    @test isconcretetype(typeof(configured))
    @test configured.data.frequency_workers isa Int
    @test configured.data.solver_threads isa Int
    @test defaults.data.mumps_ordering === nothing
    @test defaults.data.petsc_prealloc === nothing
    native = computation_options(LineCableModelsFEM,
        ComputationOptions(mumps_ordering=Int32(0),petsc_prealloc=Int32(256)))
    @test native.data.mumps_ordering === 0
    @test native.data.petsc_prealloc === 256
    for controls in ((mumps_ordering=1,), (mumps_ordering=true,),
            (mumps_ordering=0.,), (petsc_prealloc=0,), (petsc_prealloc=true,),
            (petsc_prealloc=big(typemax(Cint))+1,))
        @test_throws ArgumentError computation_options(LineCableModelsFEM,ComputationOptions(controls))
    end
    @test configured.data.getdp_executable isa String
    @test configured.data.on_result === callback
    @test configured.data.trace === Val(true)
    resources(options) = options.data.frequency_workers * options.data.solver_threads
    @test (@inferred resources(configured)) == 8
    @test formulation.options.data.physics === Symbol("quasi-fw")
    for invalid in ((ui=1,), (plot_field_maps=:yes,), (keep_run_directory=1,),
        (mesh_path="",), (getdp_executable=1,), (log_file="",),
        (frequency_workers=true,), (solver_threads=1.5,), (gmsh_verbosity=true,),
        (getdp_verbosity=6,), (resume_run_directory="",), (resume_run_directory=:bad,))
        @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions(invalid))
    end
    for invalid in ((pml_layers=0,), (pml_layers=true,), (pml_layers=1.5,),
        (pml_layers=(1,2),), (pml_layers=(1,true,3),), (pml_layers=(1,0,3),),
        (pml_layers=(1,2.,3),), (pml_layers=typemax(Cint),),
        (pml_layers=big(typemax(Int))+1,), (pml_layers=[1,2,3],),
        (pml_grading=-1,), (pml_grading=Inf,), (pml_grading=NaN,),
        (pml_grading=true,), (pml_grading=(0,1),), (pml_grading=(0,1,Inf),),
        (pml_grading=(0,true,1),), (pml_grading=[0,1,2],),
        (pml_grading=big"1e-400",), (pml_grading=big"1e400",),
        (pml_reflection=0.,), (pml_reflection=1.,), (pml_reflection=NaN,),
        (pml_thickness=0.,), (pml_thickness=(1,2),), (pml_thickness=(1,Inf,2),),
        (pml_thickness=true,), (pml_thickness_factor=0.,), (pml_thickness_factor=Inf,),
        (pml_thickness_factor=true,), (mesh_size_factor=0.,), (mesh_size_factor=Inf,),
        (volume_quadrature=5,), (volume_quadrature=12.,),
        (exterior_mesh_size_factor=0.5,), (exterior_mesh_size_factor=Inf,),
        (exterior_mesh_size_factor=true,))
        @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions(invalid))
    end
    for invalid in (0, -1, Inf, NaN, true, "2")
        @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((;domain_skin_depths=invalid)))
    end
end
