using LineCableModels, Measurements, JLD2, Random, Test
mode, path = ARGS
gauntlet_file=joinpath(pkgdir(LineCableModels),"gauntlet","Gauntlet.jl")
if mode == "legacy_write"
    # Load the checkpoint-era callable instead of redefining it after module
    # loading: redefinition itself increments Julia's generated closure names.
    case_file=joinpath(dirname(gauntlet_file),"cases.jl")
    source=read(case_file,String)
    first_method=findfirst("function _materialize_case(",source)
    next_method=findfirst("function _materialize_joint_case(",source)
    legacy=read(joinpath(@__DIR__,"gauntlet_legacy_materializer.jl"),String)
    source=source[begin:(first(first_method)-1)]*legacy*"\n"*source[first(next_method):end]
    source=replace(source,"_materialize_case(definition, sources, nominal_problem, variation)"=>
        "_materialize_case(definition, sources, nominal_problem)")
    module_source=replace(read(gauntlet_file,String),"include(\""=>"include(\"$(dirname(gauntlet_file))/")
    module_source=replace(module_source,"include($(repr(case_file)))"=>
        "Base.include_string(@__MODULE__, $(repr(source)), $(repr(case_file)))")
    Base.include_string(Main,module_source,gauntlet_file)
else
    include(gauntlet_file)
end
using .Gauntlet

if mode == "write"
    definition = Gauntlet._catalogue_uq_benchmark(:two_insulated_wires,@__FILE__;
        frequencies=[0.1,50.,1e7])
    sampled = rand(Xoshiro(998),definition.model.problem;distribution=:uniform)
    JLD2.jldsave(path;definitions_bytes=Gauntlet._execution_bytes([definition]),
        expected=numerical_input_sha256(sampled),correlation=correlation_record(definition.model))
elseif mode == "read"
    Base.include(Gauntlet,case_index()[:two_insulated_wires])
    definition = only(Gauntlet._read_execution(path,"definitions"))
    @testset "fresh-process joint declaration replay" begin
        @test definition.reference.formulation.options.distribution === :uniform
        @test correlation_record(definition.model) == JLD2.load(path,"correlation")
        sampled = rand(Xoshiro(998),definition.model.problem;distribution=:uniform)
        @test numerical_input_sha256(sampled) == JLD2.load(path,"expected")
        joint = only(filter(v->v isa JointParameterGrids,
            Gauntlet._variation_leaves(definition.model.variation)))
        p = only(joint.source)
        @test uncertainty(p.core_radius/p.insulation_thickness) < 1e-12
        @test Measurements.cov(p.core_radius,p.insulation_thickness) ≈
            0.01nominal(p.core_radius)*nominal(p.insulation_thickness)
        @test joint.record.law === :bounded_correlated_geometry_v1
    end
elseif mode == "legacy_write"
    loaded=load_case(:two_insulated_wires;variation=compose_variations(
        ExactOverrides(frequencies=[50.]),
        RelativeStandardUncertainty(10.;tags=(:geometry,:cable_layer))))
    problem=Gauntlet._materialize_case(loaded.definition,loaded.sources,loaded.nominal_problem)
    names=fieldnames(typeof(loaded))
    model=LoadedCase((name === :problem ? problem : getfield(loaded,name) for name in names)...)
    sampled=rand(Xoshiro(998),model.problem;distribution=:normal)
    JLD2.jldsave(path;model_bytes=Gauntlet._execution_bytes(model),expected=numerical_input_sha256(sampled))
elseif mode == "legacy_read"
    Base.include(Gauntlet,case_index()[:two_insulated_wires])
    model=Gauntlet._read_execution(path,"model")
    @testset "checkpoint-era independent declaration replay" begin
        @test correlation_record(model).rule === :parameter_identity
        @test all(v->!(v isa JointParameterGrids),Gauntlet._variation_leaves(model.variation))
        @test numerical_input_sha256(rand(Xoshiro(998),model.problem;distribution=:normal)) ==
            JLD2.load(path,"expected")
    end
else
    error("unknown fixture mode")
end
