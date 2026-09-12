const CATALOGUE_CASE_IDS = (
    :cable_320kv_armoured_dc_bipole,
    :cable_320kv_no_armour_dc_bipole,
    :cable_380kv_armoured_ac_flat,
    :cable_380kv_no_armour_ac_flat,
    :cable_525kv_land_no_armour_ac_flat,
    :cable_525kv_land_no_armour_dc_bipole,
    :cable_525kv_subsea_armoured_ac_flat,
    :cable_525kv_subsea_armoured_dc_bipole,
    :cable_132kv_630mm2_flathor,
    :cable_18kv_1000mm2_trefoil,
    :cable_18kv_1000mm2_trefoil_homogenized,
    :cable_30kv_na2xs2y_630mm2_trefoil,
    :cable_220kv_eaxecew_1x2500_252_trefoil,
    :cable_132kv_cigre_tb880_case0_630cu_trefoil,
    :cable_380kv_2000mm2_flatver,
    :cable_525kv_1600mm2_bipole,
    :cable_640kv_2000mm2_bipole,
    :solid_1000mm2_single,
    :two_bare_wires,
    :two_insulated_wires,
)

const PSCAD_CATALOGUE_CASE_IDS = Tuple(
    case_id for case_id in CATALOGUE_CASE_IDS if case_id !== :two_bare_wires
)

const _CATALOGUE_PHYSICAL_OPTIONS = (
    reduce_bundle = false,
    kron_reduction = false,
    ideal_transposition = false,
)

const _CATALOGUE_LINE_REPORT = (
    quantities = (:Z, :Y, :R, :L, :G, :C),
    bands = (:all, :dc, :harmonic, :narrow, :wide),
)

const _CATALOGUE_UQ_REPORT = (
    quantities = (:R, :L, :C, :G),
    statistics = (:mean, :std),
    bands = (:all,),
)

const _CATALOGUE_UQ_SEEDS = Dict{Symbol, UInt64}(
    :cable_132kv_630mm2_flathor => 0x132630,
    :cable_18kv_1000mm2_trefoil => 0x181000,
    :cable_380kv_2000mm2_flatver => 0x380200,
    :cable_525kv_1600mm2_bipole => 0x525160,
    :cable_640kv_2000mm2_bipole => 0x640200,
    :solid_1000mm2_single => 0x501000,
    :two_insulated_wires => 0x200002,
)

const _CATALOGUE_HIGH_TRIAL_CASE_IDS = (
    :cable_525kv_1600mm2_bipole,
    :two_insulated_wires,
)

function _catalogue_benchmark_id(case_id::Symbol, family::Symbol)
    case_name = replace(string(case_id), r"^cable_" => "")
    suffix = family === :uq ? "lep_montecarlo" : string(family)
    return Symbol("benchmark_", case_name, "_", suffix)
end

function _catalogue_variation(frequencies, variation::AbstractCaseVariation)
    frequencies === nothing && return variation
    return compose_variations(variation, ExactOverrides(; frequencies))
end

function _catalogue_candidate_formulations(model::LoadedCase)
    positions = model.nominal_problem.system.positions
    saad_applicable = all(
        left == right || !iszero(positions[left].x - positions[right].x)
        for left in eachindex(positions), right in eachindex(positions)
    )
    earth_impedance = saad_applicable ? (
        :default,
        :Pollaczek1926,
        :Saad1996,
        :WedepohlWilcox1973,
        :Xue2018,
    ) : (
        :default,
        :Pollaczek1926,
        :WedepohlWilcox1973,
        :Xue2018,
    )
    earth_admittance = saad_applicable ? (
        :default,
        :Pollaczek1926,
        :default,
        :default,
        :Xue2018,
    ) : (
        :default,
        :Pollaczek1926,
        :default,
        :Xue2018,
    )
    return Formulation(
        earth_impedance = Grid(earth_impedance),
        earth_admittance = Grid(earth_admittance),
        insulation_admittance = formula(:default);
        combine = :zip,
        options = _CATALOGUE_PHYSICAL_OPTIONS,
    )
end

function _catalogue_fem_benchmark(
        case_id::Symbol,
        source_file::AbstractString;
        frequencies = nothing,
        reference_options::NamedTuple = (;),
        candidate_options::NamedTuple = (;),
        variation::AbstractCaseVariation = NoVariation(),
)
    model = load_case(case_id;
        variation = _catalogue_variation(frequencies, variation))
    execution = merge((
        mesh_policy = :remesh,
        keep_run_directory = true,
        gmsh_verbosity = 0,
        getdp_verbosity = 4,
        frequency_workers = 2,
        solver_threads = 1,
    ), reference_options)
    reference = BenchmarkCalculation(
        :fem,
        model.problem,
        Formulation(
            :LineCableModelsFEM;
            options = _CATALOGUE_PHYSICAL_OPTIONS,
        );
        options = merge((trace = true, verbosity = (default = 0,)),
            execution),
    )
    return benchmark_definition(
        _catalogue_benchmark_id(case_id, :fem),
        case_id,
        :fem,
        source_file,
        model,
        reference,
        BenchmarkCalculation(
            :lcm,
            model.problem,
            _catalogue_candidate_formulations(model);
            options = candidate_options,
        ),
        _CATALOGUE_LINE_REPORT,
        (;),
    )
end

function _catalogue_pscad_benchmark(
        case_id::Symbol,
        source_file::AbstractString;
        frequencies = nothing,
        reference_options::NamedTuple = (;),
        candidate_options::NamedTuple = (;),
        variation::AbstractCaseVariation = NoVariation(),
)
    case_id === :two_bare_wires && throw(ArgumentError(
        "the two-bare-wire case is intentionally excluded from PSCAD campaigns",
    ))
    model = load_case(case_id;
        variation = _catalogue_variation(frequencies, variation))
    reference = BenchmarkCalculation(
        :pscad,
        model.problem,
        Formulation(:pscad; earth_impedance = :WedepohlWilcox1973);
        options = merge((verbosity=(default=0,PSCAD=0),),reference_options),
    )
    return benchmark_definition(
        _catalogue_benchmark_id(case_id, :pscad),
        case_id,
        :pscad,
        source_file,
        model,
        reference,
        BenchmarkCalculation(
            :lcm,
            model.problem,
            _catalogue_candidate_formulations(model);
            options = candidate_options,
        ),
        _CATALOGUE_LINE_REPORT,
        (;),
    )
end

function _catalogue_uq_seed(case_id::Symbol)
    haskey(_CATALOGUE_UQ_SEEDS, case_id) && return _CATALOGUE_UQ_SEEDS[case_id]
    index = findfirst(==(case_id), CATALOGUE_CASE_IDS)
    index === nothing && throw(ArgumentError("unknown catalogue case :$case_id"))
    return UInt64(0x5eed0000) + UInt64(index)
end

const _CATALOGUE_GEOMETRY_LAW = :bounded_correlated_geometry_v1
const _CATALOGUE_SCALE_SUPPORT = (1 - sqrt(3)*0.1, 1 + sqrt(3)*0.1)

# This certificate concerns the current rounded 60-degree sector construction,
# not a general proof for arbitrary user builders or distributions.
function _catalogue_sector_support(p)
    p.core_sectors == 6 || throw(ArgumentError(
        "the catalogue joint law requires a support certificate for non-six-sector cores"))
    lower_f, upper_f = p.core_fillet_factor .* _CATALOGUE_SCALE_SUPPORT
    0 < lower_f <= upper_f < 1/3 || throw(ArgumentError(
        "Milliken fillet support crosses a sector contact transition"))
    k = cot(pi/6) + 2cot(5pi/24) + cot(pi/12) - pi
    wire_area = pi*p.core_wire_radius^2
    lower = (p.core_outer_radius^2*(1/2-k*upper_f^2)-wire_area)/wire_area
    upper = (pi*p.core_outer_radius^2/6-wire_area)/wire_area
    courses = floor(Int,(sqrt(1+4lower/3)-1)/2)
    courses >= 1 && 3courses*(courses+1) <= lower &&
        upper < 3(courses+1)*(courses+2) || throw(ArgumentError(
        "Milliken support does not certify a fixed strand inventory: capacity in [$lower,$upper]"))
    return (fillet_support=(lower_f,upper_f), capacity_bounds=(lower,upper),
        courses, strands_per_sector=1+3courses*(courses+1))
end

function _catalogue_geometry_variation(definition::CaseDefinition,
        variation::AbstractCaseVariation)
    sources = _case_sources(definition, variation)
    names = Tuple(p.id for p in values(definition.parameters)
        if _matches(p, (:geometry,:cable_layer)))
    isempty(names) && throw(ArgumentError("catalogue case has no geometric uncertainty inputs"))
    joints = filter(v -> v isa JointParameterGrids, _variation_leaves(variation))
    owned = Set(Iterators.flatten(_joint_names(v.source) for v in joints))
    if !isempty(intersect(owned,names))
        all(name -> name in owned, names) || throw(ArgumentError(
            "a custom catalogue joint study must supply the complete geometric law"))
        return variation
    end
    roles = map(names) do name
        parameter = getproperty(definition.parameters,name)
        role = filter(tag -> tag in (:length,:area,:dimensionless), parameter.tags)
        length(role) == 1 || throw(ArgumentError(
            "selected geometric parameter :$name must declare exactly one dimensional role"))
        source = getproperty(sources,name)
        (LineCableModels.has_uncertainty(source) ||
            (source isa Real && !iszero(LineCableModels.uncertainty(source)))) &&
            throw(ArgumentError("uncertain source :$name conflicts with the default joint law; " *
                "supply a complete JointParameterGrids study instead"))
        return only(role)
    end
    base_grids = map(names) do name
        source = getproperty(sources,name)
        source isa Union{AbstractGrid,Gridspace} ? source : Grid((source,))
    end
    base = Gridspace{NamedTuple{names}}((args...)->NamedTuple{names}(args),base_grids)
    # Eager finite-base checks are owned declaration work. Sampling remains in
    # the engine, and no checks or probes are added to its calculation loops.
    base_records = NamedTuple[]
    for p in base
        all(value -> value isa Real && isfinite(value) && value > 0 &&
            iszero(LineCableModels.uncertainty(value)), values(p)) || throw(ArgumentError(
            "catalogue joint base dimensions must be positive, finite and deterministic"))
        push!(base_records,p)
    end
    isempty(base_records) && throw(ArgumentError("catalogue geometric base space is empty"))
    # Validate all finite base models, including count/design overrides. This
    # does not purport to certify unrelated user-supplied uncertain coordinates.
    base_problem = _materialize_case(definition,sources,nothing,variation)
    base_problems = base_problem isa Gridspace ? base_problem : (base_problem,)
    for problem in base_problems
        _validate_loaded_problem(definition,problem)
    end
    certificates = NamedTuple[]
    if definition.id === :cable_220kv_eaxecew_1x2500_252_trefoil
        sectors = getproperty(sources,:core_sectors)
        sectors = sectors isa Union{AbstractGrid,Gridspace} ? sectors : (sectors,)
        for p in base_records, count in sectors
            push!(certificates,_catalogue_sector_support(merge(p,(core_sectors=count,))))
        end
    end
    independent_names = Tuple(name for (name,role) in zip(names,roles) if role !== :length)
    has_lengths = :length in roles
    primitive_names = ((has_lengths ? (:length_scale,) : ())..., independent_names...)
    inputs = NamedTuple{primitive_names}(ntuple(_ -> (
        nominal=1.0,std=0.1,support=_CATALOGUE_SCALE_SUPPORT,distribution=:uniform,
        interpretation=:multiplicative_factor), length(primitive_names)))
    dependencies = NamedTuple{names}(map(zip(names,roles)) do (name,role)
        (role === :length ? :length_scale : name,)
    end |> Tuple)
    builder = function (p,factors...)
        factors_by_name = NamedTuple{primitive_names}(factors)
        return NamedTuple{names}(map(names) do name
            getproperty(p,name)*getproperty(factors_by_name,only(getproperty(dependencies,name)))
        end)
    end
    source = Gridspace{NamedTuple{names}}(builder,
        (base,ntuple(_ -> Grid(1.0,10.0),length(primitive_names))...))
    record = (law=_CATALOGUE_GEOMETRY_LAW, independent_inputs=inputs,
        dependencies, roles=NamedTuple{names}(roles), base_points=base_records,
        support_certificates=certificates,
        interpretation=:synthetic_correlated_geometry,
        note="shared length factor; independent area/lay/fillet factors; not measured manufacturing correlations")
    return compose_variations(variation,JointParameterGrids(source;record))
end

function _catalogue_uq_benchmark(
        case_id::Symbol,
        source_file::AbstractString;
        frequencies = nothing,
        reference_options::NamedTuple = (;),
        candidate_options::NamedTuple = (;),
        variation::AbstractCaseVariation = NoVariation(),
)
    definition = Base.include(@__MODULE__,case_index()[case_id])
    selected_variation = _catalogue_geometry_variation(definition,
        _catalogue_variation(frequencies, variation))
    model = load_case(case_id; variation = selected_variation)
    problem = ParametricProblem(model.problem)
    inner = uq_inner_formulation()
    trials = case_id in _CATALOGUE_HIGH_TRIAL_CASE_IDS ?
             4 * UQ_MONTE_CARLO_TRIALS : UQ_MONTE_CARLO_TRIALS
    reference = BenchmarkCalculation(
        :monte_carlo,
        problem,
        MonteCarlo(
            inner;
            trials,
            seed = _catalogue_uq_seed(case_id),
            distribution = :uniform,
            return_samples = false,
            return_histograms = false,
        );
        options = reference_options,
    )
    candidate = BenchmarkCalculation(
        :linear_error,
        problem,
        LinearError(inner);
        options = candidate_options,
    )
    return benchmark_definition(
        _catalogue_benchmark_id(case_id, :uq),
        case_id,
        :uq,
        source_file,
        model,
        reference,
        candidate,
        _CATALOGUE_UQ_REPORT,
        uq_moment_tolerances(),
    )
end
