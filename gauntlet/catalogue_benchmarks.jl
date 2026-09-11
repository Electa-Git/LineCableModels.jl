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
        fem_options::NamedTuple = (;),
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
    ), fem_options)
    reference = BenchmarkCalculation(
        :fem,
        model.problem,
        Formulation(
            :LineCableModelsFEM;
            options = _CATALOGUE_PHYSICAL_OPTIONS,
            fem_options = execution,
        );
        options = merge((trace = true, verbosity = (default = 0,)),
            reference_options),
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

function _catalogue_uq_benchmark(
        case_id::Symbol,
        source_file::AbstractString;
        frequencies = nothing,
        reference_options::NamedTuple = (;),
        candidate_options::NamedTuple = (;),
        variation::AbstractCaseVariation = NoVariation(),
)
    selected_variation = compose_variations(
        _catalogue_variation(frequencies, variation),
        RelativeStandardUncertainty(10.0; tags = (:geometry, :cable_layer)),
    )
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
            distribution = :normal,
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
