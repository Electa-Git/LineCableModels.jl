using Dates
using JLD2
using LinearAlgebra
using LineCableModels
using SHA

module FEMCatalogueCaseLoader
using SHA
import TOML
import LineCableModels
using LineCableModels: AbstractGrid, Grid, Gridspace

const GAUNTLET_ROOT = @__DIR__
include(joinpath(GAUNTLET_ROOT, "case_loader.jl"))
include(joinpath(GAUNTLET_ROOT, "provenance.jl"))
include(joinpath(GAUNTLET_ROOT, "reference_grid.jl"))
end

using .FEMCatalogueCaseLoader: ExactOverrides, case_index, load_case,
                               reference_case, numerical_input_sha256

const EARTH_IMPEDANCE = LineCableModels.Engine.EarthImpedance
const EARTH_ADMITTANCE = LineCableModels.Engine.EarthAdmittance
const SEMICON_ADMITTANCE = LineCableModels.Engine.SemiconAdmittance

const INTERNAL_FORMULA = :Schelkunoff1934
const BASE_Z = :Papadopoulos2010
const BASE_Y = :Papadopoulos2010
const CORPUS_POLICY = :faithful_literature_observation
const FORMULA_TREATMENT = :as_registered_no_regularization_no_fallback
const DIAGNOSTIC_POLICY = :isolated_candidate_reprobe_not_measured_fallback
const ARTIFACT_SCHEMA_VERSION = 7
const HOMOGENEOUS_UNDERGROUND_Z = Set((
    :Ametani2009,
    :Bridges1995,
    :Lucca1994,
    :Magalhaes2018,
    :MartinsBritto2024,
    :Papadopoulos2010,
    :Petrache2005,
    :Pollaczek1926,
    :Saad1996,
    :Theethayi2007,
    :Vance1978,
    :WedepohlWilcox1973,
    :Xue2018
))
const STRATIFIED_Z = Set((
    :Ametani1974,
    :Nakagawa1973,
    :Papadopoulos2009,
    :Papadopoulos2011,
    :Sunde1968,
    :Tsiamitros2008
))
const OVERHEAD_Z = Set((
    :AlvaradoBetancourt1983,
    :Carson1926,
    :Gary1976,
    :Noda2006,
    :Pettersson1994,
    :Theodoulidis2015,
    :Wise1934
))
const HOMOGENEOUS_UNDERGROUND_Y = Set((
    :IdealGround,
    :Magalhaes2018,
    :MartinsBritto2024,
    :Papadopoulos2010,
    :Pollaczek1926,
    :Theethayi2007,
    :Xue2018,
    :Xue2021
))
const STRATIFIED_Y = Set((:Papadopoulos2009, :Papadopoulos2011))
const OVERHEAD_Y = Set((:Ametani2021, :Pettersson1994, :Wise1948))
const SELF_ONLY_Z = Set((:Bridges1995, :Vance1978))
const HORIZONTAL_ONLY_Z = Set((
    :Petrache2005,
    :Saad1996,
    :Theethayi2007,
    :WedepohlWilcox1973
))
const HORIZONTAL_ONLY_Y = Set((:Theethayi2007, :Xue2021))

include(joinpath(@__DIR__, "formulas", "fem_lossless_semicon.jl"))

const GAUNTLET_ROOT = joinpath(
    pkgdir(LineCableModels),
    ".linecablemodels",
    "fem",
    "gauntlet"
)
const REFERENCE_ROOT = joinpath(GAUNTLET_ROOT, "corrected_fullband")
const OUTPUT_ROOT = joinpath(GAUNTLET_ROOT, "analytical_fullband_corpus")

function relative_error(actual, reference)
    norm(actual - reference) / max(norm(reference), eps(Float64))
end
rms_error(actual, reference) = norm(actual - reference) / sqrt(length(reference))
symmetry_error(matrix) = norm(matrix - transpose(matrix)) / max(norm(matrix), eps(Float64))

function coverage_record(kind, identifier)
    applicable, category, reason = if kind === :earth_impedance
        if identifier in HOMOGENEOUS_UNDERGROUND_Z
            (true, :homogeneous_underground, "")
        elseif identifier in STRATIFIED_Z
            (
                false,
                :stratified,
                identifier === :Papadopoulos2011 ?
                "registered underground route requires a two-layer earth; " *
                "multilayer excluded by benchmark scope" :
                "registered implementation has no homogeneous-underground " *
                "route in benchmark scope"
            )
        elseif identifier in OVERHEAD_Z
            (
                false,
                :overhead,
                "registered implementation exposes an overhead-only route"
            )
        else
            error("unclassified registered earth-impedance formula $identifier")
        end
    else
        if identifier in HOMOGENEOUS_UNDERGROUND_Y
            (true, :homogeneous_underground, "")
        elseif identifier in STRATIFIED_Y
            (
                false,
                :stratified,
                identifier === :Papadopoulos2011 ?
                "registered underground route requires a two-layer earth; " *
                "multilayer excluded by benchmark scope" :
                "registered implementation has no homogeneous-underground " *
                "route in benchmark scope"
            )
        elseif identifier in OVERHEAD_Y
            (
                false,
                :overhead,
                "registered implementation exposes an overhead-only route"
            )
        else
            error("unclassified registered earth-admittance formula $identifier")
        end
    end
    catalogue_module = kind === :earth_impedance ?
                       EARTH_IMPEDANCE : EARTH_ADMITTANCE
    formula_object = catalogue_module.Formula(identifier)
    family = kind === :earth_impedance ? "earthimpedance" : "earthadmittance"
    source_path = joinpath(
        pkgdir(LineCableModels),
        "src",
        "engine",
        family,
        "formulas",
        lowercase(string(identifier)) * ".jl"
    )
    isfile(source_path) || error(
        "registered $kind formula $identifier has no source file at $source_path"
    )
    route_limitations = if kind === :earth_impedance && identifier in SELF_ONLY_Z
        "mutual route intentionally rejects because the source supplies only self"
    elseif kind === :earth_impedance && identifier in HORIZONTAL_ONLY_Z ||
           kind === :earth_admittance && identifier in HORIZONTAL_ONLY_Y
        "mutual route requires nonzero horizontal cable separation"
    else
        ""
    end
    benchmark_route = if !applicable
        "not evaluated in the homogeneous-underground gauntlet"
    elseif kind === :earth_impedance && identifier in (:Ametani2009, :Lucca1994)
        "registered pair-complete recipe delegates all-underground pairs to " *
        "its Pollaczek1926 underground leaf; its distinctive mixed route is " *
        "not exercised by these cases"
    elseif kind === :earth_admittance && identifier === :Xue2018
        "registered default self/mutual route: infinite voltage reference; " *
        "surface and penetration routes remain separate unselected leaves"
    elseif kind === :earth_impedance && identifier in SELF_ONLY_Z
        "registered self route only; evaluated only for single-cable cases"
    else
        "registered default self/mutual routes"
    end
    return (;
        kind,
        identifier,
        description = LineCableModels.Engine.description(formula_object),
        propagation = string(catalogue_module.propagation(formula_object)),
        routes = join(string.(keys(catalogue_module.routes(formula_object))), ","),
        route_limitations,
        benchmark_route,
        source = relpath(source_path, pkgdir(LineCableModels)),
        source_sha256 = bytes2hex(sha256(read(source_path))),
        applicable,
        category,
        reason
    )
end

function catalogue()
    records = NamedTuple[]
    for identifier in EARTH_IMPEDANCE.formulas()
        push!(records, coverage_record(:earth_impedance, identifier))
    end
    for identifier in EARTH_ADMITTANCE.formulas()
        push!(records, coverage_record(:earth_admittance, identifier))
    end
    length(records) == 39 || error(
        "registered earth catalogue changed from 39 to $(length(records)) formulas"
    )
    return records
end

function variant(record)
    prefix = record.kind === :earth_impedance ? "earth_impedance" :
             "earth_admittance"
    return (
        id = Symbol(prefix, "_", lowercase(string(record.identifier))),
        kind = record.kind,
        identifier = record.identifier,
        earth_impedance = record.kind === :earth_impedance ?
                          record.identifier : BASE_Z,
        earth_admittance = record.kind === :earth_admittance ?
                           record.identifier : BASE_Y,
        coverage = record
    )
end

function formulation(selected)
    return Formulation(
        internal_impedance = INTERNAL_FORMULA,
        insulation_admittance = :Gustavsen2013,
        semicon_admittance = FEM_LOSSLESS_SEMICON,
        earth_impedance = selected.earth_impedance,
        earth_admittance = selected.earth_admittance,
        earth_properties = nothing,
        equivalent_earth = :default,
        options = (
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = false,
            temperature_correction = true
        )
    )
end

function candidate_provenance(model, selected)
    selected_formulation = formulation(selected)
    return (
        input_sha256 = numerical_input_sha256(model.nominal_problem),
        implementation = FEMCatalogueCaseLoader.implementation_record(
            selected_formulation;
            external_sources = ("test/gauntlet/formulas/fem_lossless_semicon.jl",)
        ),
        repository = FEMCatalogueCaseLoader.repository_provenance()
    )
end

function load_reference(case_id, model)
    path = joinpath(REFERENCE_ROOT, string(case_id), "reference.jld2")
    isfile(path) || error("full-band FEM reference is missing for $case_id: $path")
    document = JLD2.load(path)
    document["schema_version"] in (3, 4) || error(
        "$case_id FEM reference has unsupported schema"
    )
    input_matches = haskey(document, "input_sha256") ?
                    document["input_sha256"] == numerical_input_sha256(model.nominal_problem) :
                    document["case_source_sha256"] == model.source_sha256
    input_matches || error(
        "$case_id FEM reference case digest does not match the current case"
    )
    document["frequencies"] == model.problem.frequencies || error(
        "$case_id FEM reference does not use the selected reference frequency vector"
    )
    document["port_order"] == model.port_order || error(
        "$case_id FEM reference terminal order does not match the current case"
    )
    return path, document, bytes2hex(sha256(read(path)))
end

function component_matrices(z, y, frequencies)
    inductance = similar(real(z))
    capacitance = similar(real(y))
    for index in eachindex(frequencies)
        omega = 2pi * frequencies[index]
        inductance[:, :, index] = imag(@view z[:, :, index]) / omega
        capacitance[:, :, index] = imag(@view y[:, :, index]) / omega
    end
    return (
        R = real(z),
        L = inductance,
        G = real(y),
        C = capacitance
    )
end

function frequency_errors(actual, reference)
    count = size(actual, 3)
    relative = Vector{Float64}(undef, count)
    rms_absolute = Vector{Float64}(undef, count)
    max_absolute = Vector{Float64}(undef, count)
    for index in 1:count
        candidate = @view actual[:, :, index]
        truth = @view reference[:, :, index]
        relative[index] = relative_error(candidate, truth)
        rms_absolute[index] = rms_error(candidate, truth)
        max_absolute[index] = maximum(abs, candidate - truth)
    end
    return (;
        aggregate_relative = relative_error(actual, reference),
        aggregate_rms_absolute = rms_error(actual, reference),
        aggregate_max_absolute = maximum(abs, actual - reference),
        relative,
        rms_absolute,
        max_absolute
    )
end

function result_metrics(result, reference)
    frequencies = result.f
    actual_z = result.Z.values
    actual_y = result.Y.values
    reference_z = reference["Z"]
    reference_y = reference["Y"]
    all(isfinite, actual_z) || error("analytical Z contains non-finite values")
    all(isfinite, actual_y) || error("analytical Y contains non-finite values")

    z_symmetry = [symmetry_error(@view actual_z[:, :, index])
                  for index in eachindex(frequencies)]
    y_symmetry = [symmetry_error(@view actual_y[:, :, index])
                  for index in eachindex(frequencies)]
    actual_components = component_matrices(actual_z, actual_y, frequencies)
    reference_components = component_matrices(
        reference_z,
        reference_y,
        frequencies
    )
    components = NamedTuple{keys(actual_components)}(
        map(keys(actual_components)) do quantity
        frequency_errors(
            getproperty(actual_components, quantity),
            getproperty(reference_components, quantity)
        )
    end
    )
    spectra = NamedTuple{keys(actual_components)}(
        ntuple(_ -> Vector{Float64}(undef, length(frequencies)), 4)
    )
    tolerances = NamedTuple{keys(actual_components)}(
        ntuple(_ -> Vector{Float64}(undef, length(frequencies)), 4)
    )
    violations = NamedTuple[]
    for index in eachindex(frequencies)
        z_symmetry[index] <= 1.0e-12 || push!(violations,
            (
                quantity = :Z_reciprocity,
                frequency_index = index,
                frequency_hz = frequencies[index],
                minimum_eigenvalue = NaN,
                tolerance = 1.0e-12
            ))
        y_symmetry[index] <= 1.0e-12 || push!(violations,
            (
                quantity = :Y_reciprocity,
                frequency_index = index,
                frequency_hz = frequencies[index],
                minimum_eigenvalue = NaN,
                tolerance = 1.0e-12
            ))
        for quantity in keys(actual_components)
            matrix = @view getproperty(actual_components, quantity)[:, :, index]
            lower = eigmin(Symmetric(matrix))
            tolerance = 1.0e-8 * max(opnorm(matrix), eps(Float64))
            getproperty(spectra, quantity)[index] = lower
            getproperty(tolerances, quantity)[index] = tolerance
            lower >= -tolerance || push!(violations,
                (
                    quantity,
                    frequency_index = index,
                    frequency_hz = frequencies[index],
                    minimum_eigenvalue = lower,
                    tolerance
                ))
        end
    end
    return (
        Z = frequency_errors(actual_z, reference_z),
        Y = frequency_errors(actual_y, reference_y),
        components,
        z_symmetry,
        y_symmetry,
        spectra,
        tolerances,
        violations
    )
end

function artifact_path(case_id, selected)
    joinpath(
        OUTPUT_ROOT,
        "cases",
        string(case_id),
        string(selected.id) * ".jld2"
    )
end

function partition_violations(selected, metrics)
    candidate_quantities = selected.kind === :earth_impedance ?
                           Set((:R, :L, :Z_reciprocity)) :
                           Set((:G, :C, :Y_reciprocity))
    candidate = filter(
        violation -> violation.quantity in candidate_quantities,
        metrics.violations
    )
    context = filter(
        violation -> violation.quantity ∉ candidate_quantities,
        metrics.violations
    )
    return candidate, context
end

function summary_row(
        case_id,
        selected,
        metrics,
        frequencies,
        elapsed,
        path,
        reference
)
    candidate_violations, context_violations = partition_violations(
        selected,
        metrics
    )
    candidate_frequencies = getproperty.(candidate_violations, :frequency_hz)
    context_frequencies = getproperty.(context_violations, :frequency_hz)
    candidate_quantities = sort!(unique(string.(
        getproperty.(candidate_violations, :quantity)
    )))
    return (
        case = string(case_id),
        variant = string(selected.id),
        kind = string(selected.kind),
        formula = string(selected.identifier),
        earth_impedance = string(selected.earth_impedance),
        earth_admittance = string(selected.earth_admittance),
        fem_reference_status = string(reference["status"]),
        fem_reference_observations = length(get(
            reference, "observations", NamedTuple[]
        )),
        status = isempty(candidate_violations) ?
                 "observed_no_invariant_violation" :
                 "observed_invariant_violation",
        candidate_invariant_violations = length(candidate_violations),
        candidate_violation_quantities = join(candidate_quantities, ","),
        candidate_first_violation_hz = isempty(candidate_frequencies) ?
                                       NaN : minimum(candidate_frequencies),
        candidate_last_violation_hz = isempty(candidate_frequencies) ?
                                      NaN : maximum(candidate_frequencies),
        context_invariant_violations = length(context_violations),
        context_first_violation_hz = isempty(context_frequencies) ?
                                     NaN : minimum(context_frequencies),
        context_last_violation_hz = isempty(context_frequencies) ?
                                    NaN : maximum(context_frequencies),
        elapsed_seconds = elapsed,
        Z_relative_frobenius = metrics.Z.aggregate_relative,
        Y_relative_frobenius = metrics.Y.aggregate_relative,
        R_relative_frobenius = metrics.components.R.aggregate_relative,
        L_relative_frobenius = metrics.components.L.aggregate_relative,
        G_relative_frobenius = metrics.components.G.aggregate_relative,
        C_relative_frobenius = metrics.components.C.aggregate_relative,
        Z_worst_frequency_relative = maximum(metrics.Z.relative),
        Z_worst_frequency_hz = frequencies[argmax(metrics.Z.relative)],
        Y_worst_frequency_relative = maximum(metrics.Y.relative),
        Y_worst_frequency_hz = frequencies[argmax(metrics.Y.relative)],
        R_worst_frequency_relative = maximum(metrics.components.R.relative),
        R_worst_frequency_hz = frequencies[argmax(metrics.components.R.relative)],
        L_worst_frequency_relative = maximum(metrics.components.L.relative),
        L_worst_frequency_hz = frequencies[argmax(metrics.components.L.relative)],
        G_worst_frequency_relative = maximum(metrics.components.G.relative),
        G_worst_frequency_hz = frequencies[argmax(metrics.components.G.relative)],
        C_worst_frequency_relative = maximum(metrics.components.C.relative),
        C_worst_frequency_hz = frequencies[argmax(metrics.components.C.relative)],
        Z_symmetry = maximum(metrics.z_symmetry),
        Y_symmetry = maximum(metrics.y_symmetry),
        min_R_eigenvalue = minimum(metrics.spectra.R),
        min_L_eigenvalue = minimum(metrics.spectra.L),
        min_G_eigenvalue = minimum(metrics.spectra.G),
        min_C_eigenvalue = minimum(metrics.spectra.C),
        artifact = path
    )
end

function write_hash(path)
    open(path * ".sha256", "w") do io
        println(io, bytes2hex(sha256(read(path))), "  ", basename(path))
    end
    return path
end

function archive_stale_artifact(path)
    digest = bytes2hex(sha256(read(path)))
    stem, extension = splitext(path)
    archived = "$(stem)_stale_$(first(digest, 12))$(extension)"
    isfile(archived) || cp(path, archived)
    write_hash(archived)
    return archived
end

function record_result(
        case_id,
        model,
        selected,
        result,
        metrics,
        elapsed,
        reference_path,
        reference_sha256,
        reference
)
    path = artifact_path(case_id, selected)
    mkpath(dirname(path))
    temporary = path * ".new"
    candidate_violations, context_violations = partition_violations(
        selected,
        metrics
    )
    provenance = candidate_provenance(model, selected)
    JLD2.jldsave(
        temporary;
        schema_version = ARTIFACT_SCHEMA_VERSION,
        kind = :fem_analytical_fullband_comparison,
        corpus_policy = CORPUS_POLICY,
        formula_treatment = FORMULA_TREATMENT,
        diagnostic_policy = DIAGNOSTIC_POLICY,
        status = isempty(candidate_violations) ?
                 :observed_no_invariant_violation :
                 :observed_invariant_violation,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        input_sha256 = provenance.input_sha256,
        implementation = provenance.implementation,
        repository_commit = provenance.repository.commit,
        repository_dirty = provenance.repository.dirty,
        frequencies = copy(result.f),
        port_order = copy(model.port_order),
        internal_impedance = INTERNAL_FORMULA,
        insulation_admittance = :Gustavsen2013,
        semicon_admittance = :FEMLossless,
        earth_impedance = selected.earth_impedance,
        earth_admittance = selected.earth_admittance,
        earth_properties = nothing,
        equivalent_earth = :bottom_homogeneous_layer,
        variant_kind = selected.kind,
        variant_id = selected.id,
        formula = selected.identifier,
        Z = copy(result.Z.values),
        Y = copy(result.Y.values),
        metrics,
        candidate_violations,
        context_violations,
        elapsed_seconds = elapsed,
        fem_reference_path = reference_path,
        fem_reference_sha256 = reference_sha256,
        fem_reference_status = reference["status"],
        fem_reference_observations = get(
            reference, "observations", NamedTuple[]
        ),
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    return path
end

function record_skip(
        case_id,
        model,
        selected,
        reason,
        reference_path,
        reference_sha256,
        reference
)
    path = artifact_path(case_id, selected)
    provenance = candidate_provenance(model, selected)
    mkpath(dirname(path))
    temporary = path * ".new"
    JLD2.jldsave(
        temporary;
        schema_version = ARTIFACT_SCHEMA_VERSION,
        kind = :fem_analytical_fullband_comparison,
        corpus_policy = CORPUS_POLICY,
        formula_treatment = FORMULA_TREATMENT,
        diagnostic_policy = DIAGNOSTIC_POLICY,
        status = :not_applicable,
        reason,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        input_sha256 = provenance.input_sha256,
        implementation = provenance.implementation,
        repository_commit = provenance.repository.commit,
        repository_dirty = provenance.repository.dirty,
        frequencies = copy(model.problem.frequencies),
        port_order = copy(model.port_order),
        internal_impedance = INTERNAL_FORMULA,
        earth_impedance = selected.earth_impedance,
        earth_admittance = selected.earth_admittance,
        variant_kind = selected.kind,
        variant_id = selected.id,
        formula = selected.identifier,
        fem_reference_path = reference_path,
        fem_reference_sha256 = reference_sha256,
        fem_reference_status = reference["status"],
        fem_reference_observations = get(
            reference, "observations", NamedTuple[]
        ),
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    return path
end

function diagnose_execution_failure(case_id, selected, frequencies, reference)
    failures = NamedTuple[]
    successful_indices = Int[]
    successful_z = fill(
        ComplexF64(NaN, NaN),
        size(reference["Z"])
    )
    successful_y = fill(
        ComplexF64(NaN, NaN),
        size(reference["Y"])
    )
    for (index, frequency) in pairs(frequencies)
        probe = load_case(
            case_id;
            variation = ExactOverrides(frequencies = [frequency])
        )
        try
            result = compute(probe.problem, formulation(selected))
            reference_slice = Dict(
                "Z" => reference["Z"][:, :, index:index],
                "Y" => reference["Y"][:, :, index:index]
            )
            size(result.Z.values) == size(reference_slice["Z"]) || error(
                "analytical and FEM Z shapes differ"
            )
            size(result.Y.values) == size(reference_slice["Y"]) || error(
                "analytical and FEM Y shapes differ"
            )
            result_metrics(result, reference_slice)
            successful_z[:, :, index] .= result.Z.values[:, :, 1]
            successful_y[:, :, index] .= result.Y.values[:, :, 1]
            push!(successful_indices, index)
        catch exception
            isolated = selected.kind === :earth_impedance ?
                       merge(selected, (earth_admittance = :IdealGround,)) :
                       merge(selected, (earth_impedance = :Pollaczek1926,))
            attribution = try
                isolated_result = compute(
                    probe.problem,
                    formulation(isolated)
                )
                owned = selected.kind === :earth_impedance ?
                        isolated_result.Z.values : isolated_result.Y.values
                all(isfinite, owned) || error(
                    "isolated candidate returned non-finite values"
                )
                :context_baseline
            catch
                :candidate_or_shared
            end
            push!(failures,
                (
                    frequency_index = index,
                    frequency_hz = frequency,
                    attribution,
                    exception_type = string(typeof(exception)),
                    message = sprint(showerror, exception)
                ))
        end
    end
    partial_metrics = nothing
    candidate_violations = NamedTuple[]
    context_violations = NamedTuple[]
    if !isempty(successful_indices)
        partial_result = (
            f = frequencies[successful_indices],
            Z = (values = successful_z[:, :, successful_indices],),
            Y = (values = successful_y[:, :, successful_indices],)
        )
        partial_reference = Dict(
            "Z" => reference["Z"][:, :, successful_indices],
            "Y" => reference["Y"][:, :, successful_indices]
        )
        partial_metrics = result_metrics(partial_result, partial_reference)
        candidate_violations, context_violations = partition_violations(
            selected,
            partial_metrics
        )
    end
    return (;
        frequency_failures = failures,
        successful_frequency_indices = successful_indices,
        successful_frequencies = frequencies[successful_indices],
        Z = successful_z,
        Y = successful_y,
        partial_metrics,
        candidate_violations,
        context_violations
    )
end

function record_execution_failure(
        case_id,
        model,
        selected,
        exception,
        diagnostic,
        elapsed,
        reference_path,
        reference_sha256,
        reference
)
    path = artifact_path(case_id, selected)
    provenance = candidate_provenance(model, selected)
    mkpath(dirname(path))
    temporary = path * ".new"
    JLD2.jldsave(
        temporary;
        schema_version = ARTIFACT_SCHEMA_VERSION,
        kind = :fem_analytical_fullband_comparison,
        corpus_policy = CORPUS_POLICY,
        formula_treatment = FORMULA_TREATMENT,
        diagnostic_policy = DIAGNOSTIC_POLICY,
        status = :execution_failure,
        case_id = string(case_id),
        case_source_sha256 = model.source_sha256,
        input_sha256 = provenance.input_sha256,
        implementation = provenance.implementation,
        repository_commit = provenance.repository.commit,
        repository_dirty = provenance.repository.dirty,
        frequencies = copy(model.problem.frequencies),
        port_order = copy(model.port_order),
        internal_impedance = INTERNAL_FORMULA,
        earth_impedance = selected.earth_impedance,
        earth_admittance = selected.earth_admittance,
        variant_kind = selected.kind,
        variant_id = selected.id,
        formula = selected.identifier,
        exception_type = string(typeof(exception)),
        message = sprint(showerror, exception),
        frequency_failures = diagnostic.frequency_failures,
        successful_frequency_indices = diagnostic.successful_frequency_indices,
        successful_frequencies = diagnostic.successful_frequencies,
        Z = diagnostic.Z,
        Y = diagnostic.Y,
        partial_metrics = diagnostic.partial_metrics,
        candidate_violations = diagnostic.candidate_violations,
        context_violations = diagnostic.context_violations,
        elapsed_seconds = elapsed,
        fem_reference_path = reference_path,
        fem_reference_sha256 = reference_sha256,
        fem_reference_status = reference["status"],
        fem_reference_observations = get(
            reference, "observations", NamedTuple[]
        ),
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, path; force = true)
    write_hash(path)
    return path
end

function valid_existing(path, model, selected, reference_sha256)
    isfile(path) || return false
    try
        document = JLD2.load(path)
        provenance = candidate_provenance(model, selected)
        return document["schema_version"] == ARTIFACT_SCHEMA_VERSION &&
               document["input_sha256"] == provenance.input_sha256 &&
               document["implementation"] == provenance.implementation &&
               document["frequencies"] == model.problem.frequencies &&
               document["variant_id"] == selected.id &&
               document["formula"] == selected.identifier &&
               document["fem_reference_sha256"] == reference_sha256
    catch
        return false
    end
end

function execution_breakdown_summary(
        case_id,
        selected,
        frequency_failures,
        exception_type,
        message,
        path;
        successful_frequencies = Float64[],
        partial_metrics = nothing,
        candidate_violations = NamedTuple[],
        context_violations = NamedTuple[]
)
    candidate = filter(
        failure -> hasproperty(failure, :attribution) &&
                   failure.attribution === :candidate_or_shared,
        frequency_failures
    )
    context = filter(
        failure -> hasproperty(failure, :attribution) &&
                   failure.attribution === :context_baseline,
        frequency_failures
    )
    frequencies = getproperty.(frequency_failures, :frequency_hz)
    candidate_frequencies = getproperty.(candidate, :frequency_hz)
    context_frequencies = getproperty.(context, :frequency_hz)
    partial_candidate_frequencies = getproperty.(
        candidate_violations, :frequency_hz
    )
    partial_context_frequencies = getproperty.(
        context_violations, :frequency_hz
    )
    return (
        case = string(case_id),
        variant = string(selected.id),
        kind = string(selected.kind),
        formula = string(selected.identifier),
        failed_frequencies = length(frequency_failures),
        first_breakdown_hz = isempty(frequencies) ? NaN : minimum(frequencies),
        last_breakdown_hz = isempty(frequencies) ? NaN : maximum(frequencies),
        candidate_or_shared_failures = length(candidate),
        candidate_first_breakdown_hz = isempty(candidate_frequencies) ?
                                       NaN : minimum(candidate_frequencies),
        candidate_last_breakdown_hz = isempty(candidate_frequencies) ?
                                      NaN : maximum(candidate_frequencies),
        context_baseline_failures = length(context),
        context_first_breakdown_hz = isempty(context_frequencies) ?
                                     NaN : minimum(context_frequencies),
        context_last_breakdown_hz = isempty(context_frequencies) ?
                                    NaN : maximum(context_frequencies),
        successful_frequencies = length(successful_frequencies),
        first_success_hz = isempty(successful_frequencies) ?
                           NaN : minimum(successful_frequencies),
        last_success_hz = isempty(successful_frequencies) ?
                          NaN : maximum(successful_frequencies),
        partial_candidate_invariant_violations = length(candidate_violations),
        partial_candidate_first_violation_hz =
        isempty(partial_candidate_frequencies) ?
        NaN : minimum(partial_candidate_frequencies),
        partial_candidate_last_violation_hz =
        isempty(partial_candidate_frequencies) ?
        NaN : maximum(partial_candidate_frequencies),
        partial_context_invariant_violations = length(context_violations),
        partial_context_first_violation_hz =
        isempty(partial_context_frequencies) ?
        NaN : minimum(partial_context_frequencies),
        partial_context_last_violation_hz =
        isempty(partial_context_frequencies) ?
        NaN : maximum(partial_context_frequencies),
        partial_Z_relative_frobenius = partial_metrics === nothing ?
                                       NaN : partial_metrics.Z.aggregate_relative,
        partial_Y_relative_frobenius = partial_metrics === nothing ?
                                       NaN : partial_metrics.Y.aggregate_relative,
        partial_R_relative_frobenius = partial_metrics === nothing ?
                                       NaN :
                                       partial_metrics.components.R.aggregate_relative,
        partial_L_relative_frobenius = partial_metrics === nothing ?
                                       NaN :
                                       partial_metrics.components.L.aggregate_relative,
        partial_G_relative_frobenius = partial_metrics === nothing ?
                                       NaN :
                                       partial_metrics.components.G.aggregate_relative,
        partial_C_relative_frobenius = partial_metrics === nothing ?
                                       NaN :
                                       partial_metrics.components.C.aggregate_relative,
        exception_type = string(exception_type),
        message = string(message),
        artifact = path
    )
end

function recover_existing(path, selected, reference)
    document = JLD2.load(path)
    status = document["status"]
    if status === :not_applicable
        return (
            row = nothing,
            skipped = (
                case = document["case_id"],
                variant = string(selected.id),
                kind = string(selected.kind),
                formula = string(selected.identifier),
                reason = document["reason"],
                artifact = path
            ),
            failure = nothing
        )
    elseif status === :execution_failure
        return (
            row = nothing,
            skipped = nothing,
            failure = execution_breakdown_summary(
                document["case_id"],
                selected,
                document["frequency_failures"],
                document["exception_type"],
                document["message"],
                path;
                successful_frequencies = get(
                    document, "successful_frequencies", Float64[]
                ),
                partial_metrics = get(document, "partial_metrics", nothing),
                candidate_violations = get(
                    document, "candidate_violations", NamedTuple[]
                ),
                context_violations = get(
                    document, "context_violations", NamedTuple[]
                )
            )
        )
    end
    metrics = document["metrics"]
    return (
        row = summary_row(
            Symbol(document["case_id"]),
            selected,
            metrics,
            document["frequencies"],
            document["elapsed_seconds"],
            path,
            reference
        ),
        skipped = nothing,
        failure = nothing
    )
end

function case_skip_reason(model, selected)
    selected.coverage.applicable || return selected.coverage.reason
    if selected.kind === :earth_impedance &&
       selected.identifier in SELF_ONLY_Z &&
       length(model.problem.system.designs) > 1
        return "source formula supplies only a self term; a multi-cable FEM " *
               "matrix cannot be assembled without injecting a different " *
               "mutual formula"
    end
    positions = model.problem.system.positions
    has_vertical_pair = length(unique(getproperty.(positions, :x))) <
                        length(positions)
    if has_vertical_pair && (
        selected.kind === :earth_impedance &&
        selected.identifier in HORIZONTAL_ONLY_Z ||
        selected.kind === :earth_admittance &&
        selected.identifier in HORIZONTAL_ONLY_Y
    )
        return "source closed form is undefined at zero horizontal separation; " *
               "the case contains a vertical cable pair"
    end
    return nothing
end

function write_tsv(path, rows)
    isempty(rows) && error("cannot write an empty table to $path")
    temporary = path * ".new"
    names = collect(keys(first(rows)))
    open(temporary, "w") do io
        println(io, join(string.(names), '\t'))
        for record in rows
            println(io, join((getproperty(record, name) for name in names), '\t'))
        end
    end
    mv(temporary, path; force = true)
    return path
end

function write_failure_tsv(path, failures)
    temporary = path * ".new"
    names = (
        :case,
        :variant,
        :kind,
        :formula,
        :failed_frequencies,
        :first_breakdown_hz,
        :last_breakdown_hz,
        :candidate_or_shared_failures,
        :candidate_first_breakdown_hz,
        :candidate_last_breakdown_hz,
        :context_baseline_failures,
        :context_first_breakdown_hz,
        :context_last_breakdown_hz,
        :successful_frequencies,
        :first_success_hz,
        :last_success_hz,
        :partial_candidate_invariant_violations,
        :partial_candidate_first_violation_hz,
        :partial_candidate_last_violation_hz,
        :partial_context_invariant_violations,
        :partial_context_first_violation_hz,
        :partial_context_last_violation_hz,
        :partial_Z_relative_frobenius,
        :partial_Y_relative_frobenius,
        :partial_R_relative_frobenius,
        :partial_L_relative_frobenius,
        :partial_G_relative_frobenius,
        :partial_C_relative_frobenius,
        :exception_type,
        :message,
        :artifact
    )
    open(temporary, "w") do io
        println(io, join(string.(names), '\t'))
        for record in failures
            values = replace.(
                string.((getproperty(record, name) for name in names)),
                '\t' => ' ',
                '\n' => ' '
            )
            println(io, join(values, '\t'))
        end
    end
    mv(temporary, path; force = true)
    return path
end

function write_outputs(rows, skipped, failures, coverage, reference_grids)
    mkpath(OUTPUT_ROOT)
    summary_path = write_tsv(joinpath(OUTPUT_ROOT, "summary.tsv"), rows)
    skipped_path = write_tsv(joinpath(OUTPUT_ROOT, "skipped.tsv"), skipped)
    failure_path = write_failure_tsv(
        joinpath(OUTPUT_ROOT, "failures.tsv"),
        failures
    )
    coverage_path = write_tsv(joinpath(OUTPUT_ROOT, "coverage.tsv"), coverage)
    jld_path = joinpath(OUTPUT_ROOT, "summary.jld2")
    temporary = jld_path * ".new"
    JLD2.jldsave(
        temporary;
        schema_version = ARTIFACT_SCHEMA_VERSION,
        kind = :fem_analytical_fullband_gauntlet_summary,
        corpus_policy = CORPUS_POLICY,
        formula_treatment = FORMULA_TREATMENT,
        diagnostic_policy = DIAGNOSTIC_POLICY,
        reference_grids,
        internal_impedance = INTERNAL_FORMULA,
        insulation_admittance = :Gustavsen2013,
        semicon_admittance = :FEMLossless,
        registered_earth_impedance = collect(EARTH_IMPEDANCE.formulas()),
        registered_earth_admittance = collect(EARTH_ADMITTANCE.formulas()),
        rows,
        skipped,
        failures,
        coverage,
        recorded_at_utc = string(now(UTC))
    )
    mv(temporary, jld_path; force = true)
    open(jld_path * ".sha256", "w") do io
        println(io, bytes2hex(sha256(read(jld_path))), "  ", basename(jld_path))
    end
    return summary_path, skipped_path, failure_path, coverage_path, jld_path
end

function main(args = ARGS)
    coverage = catalogue()
    variants = variant.(coverage)
    length(unique(getproperty.(variants, :id))) == 39 || error(
        "catalogue variant IDs are not unique"
    )
    available = sort!(collect(keys(case_index())); by = string)
    requested = isempty(args) ? available : Symbol.(args)
    all(case_id -> case_id in available, requested) || error(
        "unknown case requested; available cases are $(join(available, ", "))"
    )

    rows = NamedTuple[]
    skipped = NamedTuple[]
    failures = NamedTuple[]
    reference_grids = NamedTuple[]
    println(
        "PLAN\tcases=", length(requested),
        "\tregistered_z=", length(EARTH_IMPEDANCE.formulas()),
        "\tregistered_y=", length(EARTH_ADMITTANCE.formulas()),
        "\tvariants_per_case=", length(variants),
        "\tfrequencies=101\trange_hz=case_reference_grid",
        "\tpolicy=", CORPUS_POLICY,
        "\ttreatment=", FORMULA_TREATMENT,
        "\tdiagnostics=", DIAGNOSTIC_POLICY
    )
    flush(stdout)
    for case_id in requested
        model = reference_case(case_id)
        push!(reference_grids,
            (;
                case = string(case_id),
                count = length(model.problem.frequencies),
                first_hz = first(model.problem.frequencies),
                last_hz = last(model.problem.frequencies)
            ))
        reference_path, reference, reference_sha256 = load_reference(case_id, model)
        println("BEGIN_CASE\t", case_id, "\tterminals=", length(model.port_order))
        flush(stdout)
        for selected in variants
            path = artifact_path(case_id, selected)
            if valid_existing(path, model, selected, reference_sha256)
                write_hash(path)
                recovered = recover_existing(path, selected, reference)
                recovered.row === nothing || push!(rows, recovered.row)
                recovered.skipped === nothing || push!(skipped, recovered.skipped)
                recovered.failure === nothing || push!(failures, recovered.failure)
                println("REUSE\t", case_id, "\t", selected.id)
                flush(stdout)
                continue
            end
            if isfile(path)
                archived = archive_stale_artifact(path)
                println(
                    "STALE_ARCHIVED\t", case_id,
                    "\t", selected.id,
                    "\t", archived
                )
                flush(stdout)
            end
            reason = case_skip_reason(model, selected)
            if reason !== nothing
                path = record_skip(
                    case_id,
                    model,
                    selected,
                    reason,
                    reference_path,
                    reference_sha256,
                    reference
                )
                push!(skipped,
                    (
                        case = string(case_id),
                        variant = string(selected.id),
                        kind = string(selected.kind),
                        formula = string(selected.identifier),
                        reason,
                        artifact = path
                    ))
                println("SKIP\t", case_id, "\t", selected.id, "\t", reason)
                flush(stdout)
                continue
            end

            started = time()
            try
                result = compute(model.problem, formulation(selected))
                size(result.Z.values) == size(reference["Z"]) || error(
                    "analytical and FEM Z shapes differ"
                )
                size(result.Y.values) == size(reference["Y"]) || error(
                    "analytical and FEM Y shapes differ"
                )
                result.f == model.problem.frequencies || error(
                    "analytical result changed the frequency vector"
                )
                metrics = result_metrics(result, reference)
                elapsed = time() - started
                path = record_result(
                    case_id,
                    model,
                    selected,
                    result,
                    metrics,
                    elapsed,
                    reference_path,
                    reference_sha256,
                    reference
                )
                push!(rows,
                    summary_row(
                        case_id,
                        selected,
                        metrics,
                        model.problem.frequencies,
                        elapsed,
                        path,
                        reference
                    ))
                candidate_violations, context_violations = partition_violations(
                    selected,
                    metrics
                )
                state = isempty(candidate_violations) ?
                        "OBSERVED" : "OBSERVED_BREAK"
                println(
                    state, "\t", case_id,
                    "\t", selected.id,
                    "\tZ=", metrics.Z.aggregate_relative,
                    "\tY=", metrics.Y.aggregate_relative,
                    "\tcandidate_violations=", length(candidate_violations),
                    "\tcontext_violations=", length(context_violations),
                    "\telapsed_s=", round(elapsed; digits = 3)
                )
                flush(stdout)
            catch exception
                diagnostic = diagnose_execution_failure(
                    case_id,
                    selected,
                    model.problem.frequencies,
                    reference
                )
                elapsed = time() - started
                path = record_execution_failure(
                    case_id,
                    model,
                    selected,
                    exception,
                    diagnostic,
                    elapsed,
                    reference_path,
                    reference_sha256,
                    reference
                )
                push!(failures,
                    execution_breakdown_summary(
                        case_id,
                        selected,
                        diagnostic.frequency_failures,
                        typeof(exception),
                        sprint(showerror, exception),
                        path;
                        successful_frequencies = diagnostic.successful_frequencies,
                        partial_metrics = diagnostic.partial_metrics,
                        candidate_violations = diagnostic.candidate_violations,
                        context_violations = diagnostic.context_violations
                    ))
                println(
                    stderr,
                    "OBSERVED_EXECUTION_BREAKDOWN\t", case_id,
                    "\t", selected.id,
                    "\tformula=", selected.identifier,
                    "\tfailed_frequencies=", length(diagnostic.frequency_failures),
                    "\tsuccessful_frequencies=", length(
                        diagnostic.successful_frequencies
                    ),
                    "\t", sprint(showerror, exception)
                )
                flush(stderr)
            end
        end
    end

    summary_path, skipped_path, failure_path, coverage_path, jld_path = write_outputs(
        rows,
        skipped,
        failures,
        coverage,
        reference_grids
    )
    println(
        "COMPLETE\tevaluated=", length(rows),
        "\tskipped=", length(skipped),
        "\texecution_breakdowns=", length(failures),
        "\tplanned=", length(requested) * length(variants),
        "\tsummary=", summary_path,
        "\tskips=", skipped_path,
        "\tfailures=", failure_path,
        "\tcoverage=", coverage_path,
        "\tjld2=", jld_path
    )
    return isempty(failures) ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main())
