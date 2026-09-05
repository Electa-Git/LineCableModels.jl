using Dates
using JLD2
using LinearAlgebra
using LineCableModels
using SHA
using Statistics

const ROOT = joinpath(
    pkgdir(LineCableModels),
    ".linecablemodels",
    "fem",
    "gauntlet",
    "analytical_fullband_corpus"
)
const SOURCE = joinpath(ROOT, "summary.jld2")

relevant_error(row) = row.kind == "earth_impedance" ?
                      row.Z_relative_frobenius : row.Y_relative_frobenius

function case_rows(rows, skipped)
    output = NamedTuple[]
    for case_id in sort!(unique(getproperty.(rows, :case)))
        selected = filter(row -> row.case == case_id, rows)
        impedance = filter(row -> row.kind == "earth_impedance", selected)
        admittance = filter(row -> row.kind == "earth_admittance", selected)
        push!(output, (;
            case = case_id,
            fem_reference_status = first(selected).fem_reference_status,
            fem_reference_observations = first(selected).fem_reference_observations,
            applicable_variants = length(selected),
            skipped_variants = count(row -> row.case == case_id, skipped),
            candidate_break_variants = count(
                row -> row.candidate_invariant_violations > 0, selected
            ),
            context_break_variants = count(
                row -> row.context_invariant_violations > 0, selected
            ),
            minimum_Z_relative = minimum(
                getproperty.(impedance, :Z_relative_frobenius)
            ),
            maximum_Z_relative = maximum(
                getproperty.(impedance, :Z_relative_frobenius)
            ),
            minimum_Y_relative = minimum(
                getproperty.(admittance, :Y_relative_frobenius)
            ),
            maximum_Y_relative = maximum(
                getproperty.(admittance, :Y_relative_frobenius)
            )
        ))
    end
    return output
end

function formula_rows(rows)
    output = NamedTuple[]
    identifiers = sort!(unique((row.kind, row.formula) for row in rows))
    for (kind, formula) in identifiers
        selected = filter(
            row -> row.kind == kind && row.formula == formula,
            rows
        )
        errors = relevant_error.(selected)
        worst = argmax(errors)
        push!(output, (;
            kind,
            formula,
            applicable_cases = length(selected),
            skipped_cases = 7 - length(selected),
            candidate_break_cases = count(
                row -> row.candidate_invariant_violations > 0, selected
            ),
            candidate_invariant_violations = sum(
                getproperty.(selected, :candidate_invariant_violations)
            ),
            context_break_cases = count(
                row -> row.context_invariant_violations > 0, selected
            ),
            minimum_relevant_relative = minimum(errors),
            median_relevant_relative = median(errors),
            maximum_relevant_relative = maximum(errors),
            maximum_error_case = selected[worst].case
        ))
    end
    return output
end

function break_rows(rows)
    output = NamedTuple[]
    selected = sort!(
        filter(row -> row.candidate_invariant_violations > 0, rows);
        by = row -> (row.case, row.kind, row.formula)
    )
    for row in selected
        push!(output, (;
            case = row.case,
            kind = row.kind,
            formula = row.formula,
            quantities = row.candidate_violation_quantities,
            violations = row.candidate_invariant_violations,
            first_violation_hz = row.candidate_first_violation_hz,
            last_violation_hz = row.candidate_last_violation_hz,
            relevant_relative_error = relevant_error(row),
            context_violations = row.context_invariant_violations,
            artifact = row.artifact
        ))
    end
    return output
end

function violation_mode_rows(rows)
    output = NamedTuple[]
    selected = sort!(
        filter(
            row -> row.candidate_invariant_violations > 0 ||
                   row.context_invariant_violations > 0,
            rows
        );
        by = row -> (row.case, row.kind, row.formula)
    )
    for row in selected
        document = JLD2.load(row.artifact)
        for (attribution, field) in (
                ("candidate", "candidate_violations"),
                ("context", "context_violations")
            )
            violations = get(document, field, NamedTuple[])
            for quantity in sort!(unique(getproperty.(violations, :quantity));
                    by = string)
                matching = filter(
                    violation -> violation.quantity == quantity,
                    violations
                )
                severities = map(matching) do violation
                    isfinite(violation.minimum_eigenvalue) ?
                    -violation.minimum_eigenvalue / violation.tolerance : NaN
                end
                finite = findall(isfinite, severities)
                worst = isempty(finite) ? 0 :
                        finite[argmax(severities[finite])]
                push!(output, (;
                    case = row.case,
                    kind = row.kind,
                    formula = row.formula,
                    attribution,
                    quantity = string(quantity),
                    violations = length(matching),
                    first_violation_hz = minimum(getproperty.(
                        matching, :frequency_hz
                    )),
                    last_violation_hz = maximum(getproperty.(
                        matching, :frequency_hz
                    )),
                    worst_frequency_hz = worst == 0 ? NaN :
                                         matching[worst].frequency_hz,
                    worst_minimum_eigenvalue = worst == 0 ? NaN :
                                               matching[worst].minimum_eigenvalue,
                    worst_tolerance = worst == 0 ? NaN :
                                      matching[worst].tolerance,
                    worst_severity = worst == 0 ? NaN : severities[worst],
                    artifact = row.artifact
                ))
            end
        end
    end
    return output
end

function relevant_matrices(rows)
    matrices = Dict{Tuple{String, String, String}, Array{ComplexF64, 3}}()
    for row in rows
        document = JLD2.load(row.artifact)
        quantity = row.kind == "earth_impedance" ? "Z" : "Y"
        matrices[(row.case, row.kind, row.formula)] = document[quantity]
    end
    return matrices
end

function near_equivalence_rows(rows, matrices)
    output = NamedTuple[]
    for kind in ("earth_impedance", "earth_admittance")
        selected = filter(row -> row.kind == kind, rows)
        formulas = sort!(unique(getproperty.(selected, :formula)))
        for first_index in eachindex(formulas)
            for second_index in (first_index + 1):length(formulas)
                first_formula = formulas[first_index]
                second_formula = formulas[second_index]
                first_cases = Set(
                    row.case for row in selected if row.formula == first_formula
                )
                second_cases = Set(
                    row.case for row in selected if row.formula == second_formula
                )
                common = sort!(collect(intersect(first_cases, second_cases)))
                isempty(common) && continue
                distances = map(common) do case_id
                    first_matrix = matrices[(case_id, kind, first_formula)]
                    second_matrix = matrices[(case_id, kind, second_formula)]
                    norm(first_matrix - second_matrix) / max(
                        norm(first_matrix), norm(second_matrix), eps(Float64)
                    )
                end
                maximum(distances) <= 1.0e-8 || continue
                push!(output, (;
                    kind,
                    first_formula,
                    second_formula,
                    common_cases = length(common),
                    exact_cases = count(iszero, distances),
                    median_relative_distance = median(distances),
                    maximum_relative_distance = maximum(distances),
                    cases = join(common, ",")
                ))
            end
        end
    end
    return output
end

function write_tsv(path, rows)
    isempty(rows) && error("cannot write empty report $path")
    temporary = path * ".new"
    names = keys(first(rows))
    open(temporary, "w") do io
        println(io, join(string.(names), '\t'))
        for row in rows
            values = replace.(
                string.((getproperty(row, name) for name in names)),
                '\t' => ' ',
                '\n' => ' '
            )
            println(io, join(values, '\t'))
        end
    end
    mv(temporary, path; force = true)
    return path
end

function write_hash(path)
    open(path * ".sha256", "w") do io
        println(io, bytes2hex(sha256(read(path))), "  ", basename(path))
    end
    return path
end

document = JLD2.load(SOURCE)
document["schema_version"] == 6 || error("unsupported corpus schema")
rows = document["rows"]
skipped = document["skipped"]
failures = document["failures"]
isempty(failures) || error("corpus contains execution breakdowns")

cases = case_rows(rows, skipped)
formulas = formula_rows(rows)
breaks = break_rows(rows)
violation_modes = violation_mode_rows(rows)
equivalences = near_equivalence_rows(rows, relevant_matrices(rows))

case_path = write_hash(write_tsv(joinpath(ROOT, "case_signatures.tsv"), cases))
formula_path = write_hash(write_tsv(
    joinpath(ROOT, "formula_signatures.tsv"), formulas
))
break_path = write_hash(write_tsv(
    joinpath(ROOT, "invariant_breaks.tsv"), breaks
))
violation_mode_path = write_hash(write_tsv(
    joinpath(ROOT, "invariant_modes.tsv"), violation_modes
))
equivalence_path = write_hash(write_tsv(
    joinpath(ROOT, "near_equivalences.tsv"), equivalences
))
jld_path = joinpath(ROOT, "catalogue_report.jld2")
temporary = jld_path * ".new"
JLD2.jldsave(
    temporary;
    schema_version = 1,
    kind = :fem_analytical_fullband_catalogue_report,
    corpus_sha256 = bytes2hex(sha256(read(SOURCE))),
    cases,
    formulas,
    breaks,
    violation_modes,
    equivalences,
    recorded_at_utc = string(now(UTC))
)
mv(temporary, jld_path; force = true)
write_hash(jld_path)
println(
    "COMPLETE\tcases=", length(cases),
    "\tformulas=", length(formulas),
    "\tbreaks=", length(breaks),
    "\tviolation_modes=", length(violation_modes),
    "\tnear_equivalences=", length(equivalences),
    "\tcase_report=", case_path,
    "\tformula_report=", formula_path,
    "\tbreak_report=", break_path,
    "\tviolation_mode_report=", violation_mode_path,
    "\tequivalence_report=", equivalence_path,
    "\tjld=", jld_path
)
