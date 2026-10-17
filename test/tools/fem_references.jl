# The FEM references against this branch. Run
# `julia --project=test test/tools/fem_references.jl [FEM_CHECKOUT]`
# from the repository root. FEM_CHECKOUT defaults to the FEM investigation worktree next to
# this repository, `LinecableModels-fem-investigation`.
#
# Each case transcribes a script of `FEM_CHECKOUT/test/manual/cases/` with the API of this
# branch, in `known_case_problem` of `test/support/scenarios.jl`. It keeps the geometry,
# materials, placements, connections, earth, temperature and frequencies of the script. The shunt model `:coaxial` of the FEM checkout is
# `:equivalent` on this branch. The newest run folder of each case,
# `FEM_CHECKOUT/.linecablemodels/fem/<case>/<run>/`, holds `matrices.csv`. That file lists
# every Z and Y entry of the analytical run and of the FEM run of that checkout.
#
# The script compares the analytical result of this branch at the shared frequencies.
# - Against the FEM columns, it applies the metric of the user. The error e_Z divides
#   max|ΔZ_ij| by max|Z_ii| of this branch. The error e_Y does the same for Y. Each must
#   stay within 1 % at every frequency.
# - Against the analytical columns of the FEM checkout, it reports the largest relative
#   difference of one entry, |ΔZ_ij| / |Z_ij|.
#
# The tables go to `.linecablemodels/checkpoint/fem/comparison.md`.
const REPOSITORY = dirname(dirname(@__DIR__))
const OUTPUT = joinpath(REPOSITORY, ".linecablemodels", "checkpoint", "fem")
const ACCEPTANCE = 0.01

using LineCableModels
using LinearAlgebra: diag
using Printf: @sprintf

# `known_case_problem` of the test scenarios builds each case, at the 101 frequencies of
# the case scripts.
Base.include(@__MODULE__, joinpath(REPOSITORY, "test", "support", "scenarios.jl"))
using .CurrentScenarios: known_case_problem
const CASES = (:two_bare_wires, :two_insulated_wires, :solid_1000mm2_single,
    :cable_525kv_land_no_armour_dc_bipole, :cable_320kv_no_armour_dc_bipole,
    :cable_380kv_no_armour_ac_flat, :cable_132kv_cigre_tb880_case0_630cu_trefoil)

const FORMULATION = Formulation(earth_impedance = :unified, earth_admittance = :unified,
    shunt_model = :equivalent,
    options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))

# The conductor labels and the matrices of `matrices.csv`, by source ("analytical" or
# "fem") and quantity ("Z" or "Y"), then by frequency.
function reference(path)
    labels = String[]
    matrices = Dict{Tuple{String, String}, Dict{Float64, Matrix{ComplexF64}}}()
    for line in eachline(path)
        if startswith(line, "# Conductors in row/column order: ")
            labels = [last(split(entry, " = "))
                      for entry in split(line[(length("# Conductors in row/column order: ") + 1):end], ", ")]
        elseif !startswith(line, "#") && !startswith(line, "source,")
            source, frequency, quantity, row, column, real_part, imaginary_part = split(line, ',')
            n = length(labels)
            byfrequency = get!(matrices, (source, quantity), Dict{Float64, Matrix{ComplexF64}}())
            matrix = get!(() -> zeros(ComplexF64, n, n), byfrequency, parse(Float64, frequency))
            matrix[parse(Int, row), parse(Int, column)] =
                complex(parse(Float64, real_part), parse(Float64, imaginary_part))
        end
    end
    return labels, matrices
end

# The newest run folder of `case` that holds `matrices.csv`. The folder names are
# timestamps. The last name in sorted order is the newest.
function newest(checkout, case)
    root = joinpath(checkout, ".linecablemodels", "fem", String(case))
    runs = sort!(filter(name -> isfile(joinpath(root, name, "matrices.csv")), readdir(root)))
    isempty(runs) && error("no run folder with matrices.csv under $root")
    return joinpath(root, last(runs))
end

# The metric of the user divides the largest entry difference by the largest diagonal
# entry of `scale`. It comes back with the entry where the difference occurs.
function scaled(a, b, scale)
    difference, worst = findmax(abs.(a .- b))
    return difference / maximum(abs, diag(scale)), Tuple(worst)
end

# The largest relative difference of one entry, and the entry where it occurs. The loop
# skips an entry that has the same value on both sides, zero included. Without any
# difference, the entry is (0, 0).
function entrywise(a, b)
    largest, worst = 0.0, (0, 0)
    for index in CartesianIndices(a)
        a[index] == b[index] && continue
        relative = abs(a[index] - b[index]) / abs(b[index])
        relative > largest && ((largest, worst) = (relative, Tuple(index)))
    end
    return largest, worst
end

percent(x) = @sprintf("%.3f", 100x)
entry((i, j)) = i == 0 ? "—" : "($i,$j)"

function compare(case, problem, checkout, lines)
    run = newest(checkout, case)
    labels, matrices = reference(joinpath(run, "matrices.csv"))
    seconds = @elapsed result = compute(problem, FORMULATION)
    coordinates = details(result).data.coordinates
    coordinates == labels || error("$case: conductor order $coordinates differs from $labels")
    Z_branch, Y_branch = observe(result, Z), observe(result, Y)
    column = Dict(f => k for (k, f) in enumerate(observe(result, LineCableModels.frequencies)))

    push!(lines, "## $case", "", "Reference run: `$(relpath(run, checkout))` in `$checkout`.",
        "Conductors: $(join(labels, ", ")). This branch computed in $(@sprintf("%.1f", seconds)) s.",
        "", "### a) Against the FEM columns", "",
        "e_Z and e_Y in %, scaled by the diagonal of this branch. The last two columns give " *
        "the same metric for the analytical columns of the FEM checkout.", "",
        "| f [Hz] | e_Z | worst Z | e_Y | worst Y | ≤ 1 % | e_Z, FEM checkout | e_Y, FEM checkout |",
        "| ---: | ---: | --- | ---: | --- | --- | ---: | ---: |")
    passed = true
    for f in sort!(collect(keys(matrices[("fem", "Z")])))
        haskey(column, f) || error("$case: FEM frequency $f is not on the analytical grid")
        k = column[f]
        Z_fem, Y_fem = matrices[("fem", "Z")][f], matrices[("fem", "Y")][f]
        Z_old, Y_old = matrices[("analytical", "Z")][f], matrices[("analytical", "Y")][f]
        e_Z, Z_worst = scaled(Z_branch[:, :, k], Z_fem, Z_branch[:, :, k])
        e_Y, Y_worst = scaled(Y_branch[:, :, k], Y_fem, Y_branch[:, :, k])
        verdict = max(e_Z, e_Y) <= ACCEPTANCE
        passed &= verdict
        push!(lines, "| $(@sprintf("%.4g", f)) | $(percent(e_Z)) | $(entry(Z_worst)) | " *
            "$(percent(e_Y)) | $(entry(Y_worst)) | $(verdict ? "yes" : "**no**") | " *
            "$(percent(first(scaled(Z_old, Z_fem, Z_old)))) | " *
            "$(percent(first(scaled(Y_old, Y_fem, Y_old)))) |")
    end

    push!(lines, "", "### b) Against the analytical columns of the FEM checkout", "",
        "The largest relative difference of one entry, |Δ| / |reference|, at each frequency.",
        "", "| f [Hz] | Z | worst Z | Y | worst Y |", "| ---: | ---: | --- | ---: | --- |")
    largest = Dict("Z" => (0.0, 0.0, (0, 0)), "Y" => (0.0, 0.0, (0, 0)))
    for f in sort!(collect(keys(matrices[("analytical", "Z")])))
        haskey(column, f) || error("$case: analytical frequency $f is missing on this branch")
        k = column[f]
        r_Z, Z_worst = entrywise(Z_branch[:, :, k], matrices[("analytical", "Z")][f])
        r_Y, Y_worst = entrywise(Y_branch[:, :, k], matrices[("analytical", "Y")][f])
        r_Z > first(largest["Z"]) && (largest["Z"] = (r_Z, f, Z_worst))
        r_Y > first(largest["Y"]) && (largest["Y"] = (r_Y, f, Y_worst))
        push!(lines, "| $(@sprintf("%.4g", f)) | $(@sprintf("%.2e", r_Z)) | $(entry(Z_worst)) | " *
            "$(@sprintf("%.2e", r_Y)) | $(entry(Y_worst)) |")
    end
    shared = length(matrices[("analytical", "Z")])
    summary = [iszero(largest[q][1]) ? "$q identical at all $shared frequencies" :
               @sprintf("%s %.2e at %.4g Hz, entry %s", q, largest[q][1], largest[q][2],
                   entry(largest[q][3])) for q in ("Z", "Y")]
    push!(lines, "", "Largest: $(join(summary, "; ")).", "")
    println(case, ": computed in ", @sprintf("%.1f", seconds), " s; FEM ",
        passed ? "within" : "outside", " 1 % at every shared frequency; ",
        join(summary, "; "))
    return seconds
end

function main(arguments)
    length(arguments) <= 1 ||
        error("Usage: julia --project=test test/tools/fem_references.jl [FEM_CHECKOUT]")
    checkout = abspath(get(arguments, 1,
        joinpath(dirname(REPOSITORY), "LinecableModels-fem-investigation")))
    isdir(checkout) || error("FEM checkout not found: $checkout")
    lines = ["# FEM references against this branch", ""]
    started = time()
    seconds = Dict(case => compare(case, known_case_problem(case), checkout, lines) for case in CASES)
    total = time() - started
    push!(lines, "## Run times", "", "| Case | Computation [s] |", "| --- | ---: |",
        ["| $case | $(@sprintf("%.1f", seconds[case])) |" for case in CASES]...,
        "", "The whole comparison took $(@sprintf("%.1f", total)) s. The first case includes compilation.")
    mkpath(OUTPUT)
    write(joinpath(OUTPUT, "comparison.md"), join(lines, '\n') * '\n')
    println("Tables: ", relpath(joinpath(OUTPUT, "comparison.md"), REPOSITORY))
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
