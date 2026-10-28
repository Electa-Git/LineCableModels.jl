# Summary tables and figures for run_two_bare_wires_placements.jl: FEM against
# the Unified earth-return formula, labelled "Proposed".
#
# Run it in the REPL that ran the runner, whose results are reused, or in a fresh
# session from the repository root. A fresh session includes the runner first,
# which reloads the kept FEM runs (resume_run_directory = :latest) without solving:
#     julia -i --project=test test/manual/calculations/summarize_two_bare_wires_placements.jl
#
# Grouping. Figures never mix placement classes: overhead pair (both in air),
# buried pair (both in earth) and overhead-buried pair (1 air, 2 earth). The
# 1 earth, 2 air placement is the same pair mirrored, so it is not plotted; it
# checks mirror consistency in Table III. Per class, one figure family pins the
# radius and overlays the resistivities, another pins the resistivity and
# overlays the radii. Each family has R, X, G and B pages in the 2x2 matrix
# layout. Color encodes the swept value (lightest = smallest), line style the
# model (FEM solid, Proposed dashed).
#
# Tables use the scaled deviation per frequency, e_Z = max|ΔZ_ij|/max|Z_ii| and
# e_Y alike, normalized by Proposed, and report its maximum within each band.
# SVG (publication theme), PDF (when Inkscape is installed), CSV and LaTeX files
# go to a new folder under .linecablemodels/reports/two-bare-wires.

using LineCableModels
using LinearAlgebra
using Printf
using Dates
using CairoMakie

@isdefined(fem_results) && @isdefined(unified_results) ||
    include(joinpath(@__DIR__, "run_two_bare_wires_placements.jl"))

pinned_radius = 4.25e-2  # [m], resistivity sweeps
pinned_rho = 100.0       # [Ω·m], radius sweeps
bands = ["0.1 Hz–1 kHz" => (0.1, 1e3), "1 kHz–1 MHz" => (1e3, 1e6),
    "1 MHz–100 MHz" => (1e6, 1e8)]
classes = [
    (key = "overhead", label = "Overhead", placement = "both in air"),
    (key = "buried", label = "Buried", placement = "both in earth"),
    (key = "mixed", label = "Overhead-buried", placement = "1 air, 2 earth")
]
mirror_placement = "1 earth, 2 air"
ramp = ["#6da7ec", "#2a78d6", "#184f95", "#0d366b"]  # ordinal blue, validated on white
summary_directory = joinpath(pwd(), ".linecablemodels", "reports", "two-bare-wires",
    Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(summary_directory)

fmt(x) = isinteger(x) ? string(Int(x)) : string(x)
# Points vary radius fastest, then placement, then resistivity (see the runner).
function point(r, placement, ρ)
    i = LinearIndices((length(radius_values), length(placements), length(rho_values)))[
        findfirst(==(r), radius_values), findfirst(==(placement), placements),
        findfirst(==(ρ), rho_values)]
    points[i] == "$placement, r = $(100r) cm, ρ = $ρ Ω·m" ||
        error("results do not follow the two-wire runner order: $(points[i])")
    return i
end

scaled(a, e) = maximum(abs, a - e) / maximum(abs, diag(e))
deviation(quantity) = [scaled(observe(f, quantity)[:, :, k], observe(p, quantity)[:, :, k])
                       for (f, p) in zip(fem_results, unified_results),
                       k in eachindex(frequencies)]
e_Z, e_Y = deviation(Z), deviation(Y)
masks = [[lo <= f && (f < hi || f == hi == last(frequencies)) for f in frequencies]
         for (_, (lo, hi)) in bands]
worst(e, i, mask = trues(length(frequencies))) = 100maximum(e[i, mask])

# Mirror consistency: 1 earth, 2 air equals 1 air, 2 earth with conductors swapped.
swap = [0 1; 1 0]
mirror(results, quantity) = maximum(Iterators.product(radius_values, rho_values)) do (r, ρ)
    a = observe(results[point(r, classes[3].placement, ρ)], quantity)
    b = observe(results[point(r, mirror_placement, ρ)], quantity)
    maximum(k -> scaled(a[:, :, k], swap * b[:, :, k] * swap), eachindex(frequencies))
end
mirror_rows = [(model, 100mirror(results, Z), 100mirror(results, Y))
               for (model, results) in ("FEM" => fem_results, "Proposed" => unified_results)]

println("\nTable I. Maximum scaled deviation of FEM from Proposed, 0.1 Hz-100 MHz, e_Z/e_Y [%]")
@printf("%-16s %-10s", "placement", "ρ [Ω·m]")
foreach(r -> @printf("%16s", "r = $(fmt(100r)) cm"), radius_values); println()
for class in classes, ρ in rho_values
    @printf("%-16s %-10s", class.label, fmt(ρ))
    foreach(radius_values) do r
        i = point(r, class.placement, ρ)
        @printf("%16s", @sprintf("%.2f/%.2f", worst(e_Z, i), worst(e_Y, i)))
    end
    println()
end
println("\nTable II. Worst scaled deviation per band over all radii and resistivities, e_Z/e_Y [%]")
@printf("%-16s", "placement"); foreach(b -> @printf("%18s", first(b)), bands)
@printf("%14s\n", "points <= 1 %")
for class in classes
    ids = [point(r, class.placement, ρ) for r in radius_values, ρ in rho_values]
    @printf("%-16s", class.label)
    foreach(masks) do mask
        @printf("%18s", @sprintf("%.2f/%.2f", maximum(i -> worst(e_Z, i, mask), ids),
            maximum(i -> worst(e_Y, i, mask), ids)))
    end
    @printf("%11d/%d\n", count(i -> max(worst(e_Z, i), worst(e_Y, i)) <= 1, ids), length(ids))
end
println("\nTable III. Mirror consistency, max scaled |M(1 air, 2 earth) - S M(1 earth, 2 air) S| [%]")
foreach(row -> @printf("%-10s e_Z %.2e   e_Y %.2e\n", row...), mirror_rows)

open(joinpath(summary_directory, "deviations.csv"), "w") do io
    println(io, "placement,radius_m,rho_ohm_m,band,e_Z_percent,e_Y_percent")
    for class in classes, r in radius_values, ρ in rho_values,
        (mask, band) in zip([[trues(length(frequencies))]; masks],
            ["0.1 Hz–100 MHz"; first.(bands)])

        i = point(r, class.placement, ρ)
        @printf(io, "%s,%.6g,%.6g,%s,%.6g,%.6g\n", class.placement, r, ρ, band,
            worst(e_Z, i, mask), worst(e_Y, i, mask))
    end
end
open(joinpath(summary_directory, "mirror.csv"), "w") do io
    println(io, "model,e_Z_percent,e_Y_percent")
    foreach(row -> @printf(io, "%s,%.6g,%.6g\n", row...), mirror_rows)
end

# IEEEtran tables: caption above, \hline rules, values in percent.
open(joinpath(summary_directory, "tables.tex"), "w") do io
    print(io, """
    \\begin{table*}[!t]
    \\renewcommand{\\arraystretch}{1.2}
    \\caption{Maximum Scaled Deviation of FEM From Proposed, 0.1~Hz--100~MHz,
    \$e_Z/e_Y\$ (\\%). \$e_Z(f)=\\max_{ij}|Z^{\\mathrm{FEM}}_{ij}-Z^{\\mathrm{P}}_{ij}|/\\max_i|Z^{\\mathrm{P}}_{ii}|\$, \$e_Y\$ Alike}
    \\label{tab:fem-proposed-max}
    \\centering
    \\begin{tabular}{ll$(repeat("c", length(radius_values)))}
    \\hline
    Placement & \$\\rho\$ (\$\\Omega\\,\$m) & \\multicolumn{$(length(radius_values))}{c}{Conductor radius \$r\$ (cm)}\\\\
     & & $(join(fmt.(100 .* radius_values), " & "))\\\\
    \\hline
    """)
    for class in classes
        for (j, ρ) in enumerate(rho_values)
            cells = [@sprintf("%.2f/%.2f", worst(e_Z, i), worst(e_Y, i))
                     for i in (point(r, class.placement, ρ) for r in radius_values)]
            println(io, j == 1 ? class.label : "", " & ", fmt(ρ), " & ", join(cells, " & "), "\\\\")
        end
        class === last(classes) || println(io, "\\hline")
    end
    print(io, """
    \\hline
    \\end{tabular}
    \\end{table*}

    \\begin{table}[!t]
    \\renewcommand{\\arraystretch}{1.2}
    \\setlength{\\tabcolsep}{4pt}
    \\caption{Worst Scaled Deviation of FEM From Proposed per Band, All Radii and Resistivities, \$e_Z/e_Y\$ (\\%)}
    \\label{tab:fem-proposed-bands}
    \\centering
    \\begin{tabular}{l$(repeat("c", length(bands)))}
    \\hline
    Placement & $(join((replace(first(b), "–" => "--", " Hz" => "~Hz", " kHz" => "~kHz", " MHz" => "~MHz") for b in bands), " & "))\\\\
    \\hline
    """)
    for class in classes
        ids = [point(r, class.placement, ρ) for r in radius_values, ρ in rho_values]
        cells = [@sprintf("%.2f/%.2f", maximum(i -> worst(e_Z, i, mask), ids),
                     maximum(i -> worst(e_Y, i, mask), ids)) for mask in masks]
        println(io, class.label, " & ", join(cells, " & "), "\\\\")
    end
    print(io, """
    \\hline
    \\end{tabular}
    \\end{table}

    \\begin{table}[!t]
    \\renewcommand{\\arraystretch}{1.2}
    \\caption{Mirror Consistency of the Overhead-Buried Pair, Maximum Scaled Deviation Between Mirrored Placements (\\%)}
    \\label{tab:mirror}
    \\centering
    \\begin{tabular}{lcc}
    \\hline
    Model & \$e_Z\$ & \$e_Y\$\\\\
    \\hline
    """)
    foreach(row -> @printf(io, "%s & %.1e & %.1e\\\\\n", row...), mirror_rows)
    print(io, """
    \\hline
    \\end{tabular}
    \\end{table}
    """)
end

# Figures: one family per placement class and sweep, R/X/G/B pages of the
# 2x2 matrix, sized for an IEEE double-column width (688 units = 7.16 in).
# Signed log axes label their ticks as signed powers of ten.
superscript(n) = map(c -> c == '-' ? '⁻' : "⁰¹²³⁴⁵⁶⁷⁸⁹"[nextind("⁰¹²³⁴⁵⁶⁷⁸⁹", 0, c - '0' + 1)], string(n))
signed_decades(values) = map(values) do v
    v == 0 && return "0"
    e = round(Int, log10(abs(v)))
    isapprox(abs(v), 10.0^e; rtol = 1e-6) || return @sprintf("%.3g", v)
    return (v < 0 ? "−" : "") * "10" * superscript(e)
end
inkscape = Sys.which("inkscape")
summary_plots = Dict{Tuple{String, Symbol}, Any}()
for class in classes, sweep in (:rho, :radius)
    values, ids, legend = if sweep === :rho
        rho_values, [point(pinned_radius, class.placement, ρ) for ρ in rho_values],
        ρ -> "ρ = $(fmt(ρ)) Ω·m"
    else
        radius_values, [point(r, class.placement, pinned_rho) for r in radius_values],
        r -> "r = $(fmt(100r)) cm"
    end
    pages = Makie.plot(fem_results[ids]..., unified_results[ids]...;
        ydata = (Z, Y),
        layout = (2, 2),
        length_unit = :base,
        series_labels = [m * ", " * legend(v) for m in ("FEM", "Proposed") for v in values],
        series_attributes = [(color = c, linestyle = s, linewidth = 1.5, marker = nothing)
                             for s in (:solid, :dash) for c in ramp],
        xscale = :log10,
        yscale = :log10,
        ytickformat = signed_decades,
        backend = :cairo,
        display_plot = false,
        controls = false,
        fig_size = (688, 600),
        legend_position = :bottom,
        legend_attributes = (; orientation = :horizontal, nbanks = 2),
        legend_overflow = :show_all
    )
    summary_plots[(class.key, sweep)] = pages
    for (page, quantity) in zip(pages, ("R", "X", "G", "B"))
        svg = export_svg(page; theme = :publication, open_file = false,
            path = joinpath(summary_directory, "$(class.key)_$(sweep)_$(quantity).svg"))
        inkscape === nothing ||
            run(pipeline(`$inkscape $svg --export-filename=$(first(splitext(svg))).pdf`;
                stdout = devnull, stderr = devnull))
    end
end

println("\nSummary written to ", summary_directory)
nothing
