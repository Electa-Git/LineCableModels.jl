# Plot the independently computed earth matrices against the original thin-wire
# FEM replay. Run with --project=. after earth_matrices_manual.jl and the replay.
using CairoMakie
using JSON3
using Printf

root = normpath(joinpath(@__DIR__, "..", ".."))
directory = joinpath(root, ".linecablemodels", "qa", "thin-wire-trail")
analytic_path = joinpath(root, ".linecablemodels", "qa", "earth-matrices-manual", "matrices.json")
analytic = JSON3.read(read(analytic_path, String))
fem = JSON3.read(read(joinpath(directory, "replay.json"), String))
historical = JSON3.read(read(joinpath(directory, "historical-sheet4.json"), String))
decode(value) = reshape(complex.(Float64.(value.real), Float64.(value.imag)), Tuple(Int.(value.shape)))
fem_entry(value, row, col) = complex(Float64(value.real[row][col]), Float64(value.imag[row][col]))
frequency = Float64.(getproperty.(analytic.rows, :frequency))
@assert frequency == Float64.(getproperty.(fem, :frequency))
dense_frequency = Float64.(getproperty.(analytic.dense, :frequency))
models = (:field_average, :full_current, :xue)
names = ["Proposed: field average · Cf (your table)", "Proposed: full manuscript · L⁻¹", "Xue"]
colors = ["#176dad", "#b12b85", "#333333"]
styles = [:solid, :dash, :dot]
fem_color = "#dd7027"
comsol_color = "#188465"
ticks = (frequency, ["0.1", "1", "10", "100", "1k", "10k", "100k", "1M"])

function historical_complex(text)
    text = strip(String(text))
    if occursin("j", text)
        return parse(ComplexF64, replace(text, "+j" => "+", "-j" => "-") * "im")
    end
    return parse(ComplexF64, replace(text, "i" => "im"))
end

curve(quantity, model, row) = [decode(sample.models[model][quantity])[row,1] for sample in analytic.dense]
sampled(quantity, model, row) = [decode(sample.models[model][quantity])[row,1] for sample in analytic.rows]
fem_values(quantity, row) = [fem_entry(sample[quantity == :Ze ? :Z : :Y], row, 1) for sample in fem]
function comsol_values(quantity, row)
    start = quantity == :Ze ? (row == 1 ? 6 : 17) : (row == 1 ? 28 : 39)
    return [historical_complex(historical["G$(start+i)"]) for i in 0:7]
end

set_theme!(Theme(fontsize=16, backgroundcolor=:white,
    Axis=(titlealign=:left, titlesize=19, xgridcolor=(:black,0.1), ygridcolor=(:black,0.1))))

function save_figure(figure, filename)
    for ext in ("png", "pdf", "svg")
        path = joinpath(directory, "$filename.$ext")
        ext == "png" ? save(path, figure; px_per_unit=1.4) : save(path, figure)
        println(path)
    end
end

function comparison_figure(quantities, filename)
    entries = [(quantity, row) for quantity in quantities for row in 1:2]
    n = length(entries)
    figure = Figure(size=(1480, 310n+340), figure_padding=26)
    Label(figure[0,1:2], "Proposal and Xue against the original FEM", fontsize=29, font=:bold, halign=:left)
    Label(figure[1,1:2], "Two bare wires · r = 42.5 mm · depth = 1 m · separation = 1 m · earth ρ = 0.1 Ω·m",
        fontsize=17, halign=:left)
    handles = Any[]
    for (i,(quantity,row)) in enumerate(entries), (col,(part,label)) in enumerate(((real,"Real"),(imag,"Imaginary")))
        symbol = quantity == :Ze ? "Z" : "Y"
        unit = quantity == :Ze ? "Ω/m" : "S/m"
        entry = row == 1 ? "₁₁ · self" : "₂₁ · mutual"
        logarithmic = row == 1 && (quantity == :Ze || col == 2)
        ax = Axis(figure[i+1,col], title="$symbol$entry · $label", ylabel=unit,
            xlabel=i == n ? "Frequency (Hz)" : "", xticklabelsvisible=i == n,
            xscale=log10, yscale=logarithmic ? log10 : identity, xticks=ticks)
        if logarithmic
            values = vcat([part.(curve(quantity,model,row)) for model in models]...,
                part.(fem_values(quantity,row)),part.(comsol_values(quantity,row)))
            powers = floor(Int,log10(minimum(values))):ceil(Int,log10(maximum(values)))
            ax.yticks = (10.0 .^ powers, [@sprintf("%.0e",10.0^p) for p in powers])
        end
        for (k,model) in enumerate(models)
            h = lines!(ax, dense_frequency, part.(curve(quantity,model,row));
                color=colors[k], linestyle=styles[k], linewidth=2.5)
            i == 1 && col == 1 && push!(handles,h)
        end
        hc = scatter!(ax, frequency, part.(comsol_values(quantity,row)); color=comsol_color,
            marker=:xcross, markersize=13, strokewidth=0)
        hf = scatter!(ax, frequency, part.(fem_values(quantity,row)); color=fem_color,
            marker=:circle, markersize=7, strokecolor=:white, strokewidth=0.5)
        i == 1 && col == 1 && append!(handles,[hf,hc])
        row == 2 && hlines!(ax,[0]; color=(:black,0.2), linewidth=0.8)
        xlims!(ax, first(frequency), last(frequency))
    end
    Legend(figure[n+2,1:2], handles, vcat(names,["Fresh original GetDP", "Archived COMSOL"]);
        orientation=:horizontal, nbanks=3, framevisible=false)
    footer = "Analytical curves: 281 quadrature samples; Ye = jω Pe⁻¹. Markers: eight frequencies; COMSOL was not rerun.\n" *
        "Cf uses isolated primary-current normalization; L⁻¹ uses the full manuscript current map."
    if :Ze in quantities
        footer *= "\nOriginal FEM Z uses a conducting active wire and a nonconducting receiver; FEM Y uses equipotential electrodes."
    else
        footer *= "\nFEM Y comes directly from the original Electric branch's electrode reactions."
    end
    Label(figure[n+3,1:2], footer; fontsize=13, halign=:left, color="#4c5560")
    rowgap!(figure.layout, 12)
    save_figure(figure,filename)
end

comparison_figure((:Ze,:Ye), "proposal-xue-fem-ZY")
comparison_figure((:Ze,), "proposal-xue-fem-Z")
comparison_figure((:Ye,), "proposal-xue-fem-Y")

# Relative complex differences show effects hidden when absolute curves overlap.
figure = Figure(size=(1400,1050), figure_padding=26)
Label(figure[0,1:2], "Difference from the fresh original GetDP results", fontsize=27, font=:bold, halign=:left)
handles = Any[]
for (i,quantity) in enumerate((:Ze,:Ye)), row in 1:2
    symbol = quantity == :Ze ? "Z" : "Y"
    entry = row == 1 ? "₁₁ · self" : "₂₁ · mutual"
    ax = Axis(figure[i,row], title="$symbol$entry", ylabel="|analytical − FEM| / |FEM| (%)",
        xlabel=i == 2 ? "Frequency (Hz)" : "", xscale=log10, yscale=log10, xticks=ticks)
    reference = fem_values(quantity,row)
    for (k,model) in enumerate(models)
        difference = 100abs.(sampled(quantity,model,row).-reference)./abs.(reference)
        h = scatterlines!(ax,frequency,difference; color=colors[k], linestyle=styles[k],
            linewidth=2, markersize=8)
        i == 1 && row == 1 && push!(handles,h)
    end
    xlims!(ax,first(frequency),last(frequency))
end
Legend(figure[3,1:2], handles, names; orientation=:horizontal, nbanks=2, framevisible=false)
Label(figure[4,1:2], "Eight solved frequencies; connecting lines guide the eye. FEM discretization error has not been converged.\nThese differences compare the stated analytical normalizations against the original Electric/Magnetic branches.",
    fontsize=13, halign=:left)
save_figure(figure,"proposal-xue-fem-relative")

open(joinpath(directory,"proposal-xue-fem.csv"),"w") do io
    println(io,"frequency_Hz,quantity,model,row,column,real,imag,fem_real,fem_imag,relative_complex_difference_percent")
    for (i,f) in enumerate(frequency), quantity in (:Ze,:Ye), model in models, row in 1:2, col in 1:2
        value = decode(analytic.rows[i].models[model][quantity])[row,col]
        reference = fem_entry(fem[i][quantity == :Ze ? :Z : :Y],row,col)
        difference = 100abs(value-reference)/abs(reference)
        @printf(io,"%.17g,%s,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            f,String(quantity),String(model),row,col,real(value),imag(value),real(reference),imag(reference),difference)
    end
end
