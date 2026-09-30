# Plot retained qualification matrices through the public observation/plot API.
# Default: interactive GLMakie. For report export: FEM_PLOT_BACKEND=cairo.
# julia --project=. plot_spectrum.jl SUMMARY/matrices.csv NEW_OUTPUT
using LineCableModels
using CairoMakie
const backend = Symbol(get(ENV, "FEM_PLOT_BACKEND", "gl"))
if backend === :gl
    using GLMakie
end

function plot_spectrum(input, output)
    ispath(output) && error("Use a fresh plot output directory")
    mkpath(output)
    cp(@__FILE__, joinpath(output, "plot_spectrum.jl"))
    rows = map(readlines(input)[2:end]) do line
        x = split(line, ',')
        (layout=x[1], rho=parse(Float64,x[2]), physics=parse(Int,x[3]),
         frequency=parse(Float64,x[4]), quantity=Symbol(x[5]),
         i=parse(Int,x[6]), j=parse(Int,x[7]),
         value=complex(parse(Float64,x[8]),parse(Float64,x[9])),
         reference=complex(parse(Float64,x[10]),parse(Float64,x[11])))
    end
    rhos = [0.1, 1., 100., 1000.]
    colors = [:dodgerblue3, :darkorange2, :forestgreen, :purple3]
    plots = UIPlot[]
    for physics in (1,0), layout in ("air_air", "air_earth", "earth_earth")
        label = physics == 1 ? "quasi-fw" : "quasi-TEM"
        placement = layout == "air_air" ? "Two aerial wires" :
            layout == "earth_earth" ? "Two buried wires" : "Mixed: conductor 1 aerial, conductor 2 buried"
        observations, labels, styles = ObservedResult[], String[], NamedTuple[]
        for reference in (false,true), (c,rho) in enumerate(rhos)
            selected = filter(row -> row.layout == layout && row.rho == rho && row.physics == physics, rows)
            frequencies = sort!(unique(getproperty.(selected, :frequency)))
            tensors = Dict(:Z=>zeros(ComplexF64,2,2,length(frequencies)),
                           :Y=>zeros(ComplexF64,2,2,length(frequencies)))
            for row in selected
                row.quantity in (:Z,:Y) || continue
                k = searchsortedfirst(frequencies,row.frequency)
                tensors[row.quantity][row.i,row.j,k] = reference ? row.reference : row.value
            end
            raw = LineParameters(tensors[:Z],tensors[:Y],frequencies)
            push!(observations, ObservedResult(raw,(R,X,G,B);clip=false,length_unit=:base))
            push!(labels, "$(reference ? "Analytical" : "FEM"): ρ = $rho Ω·m")
            push!(styles, (color=colors[c],linestyle=reference ? :dash : :solid,
                           marker=reference ? :utriangle : :circle,markersize=5))
        end
        println("PLOT ",label," ",layout); flush(stdout)
        pages = LineCableModels.plot(observations; ydata=(R,G,B), layout=(2,2),
            backend, display_plot=backend===:gl, controls=backend===:gl, open_export=false,
            xscale=:log10, yscale=:log10, quantity_units=:base,
            series_labels=labels, series_attributes=styles,
            figure_title="$placement — $label; radius 42.5 mm",
            legend_title="Unclipped values; solid FEM, dashed analytical",
            legend_overflow=:show_all, figure=(size=(1150,800),fontsize=14))
        for (quantity,page) in zip((:R,:G,:B),pages)
            name = "$(label)-$(layout)-$(quantity)"
            export_svg(page; path=joinpath(output,name*".svg"), open_file=false)
            CairoMakie.save(joinpath(output,name*".png"),page.figure;update=false)
            push!(plots,page)
        end
    end
    println("Exported ",length(plots)," API figures to ",output); flush(stdout)
    return plots
end

qualification_plots = plot_spectrum(abspath(ARGS[1]),abspath(ARGS[2]));
