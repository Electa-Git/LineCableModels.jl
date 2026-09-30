# Plot saved public-compute matrices only; this file never runs FEM.
# julia --project=. plot_bare_wire_mvp.jl OUTPUT_ROOT
using LineCableModels, CairoMakie

function plot_bare_wire_mvp(root)
    for (directory,_,files) in walkdir(root)
        all(in(files), ("complete.txt","matrices.csv","analytical.csv")) || continue
        output = joinpath(directory,"plots")
        mkpath(output)
        # Completed matrix snapshots are immutable; resume only missing plots.
        if all(isfile(joinpath(output,"$q.$ext")) for q in (:R,:G,:B) for ext in ("svg","png"))
            println("REUSE PLOTS ",relpath(output,root)); flush(stdout)
            continue
        end
        observations = ObservedResult[]
        n = 0
        for name in ("matrices.csv","analytical.csv")
            rows = split.(readlines(joinpath(directory,name))[2:end],',')
            frequencies = sort!(unique(parse(Float64,row[2]) for row in rows))
            n = maximum(parse(Int,row[3]) for row in rows)
            tensors = Dict(q=>zeros(ComplexF64,n,n,length(frequencies)) for q in ("Z","Y"))
            for row in rows
                k = searchsortedfirst(frequencies,parse(Float64,row[2]))
                tensors[row[1]][parse(Int,row[3]),parse(Int,row[4]),k] =
                    complex(parse(Float64,row[5]),parse(Float64,row[6]))
            end
            parameters = LineParameters(tensors["Z"],tensors["Y"],frequencies)
            push!(observations,ObservedResult(parameters,(R,X,G,B);clip=false,length_unit=:base))
        end
        pages = LineCableModels.plot(observations;ydata=(R,G,B),layout=(n,n),
            backend=:cairo,display_plot=false,controls=false,open_export=false,
            xscale=:log10,yscale=:log10,quantity_units=:base,
            series_labels=["FEM","Analytical (Γ = 0)"],
            series_attributes=[(color=:dodgerblue3,linestyle=:solid,marker=:circle),
                (color=:darkorange2,linestyle=:dash,marker=:utriangle)],
            figure_title=relpath(directory,root),legend_title="Unclipped per-metre values",
            figure=(size=(1200,950),fontsize=13))
        for (quantity,page) in zip((:R,:G,:B),pages)
            svg, png = joinpath(output,"$quantity.svg"), joinpath(output,"$quantity.png")
            isfile(svg) || export_svg(page;path=svg,open_file=false)
            isfile(png) || CairoMakie.save(png,page.figure;update=false)
        end
        println("PLOTS ",relpath(output,root)); flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && plot_bare_wire_mvp(abspath(only(ARGS)))
