# Render saved qualification matrices through the public plotting API.
# This script never invokes FEM and resumes missing exports only.
using LineCableModels, CairoMakie

function plot_cable_fixtures(root)
    physics = "quasi_fw"
    for name in ("screen","tube","sector")
        paths = [joinpath(root,name,level,physics,"matrices.csv") for level in ("normal","refined")]
        all(isfile,paths) || continue
        output = joinpath(root,name,"plots",physics)
        mkpath(output)
        all(isfile(joinpath(output,"$q.$ext")) for q in (:R,:G,:B) for ext in ("svg","png")) && continue
        observations = ObservedResult[]
        n = 0
        for path in paths
            rows = split.(readlines(path)[2:end],',')
            frequencies = sort!(unique(parse(Float64,row[2]) for row in rows))
            n = maximum(parse(Int,row[3]) for row in rows)
            tensors = Dict(q=>zeros(ComplexF64,n,n,length(frequencies)) for q in ("Z","Y"))
            for row in rows
                k = searchsortedfirst(frequencies,parse(Float64,row[2]))
                tensors[row[1]][parse(Int,row[3]),parse(Int,row[4]),k] =
                    complex(parse(Float64,row[5]),parse(Float64,row[6]))
            end
            push!(observations,ObservedResult(LineParameters(tensors["Z"],tensors["Y"],frequencies),
                (R,X,G,B);clip=false,length_unit=:base))
        end
        pages = LineCableModels.plot(observations;ydata=(R,G,B),layout=(n,n),
            backend=:cairo,display_plot=false,controls=false,open_export=false,
            xscale=:log10,yscale=:log10,quantity_units=:base,
            series_labels=["Prescribed mesh","Refined 1 MHz control"],
            series_attributes=[(color=:dodgerblue3,linestyle=:solid,marker=:circle),
                (color=:darkorange2,linestyle=:dash,marker=:utriangle)],
            figure_title="$name / $physics",legend_title="Unclipped per-metre values",
            figure=(size=(max(1200,400n),max(950,300n)),fontsize=13))
        for (quantity,page) in zip((:R,:G,:B),pages)
            svg,png = joinpath(output,"$quantity.svg"),joinpath(output,"$quantity.png")
            isfile(svg) || export_svg(page;path=svg,open_file=false)
            isfile(png) || CairoMakie.save(png,page.figure;update=false)
        end
        println("PLOTS $name/$physics"); flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && plot_cable_fixtures(abspath(only(ARGS)))
