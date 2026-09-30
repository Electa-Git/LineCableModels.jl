# Plot the saved results of the historical contour-averaged voltage study.
# Its preparation-based solve modes have been retired. This reader never meshes
# or solves and preserves the historical reference/averaging labels in the CSV.
#   julia --project=test test/manual/calculations/run_fem_voltage_reference.jl --plot /path/to/study
using Printf
import CairoMakie as CM

function plot_study(directory)
    rows = [split(row, ',') for row in Iterators.drop(eachline(joinpath(directory,"matrices.csv")),1)]
    names = ("Yhistorical", "Ydeep", "Ysurface", "Yanalytical")
    labels = ("Historical point / deep earth", "Matched average / deep earth",
        "Corrected / local surface", "Analytical: unified")
    colors = (:gray55, :darkorange2, :dodgerblue3, :black)
    styles = (:dot, :dash, :solid, :dashdot)
    n = count(r -> r[1]=="Ysurface" && r[3]=="1" && r[4]=="1", rows)
    fig = CM.Figure(size=(1200,850))
    CM.Label(fig[0,1:2], "$(basename(directory)) — signed conductance, $n frequencies", fontsize=21)
    handles = Any[]
    for i in 1:2, j in 1:2
        magnitudes=filter(!iszero,[abs(parse(Float64,r[5])) for r in rows
            if r[1] in names && r[3]==string(i) && r[4]==string(j)])
        signed_log=!isempty(magnitudes) && maximum(magnitudes)/minimum(magnitudes)>1000
        low=isempty(magnitudes) ? -12 : floor(Int,log10(minimum(magnitudes)))
        high=isempty(magnitudes) ? -12 : ceil(Int,log10(maximum(magnitudes)))
        decades=10.0.^(low:2:high)
        ticks=vcat(-reverse(decades),0.0,decades)
        ax = CM.Axis(fig[i,j], title="G[$i,$j]", xlabel="Frequency [Hz]", ylabel="G [S/m]",
            xscale=log10, yscale=signed_log ? CM.Makie.Symlog10(10.0^low) : identity,
            yticks=signed_log ? (ticks,[iszero(t) ? "0" : @sprintf("%.0e",t) for t in ticks]) : CM.Makie.automatic)
        CM.hlines!(ax,[0.0];color=:gray80,linewidth=1,yautolimits=false)
        for k in eachindex(names)
            selected = sort(filter(r -> r[1]==names[k] && r[3]==string(i) && r[4]==string(j),rows);
                by=r -> parse(Float64,r[2]))
            f = [parse(Float64,r[2]) for r in selected]
            g = [parse(Float64,r[5]) for r in selected]
            line = CM.lines!(ax,f,g;color=colors[k],linestyle=styles[k],linewidth=2)
            # Mark the actual samples; four-point smoke lines are not a dense sweep.
            CM.scatter!(ax,f,g;color=colors[k],marker=(:rect,:utriangle,:diamond,:circle)[k],
                markersize=n<=4 ? 8 : 3)
            i==1 && j==1 && push!(handles,line)
        end
    end
    CM.Legend(fig[3,1:2],handles,collect(labels);orientation=:horizontal,nbanks=2)
    CM.Label(fig[4,1:2],"Unclipped real(Y). Wide-range panels use signed log axes, linear near zero. Lines join computed samples.",fontsize=14)
    for extension in ("png","pdf")
        path=joinpath(directory,"conductance.$extension")
        CM.save(path,fig)
        println("Conductance plot: ",path)
    end
    return fig
end

function main(args)
    length(args) == 2 && first(args) == "--plot" ||
        throw(ArgumentError("Usage: --plot STUDY_DIRECTORY (saved matrices.csv only)"))
    return plot_study(args[2])
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
