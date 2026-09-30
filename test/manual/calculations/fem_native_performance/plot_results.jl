# Read retained native results only. No meshing, FEM, or analytical recomputation.
# Interactive: julia -i --project=test .../plot_results.jl screen
# Export all:  julia --project=test .../plot_results.jl --save --all
println("Loading plotting packages; this script does not run FEM."); flush(stdout)
using LineCableModels
using Printf
using CairoMakie
const save_only = "--save" in ARGS
if !save_only
    using GLMakie
end

const evidence = joinpath(pkgdir(LineCableModels), ".linecablemodels/fem/native-performance-20260930")
const plot_directory = joinpath(evidence, "review-plots")
const styles = Dict(
    "baseline" => (label="Reference FEM", color=:black, marker=:circle, markersize=14, linestyle=:solid),
    "analytical" => (label="Analytical", color=:dodgerblue3, marker=:rect, markersize=11, linestyle=:dash),
    "combined" => (label="Combined assembly", color=:darkorange2, marker=:xcross, markersize=11, linestyle=:dash),
    "harmonic" => (label="Harmonic assembly only", color=:seagreen, marker=:utriangle, markersize=9, linestyle=:dot),
    "physical3" => (label="Physical quadrature: 3 points", color=:purple, marker=:dtriangle, markersize=8, linestyle=:dashdot),
    "quad9" => (label="PML quadrangles: 9 points", color=:crimson, marker=:diamond, markersize=7, linestyle=:dot),
    "quad16" => (label="PML quadrangles: 16 points", color=:magenta, marker=:cross, markersize=13, linestyle=:dash))
const variants = ("baseline", "analytical", "combined", "harmonic", "physical3", "quad9", "quad16")

function native_matrix(directory, quantity)
    path = joinpath(directory, "results/f0001-quasi-fw-b0000/matrices/$quantity.tsv")
    rows = split.(readlines(path)[3:end], '\t')
    n = maximum(parse(Int,r[1]) for r in rows)
    matrix = zeros(ComplexF64,n,n)
    for r in rows
        matrix[parse(Int,r[1]),parse(Int,r[2])] = complex(parse(Float64,r[5]),parse(Float64,r[6]))
    end
    matrix
end

function analytical_matrices(directory)
    rows = split.(readlines(joinpath(directory,"components-combined.csv"))[2:end], ',')
    matrices = Dict(q=>zeros(2,2) for q in ("R","X","G","B"))
    for r in rows
        matrices[r[1]][parse(Int,r[2]),parse(Int,r[3])] = parse(Float64,r[6])
    end
    complex.(matrices["R"],matrices["X"]), complex.(matrices["G"],matrices["B"])
end

function casespec(name)
    low = endswith(name,"-low")
    base = low ? name[1:end-4] : name
    finite = endswith(base,"-finite")
    base = finite ? base[1:end-7] : base
    if base in ("air","soil","mixed")
        gamma = finite ? .99 : 0.
        fs = low ? [.1] : finite ? [.1,1e6] : [.1,1e3,1e6]
        paths = [joinpath(evidence,"$base-f$f-gamma$gamma") for f in fs]
        title = "Two bare wires: $base — ρ = 0.1 Ω·m, r = 4.25 cm, " *
            (finite ? "Γ = 0.99γearth" : "Γ = 0")
    else
        fs = low ? [.1] : [.1,1e6]
        paths = [joinpath(evidence,"fixtures","$base-f$f") for f in fs]
        title = Dict("three"=>"Three conductors in different layers",
            "screen"=>"49-wire screen — 1: core, 2: screen, 3: foil",
            "tube"=>"Tubular sheath — 1: core, 2: sheath",
            "sector"=>"Sector cable — 1–3: phases, 4: neutral")[base] * " — Γ = 0"
    end
    sampling = low ? "0.1 Hz close-up" : length(fs)==2 ? "saved endpoints only" : "3 saved frequency points"
    return fs, paths, "$title\n$sampling; markers are computed samples; no clipping"
end

function saved_observations(name)
    fs, paths, title = casespec(name)
    observed = ObservedResult[]
    labels = String[]
    attributes = NamedTuple[]
    n = 0
    for variant in variants
        frequencies = Float64[]
        zs, ys = Matrix{ComplexF64}[], Matrix{ComplexF64}[]
        for (f,path) in zip(fs,paths)
            if variant=="analytical"
                isfile(joinpath(path,"components-combined.csv")) || continue
                z,y = analytical_matrices(path)
            else
                dir = joinpath(path,variant)
                isfile(joinpath(dir,"results/f0001-quasi-fw-b0000/completed.txt")) || continue
                z,y = native_matrix(dir,"Z"),native_matrix(dir,"Y")
            end
            push!(frequencies,f); push!(zs,z); push!(ys,y)
        end
        isempty(frequencies) && continue
        n = size(first(zs),1)
        parameters = LineParameters(cat(zs...;dims=3),cat(ys...;dims=3),frequencies)
        push!(observed,ObservedResult(parameters,(R,X,G,B);clip=false,
            atol=(R=0.,X=0.,G=0.,B=0.),length_unit=:base))
        s = styles[variant]
        samples = join((@sprintf("%g",f) for f in frequencies),", ")
        push!(labels,"$(s.label) [$samples Hz]")
        push!(attributes,(color=s.color,marker=s.marker,markersize=s.markersize,
            linestyle=s.linestyle,linewidth=1.6))
    end
    return observed, labels, attributes, n, title
end

function plot_saved(name)
    println("READ/PLOT ",name); flush(stdout)
    observed, labels, attributes, n, title = saved_observations(name)
    output = joinpath(plot_directory,name)
    quantities = endswith(name,"-low") ? (G,B) : (G,B,R,X)
    pages = LineCableModels.plot(observed;ydata=quantities,layout=(n,n),
        backend=save_only ? :cairo : :gl,display_plot=!save_only,controls=!save_only,
        open_export=false,xscale=:log10,yscale=:log10,quantity_units=:base,
        series_labels=labels,series_attributes=attributes,figure_title=title,
        legend_title="Retained results — scientific acceptance is the researcher's decision",
        legend_attributes=(nbanks=2,),legend_overflow=:show_all,
        figure=(size=(max(1450,400n),max(1050,325n)),fontsize=14))
    # Native limits store origin + width: tiny opposite-sign endpoints can
    # disappear from that sum across >16 decades. Symmetric limits retain
    # both signs. UIPlot.axes and Makie.ylims! are public, caller-owned APIs.
    for (q,p) in zip(quantities,pages)
        values = [observe(o,q) for o in observed]
        for (k,axis) in enumerate(p.axes)
            i,j = divrem(k-1,n) .+ (1,1)
            lo,hi = extrema(v for a in values for v in a[i,j,:])
            if lo < 0 < hi
                bound = 1.1max(abs(lo),abs(hi))
                CairoMakie.ylims!(axis,-bound,bound)
            end
        end
    end
    if save_only
        mkpath(output)
        for (q,p) in zip(quantities,pages)
            svg,png = joinpath(output,"$q.svg"),joinpath(output,"$q.png")
            if "--replace" in ARGS
                rm(svg;force=true); rm(png;force=true)
            end
            isfile(svg) || export_svg(p;path=svg,open_file=false)
            isfile(png) || CairoMakie.save(png,p.figure;update=false)
        end
        println("PLOTS ",name," -> ",output); flush(stdout)
    end
    return pages
end

const available_cases = ("air","soil","mixed","air-finite","soil-finite","mixed-finite",
    "three","screen","tube","sector","air-low","screen-low","tube-low")
const requested_cases = "--all" in ARGS ? collect(available_cases) :
    filter(x->!startswith(x,"--"),ARGS)
isempty(requested_cases) && push!(requested_cases,"air")
review_plots = [plot_saved(name) for name in requested_cases];
