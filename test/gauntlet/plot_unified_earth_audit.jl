# Replot retained audit data from the repository root with --project=.
# The exploratory audit runner is retired.
# Produces scientific PNG/PDF/SVG figures and the exact plotted comparison data.
using CairoMakie
using JSON3
using LinearAlgebra
using Printf

decode_matrix(value) = reshape(complex.(Float64.(value.real), Float64.(value.imag)),
    Tuple(Int.(value.shape)))

function plot_earth_audit(directory)
    analytical = JSON3.read(read(joinpath(directory,"analytical.json"),String))
    dense = JSON3.read(read(joinpath(directory,"analytical-dense.json"),String))
    fem = JSON3.read(read(joinpath(directory,"fem-pec-boundary.json"),String))
    @assert analytical.candidate.internal_impedance == "PEC: Zc=0"
    @assert fem.core_material == "pec"
    @assert fem.representation == "exact PEC boundary; normalized Gamma=0 quasi-TEM"
    frequency = Float64.(analytical.frequencies)
    @assert frequency == Float64.(fem.frequencies)
    curve_frequency = Float64.(dense.frequencies)
    Zcurve = cat((decode_matrix(row.Z) for row in dense.raw)...; dims=3)
    Ycurve = cat((decode_matrix(row.Y) for row in dense.raw)...; dims=3)
    Zf, Yf = decode_matrix(fem.Z), decode_matrix(fem.Y)
    exact = filter(row -> row.method == "quad", analytical.raw)
    Za = cat((decode_matrix(row.Zearth) for row in exact)...; dims=3)
    Ya = cat((decode_matrix(row.Y) for row in exact)...; dims=3)
    output = joinpath(directory,"plots")
    mkpath(output)
    blue, orange = "#1464a0", "#cf5b22"
    set_theme!(Theme(fontsize=17, backgroundcolor=:white,
        Axis=(titlealign=:left, titlesize=19, xgridcolor=(:black,0.10),
            ygridcolor=(:black,0.10), xminorgridvisible=true,
            xminorgridcolor=(:black,0.04), xminorticks=IntervalsBetween(9))))

    figure = Figure(size=(1450,1540), figure_padding=(25,30,22,24))
    Label(figure[0,1:2], "Earth return · two bare PEC wires", fontsize=30,
        font=:bold, halign=:left)
    Label(figure[1,1:2],
        "r = 42.5 mm   ·   depth = 1 m   ·   separation = 1 m   ·   earth ρ = 0.1 Ω·m",
        fontsize=18, halign=:left)
    handles = Any[]
    entries = (("Z₁₁ · self", Zcurve, Zf, "Ω/m", 1),
        ("Z₁₂ · mutual", Zcurve, Zf, "Ω/m", 2),
        ("Y₁₁ · self", Ycurve, Yf, "S/m", 1),
        ("Y₁₂ · mutual", Ycurve, Yf, "S/m", 2))
    xticks = (10.0 .^ (-1:6), ["0.1","1","10","100","1k","10k","100k","1M"])
    for (i,(name, curves, reference, unit, column)) in enumerate(entries)
        for (j,(component, component_name)) in enumerate(((real,"Real"),(imag,"Imaginary")))
            av = component.(curves[1,column,:])
            fv = component.(reference[1,column,:])
            values = vcat(av,fv)
            # Linear signed mutual axes show the actual zero crossings without
            # compressing them into artificial-looking vertical symlog jumps.
            logarithmic = i == 1 || (i == 3 && j == 2)
            scale = logarithmic ? log10 : identity
            axis = Axis(figure[i+1,j], title="$name — $component_name",
                xlabel=i==4 ? "Frequency (Hz)" : "", ylabel=unit,
                xscale=log10, yscale=scale, xticks=xticks,
                xticklabelsvisible=i==4)
            if logarithmic
                ticks = 10.0 .^ (ceil(Int,log10(minimum(values))):floor(Int,log10(maximum(values))))
                axis.yticks = (ticks, [@sprintf("%.0e",value) for value in ticks])
            end
            analytic_handle = lines!(axis,curve_frequency,av; color=blue,linewidth=2.6)
            fem_handle = scatter!(axis,frequency,fv; color=orange,marker=:circle,
                markersize=11, strokecolor=:white, strokewidth=1.0)
            minimum(values) <= 0 && hlines!(axis,[0.0]; color=(:black,0.30),linewidth=0.8)
            xlims!(axis,first(frequency),last(frequency))
            if i==1 && j==1
                append!(handles,(analytic_handle,fem_handle))
            end
        end
    end
    Legend(figure[6,1:2],handles,["Supplied framework · earth only (Zc = 0)",
        "FEM · explicit PEC boundary"], orientation=:horizontal,framevisible=false)
    Label(figure[7,1:2],
        "PEC conductor interiors excluded; no internal impedance. FEM uses the normalized Γ = 0 quasi-TEM equations.\n" *
        "Self Z and imaginary self Y use logarithmic y axes; mutual values retain their signs. FEM: eight requested frequencies.",
        fontsize=14, halign=:left, color="#4c5560")
    rowgap!(figure.layout,12)
    for extension in ("png","pdf","svg")
        path = joinpath(output,"earth-only-self-mutual-ZY.$extension")
        extension == "png" ? save(path,figure; px_per_unit=1.5) : save(path,figure)
        println(path)
    end

    open(joinpath(output,"earth-only-self-mutual-ZY.csv"),"w") do io
        println(io,"frequency_Hz,quantity,row,column,framework_real,framework_imag,FEM_real,FEM_imag,absolute_error,relative_error")
        for (quantity,a,b) in (("Z",Za,Zf),("Y",Ya,Yf)), k in eachindex(frequency), col in 1:2, row in 1:2
            av,bv = a[row,col,k],b[row,col,k]
            absolute = abs(av-bv)
            println(io,join((frequency[k],quantity,row,col,real(av),imag(av),real(bv),imag(bv),
                absolute,absolute/max(abs(bv),floatmin(Float64))),','))
        end
    end
    println("Frequency | Z matrix relative error | Y matrix relative error")
    for (i,f) in enumerate(frequency)
        @printf("%9g | %.8g | %.8g\n", f,
            norm(Za[:,:,i]-Zf[:,:,i])/norm(Zf[:,:,i]),
            norm(Ya[:,:,i]-Yf[:,:,i])/norm(Yf[:,:,i]))
    end
    return figure
end

directory = isempty(ARGS) ? joinpath(@__DIR__,"..","..",".linecablemodels","qa",
    "unified-earth-audit","pec") : abspath(first(ARGS))
plot_earth_audit(directory)
