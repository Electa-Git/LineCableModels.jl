# Run after earth_matrices_manual.jl, with --project=.
using CairoMakie
using JSON3
using Printf

decode(value) = reshape(complex.(Float64.(value.real),Float64.(value.imag)),Tuple(Int.(value.shape)))
directory = joinpath(@__DIR__,"..","..",".linecablemodels","qa","earth-matrices-manual")
data = JSON3.read(read(joinpath(directory,"matrices.json"),String))
models = (:field_average,:point_field,:xue,:full_current)
names = ["Field average · Cf", "Point field · Cf", "Xue", "Full manuscript · L⁻¹"]
colors = ["#1464a0","#d67e14","#292929","#af3685"]
styles = [:solid,:dash,:dot,:dashdot]
frequency = Float64.(getproperty.(data.rows,:frequency))
dense_frequency = Float64.(getproperty.(data.dense,:frequency))
curve(quantity,model,column) = [decode(row.models[model][quantity])[1,column] for row in data.dense]
set_theme!(Theme(fontsize=16,backgroundcolor=:white,
    Axis=(titlealign=:left,titlesize=18,xgridcolor=(:black,0.1),ygridcolor=(:black,0.1))))
ticks = (frequency,["0.1","1","10","100","1k","10k","100k","1M"])

function comparison_figure(quantities,filename;size=(1450,1520))
    figure = Figure(size=size,figure_padding=(24,30,20,24))
    Label(figure[0,1:2],"Earth matrices · two bare wires",fontsize=29,font=:bold,halign=:left)
    Label(figure[1,1:2],"r = 42.5 mm · depth = 1 m · separation = 1 m · ρ = 0.1 Ω·m · Γ = 0",
        fontsize=17,halign=:left)
    handles = Any[]
    entries = [(quantity,column) for quantity in quantities for column in 1:2]
    for (i,(quantity,column)) in enumerate(entries), (j,(component,label)) in enumerate(((real,"Real"),(imag,"Imaginary")))
        unit = quantity==:Ze ? "Ω/m" : quantity==:Pe ? "m/F" : "S/m"
        symbol = quantity==:Ze ? "Zₑ" : quantity==:Pe ? "Pₑ" : "Yₑ"
        entry = column==1 ? "₁₁ · self" : "₁₂ · mutual"
        values = [component.(curve(quantity,model,column)) for model in models]
        # Logarithmic self axes where positive; signed mutual traces stay linear.
        logarithmic = column==1 && (quantity==:Ze || (quantity==:Ye && j==2) || (quantity==:Pe && j==2))
        axis = Axis(figure[i+1,j],title="$symbol$entry — $label",ylabel=unit,
            xlabel=i==length(entries) ? "Frequency (Hz)" : "",
            xscale=log10,yscale=logarithmic ? log10 : identity,
            xticks=ticks,xticklabelsvisible=i==length(entries))
        allvalues = vcat(values...)
        if logarithmic
            yticks = 10.0 .^ (ceil(Int,log10(minimum(allvalues))):floor(Int,log10(maximum(allvalues))))
            axis.yticks = (yticks,[@sprintf("%.0e",v) for v in yticks])
        end
        for (k,model) in enumerate(models)
            handle = lines!(axis,dense_frequency,values[k];color=colors[k],linestyle=styles[k],linewidth=2.3)
            i==1 && j==1 && push!(handles,handle)
        end
        minimum(allvalues)<0<maximum(allvalues) && hlines!(axis,[0];color=(:black,0.25),linewidth=0.8)
        xlims!(axis,first(frequency),last(frequency))
    end
    Legend(figure[length(entries)+2,1:2],handles,names;orientation=:horizontal,framevisible=false,nbanks=2)
    Label(figure[length(entries)+3,1:2],
        "Cf = 1/[κg r K₁(κg r)] · supplied table uses isolated primary-current normalization.\n" *
        "Full manuscript uses the complete L matrix. Ye = jω Pe⁻¹ in every comparison; no internal or insulation terms.",
        fontsize=13,halign=:left,color="#4c5560")
    rowgap!(figure.layout,12)
    for ext in ("png","pdf","svg")
        path = joinpath(directory,"$filename.$ext")
        ext=="png" ? save(path,figure;px_per_unit=1.5) : save(path,figure)
        println(path)
    end
end

comparison_figure((:Ze,:Ye),"earth-Ze-Ye")
comparison_figure((:Pe,),"earth-Pe";size=(1450,900))

# Relative complex differences expose small effects hidden by overlaid curves.
figure = Figure(size=(1300,940),figure_padding=25)
Label(figure[0,1:2],"Earth matrices · complex relative difference from Xue",fontsize=26,font=:bold,halign=:left)
handles = Any[]
for (i,quantity) in enumerate((:Ze,:Ye)), column in 1:2
    symbol = quantity==:Ze ? "Zₑ" : "Yₑ"
    entry = column==1 ? "₁₁ · self" : "₁₂ · mutual"
    reference = curve(quantity,:xue,column)
    axis = Axis(figure[i,column],title="$symbol$entry",ylabel="|model − Xue| / |Xue| (%)",
        xlabel=i==2 ? "Frequency (Hz)" : "",xscale=log10,yscale=log10,xticks=ticks)
    for k in (1,2,4)
        difference = 100abs.(curve(quantity,models[k],column).-reference)./abs.(reference)
        handle = lines!(axis,dense_frequency,difference;color=colors[k],linestyle=styles[k],linewidth=2.5)
        i==1 && column==1 && push!(handles,handle)
    end
    xlims!(axis,first(frequency),last(frequency))
end
Legend(figure[3,1:2],handles,names[[1,2,4]];orientation=:horizontal,framevisible=false)
for ext in ("png","pdf")
    path=joinpath(directory,"earth-relative-differences.$ext")
    ext=="png" ? save(path,figure;px_per_unit=1.5) : save(path,figure)
    println(path)
end
