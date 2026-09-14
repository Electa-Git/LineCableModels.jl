# Run with --project=. after three_wire_earth_matrices.jl.
using CairoMakie
using JSON3
using Printf

directory = normpath(joinpath(@__DIR__,"..","..",".linecablemodels","qa","three-wire-earth-matrices"))
data = JSON3.read(read(joinpath(directory,"matrices.json"),String))
decode(value) = reshape(complex.(Float64.(value.real),Float64.(value.imag)),Tuple(Int.(value.shape)))
models = (:proposal,:xue,:field_average_cf)
names = ["Proposal · full L", "Xue", "Field average · Cf"]
colors = ["#1464a0","#d46620","#747981"]
styles = [:solid,:dash,:dot]
frequency = Float64.(data.frequencies)
dense_frequency = Float64.(getproperty.(data.dense,:frequency))
curve(quantity,model,i,j) = [decode(row.models[model][quantity])[i,j] for row in data.dense]
ticks = (frequency,["0.1","1","10","100","1k","10k","100k","1M"])
set_theme!(Theme(fontsize=16,backgroundcolor=:white,
    Axis=(titlealign=:left,titlesize=18,xgridcolor=(:black,0.1),ygridcolor=(:black,0.1))))

for quantity in (:Ze,:Ye)
    entries = quantity==:Ze ? ((1,1,"Outer self · Z₁₁ = Z₃₃"),
        (2,2,"Middle self · Z₂₂"),(1,2,"Adjacent · Z₁₂ = Z₂₁"),
        (1,3,"Outer pair · Z₁₃ = Z₃₁")) :
        ((1,1,"Outer self · Y₁₁ = Y₃₃"),(2,2,"Middle self · Y₂₂"),
        (1,2,"Adjacent · Y₁₂ = Y₃₂"),(2,1,"Adjacent · Y₂₁ = Y₂₃"),
        (1,3,"Outer pair · Y₁₃ = Y₃₁"))
    unit = quantity==:Ze ? "Ω/m" : "S/m"
    figure = Figure(size=(1450,280+300length(entries)),figure_padding=(24,30,20,24))
    Label(figure[0,1:2],"Earth $(quantity==:Ze ? "impedance" : "admittance") · three bare wires",
        fontsize=28,font=:bold,halign=:left)
    Label(figure[1,1:2],"Positions = 0, 1, 2 m · depth = 1 m · radius = 42.5 mm · ρ = 0.1 Ω·m · Γ = 0",
        fontsize=17,halign=:left)
    handles = Any[]
    for (row,(i,j,label)) in enumerate(entries), (column,(component,component_label)) in enumerate(((real,"Real"),(imag,"Imaginary")))
        values = [component.(curve(quantity,model,i,j)) for model in models]
        allvalues = vcat(values...)
        logarithmic = i==j && (quantity==:Ze || component===imag)
        axis = Axis(figure[row+1,column],title="$label — $component_label",ylabel=unit,
            xlabel=row==length(entries) ? "Frequency (Hz)" : "",
            xscale=log10,yscale=logarithmic ? log10 : identity,
            xticks=ticks,xticklabelsvisible=row==length(entries))
        if logarithmic
            yticks = 10.0 .^ (ceil(Int,log10(minimum(allvalues))):floor(Int,log10(maximum(allvalues))))
            axis.yticks = (yticks,[@sprintf("%.0e",v) for v in yticks])
        end
        for k in eachindex(models)
            handle=lines!(axis,dense_frequency,values[k];color=colors[k],linestyle=styles[k],linewidth=2.4)
            row==1 && column==1 && push!(handles,handle)
        end
        minimum(allvalues)<0<maximum(allvalues) && hlines!(axis,[0];color=(:black,0.25),linewidth=0.8)
        xlims!(axis,first(frequency),last(frequency))
    end
    Legend(figure[length(entries)+2,1:2],handles,names;orientation=:horizontal,framevisible=false)
    Label(figure[length(entries)+3,1:2],
        "Complete 3 × 3 solves · Ye = jω Pe⁻¹ · deep-earth voltage reference\n" *
        "No internal impedance, insulation, or FEM. Directional Y entries are retained as calculated.",
        fontsize=13,halign=:left,color="#4c5560")
    rowgap!(figure.layout,12)
    for extension in ("png","pdf","svg")
        path=joinpath(directory,"three-wire-$quantity.$extension")
        extension=="png" ? save(path,figure;px_per_unit=1.5) : save(path,figure)
        println(path)
    end
end
