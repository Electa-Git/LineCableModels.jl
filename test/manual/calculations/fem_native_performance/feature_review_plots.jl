using CairoMakie, Statistics
if !("--save" in ARGS)
    using GLMakie
end

function plot_native_review(fs)
    fs=TOML.parsefile(joinpath(REVIEW_ROOT,"settings.toml"))["frequencies_hz"]
    review_say("PLOT BEGIN: signed components, absolute conductance differences, timings")
    output=joinpath(REVIEW_ROOT,"plots"); mkpath(output)
    names=["baseline"; getproperty.(collect(REVIEW_CHOICES),:name); "analytical"]
    labels=["Frozen FEM: split terms; triangles 12/12";
        getproperty.(collect(REVIEW_CHOICES),:label); "Analytical: unified, Γ = 0"]
    colors=[:black,:royalblue,:darkorange2,:crimson,:seagreen]
    styles=[(color=c,marker=m,markersize=10,linestyle=s,linewidth=1.6)
        for (c,m,s) in zip(colors,[:circle,:rect,:utriangle,:diamond,:cross],
            [:solid,:solid,:dash,:dot,:dashdot])]
    records=[[TOML.parsefile(joinpath(REVIEW_ROOT,@sprintf("f%02d",k),name*".toml"))
        for k in eachindex(fs)] for name in names]
    matrices=[[review_matrices(data) for data in record] for record in records]
    parameters=[LineParameters(cat(first.(pairs)...;dims=3),cat(last.(pairs)...;dims=3),fs)
        for pairs in matrices]
    observations=[ObservedResult(p,(R,X,G,B);clip=false,
        atol=(R=0.,X=0.,G=0.,B=0.),length_unit=:base) for p in parameters]
    save_only="--save" in ARGS
    function pages(obs,qs,title,series_labels,attrs)
        handles=LineCableModels.plot(obs;ydata=qs,layout=(2,2),
            backend=save_only ? :cairo : :gl,display_plot=!save_only,controls=!save_only,
            open_export=false,xscale=:log10,yscale=:log10,quantity_units=:base,
            series_labels,series_attributes=attrs,figure_title=title,
            legend_title="Explicit prescriptions — scientific acceptance remains yours",
            legend_attributes=(nbanks=2,),legend_overflow=:show_all,
            figure=(size=(1500,1100),fontsize=14))
        handles = handles isa UIPlot ? [handles] : handles
        # Avoid loss of a tiny opposite-sign endpoint in floating-point limits.
        # This is a public axis setting; neither data nor cutoffs are altered.
        for (q,p) in zip(qs,handles), (k,axis) in enumerate(p.axes)
            i,j=divrem(k-1,2).+(1,1)
            lo,hi=extrema(v for o in obs for v in observe(o,q)[i,j,:])
            if lo < 0 < hi
                # One decade also leaves room for markers on the signed-log
                # ordinate; a 10% margin nearly vanishes across many decades.
                bound=10max(abs(lo),abs(hi)); CairoMakie.ylims!(axis,-bound,bound)
            end
        end
        handles
    end
    title="Two bare wires: mixed, ρ = 0.1 Ω·m, r = 4.25 cm, Γ = 0\n" *
        "$(length(fs)) computed frequencies; mesh factors 3/8; PML resolution 72 / 0.12; no clipping"
    plots=pages(observations,(R,X,G,B),title,labels,styles)
    diffs=[LineParameters(zeros(ComplexF64,2,2,length(fs)),
        complex.(abs.(real.(observe(p,Y).-observe(first(parameters),Y)))),fs) for p in parameters[2:end]]
    delta_obs=[ObservedResult(p,(G,B);clip=false,atol=(G=0.,B=0.),length_unit=:base) for p in diffs]
    delta=only(pages(delta_obs,(G,),"Absolute conductance difference |G − G frozen FEM| [S/m]\n" *
        "Absolute differences, not absolute conductances; all $(length(fs)) frequencies",labels[2:end],styles[2:end]))
    # CSV is the unrounded numerical comparison, including exact zero differences.
    open(joinpath(REVIEW_ROOT,"components.csv"),"w") do io
        println(io,"frequency_hz,variant,quantity,row,column,baseline,candidate,analytical,signed_difference,absolute_difference,relative_difference,sign_changed")
        for (k,f) in enumerate(fs), (v,name) in enumerate(names[2:end-1]),
                q in (R,X,G,B), i in 1:2, j in 1:2
            a=observe(observations[1],q)[i,j,k]; b=observe(observations[v+1],q)[i,j,k]
            c=observe(observations[end],q)[i,j,k]
            println(io,join((f,name,q,i,j,a,b,c,b-a,abs(b-a),iszero(a) ? NaN : abs((b-a)/a),signbit(a)!=signbit(b)),','))
        end
    end
    # Native total includes preprocessing and output, separately from Julia's
    # whole-call wall time and its measured first-use compilation component.
    open(joinpath(REVIEW_ROOT,"timings.csv"),"w") do io
        println(io,"frequency_hz,variant,assembly_seconds,solve_seconds,native_seconds,managed_seconds,compile_seconds,directory")
        for v in 1:4, (k,f) in enumerate(fs)
            r=records[v][k]
            println(io,join((f,names[v],(r[key] for key in ("assembly_seconds","solve_seconds",
                "native_seconds","managed_seconds","compile_seconds"))...,r["directory"]),','))
        end
    end
    fig=CairoMakie.Figure(size=(1500,950),fontsize=16)
    time_axes=CairoMakie.Axis[]
    CairoMakie.Label(fig[0,:],"Fixed prescriptions: measured runtime (serial, one solver thread)")
    for (k,(key,label)) in enumerate((("assembly_seconds","Native assembly"),
            ("solve_seconds","Native solve"),("native_seconds","Native total: preprocess through output"),
            ("managed_seconds","Managed whole call; first-use compilation included")))
        row,col=divrem(k-1,2).+(1,1)
        ax=CairoMakie.Axis(fig[row,col];xlabel="Frequency [Hz]",ylabel="Elapsed [s]",title=label,xscale=log10)
        push!(time_axes,ax)
        for v in (key=="managed_seconds" ? (2:4) : (1:4))
            CairoMakie.scatterlines!(ax,fs,[r[key] for r in records[v]];color=colors[v],
                marker=styles[v].marker,label=labels[v],linewidth=2)
        end
    end
    CairoMakie.Legend(fig[3,:],first(time_axes);nbanks=2,tellheight=true)
    if save_only
        for (name,p) in zip(["R","X","G","B","absolute-delta-G"],[plots;delta])
            svg,png=joinpath(output,name*".svg"),joinpath(output,name*".png")
            "--replace" in ARGS && (rm(svg;force=true);rm(png;force=true))
            isfile(svg) || export_svg(p;path=svg,open_file=false)
            isfile(png) || CairoMakie.save(png,p.figure;update=false)
        end
        for ext in ("svg","png")
            CairoMakie.save(joinpath(output,"timings."*ext),fig)
        end
    else
        display(GLMakie.Screen(),fig)
    end
    open(joinpath(REVIEW_ROOT,"review.md"),"w") do io
        println(io,"# Two-wire feature comparison\n\n",title,"\n")
        println(io,"All choices execute through feature code. No numerical acceptance threshold is applied.\n")
        println(io,"| Prescription | Native total [s] | Assembly [s] | Solve [s] | Julia compilation [s] |\n|---|---:|---:|---:|---:|")
        for v in 1:4
            totals=[sum(r[key] for r in records[v]) for key in
                ("native_seconds","assembly_seconds","solve_seconds","compile_seconds")]
            println(io,"| ",labels[v]," | ",join((@sprintf("%.3f",t) for t in totals)," | ")," |")
        end
        println(io,"\nNative total is the solver process duration, including preprocessing and output. " *
            "Managed whole-call durations and compilation are separately retained in timings.csv. " *
            "The frozen baseline reuses the two saved endpoints and solves only the eight missing points on the managed triangle mesh. " *
            "These are per-point measurements, not a repeated-run speedup certification.\n")
        println(io,"| Prescription | Sign differences from frozen FEM (R/X/G/B) | Largest relative G difference | Frequency / entry | Frozen G → candidate G [S/m] | Absolute difference there [S/m] |\n|---|---:|---:|---|---|---:|")
        for v in 2:4
            flips=sum(sign(a)!=sign(b) for q in (R,X,G,B) for (a,b) in
                zip(observe(observations[1],q),observe(observations[v],q)))
            a,b=observe(observations[1],G),observe(observations[v],G)
            ratios=map((x,y)->iszero(x) ? 0. : abs((y-x)/x),a,b)
            point=argmax(ratios); i,j,k=Tuple(point)
            println(io,"| ",labels[v]," | ",flips," | ",@sprintf("%.6g%%",100ratios[point]),
                " | ",@sprintf("%.8g Hz / G[%d,%d]",fs[k],i,j)," | ",
                @sprintf("%.9g → %.9g",a[point],b[point])," | ",
                @sprintf("%.9g",abs(b[point]-a[point]))," |")
        end
        println(io,"\nSign differences count departures from the frozen FEM at matching samples; they do not certify analytical agreement. " *
            "Exact zero reference entries have no defined percentage; components.csv retains their raw and absolute differences. " *
            "Native logs, including the pre-existing embedded-entity warnings, remain under the directories in timings.csv.\n")
        for name in ("G","B","R","X","absolute-delta-G","timings")
            println(io,"[",name,"](plots/",name,".svg)\n")
        end
        println(io,"Retained limitations: previous screen/tube tiny-conductance discrepancies remain available in " *
            "../native-performance-20260930/review-plots/index.md. Those fixtures were not rerun or used to block these options.")
    end
    review_say("PLOT COMPLETE ",output,"; components.csv, timings.csv and review.md written")
    return (;components=plots,absolute_conductance_difference=delta,timings=fig)
end
