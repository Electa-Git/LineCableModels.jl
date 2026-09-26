# Included inside GoldenFixtures. These assertions inspect the actual native
# scene; the independent recipes describe routing and rendering, not physics.
native_primitives(plots)=collect(Iterators.flatten(
    plot isa Makie.PlotList ? native_primitives(plot.plots) : [plot] for plot in plots))
function check_scene(name,handles)
    expected_pages=name in ("line_rlcg","line_zy_cartesian","line_zy_polar") ? 4 : 1
    @assert length(handles)==expected_pages "missing or duplicated quantity page"
    for handle in handles
        @assert !isempty(handle.export_name) "missing export identity"
        for axis in handle.axes
            @assert !isempty(axis.title[]) "missing panel title"
            @assert all(isfinite,axis.finallimits[].origin)
            @assert all(isfinite,axis.finallimits[].widths)
        end
    end
    if startswith(name,"line_") || name=="uncertainty_intervals"
        parameters=two_conductor_results()
        arrays=name=="line_rlcg" ? (R(parameters),L(parameters),G(parameters),C(parameters)) :
            name=="line_zy_cartesian" ? (real.(Z(parameters)),imag.(Z(parameters)),real.(Y(parameters)),imag.(Y(parameters))) :
            # :base specifies an SI prefix; the current angle display unit is
            # degrees. Derive that conversion explicitly from radian phase.
            name=="line_zy_polar" ? (abs.(Z(parameters)),abs.(Y(parameters)),180/pi.*angle.(Z(parameters)),180/pi.*angle.(Y(parameters))) : (R(parameters),)
        for (handle,array) in zip(handles,arrays)
            @assert Set(keys(handle.addon_state.panel_data))==Set(((1,1),(1,2),(2,1),(2,2)))
            for ((i,j),panel) in handle.addon_state.panel_data
                curves=filter(plot->plot isa Makie.Lines,panel.axis.scene.plots)
                @assert length(curves)==(name=="uncertainty_intervals" ? 3 : 1) "wrong series count"
                for curve in curves
                    points=curve[1][]
                    @assert length(points)==3
                    expected=[eltype(points)(parameters.f[k],array[i,j,k]) for k in 1:3]
                    @assert isapprox(points,expected;rtol=8eps(Float32),atol=0) "wrong quantity, units, frequency or matrix entry"
                    @assert curve.visible[] "required series is hidden"
                end
                @assert occursin("Frequency",sprint(show,panel.axis.xlabel[]))
                @assert !isempty(sprint(show,panel.axis.ylabel[]))
                if name=="uncertainty_intervals"
                    bars=filter(plot->plot isa Makie.Errorbars,panel.axis.scene.plots)
                    @assert length(bars)==6
                    @assert count(p->p.direction[]===:x,bars)==3
                    @assert count(p->p.direction[]===:y,bars)==3
                    @assert all(p->p.visible[] && p.whiskerwidth[]>0 && p.linewidth[]>0,bars)
                    for direction in (:x,:y)
                        selected=filter(p->p.direction[]===direction,bars)
                        for (width,bar) in enumerate(selected)
                            expected_error=direction===:x ? width*.01parameters.f :
                                [width*.01(i+j+k)*array[i,j,k] for k in 1:3]
                            @assert getindex.(bar[1][],1) ≈ parameters.f rtol=8eps(Float32)
                            @assert getindex.(bar[1][],2) ≈ array[i,j,:] rtol=8eps(Float32)
                            @assert getindex.(bar[1][],3) ≈ expected_error rtol=8eps(Float32)
                            @assert getindex.(bar[1][],4) ≈ expected_error rtol=8eps(Float32)
                        end
                    end
                end
            end
        end
    elseif startswith(name,"mc_")
        result=actual_mc()
        samples=vec(only(result.sample_values).R)
        histogram=only(only(result.histogram_values).R)
        handle=only(handles);axis=only(handle.axes)
        @assert isapprox(sum(histogram.density.*diff(histogram.edges)),1;rtol=1e-10)
        if name=="mc_hist"
            native=only(filter(p->p isa Makie.Hist,axis.scene.plots))
            @assert native[1][]==samples "histogram uses another marginal"
        elseif name=="mc_pdf"
            native=only(filter(p->p isa Makie.Stairs,native_primitives(axis.scene.plots)))
            expected=[eltype(native[1][])(x,y) for (x,y) in zip(histogram.edges,[histogram.density;last(histogram.density)])]
            @assert isapprox(native[1][],expected;rtol=8eps(Float32),atol=0)
        elseif name=="mc_ecdf"
            native=only(filter(p->p isa Makie.Stairs,native_primitives(axis.scene.plots)))
            sorted=sort(samples);unique_samples=unique(sorted)
            x=[first(sorted)-eps(first(sorted));unique_samples]
            y=[count(v->v<=value,sorted)/length(sorted) for value in x]
            expected=[eltype(native[1][])(a,b) for (a,b) in zip(x,y)]
            @assert isapprox(native[1][],expected;rtol=8eps(Float32),atol=0) "ECDF uses another population"
            @assert first(y)==0 && last(y)==1
        else
            native=only(filter(p->p isa Makie.Scatter,axis.scene.plots))
            sorted=sort(samples);mass=cumsum(histogram.density.*diff(histogram.edges))
            model=map(eachindex(sorted)) do i
                probability=(i-.5)/length(sorted)
                bin=searchsortedfirst(mass,probability)
                prior=bin==1 ? 0. : mass[bin-1]
                histogram.edges[bin]+(probability-prior)/histogram.density[bin]
            end
            expected=[eltype(native[1][])(x,y) for (x,y) in zip(sorted,model)]
            @assert isapprox(native[1][],expected;rtol=8eps(Float32),atol=0) "wrong QQ probability or marginal"
        end
    elseif name=="uncertainty_intervals"
        error("unreachable uncertainty validation")
    elseif name in ("cable_preview","cable_preview_compact","system_preview")
        handle=only(handles);axis=only(handle.axes)
        @assert occursin("m",sprint(show,axis.xlabel[])) && occursin("m",sprint(show,axis.ylabel[]))
        @assert count(p->p isa Makie.Poly,axis.scene.plots)>=4 "missing material layer"
    elseif name=="material_scale"
        @assert length(only(handles).colorbars)==3 "missing material/property category"
    elseif name in ("formulation_comparison","uq_comparison")
        handle=only(handles)
        @assert length(handle.addon_state.order)>=2 "comparison lost an operand"
        @assert length(unique(values(handle.addon_state.labels)))>=2 "comparison lost candidate identity"
        @assert length(handle.axes)==4
    elseif name=="custom_layout"
        @assert length(only(handles).axes)==2
        @assert [axis.title[] for axis in only(handles).axes] == ["Quadratic","Discrete samples"]
    end
    return true
end

function defect!(name,handle)
    if name=="material_scale"
        bar=first(handle.colorbars)
        bar.colormap[]=Makie.to_colormap([:red,:red])
        return "wrong material color"
    elseif name=="custom_layout"
        xlims!(last(handle.axes),10.,20.)
        return "required samples clipped"
    end
    plots=native_primitives(collect(Iterators.flatten(axis.scene.plots for axis in handle.axes)))
    selected=name=="uncertainty_intervals" ? filter(p->p isa Makie.Errorbars,plots) :
        occursin("preview",name) ? filter(p->p isa Makie.Poly,plots) :
        filter(p->p isa Union{Makie.Lines,Makie.Hist,Makie.Stairs,Makie.ECDFPlot,Makie.Scatter},plots)
    isempty(selected) && error("no applicable visible primitive for $name")
    if name=="uncertainty_intervals"
        foreach(p->p.visible[]=false,selected)
        return "uncertainty bars hidden"
    end
    first(selected).visible[]=false
    return occursin("preview",name) ? "material layer removed" : "required series omitted"
end
