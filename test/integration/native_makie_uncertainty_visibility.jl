@testitem "Makie addons / uncertainty legend actions preserve complete series" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0, 10.0, 100.0]
    sources = map((1.0, 10.0)) do offset
        z = [measurement(offset + i + 2j + k, 0.1offset + 0.01i + 0.001k) +
             im * measurement(2offset + i + j + k, 0.2offset)
             for i in 1:2, j in 1:2, k in eachindex(f)]
        LineParameters(z, z .* 1e-6, f)
    end
    page = LineCableModels.plot(sources...; ydata=(R,),
        backend=:cairo, display_plot=false, controls=true, open_export=false,
        series_labels=("first", "second"), legend_position=:bottom,
        series_attributes=((marker=:circle,), (marker=:utriangle,)),
        length_unit=:base, quantity_units=:base, clip=false)
    entries = last(only(page.legend.entrygroups[]))
    for entry in entries
        # Exercise the native legend action, not a direct assignment to the mean
        # line: the latter missed the exception in the error-bar child plots.
        Makie.toggle_visibility!(entry)
    end
    for axis in page.axes
        @test all(!plot.visible[] for plot in axis.scene.plots
            if plot isa Union{Makie.Lines,Makie.Errorbars,Makie.Scatter})
    end
    @test all(entry -> first(Makie.get_n_visible(entry)) == 0, entries)
    for entry in entries
        Makie.toggle_visibility!(entry, true)
    end
    for axis in page.axes
        @test all(plot.visible[] for plot in axis.scene.plots
            if plot isa Union{Makie.Lines,Makie.Errorbars,Makie.Scatter})
    end
    @test all(entry -> ==(Makie.get_n_visible(entry)...), entries)
end

@testitem "Makie addons / uncertainty coordinates and visibility survive every axis mode" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0, 10.0, 100.0]
    options = (backend=:cairo, display_plot=false, controls=true, open_export=false,
        length_unit=:base, quantity_units=:base, freq_unit=:base, clip=false)
    for (xspread,yspread,markers) in ((0,0,false),(0,1,false),(1,0,true),(1,1,true))
        frequency = measurement.(f, xspread .* 0.01f)
        sources = map((1.0,10.0)) do offset
            z = [measurement(offset+i+2j+k, yspread*(0.1offset+0.01i+0.001k)) +
                 im*measurement(2offset+i+j+k, yspread*0.2offset)
                 for i in 1:2, j in 1:2, k in eachindex(f)]
            LineParameters(z, z.*1e-6, frequency)
        end
        originals = map(p -> (copy(Z(p)),copy(Y(p)),copy(frequencies(p))),sources)
        attributes = markers ? ((marker=:circle,), (marker=:utriangle,)) : nothing
        page = LineCableModels.plot(sources...; options..., ydata=(R,),
            series_labels=("first","second"), series_attributes=attributes,
            legend_position=:bottom, panel_legends=(1,1)=>:right)
        for (index,axis) in enumerate(page.axes)
            i,j = cld(index,2),mod1(index,2)
            curves = filter(p -> p isa Makie.Lines,axis.scene.plots)
            bars = filter(p -> p isa Makie.Errorbars,axis.scene.plots)
            @test length(bars) == 2*(xspread+yspread)
            for (source,curve) in zip(sources,curves)
                expected = observe(source,R)[i,j,:]
                @test first.(curve[1][]) ≈ f
                @test last.(curve[1][]) ≈ nominal.(expected)
                matching = filter(bar -> Makie.to_color(bar.color[]) == Makie.to_color(curve.color[]),bars)
                @test length(matching) == xspread+yspread
                for bar in matching
                    @test getindex.(bar[1][],1) ≈ f
                    @test getindex.(bar[1][],2) ≈ nominal.(expected)
                    spread = bar.direction[] === :x ? uncertainty.(frequency) : uncertainty.(expected)
                    @test getindex.(bar[1][],3) ≈ spread
                    @test getindex.(bar[1][],4) ≈ spread
                end
            end
        end
        # Panel actions are local; figure actions recover a genuinely mixed
        # panel state. No assumption about private dependency storage is needed.
        panel_entry = first(last(only(page.panel_legends[(1,1)].entrygroups[])))
        Makie.toggle_visibility!(panel_entry)
        first_color = Makie.to_color(first(filter(p -> p isa Makie.Lines,first(page.axes).scene.plots)).color[])
        for (index,axis) in enumerate(page.axes), plot in axis.scene.plots
            plot isa Union{Makie.Lines,Makie.Errorbars,Makie.Scatter} || continue
            @test plot.visible[] == !(index == 1 && Makie.to_color(plot.color[]) == first_color)
        end
        entry = first(last(only(page.legend.entrygroups[])))
        visible,total = Makie.get_n_visible(entry)
        @test 0 < visible < total
        Makie.toggle_visibility!(entry,visible!=total)
        for cycle in 1:3
            Makie.toggle_visibility!(entry)
            for axis in page.axes, plot in axis.scene.plots
                plot isa Union{Makie.Lines,Makie.Errorbars,Makie.Scatter} || continue
                @test plot.visible[] == (Makie.to_color(plot.color[]) != first_color)
            end
            Makie.toggle_visibility!(entry,true)
        end
        for key in (:xlog,:ylog)
            page.controls[key].active[] = true
            @test !isempty(Makie.colorbuffer(page.figure))
            page.controls[key].active[] = false
        end
        for (source,original) in zip(sources,originals)
            @test isequal(Z(source),original[1])
            @test isequal(Y(source),original[2])
            @test isequal(frequencies(source),original[3])
        end
    end
end

@testitem "Makie addons / uncertainty legends survive recreation and publication" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0,10.0,100.0]
    z = [measurement(i+2j+k,0.1i+0.01j+0.001k) + im*measurement(i+j+k,0.2)
         for i in 1:3, j in 1:3, k in eachindex(f)]
    parameters = LineParameters(z,z.*1e-6,measurement.(f,0.01f))
    options = (backend=:cairo,display_plot=false,controls=true,open_export=false,
        length_unit=:base,quantity_units=:base,clip=false,fig_size=(1000,700))
    pages = LineCableModels.plot(parameters; options...,ydata=(R,),blocks=(2,2),
        series_labels=("uncertain",),series_attributes=(marker=:circle,),legend_position=:bottom)
    @test length.(getproperty.(pages,:axes)) == [4,2,2,1]
    pages = copy(pages) # The block-layout owner retains the original page set.
    push!(pages,Makie.plot(parameters.Z,frequencies(parameters),(R,1,1,:);
        options...,series_labels=("uncertain",),legend_position=:bottom))
    # This exercises the publication renderer, not the matrix renderer used by
    # the benchmark inspector. Both must honour the same uncertainty contract.
    artifact = report(TableReportDefinition(((R,:,:,:),); illustration=true,
        clip=false,
        plot_options=(backend=:cairo,display_plot=false,controls=true,open_export=false,
            fig_size=(1000,700),display_legend=true,legend_position=:bottom)),parameters)
    push!(pages,artifact.illustration)
    single = Makie.plot(parameters.Y,frequencies(parameters),(G,1,1,:);options...)
    @test single.legend === nothing
    single_line = only(filter(p->p isa Makie.Lines,only(single.axes).scene.plots))
    single_line.visible[] = false
    @test all(!p.visible[] for p in only(single.axes).scene.plots if p isa Makie.Errorbars)
    single_line.visible[] = true
    @test all(p.visible[] for p in only(single.axes).scene.plots if p isa Makie.Errorbars)
    for page in pages
        for position in (:right,:inside,:bottom)
            legend = figurelegend!(page;position,overflow=:show_all)
            entries = last(only(legend.entrygroups[]))
            foreach(Makie.toggle_visibility!,entries)
            @test all(!plot.visible[] for axis in page.axes for plot in axis.scene.plots
                if plot isa Union{Makie.Lines,Makie.Errorbars,Makie.Scatter})
            # Recreating a hidden legend must not resurrect its series or lose
            # its visibility state. Glyph count is deliberately not a contract.
            legend = figurelegend!(page;position,overflow=:show_all)
            @test all(entry -> first(Makie.get_n_visible(entry)) == 0,
                last(only(legend.entrygroups[])))
            foreach(entry -> Makie.toggle_visibility!(entry,true),last(only(legend.entrygroups[])))
            @test all(plot.visible[] for axis in page.axes for plot in axis.scene.plots
                if plot isa Union{Makie.Lines,Makie.Errorbars,Makie.Scatter})
        end
        @test !isempty(Makie.colorbuffer(page.figure))
    end
    page = last(pages)
    axis = only(page.axes)
    line = first(filter(p -> p isa Makie.Lines,axis.scene.plots))
    bar = first(filter(p -> p isa Makie.Errorbars,axis.scene.plots))
    bar.visible[] = false
    @test line.visible[] # Independent native edits are still allowed.
    figurelegend!(page;position=:bottom,overflow=:show_all)
    @test !bar.visible[] # A legend rebuild is not a series hide/show action.
    line.visible[] = false
    @test !bar.visible[]
    line.visible[] = true
    @test bar.visible[] # A subsequent series action restores its components.
    mktempdir() do directory
        path = export_svg(page;path=joinpath(directory,"uncertainty.svg"),open_file=false)
        @test filesize(path)>0
    end
end

@testitem "Makie addons / uncertainty legend shading follows owners after relayout" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0,10.0,100.0]
    z = reshape(complex.(measurement.([1.,2.,3.],0.4),1.0),1,1,:)
    source = LineParameters(z,z.*1e-6,f)
    page = LineCableModels.plot(source; ydata=(R,),backend=:cairo,display_plot=false,
        controls=false,open_export=false,series_labels=("uncertain",),
        series_attributes=(marker=:circle,),legend_position=:bottom,fig_size=(900,600))
    function legend_pixels(page)
        pixels = Makie.colorbuffer(page.figure)
        bounds = page.legend.layoutobservables.computedbbox[]
        viewport = page.figure.scene.viewport[]
        sx,sy = size(pixels,2)/viewport.widths[1],size(pixels,1)/viewport.widths[2]
        x1 = clamp(floor(Int,bounds.origin[1]*sx)+1,1,size(pixels,2))
        x2 = clamp(ceil(Int,(bounds.origin[1]+bounds.widths[1])*sx),x1,size(pixels,2))
        y1 = clamp(floor(Int,(viewport.widths[2]-bounds.origin[2]-bounds.widths[2])*sy)+1,1,size(pixels,1))
        y2 = clamp(ceil(Int,(viewport.widths[2]-bounds.origin[2])*sy),y1,size(pixels,1))
        return pixels[y1:y2,x1:x2]
    end
    shown = legend_pixels(page)
    Makie.toggle_visibility!(first(last(only(page.legend.entrygroups[]))))
    hidden = legend_pixels(page)
    @test shown != hidden
    figurelegend!(page;position=:bottom,overflow=:show_all)
    # Pixel comparison checks the actual shading listeners, not just the click
    # target list. Replacing a hidden legend must preserve its rendered state.
    @test legend_pixels(page) == hidden
    Makie.toggle_visibility!(first(last(only(page.legend.entrygroups[]))),true)
    @test legend_pixels(page) == shown
    figurelegend!(page;position=:bottom,overflow=:show_all,joinstyle=:round)
    @test !isempty(legend_pixels(page))
    axis = only(page.axes)
    view = Makie.Rect2d(2.0,1.0,30.0,1.0)
    axis.targetlimits[] = view
    figurelegend!(page;position=:bottom,overflow=:show_all)
    # Reinitialising native legend shading is not a visibility change and
    # must not reset a caller's current zoom/pan.
    @test axis.targetlimits[] == view
end

@testitem "Makie addons / overflowing uncertainty entries retain their owners" tags=[:visual] begin
    using CairoMakie, Measurements
    f = [1.0,10.0,100.0]
    sources = Tuple(LineParameters(fill(measurement(Float64(i),0.1)+im,1,1,3),
        fill(measurement(1e-6,1e-8)+1e-6im,1,1,3),f) for i in 1:24)
    legend_bounds = Observable(Rect2f(0,0,300,100))
    page = LineCableModels.plot(sources...;ydata=(R,),backend=:cairo,display_plot=false,
        controls=false,open_export=false,fig_size=(620,340),legend_position=:inside,
        legend_attributes=(bbox=legend_bounds,),
        series_labels=Tuple("Uncertain result $i" for i in 1:24))
    Makie.colorbuffer(page.figure)
    entries = last(only(page.legend.entrygroups[]))
    @test last(entries).label[] == "(...)"
    resize!(page.figure,1200,950)
    legend_bounds[] = Rect2f(0,0,300,900)
    Makie.colorbuffer(page.figure)
    entries = last(only(page.legend.entrygroups[]))
    @test last(entries).label[] == "Uncertain result 24"
    # This entry was absent from the compact legend. Normalising only visible
    # entries would allow its old derived-child click targets to reappear.
    Makie.toggle_visibility!(last(entries))
    axis = only(page.axes)
    @test !last(filter(p->p isa Makie.Lines,axis.scene.plots)).visible[]
    @test !last(filter(p->p isa Makie.Errorbars,axis.scene.plots)).visible[]
    @test first(filter(p->p isa Makie.Lines,axis.scene.plots)).visible[]
    Makie.toggle_visibility!(last(entries),true)
    @test all(p.visible[] for p in axis.scene.plots if p isa Union{Makie.Lines,Makie.Errorbars})
end
