@testitem "Makie addons / preview material appearance and switches" tags=[:visual] begin
    using CairoMakie
    using CairoMakie: RGB, red, green, blue
    using Measurements

    extension = Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    copper = Material(kind = :conductor, rho = 1.72e-8)
    dielectric = Material(kind = :insulator, rho = Inf, eps_r = 2.3)
    semicon = Material(kind = :semicon, rho = 0.1, eps_r = 20.0)
    design = build(CableDesign,
        "pattern-controls",
        terminal(:core,
            core(copper; r = 0.01), insulation(dielectric; t = 0.004),
            screen(semicon; t = 0.0002)))
    options = (; backend = :cairo, display_plot = false, controls = false,
        display_colorbars = false, open_export = false)
    plain = preview(design; options..., display_dielectric_pattern = false)
    patterned = preview(design; options...)
    plain_plots = only(plain.axes).scene.plots
    pattern_plots = only(patterned.axes).scene.plots
    @test length(plain_plots) == length(pattern_plots)
    @test count(p -> p.color[] isa Makie.AbstractPattern, pattern_plots) == 1
    @test all(p -> p.color[] isa CairoMakie.Colorant, plain_plots)
    for (before, after) in zip(plain_plots, pattern_plots)
        @test before[1][] == after[1][]
        if !(after.color[] isa Makie.AbstractPattern)
            @test before.color[] == after.color[]
        end
    end
    @test length(patterned.addon_state.groups) == length(plain.addon_state.groups)
    @test !isempty(Makie.colorbuffer(patterned.figure))
    figurelegend!(patterned; position = :top)
    @test length(only(patterned.axes).scene.plots) == length(pattern_plots)

    for enabled in (true, false)
        collection = preview([design, design]; options...,
            display_dielectric_pattern = enabled)
        @test all(collection.axes) do axis
            count(p -> p.color[] isa Makie.AbstractPattern, axis.scene.plots) ==
            Int(enabled)
        end
        system = build(LineCableSystem, design, (0.0, -0.1); connections = Dict(:core=>1))
        plotted = preview(system; options..., display_dielectric_pattern = enabled)
        @test count(p -> p isa Makie.Poly && p.color[] isa Makie.AbstractPattern,
            only(plotted.axes).scene.plots) == Int(enabled)
    end

    # Physical transfer endpoints and hue progression, independent of palette helpers.
    @test materialcolors(:rho).colormap == materialcolors(:rho, (1.72e-8, 2.5e-7)).colormap
    copper_color = extension._material_color(copper)
    @test [red(copper_color), green(copper_color), blue(copper_color)] ≈
          [221/255, 228/255, 236/255] atol=1e-6
    magnetic = [extension._material_color(Material(kind = :conductor, rho = 1.72e-8,
                    mu_r = mu)) for mu in (1.0, 2.0, 5.0, 10.0, 30.0, 100.0, 300.0)]
    @test first(magnetic) == extension._material_color(copper)
    @test all(c -> blue(c) > red(c) && blue(c) > green(c), magnetic[2:5])
    @test red(last(magnetic)) > red(magnetic[5])
    @test length(unique(magnetic)) == length(magnetic)
    @test extension._material_color(Material(kind = :conductor, rho = 1.72e-8,
        mu_r = measurement(10.0, 0.1))) == magnetic[4]
    anchors = ((0.1, RGB(94/255, 114/255, 131/255)),
        (100.0, RGB(139/255, 120/255, 96/255)),
        (1e4, RGB(184/255, 154/255, 99/255)))
    for (rho, expected) in anchors
        actual = extension._material_color(layer(; rho))
        @test maximum(abs.((red(actual)-red(expected), green(actual)-green(expected),
            blue(actual)-blue(expected)))) < 1e-6
    end
end

@testitem "Makie addons / dielectric pattern clips holes and concavities in rendered output" tags=[:visual] begin
    using CairoMakie
    extension = Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    material = Material(kind = :insulator, rho = Inf, eps_r = 2.3)
    outer = Point2f[(-2, -2), (2, -2), (2, 2), (0.5, 2), (
        0.5, 0.7), (-0.5, 0.7), (-0.5, 2), (-2, 2)]
    hole = Point2f[(-0.6, -1), (-0.6, -0.3), (0.6, -0.3), (
        0.6, -1)]
    figure = Figure(size = (400, 400))
    axis = Axis(figure[1, 1]; aspect = DataAspect())
    hidedecorations!(axis)
    geometry = Makie.GeometryBasics.Polygon(outer, [hole])
    poly!(axis, geometry;
        color = extension._material_color(material), strokewidth = 0)
    pattern = poly!(axis, geometry; visible = false,
        color = extension._material_color(material; pattern = true), strokewidth = 0)
    for limits in ((-2.2, 2.2, -2.2, 2.2), (-0.8, 0.8, -1.2, 0.9))
        limits!(axis, limits...)
        pattern.visible = false
        plain = copy(Makie.colorbuffer(axis; include_decorations = false))
        pattern.visible = true
        patterned = Makie.colorbuffer(axis; include_decorations = false)
        @test size(plain) == size(patterned)
        view = axis.finallimits[]
        clear_pixels = 0
        expected_clear = 0
        for row in axes(plain, 1), column in axes(plain, 2)

            x = view.origin[1] + (column-0.5)/size(plain, 2)*view.widths[1]
            y = view.origin[2] + (1-(row-0.5)/size(plain, 1))*view.widths[2]
            in_hole = -0.55 < x < 0.55 && -0.95 < y < -0.35
            in_notch = -0.45 < x < 0.45 && 0.75 < y < 1.9
            if in_hole || in_notch
                expected_clear += 1
                clear_pixels += plain[row, column] == patterned[row, column]
            end
        end
        # Compare complete hole/notch interiors, independent of the tile layout.
        @test clear_pixels == expected_clear > 100
        # Native image tiles and solid fills can round one color level apart.
        # Count visible ink rather than that 8-bit raster conversion difference.
        changed = count(eachindex(plain)) do index
            a, b = plain[index], patterned[index]
            max(abs(Makie.red(a)-Makie.red(b)), abs(Makie.green(a)-Makie.green(b)),
                abs(Makie.blue(a)-Makie.blue(b))) > 2/255
        end
        @test 0 < changed < 0.08length(plain)
    end
end

@testitem "Makie addons / earth spans follow the viewport and retain physical depths" tags=[:visual] begin
    using CairoMakie
    copper = Material(kind = :conductor, rho = 1.72e-8)
    design = build(CableDesign, "earth-viewport", terminal(:core, core(copper; r = 0.01)))
    system = build(LineCableSystem, design, (0.0, -0.2); connections = Dict(:core=>1))
    earth = build(EarthModel,
        (layer(rho = 1.0, thickness = 0.1),
            layer(rho = 100.0, thickness = 0.25), layer(rho = 1e4)))
    options = (;
        backend = :cairo, display_plot = false, controls = false, open_export = false)
    plotted = preview(system; earth_model = earth, options...)
    axis = only(plotted.axes)
    layers = [only(plotted.addon_state.groups[Symbol("earth_$i")]) for i in 1:3]
    scene_plots = copy(axis.scene.plots)
    washes = filter(p -> p isa Makie.HSpan && p ∉ layers, scene_plots)
    @test length(washes) == 1
    listener_count = length(axis.finallimits.listeners)
    for (xrange, yrange) in (((-20.0, 20.0), (-10.0, 10.0)),
        ((5.0, 6.0), (-1.0, 0.0)), ((-3.0, 3.0), (-100.0, -94.0)),
        ((-1.0, 1.0), (2.0, 4.0)), ((-0.04, 0.04), (-0.24, -0.16)))
        limits!(axis, xrange..., yrange...)
        Makie.colorbuffer(plotted.figure)
        limits = axis.finallimits[]
        bottom, top = limits.origin[2], limits.origin[2]+limits.widths[2]
        for (i, (physical_bottom, physical_top)) in
            enumerate(((-0.1, 0.0), (-0.35, -0.1), (-Inf, -0.35)))
            @test layers[i][1][] ≈ clamp(physical_bottom, bottom, top)
            @test layers[i][2][] ≈ clamp(physical_top, bottom, top)
            @test !layers[i].xautolimits[] && !layers[i].yautolimits[]
            bounds = Makie.boundingbox(only(layers[i].plots))
            @test bounds.origin[1] ≈ limits.origin[1]
            @test bounds.widths[1] ≈ limits.widths[1]
        end
        for wash in washes
            lows, highs = wash[1][], wash[2][]
            @test all(lows .<= highs)
            @test all(lows[lows .< highs] .>= 0)
            @test all(highs[lows .< highs] .<= top)
            @test wash.visible[] == (top > 0)
            if top <= 0
                @test lows == highs
            else
                @test first(lows) ≈ max(bottom, 0)
                @test last(highs) ≈ top
            end
        end
        @test axis.scene.plots == scene_plots
        @test length(axis.finallimits.listeners) == listener_count
    end
    autolimits!(axis)
    Makie.colorbuffer(plotted.figure)
    @test axis.finallimits[].widths[1] < 0.1
    @test axis.finallimits[].widths[2] < 1
    @test earth.layers[2].thickness == 0.1
    @test earth.layers[3].thickness == 0.25
    @test isinf(earth.layers[4].thickness)

    finite = preview(system; earth_model = EarthModel(100.0; thickness = 0.5), options...)
    finite_axis = only(finite.axes)
    limits!(finite_axis, -1, 1, -1, 1)
    Makie.colorbuffer(finite.figure)
    finite_layer = only(finite.addon_state.groups[:earth_1])
    @test finite_layer[1][] == -0.5
    @test finite_layer[2][] == 0.0
    flat = preview(system; earth_model = earth, display_surface_gradient = false, options...)
    @test length(only(flat.axes).scene.plots) == length(axis.scene.plots) - 1

    # The sky tint grows from transparency near z=0 to blue at the top.
    flat_axis = only(flat.axes)
    for candidate in (axis, flat_axis)
        limits!(candidate, -1, 1, -1, 1)
        hidedecorations!(candidate)
    end
    with_sky = Makie.colorbuffer(axis; include_decorations = false)
    without_sky = Makie.colorbuffer(flat_axis; include_decorations = false)
    @test size(with_sky) == size(without_sky)
    view = axis.finallimits[]
    column = round(Int, 0.25size(with_sky, 2))
    rows = [clamp(round(Int, (1 - (z - view.origin[2]) / view.widths[2]) *
                            size(with_sky, 1)), 1, size(with_sky, 1))
            for z in (0.01, 0.50, 0.99, -0.05)]
    near, middle, top, soil = with_sky[rows, column]
    @test near == without_sky[rows[1], column]
    @test Makie.red(near) > Makie.red(middle) > Makie.red(top)
    @test Makie.blue(middle) > Makie.red(middle) + 0.03
    @test Makie.blue(top) > Makie.red(top) + 0.10
    @test soil == without_sky[rows[4], column]
    for candidate in (axis, flat_axis)
        limits!(candidate, -1, 1, -0.5, -0.02)
    end
    @test Makie.colorbuffer(axis; include_decorations = false) ==
          Makie.colorbuffer(flat_axis; include_decorations = false)

    controlled = preview(system; earth_model = earth, options..., controls = true)
    controlled_axis = only(controlled.axes)
    Makie.colorbuffer(controlled.figure)
    initial = controlled_axis.finallimits[]
    limits!(controlled_axis, -10, 10, -10, 10)
    controlled.controls[:reset].clicks[] += 1
    Makie.colorbuffer(controlled.figure)
    @test controlled_axis.finallimits[] == initial
    resize!(controlled.figure.scene, 1000, 500)
    Makie.colorbuffer(controlled.figure)
    resized = controlled_axis.finallimits[]
    for index in 1:3
        span = only(controlled.addon_state.groups[Symbol("earth_$index")])
        bounds = Makie.boundingbox(only(span.plots))
        @test bounds.origin[1] ≈ resized.origin[1]
        @test bounds.widths[1] ≈ resized.widths[1]
    end

    mktempdir() do directory
        for suffix in ("png", "pdf")
            path = joinpath(directory, "earth.$suffix")
            CairoMakie.save(path, plotted.figure)
            @test filesize(path) > 1000
        end
        saved_limits = axis.finallimits[]
        path = export_svg(plotted; path = joinpath(directory, "earth.svg"), open_file = false)
        @test filesize(path) > 1000
        @test axis.finallimits[] == saved_limits
        @test axis.scene.plots == scene_plots
    end
end
