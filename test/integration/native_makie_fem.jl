@testitem "Makie FEM / independent mesh and field inspection" tags=[:visual] begin
    using CairoMakie, Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    fixtures = joinpath(pkgdir(LineCableModels), "test", "fixtures", "data", "fem")
    meshfile = joinpath(fixtures, "sparse.msh")
    fieldfile = joinpath(fixtures, "discontinuous.pos")
    settings = (; backend = :cairo, display_plot = false, open_export = false)
    withenv("LINECABLEMODELS_GETDP"=>"/unavailable/getdp", "DISPLAY"=>"") do
        mesh = import_data(:msh, meshfile)
        fields = import_data(:pos, fieldfile)
        @test !Bool(Gmsh.gmsh.is_initialized())
        p = Makie.plot(mesh; settings...)
        @test length(p.axes) == 1
        @test isempty(p.colorbars)
        meshplot = only(p.addon_state.groups[:mesh])
        @test length(meshplot[1][]) == 10 # Five distinct triangle edges.
        @test LineCableModels.plot(meshfile; settings...) isa UIPlot
        scalar = first(fields)
        q = LineCableModels.plot(scalar; part = :real, mesh, settings...)
        Makie.colorbuffer(q.figure)
        @test length(q.colorbars) == 1
        colors = only(q.addon_state.groups[:field]).color[]
        @test colors == [1, 2, 3, 9, 10, 11] # Shared corners keep both element-side values.
        @test q.addon_state.spatial_data === scalar
        @test only(q.colorbars).label[] == "Real part"
        imag = LineCableModels.plot(fieldfile; view = 1, part = :imag, settings...)
        @test only(imag.addon_state.groups[:field]).color[] == [4, 5, 6, 12, 13, 14]
        magnitude = LineCableModels.plot(last(fields); part = :magnitude, settings...)
        @test only(magnitude.addon_state.groups[:field]).color[] == [5, 5, 5]
        phase = LineCableModels.plot(scalar; part = :phase, settings...)
        @test only(phase.addon_state.groups[:field]).color[] ≈
              angle.(complex.([1, 2, 3, 9, 10, 11], [4, 5, 6, 12, 13, 14]))
        vectors = LineCableModels.plot(
            last(fields); component = 1, part = :real, arrows = true,
            arrow_attributes = (; normalize = true, lengthscale = 0.1), settings...)
        @test length(vectors.addon_state.groups[:field]) == 2
        arrows = last(vectors.addon_state.groups[:field])
        @test arrows.normalize[]
        @test arrows.lengthscale[] ≈ 0.1
        Makie.colorbuffer(vectors.figure)
        logplot = LineCableModels.plot(scalar; colorscale = log10, settings...)
        @test only(logplot.colorbars).scale[] === log10
        @test only(logplot.addon_state.groups[:field]).color[] == colors
        @test_throws ArgumentError LineCableModels.plot(
            last(fields); part = :real, component = 1,
            step = 2, colorscale = log10, settings...)
        @test haskey(q.addon_state.groups, :mesh_overlay)
        @test_throws ArgumentError LineCableModels.plot(fieldfile; settings...)
        @test_throws ArgumentError LineCableModels.plot(last(fields); part = :phase, settings...)
        @test_throws ArgumentError LineCableModels.plot(last(fields); part = :real, settings...)

        copper=Material(kind = :conductor, rho = 1.72e-8)
        wire=build(CableDesign, "overlay", terminal(:core, core(copper; r = 0.1)))
        system=build(LineCableSystem, wire, (0.25, 0.25); connections = Dict(:core=>1))
        geometry=preview(system; mesh = meshfile, settings...)
        @test haskey(geometry.addon_state.groups, :mesh)
        @test occursin("y [m]", repr(geometry.axes[1].xlabel[]))
        @test occursin("z [m]", repr(geometry.axes[1].ylabel[]))
        @test !Bool(Gmsh.gmsh.is_initialized())
        mktempdir() do directory
            file=export_svg(q; path = joinpath(directory, "field.svg"), open_file = false)
            @test isfile(file) && filesize(file)>0
        end
    end
end

@testitem "Makie FEM / Copy details sends the current complete record" tags=[:visual] begin
    using CairoMakie, Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    if Sys.islinux()
        # Exercise the native button and Julia's clipboard API, replacing only
        # the OS utilities so this test needs no display or user clipboard access.
        mktempdir() do directory
            capture = joinpath(directory, "copied.txt")
            for name in ("xclip", "xsel", "wl-copy")
                command = joinpath(directory, name)
                write(command, "#!/bin/sh\n/bin/cat > \"\$LCM_CLIPBOARD_CAPTURE\"\n")
                chmod(command, 0o755)
            end
            withenv("PATH"=>directory*":"*get(ENV,"PATH",""),
                    "LCM_CLIPBOARD_CAPTURE"=>capture, "DISPLAY"=>"") do
                fixture = joinpath(pkgdir(LineCableModels), "test", "fixtures", "data", "fem", "inspection.msh")
                mesh = import_data(:msh, fixture)
                p = LineCableModels.plot(mesh; color_by=:physical, inspect=:element,
                    backend=:cairo, display_plot=false, open_export=false)
                ext = Base.get_extension(LineCableModels, :LineCableModelsGmshMakieExt)
                view = p.addon_state.mesh_inspection.view
                quad = findfirst(view.drawing.elements) do (b,c)
                    mesh.blocks[b].element_tags[c] == 305
                end
                ext._fem_mesh_select!(view, (:element,quad))
                p.controls[:mesh_copy].clicks[] += 1
                @test read(capture,String) == view.details[]
                @test occursin("Quadrangle #305", read(capture,String))
                @test occursin("Area [m²]: 1", read(capture,String))
                @test p.status[] == "Mesh details copied"
                view.mode[] = :node
                ext._fem_mesh_select!(view, (:node,findfirst(==(7),mesh.node_tags)))
                p.controls[:mesh_copy].clicks[] += 1
                @test read(capture,String) == view.details[]
                @test startswith(read(capture,String), "Node #7\n")
            end
        end
    else
        @test_skip false # The isolated clipboard-utility fixture targets Linux.
    end
end

@testitem "Makie FEM / physical groups and native diagnostics" tags=[:visual] begin
    using CairoMakie, Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    ext = Base.get_extension(LineCableModels, :LineCableModelsGmshMakieExt)
    fixtures = joinpath(pkgdir(LineCableModels), "test", "fixtures", "data", "fem")
    mesh = import_data(:msh, joinpath(fixtures, "inspection.msh"))
    original = deepcopy(mesh)
    settings = (; backend=:cairo, display_plot=false, open_export=false)
    p = LineCableModels.plot(mesh; color_by=:physical, inspect=:element, settings...)
    view = p.addon_state.mesh_inspection.view
    @test p.controls[:mesh_color_by].selection[] === :physical
    @test p.controls[:mesh_inspect].selection[] === :element
    @test length(view.categories[].labels) == 4
    @test "2-D #8: soil; 2-D #9: computational domain" in view.categories[].labels
    @test "2-D #10: Unnamed" in view.categories[].labels
    @test "Ungrouped" in view.categories[].labels
    @test length(unique(view.faces[].color[])) == 4
    elementindex(tag) = findfirst(view.drawing.elements) do (b,c)
        mesh.blocks[b].element_tags[c] == tag
    end
    quad = elementindex(305)
    @test count(==(quad), view.face_owners) == 4
    details = ext._fem_mesh_details(view, (:element, quad))
    @test occursin("Quadrangle #305", details)
    @test occursin("Nodes: 7, 23, 64, 105", details)
    @test occursin("Area [m²]: 1\n", details)
    @test occursin("Perimeter [m]: 4\n", details)
    @test occursin("Min / max edge [m]: 1 / 1", details)
    boundary = elementindex(101)
    @test occursin("1-D #8: bottom boundary", ext._fem_mesh_details(view, (:element,boundary)))
    @test occursin("Length [m]: 1", ext._fem_mesh_details(view, (:element,boundary)))
    node = findfirst(==(7), mesh.node_tags)
    node_details = ext._fem_mesh_details(view, (:node,node))
    @test occursin("Incident elements (2): 101, 305", node_details)
    @test occursin("1-D #8: bottom boundary", node_details)
    @test occursin("2-D #8: soil", node_details)
    @test occursin("2-D #9: computational domain", node_details)
    ext._fem_mesh_select!(view, (:element,quad))
    Makie.colorbuffer(p.figure)
    limits = p.axes[1].finallimits[]
    ext._fem_mesh_select!(view, (:node,node))
    @test p.axes[1].finallimits[] == limits
    p.controls[:mesh_inspect].i_selected[] = 3
    @test view.mode[] === :node
    @test view.selection[] === nothing
    @test length(view.node_owners) == 14 # Includes a distinct coincident native node.
    p.controls[:mesh_color_by].i_selected[] = 3
    @test view.color_by[] === :entity
    @test "2-D entity #11" in view.categories[].labels
    @test_throws ArgumentError LineCableModels.plot(mesh; color_by=:material, settings...)
    @test_throws ArgumentError LineCableModels.plot(mesh; inspect=:equations, settings...)
    other = LineCableModels.plot(mesh; settings...)
    @test other.addon_state.mesh_inspection.view.mode[] === :none
    @test other.addon_state.mesh_inspection.view.selection[] === nothing

    # Lower-dimensional mesh files use the same native metadata/color contract.
    Engine = LineCableModels.Engine
    lines = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).FEMMesh(mesh.source, mesh.node_tags, mesh.coordinates,
        filter(b -> b.dimension==1, mesh.blocks), mesh.physical_names, mesh.provenance)
    lineplot = LineCableModels.plot(lines; color_by=:physical, inspect=:element, settings...)
    lineview = lineplot.addon_state.mesh_inspection.view
    @test lineview.categories[].labels == ["1-D #8: bottom boundary"]
    lineplot.controls[:mesh_color_by].i_selected[] = 1
    lineplot.controls[:mesh_color_by].i_selected[] = 2
    Makie.colorbuffer(lineplot.figure)
    points = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).FEMMesh(mesh.source, mesh.node_tags, mesh.coordinates,
        [Base.get_extension(LineCableModels,:LineCableModelsGmshExt).FEMElementBlock(15,0,0,1,71,[8],UInt64[808],reshape([1],1,1))],
        Dict((0,8)=>"terminal point"), mesh.provenance)
    pointplot = LineCableModels.plot(points; color_by=:physical, inspect=:element, settings...)
    pointview = pointplot.addon_state.mesh_inspection.view
    ext._fem_mesh_select!(pointview, (:element,1))
    @test occursin("Point #808", pointview.details[])
    @test pointview.marker.visible[]
    Makie.colorbuffer(pointplot.figure)

    blocks = [Base.get_extension(LineCableModels,:LineCableModelsGmshExt).FEMElementBlock(2,2,1,3,i,[i],UInt64[10000+i],reshape([1,2,3],3,1)) for i in 1:25]
    names = Dict((2,i)=>repeat("long name $i ",12) for i in 1:25)
    many = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).FEMMesh(mesh.source, mesh.node_tags, mesh.coordinates, blocks, names, mesh.provenance)
    paged = LineCableModels.plot(many; color_by=:physical, settings...)
    pv = paged.addon_state.mesh_inspection.view
    sidebar = paged.addon_state.mesh_inspection.sidebar
    @test length(unique(pv.categories[].palette)) == 25
    sidebar.following.clicks[] += 1
    @test sidebar.page[] == 2
    @test occursin("9–16 of 25", first(only(sidebar.legend[].entrygroups[])))
    sidebar.previous.clicks[] += 1
    @test sidebar.page[] == 1
    ext._fem_mesh_select!(pv,(:element,1))
    @test occursin("Copy details", sidebar.readout.text[])
    @test occursin(names[(2,1)],pv.details[])
    Makie.colorbuffer(paged.figure)

    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator,rho=1e14,eps_r=3.5)
    wire = build(CableDesign, "overlay", terminal(:core, core(copper; r=.15),insulation(dielectric;t=.05)))
    system = build(LineCableSystem, wire, (.5,.5); connections=Dict(:core=>1))
    geometry = preview(system; mesh, mesh_color_by=:physical, mesh_inspect=:element, settings...)
    @test geometry.addon_state.mesh_inspection.view.color_by[] === :physical
    Makie.colorbuffer(geometry.figure)
    geometry.controls[:mesh_color_by].i_selected[] = 1
    geometry.controls[:mesh_inspect].i_selected[] = 1
    Makie.colorbuffer(geometry.figure)
    geometry.controls[:mesh_color_by].i_selected[] = 2
    Makie.colorbuffer(geometry.figure)
    @test mesh.node_tags == original.node_tags
    @test mesh.coordinates == original.coordinates
    @test mesh.physical_names == original.physical_names
    @test [b.connectivity for b in mesh.blocks] == [b.connectivity for b in original.blocks]
    @test !Bool(Gmsh.gmsh.is_initialized())
    static = LineCableModels.plot(mesh; color_by=:physical, controls=false, settings...)
    @test isempty(static.controls)
    @test static.addon_state.mesh_inspection.view.mode[] === :none
    Makie.colorbuffer(static.figure)
    mktempdir() do directory
        for (name, plot) in (("mesh", p), ("preview", geometry), ("static", static))
            file = export_svg(plot; path=joinpath(directory,name*".svg"), open_file=false)
            @test isfile(file) && filesize(file)>0
            Base.get_extension(LineCableModels,:LineCableModelsMakieExt)._addon_export_presentation!(plot, :default) do
                Makie.colorbuffer(plot.figure)
                viewport = plot.figure.scene.viewport[]
                sidebar = plot.addon_state.mesh_inspection.sidebar
                for block in (sidebar.legend[], sidebar.readout)
                    bounds = block.layoutobservables.computedbbox[]
                    @test all(bounds.origin .>= viewport.origin .- 1)
                    @test all(bounds.origin+ bounds.widths .<= viewport.origin+viewport.widths .+ 1)
                end
            end
        end
    end
end

@testitem "Makie FEM / saved mesh frequency remains visible in plot and preview" tags=[:visual] begin
    using CairoMakie, Gmsh
    E=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    mesh=import_data(:msh,joinpath(pkgdir(LineCableModels),"test/fixtures/data/fem/sparse.msh"))
    provenance=(run_directory="/retained/run-example",frequency_index=3,frequency_hz=1000.,
        terminal_ids=["core"],earth_inputs=(rho=100.,eps_r=1.,mu_r=1.),gamma=0im)
    saved=E.FEMMesh(mesh.source,mesh.node_tags,mesh.coordinates,mesh.blocks,mesh.physical_names,provenance)
    copper=Material(kind=:conductor,rho=1.72e-8)
    design=build(CableDesign,"wire",terminal(:core,Region(:metal,Disk(.01),copper)))
    system=build(LineCableSystem,[design],[Pose2(0.,1.)];connections=[Dict(:core=>1)])
    settings=(backend=:cairo,display_plot=false,open_export=false)
    for p in (LineCableModels.plot(saved;settings...),preview(system;mesh=saved,settings...))
        view=p.addon_state.mesh_inspection.view
        for text in (view.details[],p.axes[1].title[])
            @test occursin("run-example",text)
            @test occursin("frequency 3",text)
            @test occursin("1000.0 Hz",text)
        end
        Base.get_extension(LineCableModels,:LineCableModelsGmshMakieExt)._fem_mesh_select!(view,(:node,1))
        @test occursin("1000.0 Hz",p.axes[1].title[])
        Makie.colorbuffer(p.figure)
    end
end
