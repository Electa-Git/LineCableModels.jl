# Run with the test project on a machine with an OpenGL display:
# julia --project=test test/manual/plotting/mesh_inspection.jl [saved-study.msh]
# Tests use hidden native windows and never mesh or invoke GetDP.
using LineCableModels, Gmsh, GLMakie, GeometryBasics, Test

function validate_mesh_inspection(; output=mktempdir(), meshfile=nothing)
    mkpath(output)
    ext = Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    fixture = joinpath(pkgdir(LineCableModels), "test", "fixtures", "data", "fem", "inspection.msh")
    mesh = import_data(:msh, fixture)
    settings = (; backend=:gl, display_plot=false, open_export=false, size=(1100,700))
    pixel(axis, point) = Point2f(Makie.project(axis.scene,:data,:pixel,Point2d(point))[1:2]) + axis.scene.viewport[].origin
    function click(view, screen, point; raw=false)
        # Hidden screens have no render loop; flush native GPU updates as that loop does.
        GLMakie.poll_updates(screen)
        GLMakie.render_frame(screen)
        xy = pixel(view.axis,point)
        if raw
            events = Makie.events(view.axis.scene)
            events.mouseposition[] = Tuple(xy)
            events.mousebutton[] = Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.press)
            events.mousebutton[] = Makie.MouseButtonEvent(Makie.Mouse.left, Makie.Mouse.release)
        else
            local_xy = xy - view.axis.scene.viewport[].origin
            event = Makie.MouseEvent(Makie.MouseEventTypes.leftclick,1.,Point2d(point),local_xy,0.,Point2d(point),local_xy)
            Makie.process_interaction(last(view.axis.interactions[view.interaction]),event,view.axis)
        end
        return view.selection[]
    end
    function tag(view)
        kind,i = view.selection[]
        kind === :node && return view.mesh.node_tags[i]
        b,c = view.drawing.elements[i]
        return view.mesh.blocks[b].element_tags[c]
    end
    screen = GLMakie.Screen(;visible=false,start_renderloop=false)
    second = GLMakie.Screen(;visible=false,start_renderloop=false)
    try
        @testset "Native mesh clicks" begin
            first_use = @timed LineCableModels.plot(mesh; color_by=:physical,inspect=:element,settings...)
            p = first_use.value
            v = p.addon_state.mesh_inspection.view
            display(screen,p.figure)
            GLMakie.render_frame(screen)
            @test click(v,screen,(.7,.2);raw=true) !== nothing
            @test tag(v) == 305
            @test click(v,screen,(.2,.7)) !== nothing
            @test tag(v) == 305 # Both display triangles refer to the original quad.
            @test click(v,screen,(1.4,.2)) !== nothing
            @test tag(v) == 901
            @test click(v,screen,(.5,0)) !== nothing
            @test tag(v) == 101 # Original boundary line, not a surface display edge.
            click(v,screen,(2.,.9))
            @test v.selection[] === nothing
            p.controls[:mesh_inspect].i_selected[] = 3
            @test click(v,screen,(0,0)) !== nothing
            @test tag(v) == 7 # Coincident node #2007 remains a distinct record.
            click(v,screen,(1,1))
            @test tag(v) == 64
            limits = p.axes[1].finallimits[]
            for _ in 1:5
                click(v,screen,(0,0))
                @test tag(v) == 7
            end
            @test p.axes[1].finallimits[] == limits
            selected = v.selection[]
            event = Makie.MouseEvent(Makie.MouseEventTypes.leftdrag,1.,Point2d(0),Point2f(10),0.,Point2d(0),Point2f(0))
            Makie.process_interaction(last(v.axis.interactions[v.interaction]),event,v.axis)
            @test v.selection[] == selected
            # Native scrolling and panning must not invalidate display/native maps.
            events = Makie.events(v.axis.scene)
            xy = pixel(v.axis,(.5,.5))
            events.mouseposition[] = Tuple(xy)
            Makie.process_interaction(last(v.axis.interactions[:scrollzoom]),Makie.ScrollEvent(0,-1),v.axis)
            Makie.process_interaction(last(v.axis.interactions[:dragpan]),
                Makie.MouseEvent(Makie.MouseEventTypes.rightdrag,1.,Point2d(0),xy.+Point2f(5,5),0.,Point2d(0),xy),v.axis)
            @test click(v,screen,(1,1)) !== nothing
            @test tag(v) == 64
            resize!(p.figure,1200,800)
            @test click(v,screen,(1,1)) !== nothing
            @test tag(v) == 64
            p.controls[:reset].clicks[] += 1
            @test click(v,screen,(1,1)) !== nothing
            @test tag(v) == 64
            for _ in 1:4
                p.controls[:mesh_inspect].i_selected[] = 1
                p.controls[:mesh_inspect].i_selected[] = 2
                p.controls[:mesh_color_by].i_selected[] = 3
                p.controls[:mesh_color_by].i_selected[] = 2
            end
            @test click(v,screen,(.7,.2)) !== nothing
            @test tag(v) == 305
            v.wireframe.visible = false
            GLMakie.poll_updates(screen)
            @test !v.faces[].visible[]
            @test v.selection[] === nothing
            v.wireframe.visible = true
            @test click(v,screen,(.7,.2)) !== nothing
            @test tag(v) == 305
            Makie.save(joinpath(output,"mesh-gl.png"),p.figure)
            warmed = @timed LineCableModels.plot(mesh; color_by=:physical,inspect=:element,settings...)
            q = warmed.value
            display(second,q.figure)
            @test q.addon_state.mesh_inspection.view.selection[] === nothing
            @test v.selection[] !== nothing
            @test click(q.addon_state.mesh_inspection.view,second,(3.4,.2)) !== nothing
            @test tag(q.addon_state.mesh_inspection.view) == 9001
            @test tag(v) == 305
            println("Fixture creation first/warm seconds: ", (first_use.time,warmed.time))
            println("Fixture creation first/warm allocated bytes: ", (first_use.bytes,warmed.bytes))
            copper = Material(kind=:conductor,rho=1.72e-8)
            dielectric = Material(kind=:insulator,rho=1e14,eps_r=3.5)
            wire = build(CableDesign,"overlay",terminal(:core,core(copper;r=.2),insulation(dielectric;t=.05)))
            system = build(LineCableSystem,wire,(.5,.5);connections=Dict(:core=>1))
            g = preview(system;mesh,mesh_color_by=:physical,mesh_inspect=:element,settings...)
            display(second,g.figure)
            gv = g.addon_state.mesh_inspection.view
            @test click(gv,second,(.55,.5)) !== nothing # Inside the geometry outline.
            @test tag(gv) == 305
            @test click(gv,second,(.725,.5)) !== nothing # Inside patterned insulation.
            @test tag(gv) == 305
            g.controls[:mesh_color_by].i_selected[] = 1
            g.controls[:mesh_inspect].i_selected[] = 1
            GLMakie.render_frame(second)
            g.controls[:mesh_color_by].i_selected[] = 2
            g.controls[:mesh_inspect].i_selected[] = 2
            @test click(gv,second,(.725,.5)) !== nothing
            @test tag(gv) == 305
            Makie.save(joinpath(output,"preview-gl.png"),g.figure)
            @test !Bool(Gmsh.gmsh.is_initialized())
            empty!(q.figure)
            @test isempty(q.figure.content)
            @test click(v,screen,(.7,.2)) !== nothing
            @test tag(v) == 305
        end
        if meshfile !== nothing
            @testset "Saved study mesh" begin
                imported = @timed import_data(:msh,meshfile)
                large = imported.value
                created = @timed LineCableModels.plot(large;color_by=:physical,inspect=:element,settings...)
                p = created.value
                v = p.addon_state.mesh_inspection.view
                display(screen,p.figure)
                # Zoom to a known original face; full-domain dense meshes are pixel ambiguous.
                b = findfirst(b -> b.dimension==2,large.blocks)
                nodes = large.blocks[b].connectivity[:,1]
                xy = large.coordinates[1:2,nodes]
                center = vec(sum(xy;dims=2)/size(xy,2))
                span = maximum(maximum(xy;dims=2)-minimum(xy;dims=2))
                limits!(v.axis,center[1]-2span,center[1]+2span,center[2]-2span,center[2]+2span)
                @test click(v,screen,center) !== nothing
                @test tag(v) == large.blocks[b].element_tags[1]
                times = [(@timed click(v,screen,center)) for _ in 1:10]
                println("Saved mesh nodes/elements: ",(length(large.node_tags),sum(length(b.element_tags) for b in large.blocks)))
                println("Import/create seconds and allocated bytes: ",((imported.time,imported.bytes),(created.time,created.bytes)))
                println("Warm click seconds: ",[t.time for t in times])
                println("Warm click allocated bytes: ",[t.bytes for t in times])
                println("Retained mesh/drawing-map bytes (excluding native render resources): ",
                    (Base.summarysize(large),Base.summarysize(v.drawing)+Base.summarysize(v.face_owners)))
                autolimits!(v.axis)
                Makie.save(joinpath(output,"study-gl.png"),p.figure)
                @test !Bool(Gmsh.gmsh.is_initialized())
            end
        end
    finally
        close(screen)
        close(second)
    end
    println("Mesh inspection artifacts: ",output)
    return output
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_mesh_inspection(;output=get(ENV,"LINECABLEMODELS_MESH_CHECK_OUTPUT",mktempdir()),
        meshfile=isempty(ARGS) ? nothing : first(ARGS))
end
