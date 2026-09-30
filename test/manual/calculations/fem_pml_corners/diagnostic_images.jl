# Read saved meshes only; no Gmsh generation or GetDP solve.
using LineCableModels, Gmsh, CairoMakie, TOML

function soil_mesh_images(root)
    meshes = [import_data(:msh,joinpath(root,variant,"study.msh"))
        for variant in ("baseline","localized")]
    function edges(mesh, tag)
        segments = Set{Tuple{Int,Int}}()
        for block in mesh.blocks
            block.dimension==2 && tag in block.physical_tags || continue
            for nodes in eachcol(block.connectivity), i in eachindex(nodes)
                push!(segments,minmax(nodes[i],nodes[mod1(i+1,length(nodes))]))
            end
        end
        [Point2f(mesh.coordinates[1,n]/1000,mesh.coordinates[2,n]/1000)
            for pair in segments for n in pair]
    end
    points = [[edges(mesh,tag) for tag in (1001,1002,6001)] for mesh in meshes]
    xmin,xmax = extrema(first.(vcat(points[1][1],points[1][2])))
    height = (xmax-xmin)/10
    fig = Figure(size=(1600,850))
    for i in 1:2
        stats = TOML.parsefile(joinpath(root,i==1 ? "baseline" : "localized","mesh.toml"))
        title = (i==1 ? "Air localized; old full-width soil refinement" : "Fresh export: both air and soil localized") *
            "\nAir: $(stats["air_triangles"]) triangles; soil: $(stats["soil_triangles"]) triangles"
        ax = Axis(fig[i,1];title,xlabel="Horizontal coordinate [km]",ylabel="z [km]")
        for (vertices,color) in zip(points[i],(:teal,:darkorange,:black))
            linesegments!(ax,vertices;color,linewidth=.5)
        end
        xlims!(ax,xmin,xmax); ylims!(ax,-height,height)
    end
    Label(fig[0,1],"Mixed wires, 0.1 Hz: same physical bounds, conductor grading, PML and voltage paths",fontsize=19)
    for ext in ("png","svg")
        CairoMakie.save(joinpath(root,"soil-interface-comparison.$ext"),fig)
    end
    nothing
end

if abspath(PROGRAM_FILE)==@__FILE__
    root = joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-corner-mesh/diagonals-and-soil/soil-localization")
    soil_mesh_images(root)
end
