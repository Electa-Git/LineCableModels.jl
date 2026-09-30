# Static qualification figures from saved meshes. No meshing or solving.
using LineCableModels, Gmsh, CairoMakie, TOML
root = normpath(joinpath(@__DIR__,"../../../../.linecablemodels/fem/local-interface-mesh"))
case = isempty(ARGS) ? "air-f0.1-gamma0.0" : ARGS[1]
candidate = length(ARGS)>1 ? ARGS[2] : "localized"
output = mkpath(joinpath(root,"both-media-figures",case,candidate))
meshes = [import_data(:msh,joinpath(root,case,variant,"study.msh"))
    for variant in ("baseline",candidate)]

function medium_edges(mesh, tag)
    edges = Set{Tuple{Int,Int}}()
    for block in mesh.blocks
        block.dimension==2 && tag in block.physical_tags || continue
        for nodes in eachcol(block.connectivity), i in 1:length(nodes)
            push!(edges,minmax(nodes[i],nodes[mod1(i+1,length(nodes))]))
        end
    end
    return [Point2f(mesh.coordinates[1,node],mesh.coordinates[2,node]) for edge in edges for node in edge]
end
edges = [[medium_edges(mesh,tag) for tag in (1001,1002,6001)] for mesh in meshes]
limits = extrema(first.(vcat(edges[1][1],edges[1][2])))
radius = (limits[2]-limits[1])/2
for (name,height) in (("physical-domain",radius),("interface-band",radius/5))
    fig = Figure(size=(1800,900))
    for i in 1:2
        stats = TOML.parsefile(joinpath(root,case,i==1 ? "baseline" : candidate,"mesh.toml"))
        title = (i==1 ? "Frozen FEM reference" : "Both-media localization") *
            " — $(stats["nodes"]) nodes\nAir/soil triangles: $(stats["air_triangles"])/$(stats["soil_triangles"])"
        ax = Axis(fig[i,1];title,xlabel="Horizontal coordinate [km]",ylabel="z [km]")
        for (points,color) in zip(edges[i],(:teal,:darkorange,:black))
            linesegments!(ax,points./1000;color,linewidth=.4)
        end
        xlims!(ax,limits[1]/1000,limits[2]/1000)
        ylims!(ax,-height/1000,height/1000)
    end
    Label(fig[0,1],"$case — identical bounds and conductor/PML/path prescriptions",fontsize=20)
    CairoMakie.save(joinpath(output,"$name.png"),fig)
    CairoMakie.save(joinpath(output,"$name.svg"),fig)
end
println("Mesh comparison figures: ",output)
