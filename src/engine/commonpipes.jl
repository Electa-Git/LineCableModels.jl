_enclosure_boundary(region)=any(
    entry->entry.pattern isa DataModel.EnclosureBoundary,region.placement.patterns
)

function _assembly_radius(conductors,dielectrics,dielectric_ranges,range)
    index=last(range)
    layers=dielectric_ranges[index]
    return isempty(layers) ? conductors[index].r_ex : dielectrics[last(layers)].r_ex
end

# Only explicit physical Enclosure boundaries establish a common pipe.
# Containment chooses the direct parent when circular enclosures are nested.
function _pipe_assemblies(design,conductors,dielectrics,dielectric_ranges,ranges)
    T=eltype(first(conductors).material)
    raw=Tuple{Int,Material{T}}[]
    for region in design.geometry.regions
        shape=region.primitive
        shape isa DataModel.DifferenceShape || continue
        _enclosure_boundary(region) || continue
        outer=shape.outer
        outer isa DataModel.Disk || continue
        centre=(T(outer.at.x),T(outer.at.y))
        matches=filter(design.geometry.regions) do wall
            a=wall.primitive
            wall.terminal!==nothing && a isa DataModel.Annulus &&
                _enclosure_boundary(wall) && isapprox(a.ri,outer.r) &&
                DataModel.same_radial_position((T(a.at.x),T(a.at.y)),centre)
        end
        isempty(matches) && continue
        length(matches)==1 || throw(ArgumentError("common pipe must have one inner wall boundary"))
        wall=only(matches)
        index=findfirst(c->c.terminal===wall.terminal,conductors)
        index===nothing && throw(ArgumentError("common pipe wall must retain a terminal"))
        c=conductors[index]; a=wall.primitive
        isapprox(c.r_in,a.ri) && isapprox(c.r_ex,a.ro) ||
            throw(ArgumentError("common pipe wall requires one homogeneous circular annulus"))
        index==first(ranges[c.assembly]) || throw(ArgumentError(
            "a nonconcentric common pipe must begin its own radial assembly"
        ))
        fill=convert(Material{T},region.source.material)
        fill.kind===:insulator && isapprox(fill.mu_r,one(T)) ||
            throw(ArgumentError("common pipe cavity requires nonmagnetic insulation"))
        push!(raw,(index,fill))
    end
    allunique(first.(raw)) || throw(ArgumentError("common pipe boundaries must be distinct"))
    parent=zeros(Int,length(ranges))
    for (assembly,range) in pairs(ranges)
        c=conductors[first(range)]
        radius=_assembly_radius(conductors,dielectrics,dielectric_ranges,range)
        candidates=Int[]
        for (pipe_index,(index,_)) in pairs(raw)
            p=conductors[index]
            p.assembly==assembly && continue
            distance=hypot(c.position[1]-p.position[1],c.position[2]-p.position[2])
            distance+radius<p.r_in && push!(candidates,pipe_index)
        end
        if !isempty(candidates)
            chosen=argmin(k->conductors[raw[k][1]].r_in,candidates)
            parent[assembly]=chosen
        end
    end
    pipes=PipeAssembly{T}[]
    for (pipe_index,(index,material)) in pairs(raw)
        children=findall(==(pipe_index),parent)
        isempty(children) && throw(ArgumentError("common pipe has no represented interior assembly"))
        push!(pipes,PipeAssembly{T}(index,children,material))
    end
    sort!(pipes;by=p->conductors[p.conductor].r_in)
    return pipes
end

function _assembly_groups(cable::LocalCableData)
    groups=[collect(range) for range in cable.assemblies]
    parent=zeros(Int,length(groups))
    for pipe in cable.pipes
        assembly=findfirst(range->pipe.conductor in range,cable.assemblies)
        assembly===nothing && throw(ArgumentError("common pipe conductor has no assembly"))
        for child in pipe.children
            parent[child]==0 || throw(ArgumentError("an assembly cannot have two direct pipes"))
            parent[child]=assembly
            append!(groups[assembly],groups[child])
        end
        sort!(groups[assembly])
    end
    return groups,parent
end

function _external_groups(cable::LocalCableData)
    groups,parent=_assembly_groups(cable)
    roots=findall(iszero,parent)
    return (indices=groups[roots],representatives=[first(cable.assemblies[i]) for i in roots])
end
