# Read-only arithmetic/mesh inspection of the four saved algebraic executions.
# No meshing, assembly, solve or production option is invoked.
include("assess_algebra.jl")
using Gmsh
const gmsh = Gmsh.gmsh
const LOCALIZATION = joinpath(OUT, "residual-localization")
const FAMILIES = Dict(1=>"axial_a", 2=>"terminal_U", 5=>"transverse_bt",
    6=>"scalar_v", 7=>"terminal_V")
const MEDIA = (1001,1002,1003,1004)
const PRIMARY = (MEDIA...,10001)

function preprocessing(path)
    lines=readlines(path)
    start=findfirst(startswith(raw"$DofData"),lines)
    records,n=parse.(Int,split(lines[start+5]))
    rows=fill((family=0,entity=0,harmonic=0),n)
    fixed=Dict{Tuple{Int,Int},Vector{Int}}()
    for line in lines[start+6:start+5+records]
        fields=split(line)
        family,entity,harmonic,kind=parse.(Int,fields[1:4])
        harmonic==0 || error("Unexpected harmonic index")
        if kind==1
            row=parse(Int,fields[5])
            rows[row].family==0 || error("Duplicate row")
            rows[row]=(;family,entity,harmonic)
        elseif kind==2
            push!(get!(fixed,(family,kind),Int[]),entity)
        else
            error("Unsupported saved constraint type $kind")
        end
    end
    all(r->haskey(FAMILIES,r.family),rows) || error("Unmapped equation family")
    start=findfirst(==(raw"$ElementsXEdges"),lines)
    _,count=parse.(Int,split(lines[start+1]))
    incidence=[parse.(Int,split(line)) for line in lines[start+2:start+1+count]]
    return (;rows,fixed,incidence)
end

function import_mesh(path,pre)
    gmsh.clear()
    gmsh.open(path)
    nt,xyz,_=gmsh.model.mesh.get_nodes()
    coordinates=zeros(3,maximum(nt))
    coordinates[:,Int.(nt)]=reshape(xyz,3,:)
    _,alltags,_=gmsh.model.mesh.get_elements()
    maxtag=maximum(maximum(t) for t in alltags)
    elements=NamedTuple{(:num,:tag,:kind,:region,:entity,:nodes),Tuple{Int,Int,Int,Int,Int,Vector{Int}}}[]
    terminalnodes=Dict(3001=>Set{Int}(),3002=>Set{Int}())
    boundarypairs=Dict{Tuple{Int,Int},Set{Int}}()
    # Reproduce GetDP 3.5.0 Geo_ReadFileWithGmsh: entity -> physical group ->
    # element type -> element; additional memberships get new tags, then sort.
    for (dim,entity) in gmsh.model.get_entities()
        kinds,tags,conn=gmsh.model.mesh.get_elements(dim,entity)
        physical=Int.(gmsh.model.get_physical_groups_for_entity(dim,entity))
        for (p,region) in enumerate(physical), k in eachindex(kinds)
            kind=kinds[k]
            kind in (1,2,15) || error("Only retained linear triangles/lines/points expected")
            width=kind==2 ? 3 : kind==1 ? 2 : 1
            nodes=reshape(conn[k],width,:)
            for j in eachindex(tags[k])
                tag=Int(tags[k][j]); ns=Int.(nodes[:,j])
                if p==1
                    num=tag
                else
                    maxtag+=1; num=maxtag
                end
                push!(elements,(;num,tag,kind,region,entity,nodes=ns))
                region in (3001,3002) && union!(terminalnodes[region],ns)
                if kind==1
                    union!(get!(boundarypairs,minmax(ns...),Set{Int}()),(region,))
                end
            end
        end
    end
    sort!(elements;by=e->e.num)
    edgepairs=Dict{Int,Tuple{Int,Int}}()
    pairedges=Dict{Tuple{Int,Int},Int}()
    checked=0
    for row in pre.incidence
        index,ne=row[1:2]
        e=elements[index+1] # .pre stores a zero-based sorted imported-element index
        pairs=e.kind==2 ? ((1,2),(1,3),(2,3)) : e.kind==1 ? ((1,2),) : ()
        length(pairs)==ne || error("Saved edge count disagrees with imported element")
        for (s,(a,b)) in enumerate(pairs)
            id=row[s+2]; pair=minmax(e.nodes[a],e.nodes[b]); edge=abs(id)
            sign(id)==sign(e.nodes[b]-e.nodes[a]) || error("Edge orientation mismatch")
            get(edgepairs,edge,pair)==pair || error("Edge maps to different node pairs")
            get(pairedges,pair,edge)==edge || error("Node pair maps to different edges")
            edgepairs[edge]=pair; pairedges[pair]=edge
            checked+=1
        end
    end
    nodecells=[Int[] for _ in axes(coordinates,2)]
    edgecells=Dict{Int,Vector{Int}}()
    primary=filter(e->e.kind==2 && e.region in PRIMARY,elements)
    for (idx,e) in enumerate(primary)
        for node in e.nodes
            push!(nodecells[node],idx)
        end
        for (a,b) in ((1,2),(1,3),(2,3))
            pair=minmax(e.nodes[a],e.nodes[b])
            # Edges wholly in metal need not have a transverse basis or .pre entry.
            haskey(pairedges,pair) || continue
            push!(get!(edgecells,pairedges[pair],Int[]),idx)
        end
    end
    for r in pre.rows
        r.family==5 && !haskey(edgepairs,r.entity) && error("Unmapped transverse edge")
    end
    return (;coordinates,elements,primary,nodecells,edgecells,edgepairs,terminalnodes,boundarypairs,
        checked,edge_count=length(edgepairs))
end

function support_maps(mesh,pre,data)
    halfwidth=parse(Float64,match(r"DomainHalfwidths\(\) = \{([^}]+)\}",data)[1])
    xc=parse(Float64,match(r"Xcenter = ([^;]+)",data)[1])
    function region_name(e)
        e.region==10001 && return "conductor"
        medium=e.region in (1001,1003) ? "air" : "earth"
        e.region in (1001,1002) && return medium*"/physical"
        center=sum(mesh.coordinates[:,node] for node in e.nodes)/3
        pieces=String[]
        center[1]<xc-halfwidth && push!(pieces,"left")
        center[1]>xc+halfwidth && push!(pieces,"right")
        center[2]>halfwidth && push!(pieces,"top")
        center[2]<-halfwidth && push!(pieces,"bottom")
        isempty(pieces) && error("PML element in physical domain")
        return medium*"/PML/"*join(pieces,"+")
    end
    names=region_name.(mesh.primary)
    function row_support(r)
        if r.family==5
            indices=get(mesh.edgecells,r.entity,Int[])
        elseif r.family in (1,6)
            indices=mesh.nodecells[r.entity]
        else
            indices=sort!(unique!(reduce(vcat,(mesh.nodecells[node] for node in mesh.terminalnodes[r.entity]);init=Int[])))
        end
        if r.family in (5,6,7)
            indices=filter(i->mesh.primary[i].region in MEDIA,indices)
        elseif r.family==2
            indices=filter(i->mesh.primary[i].region==10001,indices)
        end
        isempty(indices) && error("No active weak-form support: $r")
        # Cross-region support is one explicit union category, never apportioned.
        return join(sort!(unique(names[indices]))," | ")
    end
    labels=row_support.(pre.rows)
    return (;labels,names,halfwidth,xc)
end

function residual_vector(A,b,x)
    setprecision(BigFloat,128) do
        xb=Complex{BigFloat}.(x)
        r=Vector{ComplexF64}(undef,A.n)
        for i in 1:A.n
            ri=Complex{BigFloat}(b[i])
            for k in A.rowptr[i]:A.rowptr[i+1]-1
                ri-=Complex{BigFloat}(A.values[k])*xb[A.cols[k]]
            end
            r[i]=ri
        end
        return r
    end
end

function group_energy(r,keys)
    sums=Dict{String,Float64}(); counts=Dict{String,Int}()
    for (ri,key) in zip(r,keys)
        sums[key]=get(sums,key,0.)+abs2(ri)
        counts[key]=get(counts,key,0)+1
    end
    total=sum(abs2,r)
    return [Dict("group"=>key,"rows"=>counts[key],"squared_residual"=>energy,"fraction"=>energy/total)
        for (key,energy) in sort!(collect(sums);by=p->-last(p))]
end

function row_metadata(i,pre,mesh,support)
    row=pre.rows[i]
    nodes=row.family==5 ? collect(mesh.edgepairs[row.entity]) : row.family in (1,6) ? [row.entity] : sort!(collect(mesh.terminalnodes[row.entity]))
    center=sum(mesh.coordinates[:,n] for n in nodes)/length(nodes)
    record=Dict{String,Any}("row"=>i,"basis_code"=>row.family,"entity"=>row.entity,
        "equation_family"=>FAMILIES[row.family],"support_union"=>support.labels[i],
        "nodes"=>nodes,"x"=>center[1],"y"=>center[2])
    if row.family==5
        record["line_physical_groups"]=sort!(collect(get(mesh.boundarypairs,Tuple(nodes),Set{Int}())))
        cells=filter(j->mesh.primary[j].region in MEDIA,get(mesh.edgecells,row.entity,Int[]))
    elseif row.family in (1,6)
        cells=mesh.nodecells[row.entity]
        row.family==6 && (cells=filter(j->mesh.primary[j].region in MEDIA,cells))
    else
        cells=Int[]
    end
    record["support_element_tags"]=[mesh.primary[j].tag for j in cells]
    record["support_geometry_entities"]=sort!(unique([mesh.primary[j].entity for j in cells]))
    return record
end

function row_terms(A,b,x,i,pre)
    setprecision(BigFloat,128) do
        sums=Dict(name=>zero(Complex{BigFloat}) for name in values(FAMILIES))
        terms=Dict{String,Any}[]
        for k in A.rowptr[i]:A.rowptr[i+1]-1
            j=A.cols[k]; family=FAMILIES[pre.rows[j].family]
            value=Complex{BigFloat}(A.values[k])*Complex{BigFloat}(x[j])
            sums[family]+=value
            iszero(value) && continue
            push!(terms,Dict("column"=>j,"family"=>family,"entity"=>pre.rows[j].entity,
                "A_real"=>real(A.values[k]),"A_imag"=>imag(A.values[k]),
                "x_real"=>real(x[j]),"x_imag"=>imag(x[j]),
                "Ax_real"=>Float64(real(value)),"Ax_imag"=>Float64(imag(value))))
        end
        residual=Complex{BigFloat}(b[i])-sum(values(sums))
        return Dict("rhs_real"=>real(b[i]),"rhs_imag"=>imag(b[i]),
            "residual_real"=>Float64(real(residual)),"residual_imag"=>Float64(imag(residual)),
            "block_Ax"=>Dict(k=>[Float64(real(v)),Float64(imag(v))] for (k,v) in sums),
            "individual_terms"=>sort!(terms;by=t->-hypot(t["Ax_real"],t["Ax_imag"])))
    end
end

function compare_strip_rows(reference,current,selected)
    output=Dict{String,Any}()
    refnode=Dict(Tuple(reference.mesh.coordinates[:,r.entity])=>i
        for (i,r) in enumerate(reference.pre.rows) if r.family==1)
    for i in selected
        point=Tuple(current.mesh.coordinates[:,current.pre.rows[i].entity])
        j=refnode[point]
        item=Dict{String,Any}("candidate_row"=>i,"baseline_row"=>j)
        stencils=Dict()
        for (name,state,row) in (("baseline",reference,j),("metric-1.35",current,i))
            A,pre=state.A,state.pre
            # These are axial-only rows in the Gamma=0 physical/strip stencil.
            entries=Dict{Tuple{Int,NTuple{3,Float64}},ComplexF64}()
            for k in A.rowptr[row]:A.rowptr[row+1]-1
                iszero(A.values[k]) && continue
                col=pre.rows[A.cols[k]]
                col.family==1 || error("Expected axial-only stencil")
                entries[(col.family,Tuple(state.mesh.coordinates[:,col.entity]))]=A.values[k]
            end
            stencils[name]=entries
            for config in CONFIGURATIONS,basis in 1:2
                vs=state.Vs[config]
                item["$name-$config-source$basis"]=row_terms(A,vs[2basis-1],vs[2basis],row,pre)
            end
        end
        item["nonzero_stencil_keys_equal"]=keys(stencils["baseline"])==keys(stencils["metric-1.35"])
        item["nonzero_stencil_keys_equal"] || error("Unmatched strip stencil")
        difference=maximum(abs(v-stencils["baseline"][key]) for (key,v) in stencils["metric-1.35"])
        item["maximum_coefficient_difference"]=difference
        item["max_difference_over_max_coefficient"]=difference/maximum(abs,values(stencils["baseline"]))
        output["row$i"]=item
    end
    return output
end

function localize()
    mkpath(LOCALIZATION)
    gmsh.initialize(); gmsh.option.set_number("General.Terminal",0)
    reference=nothing
    try
        for name in MESHES
            say("MAP saved preprocessing and mesh ",name)
            dir=joinpath(OUT,name,"original")
            pre=preprocessing(joinpath(dir,"study.pre"))
            mesh=import_mesh(joinpath(dir,"study.msh"),pre)
            support=support_maps(mesh,pre,read(joinpath(dir,"study_data.pro"),String))
            record(joinpath(LOCALIZATION,name*"-mapping.toml"),Dict("all_rows_mapped"=>true,
                "rows"=>length(pre.rows),"checked_signed_edge_incidences"=>mesh.checked,
                "unique_edges"=>mesh.edge_count,"imported_elements"=>length(mesh.elements),
                "fixed_dof_counts"=>Dict("basis$(k[1])"=>length(v) for (k,v) in pre.fixed),
                "primary_triangles"=>length(mesh.primary)))
            A=petsc_matrix(joinpath(dir,"file_mat_before1.m.bin"))
            Rs=Dict{Tuple{String,Int},Vector{ComplexF64}}()
            Vs=Dict{String,Vector{Vector{ComplexF64}}}()
            dominant=Set{Int}()
            for config in CONFIGURATIONS
                Vs[config]=binary_vectors(joinpath(OUT,name,config,"study.res"),A.n)
                for basis in 1:2
                    say("LOCALIZE original-coordinate residual ",name,"/",config," source=",basis)
                    vs=Vs[config]; r=residual_vector(A,vs[2basis-1],vs[2basis]); Rs[(config,basis)]=r
                    order=partialsortperm(abs2.(r),1:12;rev=true)
                    union!(dominant,order)
                    total=sum(abs2,r)
                    top=[merge(row_metadata(i,pre,mesh,support),Dict("fraction"=>abs2(r[i])/total,
                        "residual_real"=>real(r[i]),"residual_imag"=>imag(r[i]))) for i in order]
                    record(joinpath(LOCALIZATION,"$name-$config-source$basis.toml"),Dict(
                        "squared_residual"=>total,"equation_families"=>group_energy(r,[FAMILIES[row.family] for row in pre.rows]),
                        "support_regions"=>group_energy(r,support.labels),
                        "family_and_support"=>group_energy(r,[FAMILIES[row.family]*" / "*label for (row,label) in zip(pre.rows,support.labels)]),
                        "dominant_rows"=>top))
                end
            end
            traces=Dict{String,Any}()
            # The union of the top rows is inspected under ALL four saved solutions.
            for i in sort!(collect(dominant))
                trace=row_metadata(i,pre,mesh,support)
                for config in CONFIGURATIONS,basis in 1:2
                    vs=Vs[config]
                    trace["$config-source$basis"]=row_terms(A,vs[2basis-1],vs[2basis],i,pre)
                end
                traces["row$i"]=trace
            end
            record(joinpath(LOCALIZATION,name*"-row-terms.toml"),traces)
            if name=="baseline"
                reference=(;A,pre,mesh,Vs)
            else
                indices=findall(i->pre.rows[i].family==1 && support.labels[i]=="earth/PML/left",eachindex(pre.rows))
                sort!(indices;by=i->-abs2(Rs[("original",1)][i]))
                record(joinpath(LOCALIZATION,"unchanged-strip-row-comparison.toml"),
                    compare_strip_rows(reference,(;A,pre,mesh,Vs),indices[1:2]))
            end
            say("MAPPED ",name," rows=",A.n," signed edge incidences checked=",mesh.checked)
        end
    finally
        gmsh.finalize()
    end
    say("COMPLETE residual localization arithmetic; no FEM execution")
end

if abspath(PROGRAM_FILE)==@__FILE__
    localize()
end
