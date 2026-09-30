# Historical pre-feature qualification override. This requires the original
# geometry owner recorded in production-before.sha256; it deliberately stops
# if that source no longer matches. Current reproduction uses qualify.jl and
# the public pml_resolution control instead.
# The production file is read but never
# edited. Existing conductor topology, fields and native export remain owned
# by the extension. Only the Cartesian PML patch construction is substituted.
include("strips.jl")
const PML_DIAGNOSTIC_DESIGN=Ref(:interpolation)
const PML_QUALIFICATION_LABEL=Ref("diagnostic")

function proposed_strips(plan,model)
    mu0=4pi*1e-7; eps0=8.8541878128e-12; omega=2pi*plan.frequency
    earth=model.earth_materials[plan.frequency_index]
    air=model.problem.earth_props.layers[1]
    gamma_air=sqrt(complex(-omega^2*mu0*air.mu_r*eps0*air.eps_r))
    gamma_earth=sqrt(complex(-omega^2*mu0*earth.mu_r*eps0*earth.eps_r,omega*mu0*earth.mu_r/earth.rho))
    occupied=maximum(zip(model.problem.system.designs,model.problem.system.positions)) do (design,position)
        max(abs(position.x-first(model.centre)),abs(position.y))+LineCableModels.outer_radius(design)
    end
    designs=map(enumerate((:side,:top,:bottom))) do (i,direction)
        legacy = (PML_DIAGNOSTIC_DESIGN[] == :side_only && direction != :side) ||
            (PML_DIAGNOSTIC_DESIGN[] == :top_only && direction != :top)
        if legacy
            count = direction == :bottom ? 96 : 144
            return [(start=0.,stop=1.,count,ratio=exp(((192/191)*log(1536))/count))]
        end
        case=(;air=gamma_air,earth=gamma_earth,Gamma=0im,L=plan.pml_thickness[i],
            clearance=plan.domain_halfwidth-occupied,
            side_strength=plan.pml_strength[i],bottom_strength=plan.pml_strength[i])
        grid=density_grid(case,direction)
        PML_DIAGNOSTIC_DESIGN[] == :guarded && return guarded_strips(grid).strips
        PML_DIAGNOSTIC_DESIGN[] == :strict_guard && return guarded_strips(grid;relative_variation=.12).strips
        strip_grid(grid,72).strips
    end
    result=(;side=designs[1],top=designs[2],bottom=designs[3])
    path=joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-physical-mesh/resolved",PML_QUALIFICATION_LABEL[],
        string(PML_DIAGNOSTIC_DESIGN[],"-",plan.frequency,".toml"))
    mkpath(dirname(path))
    open(path,"w") do io
        TOML.print(io,Dict(string(d)=>[Dict(string(k)=>v for (k,v) in pairs(s)) for s in strips]
            for (d,strips) in pairs(result)))
    end
    result
end

function install_prototype!()
    extension=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    source=read(joinpath(pkgdir(LineCableModels),"ext/LineCableModelsGmshExt/geometry.jl"),String)
    source=source[first(findfirst("function _build_geometry!(",source)):end]
    function substitute(old,new)
        occursin(old,source) || error("geometry owner changed; cannot install qualification prototype")
        source=replace(source,old=>new;count=1)
    end
    substitute("    xs = sort!(unique([left-side,left,buried_x...,right,right+side]))", """
    prescription = Main.proposed_strips(mesh_plan,model)
    ns, nt, nb = length(prescription.side), length(prescription.top), length(prescription.bottom)
    inner_xs = sort!(unique([left,buried_x...,right]))
    xs = [left .- side .* reverse([s.stop for s in prescription.side]);
        inner_xs; right .+ side .* [s.stop for s in prescription.side]]
    """)
    substitute("reverse(xs[2:end-1])","reverse(inner_xs)")
    substitute("xs[2:end-1]","inner_xs")
    start=first(findfirst("    side_layers, top_layers, bottom_layers = mesh_plan.pml_layers",source))
    stop=first(findfirst("    vertices = ",source))
    source=source[1:start-1]*"""
    ys = [-halfwidth .- bottom .* reverse([s.stop for s in prescription.bottom]);
        -halfwidth; 0.0; halfwidth; halfwidth .+ top .* [s.stop for s in prescription.top]]
    """*source[stop:end]
    substitute("    for i in 1:length(xs)-1, j in 1:4", "    for i in 1:length(xs)-1, j in 1:length(ys)-1")
    start=first(findfirst("        interior = 1 < i < length(xs)-1",source))
    stop=first(findfirst("        a, b, c, d = vertices",source))
    source=source[1:start-1]*"""
        interior = ns < i < length(xs)-ns
        physical_row = nb < j <= nb+2
        interior && physical_row && continue
        medium = ys[j] >= 0 ? 1 : 2
        nx, rx = if interior
            (max(2,ceil(Int,(xs[i+1]-xs[i])/mesh_plan.exterior_mesh_sizes[medium])),1.0)
        else
            strip = i <= ns ? prescription.side[ns+1-i] : prescription.side[i-(length(xs)-ns-1)]
            (strip.count, i <= ns ? inv(strip.ratio) : strip.ratio)
        end
        ny, ry = if physical_row
            count, ratio = vertical_grading[medium]
            (count, medium == 2 ? inv(ratio) : ratio)
        else
            strip = j <= nb ? prescription.bottom[nb+1-j] : prescription.top[j-nb-2]
            (strip.count, j <= nb ? inv(strip.ratio) : strip.ratio)
        end
    """*source[stop:end]
    substitute("push!(j >= 3 ? air_pml_surfaces : earth_pml_surfaces, surface)","push!(medium == 1 ? air_pml_surfaces : earth_pml_surfaces, surface)")
    substitute("outer = j >= 3 ? outer_air_curves : outer_earth_curves","outer = medium == 1 ? outer_air_curves : outer_earth_curves")
    substitute("j == 4 && push!(outer,upper)","j == length(ys)-1 && push!(outer,upper)")
    substitute("j == 2 && !interior && push!(interface_curves,upper)","j == nb+1 && !interior && push!(interface_curves,upper)")
    substitute("            push!(curves, mesh_line(reference, inner, bottom_layers, inv(bottom_ratio)))", """
            previous = reference
            for strip in reverse(prescription.bottom)
                next = _point!(registry, (endpoint[1], -halfwidth-bottom*strip.start);
                    mesh_size=mesh_plan.exterior_mesh_sizes[2])
                push!(curves, mesh_line(previous, next, strip.count, inv(strip.ratio)))
                previous = next
            end
    """)
    Base.include_string(extension,source,"qualification_pml_geometry.jl")
    nothing
end
