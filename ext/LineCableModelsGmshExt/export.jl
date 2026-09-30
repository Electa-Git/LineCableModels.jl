const FEM_EXPORT_MESH_OPTIONS = (
    :domain_skin_depths, :pml_thickness, :pml_thickness_factor, :pml_layers, :pml_grading,
    :pml_resolution, :pml_reflection,
    :mesh_size_factor, :exterior_mesh_size_factor,
    :interface_refinement_factor, :volume_quadrature,
    :physical_volume_quadrature, :pml_element_family, :pml_quadrature,
    :conductor_geometry_tolerance, :conductor_skin_depth_elements, :conductor_mesh_growth,
    :conductor_skin_depths, :conductor_thickness_elements)
const FEM_EXPORT_MESH_GLOBALS = (
    "Mesh.MshFileVersion", "Mesh.Binary", "Mesh.SaveAll", "Mesh.MeshSizeMin", "Mesh.MeshSizeMax",
    "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeExtendFromBoundary")
const FEM_EXPORT_SOLVER_OPTIONS = (:mumps_ordering, :petsc_prealloc)

_export_identifier(value) = replace(String(value), r"[^a-zA-Z0-9_]" => "_")
_export_material_name(index, material) = "Material_$(index)_$(_export_identifier(material.kind))"

function _write_export_data(path, model, formulation, controls)
    # Retain the shared input contract, with named native coefficients as its
    # only values. The arrays below are GetDP's frequency/region indexing.
    _write_model_data(path, model)
    original = readlines(path)
    material_lines = ("MaterialRegionTags(", "MaterialIsConductor(",
        "MaterialHasLoss(", "MaterialMu(", "MaterialSigma_", "MaterialEpsilon_")
    open(path, "w") do io
        println(io, "// Authoritative native model data. SI units; exp(+j omega t).")
        println(io, "// Edit these coefficients directly. Native execution never rewrites this file.")
        for line in original
            any(prefix -> startswith(line, prefix), material_lines) || println(io, line)
        end
        println(io, "Frequencies() = ", _pro_array(model.problem.frequencies), ";")
        println(io, "PmlSlopeValues() = ", _pro_array([plan.pml_slope for plan in model.mesh_plans]), ";")
        println(io,"FieldMapsFw() = Str[",
            join(_pro_string.(_field_quantities(formulation.options.data.physics)),", "),"];")
        println(io, "// Terminal encounter order and connections: positive phase IDs; 0 grounded.")
        for (i, name) in enumerate(model.terminal_ids)
            println(io, "// Terminal ", i, ": ", replace(name, '\n'=>' '))
            println(io, "// Terminal~{", i, "}: surface tag ", model.tags.terminal_base+i,
                "; TerminalContour~{", i, "}: curve tag ", model.tags.terminal_contour_base+i)
            println(io, "Connection_", i, " = ", model.problem.system.connection_order[i], ";")
        end
        println(io, "Connections() = {", join(["Connection_$i" for i in eachindex(model.terminal_ids)], ", "), "};")
        for (i, m) in enumerate(model.material_plans)
            name = _export_material_name(i, m)
            println(io, "\n// ", m.object_id, "; ", m.physical_name, "; selected law: ", m.field)
            println(io, name, "_Tag = ", m.physical_tag, ";")
            println(io, name, "_Mu = ", _pro_number(m.mu_r*4π*1e-7), "; // H/m")
            println(io, name, "_Sigma() = ", _pro_array(real.(m.admittivity)), "; // S/m, by frequency")
            println(io, name, "_Epsilon() = ", _pro_array(imag.(m.admittivity)./(2π.*model.problem.frequencies)), "; // F/m")
            println(io, "MaterialSigma_", i, "() = ", name, "_Sigma();")
            println(io, "MaterialEpsilon_", i, "() = ", name, "_Epsilon();")
            println(io, name, "_HasLoss = 0;")
            println(io, "For MaterialCase In {0:FrequencyCount-1}")
            println(io, "  If(",name,"_Sigma(MaterialCase) != 0) ",name,"_HasLoss = 1; EndIf")
            println(io, "EndFor")
        end
        names = [_export_material_name(i,m) for (i,m) in enumerate(model.material_plans)]
        println(io, "MaterialRegionTags() = {", join(names .* "_Tag", ", "), "};")
        println(io, "MaterialMu() = {", join(names .* "_Mu", ", "), "};")
        println(io, "MaterialIsConductor() = ", _pro_array([m.kind===:conductor for m in model.material_plans]), ";")
        println(io, "MaterialHasLoss() = {", join(names .* "_HasLoss", ", "), "};")
        for (name, values) in (
            ("DomainHalfwidths", getproperty.(model.mesh_plans, :domain_halfwidth)),
            ("VolumeQuadratures", getproperty.(model.mesh_plans, :volume_quadrature)))
            println(io, name, "() = ", _pro_array(values), ";")
        end
        for (i, side) in enumerate(("Side", "Top", "Bottom")), (suffix, field) in (("Thickness", :pml_thickness), ("Strength", :pml_strength))
            println(io, "Pml", side, suffix, "Values() = ", _pro_array([getproperty(p,field)[i] for p in model.mesh_plans]), ";")
        end
        println(io, "\n// Current amplitudes AND coefficient normalization used by the native equations.")
        println(io, "UnitSource = 1.; // axial A; constraints and normalized Z/P use this value")
        println(io, "LineLength = ", _pro_number(model.problem.system.line_length), "; // m")
        for (name, key) in (("ReduceBundle", :reduce_bundle), ("KronReduction", :kron_reduction), ("IdealTransposition", :ideal_transposition))
            println(io, name, " = ", Int(getproperty(formulation.options.data,key)), ";")
        end
        println(io, "DefineConstant[")
        frequency_choices = join(["$i=" * _pro_string(@sprintf("%.8g Hz", f))
            for (i,f) in enumerate(model.problem.frequencies)], ",")
        controls = (
            "RunFrequencyScan = {0, Choices{0,1}, Name \"Inputs/00Run frequency scan\"}",
            "ScanFrequencyIndex = {1, Choices{1:FrequencyCount}, Min 1, Max FrequencyCount, Step 1, Loop RunFrequencyScan, Visible 0, Name \"Inputs/00Scan frequency index\"}",
            "FrequencyIndex = {1, Choices{$frequency_choices}, Visible !RunFrequencyScan, Name \"Inputs/01Frequency case\"}",
            "Physics = {$(_fem_physics_code(formulation)), Choices{1=\"quasi-fw\"}, Name \"Inputs/02Physics\"}",
            "BasisTerminal = {0, Min 0, Max NumTerminals, Step 1, Name \"Inputs/03Basis (0 = full matrix)\"}",
            "MeshSizeFactor = {$(_pro_number(controls.mesh_size_factor)), Min 0.05, Name \"Mesh/01Physical size factor\"}",
            "ExteriorMeshSizeFactor = {$(_pro_number(controls.exterior_mesh_size_factor)), Min 1., Name \"Mesh/02Exterior size factor\"}",
            "InterfaceRefinementFactor = {$(_pro_number(model.interface_refinement_factor)), Min 1., Name \"Mesh/09Interface footprint factor\"}",
            "MeshRefinements = {0, Min 0, Max 4, Step 1, Name \"Mesh/08Uniform refinements\"}",
            "ConductorGeometryTolerance = {$(_pro_number(model.conductor_mesh.geometry_tolerance)), Name \"Mesh/03Conductor geometry tolerance\"}",
            "ConductorSkinDepthElements = {$(_pro_number(model.conductor_mesh.skin_depth_elements)), Name \"Mesh/04Conductor skin depth elements\"}",
            "ConductorMeshGrowth = {$(_pro_number(model.conductor_mesh.growth)), Name \"Mesh/05Conductor normal growth\"}",
            "ConductorSkinDepths = {$(_pro_number(model.conductor_mesh.skin_depths)), Name \"Mesh/06Conductor graded skin depths\"}",
            "ConductorThicknessElements = {$(model.conductor_mesh.thickness_elements), Min 1, Step 1, Name \"Mesh/07Conductor wall elements\"}",
            "RunAction = {1, Choices{0=\"Mesh only\",1=\"Mesh and solve\"}, Name \"Inputs/04Run action\"}",
            "PlotFieldMaps = {1, Choices{0,1}, Name \"Outputs/01Write field maps\"}",
            "ReuseFactorization = {1, Choices{0,1}, Name \"Numerics/01Reuse factorization\"}",
            "GetDPThreads = {1, Min 1, Max 64, Step 1, Name \"Numerics/02GetDP threads\"}",
            "MumpsOrdering = {$(something(controls.mumps_ordering,-1)), Choices{-1=\"solver default\",0=\"AMD\",2=\"AMF\",3=\"Scotch\",4=\"PORD\",5=\"METIS\",6=\"QAMD\",7=\"automatic\"}, Name \"Numerics/03MUMPS ordering\"}",
            "PetscPrealloc = {$(something(controls.petsc_prealloc,0)), Min 0, Step 1, Name \"Numerics/04Sparse row allocation (0 = default)\"}",
            "VolumeQuadrature = {$(controls.volume_quadrature), Choices{4,7,12,13}, Name \"Numerics/05Triangle volume points\"}",
            "PhysicalVolumeQuadrature = {$(something(controls.physical_volume_quadrature,0)), Choices{0=\"inherit triangle rule\",3=\"3\",4=\"4\",7=\"7\",12=\"12\",13=\"13\"}, Name \"Numerics/06Physical triangle points\"}",
            "PmlQuadrangles = {$(Int(controls.pml_element_family === :quadrangle)), Choices{0=\"triangles\",1=\"quadrangles\"}, Name \"Mesh/10PML element family\"}",
            "PmlQuadrature = {$(controls.pml_quadrature), Choices{4,9,16}, Name \"Numerics/07PML quadrangle points\"}",
            "OutputTotal = {0, Choices{0=\"per metre\",1=\"total length\"}, Name \"Outputs/02Matrix basis\"}")
        for (i,line) in enumerate(controls)
            println(io, "  ", line, i < length(controls) ? "," : "")
        end
        println(io, "];")
        println(io, "// ONELAB loops the complete Gmsh -> GetDP action over the hidden index.")
        println(io, "// Keep the manual selection intact, including after the loop resets its index.")
        println(io, "If(RunFrequencyScan && !StrCmp(OnelabAction, \"compute\"))")
        println(io, "  FrequencyIndex = ScanFrequencyIndex;")
        println(io, "EndIf")
    end
end

function _write_export_geometry(path, model, geometry, plan, controls)
    groups = [(dim, tag, gmsh.model.get_physical_name(dim,tag),
        gmsh.model.get_entities_for_physical_group(dim,tag)) for (dim,tag) in gmsh.model.get_physical_groups()]
    paths = Dict(curve => i for i in eachindex(model.terminal_ids)
        for curve in gmsh.model.get_entities_for_physical_group(1,model.tags.voltage_path_base+i))
    # Gmsh's native writer omits numeric physical tags. Emit these explicitly
    # from the same owned model, preserving overlaps and native primitives.
    gmsh.model.remove_physical_groups()
    raw = path * "_unrolled"
    gmsh.write(raw)
    labels = Dict{Int,Vector{String}}()
    for (dim, _, name, entities) in groups
        dim == 2 || continue
        for entity in entities
            push!(get!(labels,entity,String[]),name)
        end
    end
    open(path,"w") do io
        println(io, "// Fixed native geometry for ", @sprintf("%.8g",plan.frequency), " Hz; coordinates in m.")
        println(io, "// Points -> curves -> loops -> surfaces; memberships are explicit below.")
        interface_sizes = _write_export_mesh_sizes(io,model,plan,controls)
        conductor_fields = Set(Iterators.flatten(values(geometry.conductor_fields)))
        for line in eachline(raw)
            # The native writer drops the triangle arrangement and rounds
            # grading coefficients. The geometry owner retains both exactly.
            startswith(line,"Transfinite Curve") && continue
            startswith(line,"Transfinite Surface") && continue
            startswith(line,"Recombine Surface") && continue
            # Empty field lists are defaults; older native parsers reject the
            # empty braces emitted by Gmsh 4.15's serializer.
            occursin(r"^Field\[\d+\]\.\w+ = \{\};$",line) && continue
            # The native serializer rounds geometry numbers. Keep the native
            # topology it writes, but emit owned point/field values at Float64
            # round-trip precision through the public API.
            occursin(r"^cl__\d+ =",line) && continue
            point = match(r"^Point\((\d+)\)",line)
            if point !== nothing
                tag = parse(Int,point[1])
                coordinates = gmsh.model.get_value(0,tag,Float64[])
                size = only(gmsh.model.mesh.get_sizes([(0,tag)]))
                expression = coordinates[2] == 0 ? get(interface_sizes,
                    _coordinate_key(coordinates[1]),"MeshScale*$(_pro_number(size))") :
                    "MeshScale*$(_pro_number(size))"
                println(io,"Point(",tag,") = {",join(_pro_number.(coordinates),", "),
                    ", ",expression,"};")
                continue
            end
            # Rebuild surrounding-medium fields from the editable size controls.
            # Keep conductor IDs so their native material-dependent constraints
            # retain the same ownership as the Julia-managed construction.
            definition = match(r"^Field\[(\d+)\]",line)
            definition !== nothing && !(parse(Int,definition[1]) in conductor_fields) && continue
            startswith(line,"Background Field") && continue
            field = match(r"^Field\[(\d+)\]\.(\w+) = (.*);$",line)
            if field !== nothing
                tag = parse(Int,field[1])
                value = startswith(field[3],"\"") ?
                    _pro_string(gmsh.model.mesh.field.get_string(tag,String(field[2]))) :
                    startswith(field[3],"{") ?
                    _pro_array(gmsh.model.mesh.field.get_numbers(tag,String(field[2]))) :
                    _pro_number(gmsh.model.mesh.field.get_number(tag,String(field[2])))
                println(io,"Field[",tag,"].",field[2]," = ",value,";")
                continue
            end
            matched = match(r"^Plane Surface\((\d+)\)",line)
            if matched !== nothing
                for label in get(labels,parse(Int,matched[1]),String[])
                    println(io,"// ",label)
                end
            end
            println(io,line)
        end
        for (tag,(count,ratio)) in sort!(collect(geometry.transfinite_curves);by=first)
            println(io,"Transfinite Curve {",tag,"} = ",count," Using Progression ",_pro_number(ratio),";")
        end
        for (tag,(arrangement,corners)) in sort!(collect(geometry.transfinite_surfaces);by=first)
            println(io,"Transfinite Surface {",tag,"} = {",join(corners,", "),"} ",arrangement,";")
        end
        println(io,"If(PmlQuadrangles) Recombine Surface {",join(sort(geometry.pml_surfaces),", "),"}; EndIf")
        println(io,"\n// Physical memberships: these exact tags are consumed by GetDP.")
        for (dim,tag,name,entities) in groups
            println(io, "Physical ", ("Point", "Curve", "Surface")[dim+1], "(", _pro_string(name), ", ", tag, ") = {", join(entities,", "), "};")
        end
        println(io,"Mesh.MshFileVersion = 4.1; Mesh.Binary = 1; Mesh.SaveAll = 1;")
        println(io,"Mesh.MeshSizeMin = ",isempty(geometry.conductor_fields) ? "Min(MeshFine,Min(MeshWaveAir,MeshWaveEarth))" : "0", ";")
        println(io,"Mesh.MeshSizeMax = MeshRemoteMax;")
        println(io,"Mesh.MeshSizeFromPoints = 1; Mesh.MeshSizeExtendFromBoundary = 0;")
        println(io,"Mesh.MeshSizeFactor = 1; Mesh.ElementOrder = 1;")
        _write_export_mesh_fields(io,model,geometry,plan)
        _write_export_exterior_curves(io,model,geometry,plan,paths)
        _write_export_conductor_mesh(io,model,geometry,plan)
    end
    rm(raw)
end

function _write_export_mesh_sizes(io, model, plan, controls)
    println(io,"If(MeshSizeFactor <= 0 || ExteriorMeshSizeFactor < 1 || InterfaceRefinementFactor < 1)")
    println(io,"  Error(\"Physical mesh factor must be positive; exterior and interface factors must be at least one\"); EndIf")
    println(io,"MeshScale = MeshSizeFactor/",_pro_number(controls.mesh_size_factor),";")
    println(io,"MeshBulk = MeshScale*",_pro_number(plan.domain_mesh_size),";")
    println(io,"MeshFine = MeshScale*",_pro_number(model.fine_mesh_size),";")
    growth = _pro_number(model.mesh_growth_factor-1)
    println(io,"MeshRemote = Min(ExteriorMeshSizeFactor*MeshBulk,MeshBulk+",growth,"*",
        _pro_number(max(0,plan.domain_halfwidth-plan.exterior_start_radius)),");")
    for (i,medium) in enumerate(("Air","Earth"))
        println(io,"MeshWave",medium," = MeshScale*",_pro_number(plan.wave_mesh_sizes[i]),";")
    end
    air_limit = isfinite(plan.wave_size_limits[1]) ?
        "Min(MeshRemote,MeshScale*$(_pro_number(plan.wave_size_limits[1])))" : "MeshRemote"
    println(io,"MeshRemoteAir = ExteriorMeshSizeFactor == 1 ? MeshBulk : ",air_limit,";")
    println(io,"MeshRemoteEarth = MeshRemote;")
    println(io,"MeshRemoteMax = Max(MeshBulk,Max(MeshRemoteAir,MeshRemoteEarth));")
    println(io,"MeshInterface = MeshBulk;")
    sizes = Dict{Float64,String}()
    register(x,size) = (key = _coordinate_key(x);
        sizes[key] = haskey(sizes,key) ? "Min($(sizes[key]),$size)" : size)
    register(model.centre[1],"MeshInterface")
    for offset in (-2.,2.)
        register(model.centre[1]+offset,"MeshBulk")
    end
    for (i,(design,position)) in enumerate(zip(model.problem.system.designs,model.problem.system.positions))
        clearance = max(0,abs(position.y)-LineCableModels.outer_radius(design))
        println(io,"MeshCableInterface",i," = Min(MeshBulk,MeshScale*",
            _pro_number(model.cable_outer_mesh_sizes[i]),"+",growth,"*",_pro_number(clearance),");")
        println(io,"MeshInterface = Min(MeshInterface,MeshCableInterface",i,");")
        register(position.x,"MeshCableInterface$i")
    end
    for terminal in eachindex(model.terminal_ids)
        endpoint = argmin(p -> (p[2],p[1]),[_voltage_endpoint(region.shape)
            for region in model.region_plans if region.terminal_index == terminal])
        register(endpoint[1],"MeshInterface")
    end
    return sizes
end

function _write_export_mesh_fields(io, model, geometry, plan)
    println(io,"\n// Size fields use the same local targets, growth and medium restrictions as compute.")
    next_tag = maximum(Iterators.flatten(values(geometry.conductor_fields));init=0)
    function field(kind, properties...)
        tag = (next_tag += 1)
        println(io,"Field[",tag,"] = ",kind,";")
        for (name,value) in properties
            println(io,"Field[",tag,"].",name," = ",value,";")
        end
        return tag
    end
    background = Int[]
    interface_fields = Tuple{Int,Int,Int,String}[]
    growth = _pro_number(model.mesh_growth_factor-1)
    for (i,curves) in enumerate(geometry.cable_curves)
        isempty(curves) && continue
        distance = field("Distance", "CurvesList"=>_pro_array(sort!(unique(curves))),"Sampling"=>100)
        size = "MeshScale*$(_pro_number(model.cable_outer_mesh_sizes[i]))"
        push!(background,field("Threshold","InField"=>distance,"SizeMin"=>size,
            "SizeMax"=>"MeshRemoteMax","DistMin"=>0,
            "DistMax"=>"Max($size,(MeshRemoteMax-($size))/$growth)"))
    end
    # Explicit local cap is needed even if the export originally used factor 1:
    # raising the remote maximum later must not coarsen the cable interiors.
    bulk = field("MathEval","F"=>"Sprintf(\"%.17g\",MeshBulk)")
    push!(background,field("Restrict","InField"=>bulk,"IncludeBoundary"=>1,
        "SurfacesList"=>_pro_array(reduce(vcat,geometry.material_surfaces;init=Int[]))))
    for (i,(medium,surfaces)) in enumerate(zip(("Air","Earth"),
            (geometry.air_surfaces,geometry.earth_surfaces)))
        expression = "Min(%.17g,%.17g+$growth*Max(0,Sqrt((x-" *
            _pro_number(model.centre[1]) * ")^2+y^2)-$(_pro_number(plan.exterior_start_radius))))"
        radial = field("MathEval","F"=>"Sprintf(\"$expression\",MeshRemote$medium,MeshBulk)")
        push!(background,field("Restrict","InField"=>radial,"IncludeBoundary"=>1,
            "SurfacesList"=>_pro_array(surfaces)))
        curves = unique([geometry.interface_curves;reduce(vcat,geometry.cable_curves;init=Int[])])
        distance = field("Distance","CurvesList"=>_pro_array(curves),"Sampling"=>200)
        wave = field("Threshold","InField"=>distance,
            "SizeMin"=>"Min(MeshWave$medium,MeshRemote$medium)","SizeMax"=>"MeshRemote$medium",
            "DistMin"=>_pro_number(plan.wave_decay_radii[i]),
            "DistMax"=>_pro_number(2plan.wave_decay_radii[i]))
        push!(interface_fields,(i,distance,wave,medium))
        push!(background,field("Restrict","InField"=>wave,"IncludeBoundary"=>1,
            "SurfacesList"=>_pro_array(surfaces)))
    end
    for (index,(_,bulk)) in sort!(collect(geometry.conductor_fields);by=first)
        partition = get(geometry.sector_partitions,index,nothing)
        surfaces = partition === nothing ? geometry.region_surfaces[index] : [partition.core]
        push!(background,field("Restrict","InField"=>bulk,"SurfacesList"=>_pro_array(surfaces)))
    end
    for i in eachindex(model.cable_boundaries)
        surfaces = reduce(vcat,(geometry.region_surfaces[index] for (index,region) in enumerate(model.region_plans)
            if region.cable_index == i && model.material_plans[region.material_index].kind !== :conductor);init=Int[])
        isempty(surfaces) && continue
        size = "MeshScale*$(_pro_number(model.cable_outer_mesh_sizes[i]))"
        extension = field("Extend","CurvesList"=>_pro_array(_entity_boundary(surfaces)),
            "SizeMax"=>size,"DistMax"=>"($size)/$growth")
        push!(background,field("Restrict","InField"=>extension,"SurfacesList"=>_pro_array(surfaces)))
    end
    combined = field("Min","FieldsList"=>_pro_array(background))
    println(io,"Background Field = ",combined,";")
    # Leave constant fields and their existing distance sources untouched.
    # The condition follows native edits of the physical/exterior size factors.
    for (i,distance,wave,medium) in interface_fields
        condition = "MeshWave$medium < MeshRemote$medium"
        println(io,"If(",condition,")")
        println(io,"Field[",distance,"].CurvesList = ",
            _pro_array(unique(reduce(vcat,geometry.cable_curves;init=Int[]))),";")
        sources = [distance]
        for (design,position) in zip(model.problem.system.designs,model.problem.system.positions)
            x,width = _interface_footprint(design,position)
            expression = "Sqrt(y^2+Max(Abs(x-($(_pro_number(x))))-(%.17g),0)^2)"
            push!(sources,field("MathEval","F"=>
                "Sprintf(\"$expression\",InterfaceRefinementFactor*$(_pro_number(width)))"))
        end
        local_distance = field("Min","FieldsList"=>_pro_array(sources))
        println(io,"Field[",wave,"].InField = ",local_distance,";")
        println(io,"EndIf")
    end
end

function _write_export_exterior_curves(io, model, geometry, plan, paths)
    println(io,"\n// Recompute PML tangential and physical voltage-path divisions after a size edit.")
    # Normal PML strips retain their exported physical prescription. Their
    # tangential edges and the conforming measurement paths follow local sizes.
    halfwidth = plan.domain_halfwidth
    tolerance = 64eps(halfwidth)
    function graded(curve,length,first,last,reverse)
        println(io,"EdgeLast = ",last,"; EdgeFirst = Min(",first,",EdgeLast);")
        println(io,"If(EdgeFirst == EdgeLast)")
        println(io,"  EdgeCount = Max(2,Ceil(",_pro_number(length),"/EdgeLast)); EdgeRatio = 1;")
        println(io,"Else")
        println(io,"  EdgeGap = (EdgeLast-EdgeFirst)/EdgeFirst;")
        println(io,"  EdgeLog = EdgeGap < 1e-6 ? EdgeGap*(1-EdgeGap/2+EdgeGap^2/3) : Log(EdgeLast/EdgeFirst);")
        println(io,"  EdgeCount = Max(2,Ceil(",_pro_number(length),"*EdgeLog/(EdgeLast-EdgeFirst)));")
        println(io,"  EdgeRatio = Exp(EdgeLog/(EdgeCount-1));")
        println(io,"EndIf")
        println(io,"Transfinite Curve {",curve,"} = EdgeCount+1 Using Progression ",reverse ? "1/EdgeRatio" : "EdgeRatio",";")
    end
    for curve in sort!(collect(keys(geometry.transfinite_curves)))
        gmsh.model.get_type(1,curve) == "Line" || continue
        lower,upper = gmsh.model.get_parametrization_bounds(1,curve)
        a,b = gmsh.model.get_value(1,curve,lower),gmsh.model.get_value(1,curve,upper)
        medium = (a[2]+b[2])/2 >= 0 ? "Air" : "Earth"
        if abs(a[2]-b[2]) <= tolerance && abs(a[2]) >= halfwidth-tolerance &&
                max(abs(a[1]-model.centre[1]),abs(b[1]-model.centre[1])) <= halfwidth+tolerance
            println(io,"Transfinite Curve {",curve,"} = Max(2,Ceil(",_pro_number(abs(b[1]-a[1])),"/MeshRemote",medium,"))+1;")
        elseif abs(a[1]-b[1]) <= tolerance && abs(a[1]-model.centre[1]) >= halfwidth-tolerance &&
                max(abs(a[2]),abs(b[2])) <= halfwidth+tolerance
            first = "(MeshRemote$medium == MeshBulk ? MeshRemote$medium : Min(MeshBulk,MeshWave$medium))"
            graded(curve,abs(b[2]-a[2]),first,"MeshRemote$medium",abs(a[2])>abs(b[2]))
        elseif haskey(paths,curve) && min(a[2],b[2]) >= -halfwidth-tolerance &&
                max(a[2],b[2]) <= halfwidth+tolerance
            midpoint = (a.+b)./2
            physical = any(s -> gmsh.model.is_inside(2,s,midpoint)>0,
                [geometry.air_surfaces;geometry.earth_surfaces])
            cable = model.problem.system.terminal_order[paths[curve]].cable
            last = physical ? (model.cable_hosts[cable] === :air ? "MeshInterface" : "MeshRemoteEarth") : "MeshFine"
            graded(curve,hypot(b[1]-a[1],b[2]-a[2]),"MeshFine",last,a[2]<b[2])
        end
    end
end

function _write_export_conductor_mesh(io, model, geometry, plan)
    isempty(geometry.conductor_fields) && return nothing
    println(io,"\n// Native conductor constraints; prescribed sizes, no solve/refinement loop.")
    println(io,"If(ConductorGeometryTolerance <= 0 || ConductorSkinDepthElements <= 0 || ConductorMeshGrowth < 1 || ConductorSkinDepths <= 0 || ConductorThicknessElements < 1 || ConductorThicknessElements != Floor(ConductorThicknessElements))")
    println(io,"  Error(\"Conductor mesh controls must be positive; growth must be at least one\"); EndIf")
    println(io,"ConductorCircleSegments = 12*2^Max(0,Ceil(Log(Sqrt(2*Pi^2/(3*ConductorGeometryTolerance))/12)/Log(2)));")
    for (index,(layer,bulk)) in sort!(collect(geometry.conductor_fields);by=first)
        region = model.region_plans[index]
        section = _conductor_section(region.shape)
        material = region.material_index
        prefix = "ConductorRegion$index"
        assign(name,value) = println(io,prefix,name," = ",value,";")
        assign("Omega","2*Pi*Frequencies(FrequencyIndex-1)")
        assign("Mu","MaterialMu($(material-1))")
        assign("Sigma","MaterialSigma_$material(FrequencyIndex-1)")
        assign("B","$(prefix)Omega*MaterialEpsilon_$material(FrequencyIndex-1)")
        assign("H","Sqrt($(prefix)Sigma^2+$(prefix)B^2)")
        println(io,"If($(prefix)B >= 0 && $(prefix)H+$(prefix)B > 0)")
        assign("Decay","Abs($(prefix)Sigma)*Sqrt($(prefix)Omega*$(prefix)Mu/(2*($(prefix)H+$(prefix)B)))")
        println(io,"Else")
        assign("Decay","Sqrt($(prefix)Omega*$(prefix)Mu*($(prefix)H-$(prefix)B)/2)")
        println(io,"EndIf")
        assign("Phase","Sqrt($(prefix)Omega*$(prefix)Mu*($(prefix)H+$(prefix)B)/2)")
        assign("Delta","1e300")
        println(io,"If($(prefix)Decay > 0)")
        assign("Delta","1/$(prefix)Decay")
        println(io,"EndIf")
        width = _pro_number(section.width)
        divisions = iszero(section.divisions) ? "ConductorThicknessElements" : string(section.divisions)
        assign("Cap","$width/$divisions")
        assign("Bulk","Min(MeshScale*$(_pro_number(region.mesh_size)),$(prefix)Cap)")
        println(io,"If($(prefix)Decay*$width <= ConductorSkinDepths && $(prefix)Phase > 0)")
        assign("Bulk","Min($(prefix)Bulk,2*Pi/(12*$(prefix)Phase))")
        println(io,"EndIf")
        assign("Extent","Min(ConductorSkinDepths*$(prefix)Delta,$(_pro_number(section.fraction*section.width)))")
        assign("First","Min($(prefix)Delta/ConductorSkinDepthElements,$(prefix)Cap)")
        partition = get(geometry.sector_partitions,index,nothing)
        if partition !== nothing
            assign("Bulk","Min($(prefix)Bulk,$(_pro_number(section.width/30)))")
            assign("Segments","2*ConductorCircleSegments")
            depth = _pro_number(partition.depth)
            println(io,"If(ConductorMeshGrowth == 1)")
            assign("Layers","Ceil($depth/$(prefix)First)")
            println(io,"Else")
            assign("Layers","Ceil(Log(1+(ConductorMeshGrowth-1)*$depth/$(prefix)First)/Log(ConductorMeshGrowth))")
            println(io,"EndIf")
            for (outer,inner,length,turn) in partition.curve_pairs
                angular = "$(_pro_number(turn))*$(prefix)Segments/(2*Pi)"
                tangential = "$(_pro_number(length))/(($width/5)*96/$(prefix)Segments)"
                println(io,"Transfinite Curve {$outer,$inner} = Max(2,Max(Ceil($angular),Ceil($tangential)))+1;")
            end
            for spoke in partition.spokes
                ratio = spoke > 0 ? "ConductorMeshGrowth" : "1/ConductorMeshGrowth"
                println(io,"Transfinite Curve {$(abs(spoke))} = $(prefix)Layers+1 Using Progression $ratio;")
            end
            println(io,"Field[$bulk].F = Sprintf(\"%.17g\",$(prefix)Bulk);")
            continue
        end
        println(io,"If(ConductorMeshGrowth == 1)")
        assign("Layers","Ceil($(prefix)Extent/$(prefix)First)")
        assign("First","$(prefix)Extent/$(prefix)Layers")
        println(io,"Else")
        assign("Layers","Ceil(Log(1+(ConductorMeshGrowth-1)*$(prefix)Extent/$(prefix)First)/Log(ConductorMeshGrowth))")
        assign("First","$(prefix)Extent*(ConductorMeshGrowth-1)/(ConductorMeshGrowth^$(prefix)Layers-1)")
        println(io,"EndIf")
        for curve in _entity_boundary(geometry.region_surfaces[index])
            arc = _conductor_curve_geometry(region.shape,curve)
            angular = "Ceil(ConductorCircleSegments*$(_pro_number(arc.fraction))-1e-12)+1"
            local_size = "Ceil($(_pro_number(arc.length))/(MeshScale*$(_pro_number(region.mesh_size))))+1"
            println(io,"Transfinite Curve {",curve,"} = Max(2,Max($angular,$local_size));")
        end
        println(io,"Field[$layer].Size = $(prefix)First; Field[$layer].Ratio = ConductorMeshGrowth;")
        println(io,"Field[$layer].Thickness = $(prefix)Extent*(1+1e-8);")
        println(io,"Field[$bulk].F = Sprintf(\"%.17g\",$(prefix)Bulk);")
        all_surfaces = last.(gmsh.model.get_entities(2))
        excluded = setdiff(all_surfaces,geometry.region_surfaces[index])
        println(io,"If($(prefix)Delta < $width)")
        println(io,"Field[$layer].ExcludedSurfacesList = ",_pro_array(excluded),";")
        println(io,"BoundaryLayer Field = $layer;")
        println(io,"Else")
        println(io,"Field[$layer].ExcludedSurfacesList = ",_pro_array(all_surfaces),";")
        println(io,"EndIf")
    end
    return nothing
end

function _write_export_bundle(root, stem, model, formulation, controls)
    mkpath(joinpath(root,"geometry"))
    mkpath(joinpath(root,"formulations"))
    for (key,path) in pairs(_getdp_assets(joinpath(root,"formulations")))
        write(path,getproperty(FEM_GETDP_SOURCES,key))
    end
    _write_export_data(joinpath(root,stem*"_data.pro"),model,formulation,controls)
    lock(FEM_SESSION_LOCK) do
        session = _start_gmsh(0)
        saved = Dict(option=>gmsh.option.get_number(option) for option in FEM_EXPORT_MESH_GLOBALS)
        try
            for plan in model.mesh_plans
                geo = _build_geometry!(model,"export-$(basename(root))-$(plan.frequency_index)",plan)
                _configure_mesh!(model,geo,plan)
                _write_export_geometry(joinpath(root,"geometry",@sprintf("case-%04d.geo",plan.frequency_index)),model,geo,plan,controls)
                gmsh.model.remove()
            end
        finally
            for (option,value) in saved
                gmsh.option.set_number(option,value)
            end
            _finish_gmsh(session)
        end
    end
    write(joinpath(root,stem*".geo"), """
    Include "$(stem)_data.pro";
    Include Sprintf("geometry/case-%04g.geo", FrequencyIndex);
    // Remove only this project's derived views when inputs are checked/rerun.
    If(!StrCmp(OnelabAction, "check") || !StrCmp(OnelabAction, "compute"))
      Include "views.geo";
    EndIf
    // ONELAB Run and explicit CLI BuildMesh use the same native mesh operations.
    Solver.AutoMesh = -1;
    If(!StrCmp(OnelabAction, "compute") || Exists(BuildMesh))
      Mesh 2;
      If(MeshRefinements > 0)
        For refinement In {1:MeshRefinements}
          RefineMesh;
        EndFor
      EndIf
      Save StrCat(CurrentDirectory, "$(stem).msh");
    EndIf
    """)
    write(joinpath(root,stem*".pro"), """
    // Open this file in Gmsh/ONELAB. All numerical operations use native GetDP.
    // Geometry: $(stem).geo; coefficients/connections: $(stem)_data.pro.
    // Regions, constraints and equations: formulations/quasi-*.pro.
    ProjectDirectory = CurrentDirectory;
    ModelDataPath = StrCat[ProjectDirectory,"$(stem)_data.pro"];
    Include ModelDataPath;
    Include "formulations/onelab.pro";
    """)
    cp(joinpath(@__DIR__,"onelab_export","README.md"),joinpath(root,"README.md"))
    cp(joinpath(@__DIR__,"onelab_export","views.geo"),joinpath(root,"views.geo"))
    files = sort!([relpath(joinpath(dir,name),root) for (dir,_,names) in walkdir(root) for name in names])
    write(joinpath(root,".onelab-export-files"),join(files,"\n")*"\n")
end

"""
    export_data(:onelab, problem, formulation; file_name, mesh_options=(;), overwrite=false)

Export a detached Gmsh/ONELAB project with readable native geometry and GetDP
model files. Returns the absolute entry `.pro` path. Load `Gmsh` before calling.

`file_name` selects the entry file in a dedicated bundle directory.
`mesh_options` accepts `domain_skin_depths`, `pml_thickness`,
`pml_thickness_factor`, `pml_layers`, `pml_grading`, `pml_resolution`, `pml_reflection`,
`mesh_size_factor`,
`exterior_mesh_size_factor`, `interface_refinement_factor`, `volume_quadrature`,
`physical_volume_quadrature`, `pml_element_family`, and `pml_quadrature`, with their existing FEM
meanings and defaults.
`solver_options` accepts `mumps_ordering` and `petsc_prealloc`
with the same defaults as FEM computation options. These are exported as editable
ONELAB controls; ordering availability depends on the native solver build.
Geometry and evaluated material frequency cases are fixed at export.
The detached runtime uses only native Gmsh/GetDP with ONELAB; it requires neither
Julia nor Python. Export does not execute GetDP or open a GUI. `overwrite=true`
replaces recorded bundle files and removes obsolete owned files, preserving
unrelated files in the destination.
"""
function ImportExport.export_data(::Val{:onelab}, problem::LineParametersProblem,
        formulation::LineCableModelsFEM; file_name::AbstractString,
        mesh_options::NamedTuple=(;), solver_options::NamedTuple=(;), overwrite::Bool=false)
    unknown = setdiff(keys(mesh_options),FEM_EXPORT_MESH_OPTIONS)
    isempty(unknown) || throw(ArgumentError("Unsupported export mesh options: $(join(unknown, ", "))"))
    unknown_solver = setdiff(keys(solver_options),FEM_EXPORT_SOLVER_OPTIONS)
    isempty(unknown_solver) || throw(ArgumentError("Unsupported export solver options: $(join(unknown_solver, ", "))"))
    entry = abspath(file_name)
    stem, extension = splitext(basename(entry))
    extension == ".pro" && !isempty(stem) || throw(ArgumentError("file_name must end in .pro"))
    occursin(r"[\r\n\"]",entry) && throw(ArgumentError("export path cannot contain quotes or newlines"))
    root = dirname(entry)
    existing = isdir(root) ? readdir(root) : String[]
    marker = joinpath(root,".onelab-export-files")
    !isempty(existing) && (!overwrite || !isfile(marker)) && throw(ArgumentError("destination is not empty; overwrite requires an existing exported bundle"))
    opts = computation_options(LineCableModelsFEM,ComputationOptions(;mesh_options...,solver_options...))
    model = _resolved_fem_model(_preflight_fem_problem(problem),formulation,opts)
    mkpath(dirname(root))
    stage = mktempdir(dirname(root);prefix=".onelab-export-")
    try
        _write_export_bundle(stage,stem,model,formulation,opts.data)
        owned = isfile(marker) ? Set(readlines(marker)) : Set{String}()
        for file in owned
            normalized = normpath(file)
            (isabspath(file) || normalized == ".." || startswith(normalized,".." * string(Base.Filesystem.path_separator)) || normalized == ".") &&
                throw(ArgumentError("invalid path in export ownership manifest: $file"))
            parent = dirname(joinpath(root,normalized))
            while parent != root
                islink(parent) && throw(ArgumentError("export ownership path traverses a symbolic link: $file"))
                parent = dirname(parent)
            end
        end
        files = [readlines(joinpath(stage,".onelab-export-files")); ".onelab-export-files"]

        for file in files
            target = joinpath(root,file)
            ispath(target) && !(file in owned || file==".onelab-export-files") && throw(ArgumentError("refusing to replace unrelated destination file $target"))
        end
        mkpath(root)
        for file in setdiff(owned,Set(files))
            target = joinpath(root,file)
            (isfile(target) || islink(target)) && rm(target)
        end
        for file in files
            mkpath(dirname(joinpath(root,file)))
            cp(joinpath(stage,file),joinpath(root,file);force=true)
        end
    finally
        rm(stage;recursive=true,force=true)
    end
    return entry
end

function ImportExport.export_data(format::Val{:onelab}, system::DataModel.LineCableSystem,
        formulation::LineCableModelsFEM; earth_props, frequencies,
        temperature=20.0, kwargs...)
    problem = LineParametersProblem(system;earth_props,frequencies,temperature)
    return ImportExport.export_data(format,problem,formulation;kwargs...)
end
