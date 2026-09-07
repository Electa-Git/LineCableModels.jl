_physical_surface_state(input,formula,insulation,s,temperature)=nothing

function _physical_surface_state(input::LocalCableData{T},
        formula::InternalImpedance.Formula{:Ametani1992},insulation,s,temperature) where {T}
    all(indices->length(indices)==1,input.assemblies) ||
        throw(ArgumentError("Ametani1992 supplies an isolated scalar impedance, not annular transfer terms"))
    Surface=NamedTuple{(:inner,:outer,:mutual),NTuple{3,Complex{T}}}
    surfaces=Vector{Surface}(undef,length(input.terminals))
    gap_terms=zeros(Complex{T},length(input.terminals))
    for index in eachindex(input.terminals)
        section=input.sections[index]
        section===nothing && throw(ArgumentError(
            "Ametani1992 requires one physical homogeneous cross-section, not a strand or composite reduction"))
        material=section.material
        rho=temperature===nothing ? material.rho :
            material.rho*(one(T)+material.alpha*(temperature-material.T0))
        interaction=formula(Val(:section),section.area,section.perimeter,
            rho,material.mu_r,s)
        surfaces[index]=(inner=zero(s),outer=interaction(Val(:outer)),mutual=zero(s))
        gap_terms[index]=insulation(input.r_ext[index],input.r_ins_ext[index],
            input.mu_ins[index],s)
    end
    return (;surfaces,gap_terms)
end

function _physical_surface_state(input::LocalCableData{T},
        formula::InternalImpedance.Formula{:Ametani2004},insulation,s,temperature) where {T}
    profiles=[copy(layers) for layers in input.conductor_layers]
    all(!isempty,profiles) || throw(ArgumentError(
        "bonded surface impedance requires physical concentric circular metal layers"))
    gap_terms=zeros(Complex{T},length(profiles))
    # Series attachment does not alter the separate radial shunt-layer records.
    for assembly in input.assemblies
        for (position,index) in pairs(assembly)
            interval=input.dielectric_ranges[index]
            low=first(interval); high=last(interval)
            next=position<length(assembly) ? assembly[position+1] : nothing
            if !isempty(interval) && next!==nothing &&
                    all(k->input.dielectric_materials[k].kind===:semicon,interval)
                throw(ArgumentError(
                    "distinct retained conductors must have an insulating separation"))
            end
            while low<=high && input.dielectric_materials[low].kind===:semicon
                push!(profiles[index],BlueprintConductorLayer{T}(
                    input.r_layer_in[low],input.r_layer_ext[low],
                    input.dielectric_materials[low]))
                low+=1
            end
            while low<=high && input.dielectric_materials[high].kind===:semicon
                next===nothing && throw(ArgumentError(
                    "a semiconductor behind insulation needs a bonded outer conductor"))
                pushfirst!(profiles[next],BlueprintConductorLayer{T}(
                    input.r_layer_in[high],input.r_layer_ext[high],
                    input.dielectric_materials[high]))
                high-=1
            end
            for layer in low:high
                material=input.dielectric_materials[layer]
                material.kind===:insulator || throw(ArgumentError(
                    "a semiconductor separated from both metals is outside the bonded model"))
                gap_terms[index]+=insulation(input.r_layer_in[layer],
                    input.r_layer_ext[layer],material.mu_r,s)
            end
        end
    end
    Surface=NamedTuple{(:inner,:outer,:mutual),NTuple{3,Complex{T}}}
    surfaces=Vector{Surface}(undef,length(profiles))
    for index in eachindex(profiles)
        layers=map(profiles[index]) do layer
            material=layer.material
            rho=material.rho
            if temperature!==nothing
                rho*=one(T)+material.alpha*(temperature-material.T0)
            end
            (r_in=layer.r_in,r_ex=layer.r_ex,rho=rho,mu_r=material.mu_r)
        end
        functor=formula(layers,s)
        surfaces[index]=(inner=functor(Val(:inner)),outer=functor(Val(:outer)),
            mutual=functor(Val(:mutual)))
    end
    return (;surfaces,gap_terms)
end

function _physical_surface_state(input::LocalCableData{T},
        formula::InternalImpedance.Formula{:Merkushev2015},insulation,s,temperature) where {T}
    all(indices->length(indices)==1,input.assemblies) ||
        throw(ArgumentError("Merkushev supplies only an isolated scalar conductor impedance"))
    Surface=NamedTuple{(:inner,:outer,:mutual),NTuple{3,Complex{T}}}
    surfaces=Vector{Surface}(undef,length(input.terminals))
    gap_terms=zeros(Complex{T},length(input.terminals))
    for index in eachindex(input.terminals)
        profile=input.acsr[index]
        profile===nothing && throw(ArgumentError(
            "Merkushev requires one central wire and six equal-radius, equally spaced strands with one common pitch"))
        resistivity(material)=temperature===nothing ? material.rho :
            material.rho*(one(T)+material.alpha*(temperature-material.T0))
        interaction=formula(profile.radius,profile.pitch,resistivity(profile.core),
            resistivity(profile.strands),profile.core.mu_r,s)
        surfaces[index]=(inner=interaction(Val(:inner)),outer=interaction(Val(:outer)),
            mutual=interaction(Val(:mutual)))
        gap_terms[index]=insulation(input.r_ext[index],input.r_ins_ext[index],
            input.mu_ins[index],s)
    end
    return (;surfaces,gap_terms)
end
