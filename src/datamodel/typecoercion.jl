import ..Utils: coerce_to_T
using ..Validation: required_fields, keyword_fields, coercive_fields

@inline function _rebuild_part_typed_core(p, ::Type{T}) where {T}
    C0 = typeof(p).name.wrapper                       # concrete parametric type (e.g., WireArray)
    order = (required_fields(C0)..., keyword_fields(C0)...)  # positional order for tight kernel
    coer = coercive_fields(C0)                          # only these get coerced to T

    argsT = ntuple(i -> begin
            s = order[i]
            v = getfield(p, s)
            (s in coer) ? coerce_to_T(v, T) : v              # preserve Int/categorical fields
        end, length(order))

    return C0(argsT...)                                  # call the tight numeric constructor
end

# Identity when already at T (no rebuild, preserves ===)
coerce_to_T(p::AbstractConductorPart{T}, ::Type{T}) where {T} = p
# Cross-T rebuild via your existing tight numeric-core helper
coerce_to_T(p::AbstractConductorPart{S}, ::Type{T}) where {S,T} =
    _rebuild_part_typed_core(p, T)
coerce_to_T(g::ConductorGroup{T}, ::Type{T}) where {T} = g
# Cross-T: fieldwise coerce + layer coercion (no recompute)
@inline function coerce_to_T(g::ConductorGroup{S}, ::Type{T}) where {S,T}
    n = length(g.layers)
    layersT = Vector{AbstractConductorPart{T}}(undef, n)
    @inbounds for i in 1:n
        layersT[i] = coerce_to_T(g.layers[i], T)  # uses your part-level coercers
    end
    return ConductorGroup{T}(
        coerce_to_T(g.radius_in, T),
        coerce_to_T(g.radius_ext, T),
        coerce_to_T(g.cross_section, T),
        g.num_wires,                                  # keep Int as-is
        coerce_to_T(g.num_turns, T),
        coerce_to_T(g.resistance, T),
        coerce_to_T(g.alpha, T),
        coerce_to_T(g.gmr, T),
        layersT,
    )
end

@inline coerce_to_T(p::AbstractInsulatorPart{T}, ::Type{T}) where {T} = p
@inline coerce_to_T(p::AbstractInsulatorPart{S}, ::Type{T}) where {S,T} =
    _rebuild_part_typed_core(p, T)
@inline coerce_to_T(g::InsulatorGroup{T}, ::Type{T}) where {T} = g
@inline function coerce_to_T(g::InsulatorGroup{S}, ::Type{T}) where {S,T}
    n = length(g.layers)
    layersT = Vector{AbstractInsulatorPart{T}}(undef, n)
    @inbounds for i in 1:n
        layersT[i] = coerce_to_T(g.layers[i], T)   # uses the part-level coercers above
    end
    return InsulatorGroup{T}(
        coerce_to_T(g.radius_in, T),
        coerce_to_T(g.radius_ext, T),
        coerce_to_T(g.cross_section, T),
        coerce_to_T(g.shunt_capacitance, T),
        coerce_to_T(g.shunt_conductance, T),
        layersT,
    )
end

@inline coerce_to_T(c::CableComponent{T}, ::Type{T}) where {T} = c
@inline function coerce_to_T(c::CableComponent{S}, ::Type{T}) where {S,T}
    CableComponent{T}(
        c.id,
        coerce_to_T(c.conductor_group, T),
        coerce_to_T(c.insulator_group, T),
    )
end

"Identity: no allocation when already at `T`."
@inline coerce_to_T(n::NominalData{T}, ::Type{T}) where {T} = n
# Cross-T rebuild: fieldwise coercion, preserving `nothing`
@inline function coerce_to_T(n::NominalData{S}, ::Type{T}) where {S,T}
    names = fieldnames(typeof(n))                  # e.g. (:designation_code, :U0, :U, ...)
    vals = map(names) do k                        # map over tuple of names → returns a tuple
        v = getfield(n, k)
        v === nothing ? nothing : coerce_to_T(v, T)
    end
    NT = NamedTuple{names}(vals)                   # correct: pass a SINGLE tuple, not varargs
    return NominalData{T}(; NT...)                 # call typed kernel via keyword splat
end

@inline coerce_to_T(d::CableDesign{T}, ::Type{T}) where {T} = d
@inline function coerce_to_T(d::CableDesign{S}, ::Type{T}) where {S,T}
    compsT = Vector{CableComponent{T}}(undef, length(d.components))
    @inbounds for i in eachindex(d.components)
        compsT[i] = coerce_to_T(d.components[i], T)
    end
    ndT = isnothing(d.nominal_data) ? nothing : coerce_to_T(d.nominal_data, T)
    CableDesign{T}(d.cable_id, compsT; nominal_data=ndT)
end

@inline coerce_to_T(p::CablePosition{T}, ::Type{T}) where {T} = p
@inline function coerce_to_T(p::CablePosition{S}, ::Type{T}) where {S,T}
    CablePosition{T}(
        coerce_to_T(p.design_data, T),
        coerce_to_T(p.horz, T),
        coerce_to_T(p.vert, T),
        p.conn,                      # keep Int mapping as-is
    )
end

@inline coerce_to_T(sys::LineCableSystem{T}, ::Type{T}) where {T} = sys
@inline function coerce_to_T(sys::LineCableSystem{S}, ::Type{T}) where {S,T}
    cablesT = Vector{CablePosition{T}}(undef, length(sys.cables))
    @inbounds for i in eachindex(sys.cables)
        cablesT[i] = coerce_to_T(sys.cables[i], T)
    end
    # counts will be recomputed once positions are populated; preserve them now
    LineCableSystem{T}(sys.system_id, coerce_to_T(sys.line_length, T), cablesT)
end