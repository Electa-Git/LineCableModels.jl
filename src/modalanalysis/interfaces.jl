"""
$(TYPEDEF)

Select the LineCableModels modal-transformation backend.
"""
struct LineCableModelsModal end

"""
$(TYPEDEF)

Store frequency-dependent modal-to-phase voltage and current bases.

For every frequency sample, `Tv` and `Ti` satisfy `Vₚ = Tv Vₘ` and
`Iₚ = Ti Iₘ`.

$(TYPEDFIELDS)
"""
struct ModalOperators{V <: AbstractArray, I <: AbstractArray}
    "Modal-to-phase voltage tensor."
    Tv::V
    "Modal-to-phase current tensor."
    Ti::I

    function ModalOperators(Tv::V, Ti::I) where {
            V <: AbstractArray,
            I <: AbstractArray
    }
        ndims(Tv) == 3 || throw(
            DimensionMismatch("voltage operators must be an n×n×nfreq tensor")
        )
        ndims(Ti) == 3 || throw(
            DimensionMismatch("current operators must be an n×n×nfreq tensor")
        )
        size(Tv) == size(Ti) || throw(DimensionMismatch(
            "voltage and current operators must have equal n×n×nfreq dimensions"
        ))
        size(Tv, 1) == size(Tv, 2) || throw(
            DimensionMismatch("modal operators must be square")
        )
        return new{V, I}(Tv, Ti)
    end
end

Engine.validate_modal_operators(maps::ModalOperators) = size(maps.Tv)

"""
Return the modal operators carried by modal-domain line parameters.
"""
function operators(parameters::LineParameters{T, U, D}) where {T, U, D <: ModalDomain}
    parameters.domain.operators
end

function selectdomain(domain::ModalDomain, selected)
    maps = domain.operators
    Tv = Array(view(maps.Tv, :, :, selected))
    Ti = Array(view(maps.Ti, :, :, selected))
    selected_maps = ModalOperators(Tv, Ti)
    roots = Array(view(domain.gamma, :, selected))
    return ModalDomain(selected_maps, roots)
end

function selectdetails(retained::ComputationDetails, ::ModalDomain, selected)
    record=retained.data
    haskey(record,:modal) || return retained
    modal=record.modal
    haskey(modal,:diagnostics) || return retained
    diagnostics=modal.diagnostics
    original=selected isa Colon ? collect(eachindex(diagnostics.z_coupling)) : collect(selected)
    fallback=findall(in(diagnostics.fallback_frequencies),original)
    missed=findall(in(diagnostics.missed_frequencies),original)
    sliced=merge(diagnostics,(fallback_frequencies=fallback,
        missed_frequencies=missed,
        z_coupling=diagnostics.z_coupling[selected],
        y_coupling=diagnostics.y_coupling[selected],
        eigen_residual=diagnostics.eigen_residual[:,selected],
        iterations=diagnostics.iterations[:,selected],
        converged=diagnostics.converged[:,selected]))
    return ComputationDetails(merge(record,(modal=merge(modal,(diagnostics=sliced,)),)))
end
