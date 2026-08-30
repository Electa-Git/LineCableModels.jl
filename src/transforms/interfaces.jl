"""
$(TYPEDEF)

Select the LineCableModels modal-transformation backend.
"""
struct LineCableModelsModal end

"""
$(TYPEDEF)

Store frequency-dependent phase-to-modal voltage and current operators.

For every frequency sample, `voltage` stores `A` and `current` stores `B` in
`Vₘ = A Vₚ` and `Iₘ = B Iₚ`.

$(TYPEDFIELDS)
"""
struct ModalOperators{V <: AbstractArray, I <: AbstractArray}
    "Phase-to-modal voltage tensor."
    voltage::V
    "Phase-to-modal current tensor."
    current::I

    function ModalOperators(voltage::V, current::I) where {
            V <: AbstractArray,
            I <: AbstractArray
    }
        ndims(voltage) == 3 || throw(
            DimensionMismatch("voltage operators must be an n×n×nfreq tensor")
        )
        ndims(current) == 3 || throw(
            DimensionMismatch("current operators must be an n×n×nfreq tensor")
        )
        size(voltage) == size(current) || throw(DimensionMismatch(
            "voltage and current operators must have equal n×n×nfreq dimensions"
        ))
        size(voltage, 1) == size(voltage, 2) || throw(
            DimensionMismatch("modal operators must be square")
        )
        return new{V, I}(voltage, current)
    end
end

"Return the modal operators carried by modal-domain line parameters."
function operators(parameters::LineParameters{T, U, D}) where {T, U, D <: ModalDomain}
    parameters.domain.operators
end

function selectdomain(domain::ModalDomain, selected)
    maps = domain.operators
    voltage = Array(view(maps.voltage,:,:,selected))
    current = Array(view(maps.current,:,:,selected))
    return ModalDomain(ModalOperators(voltage, current), domain.formula)
end
