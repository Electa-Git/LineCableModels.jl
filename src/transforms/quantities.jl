"""
$(TYPEDSIGNATURES)

Return the diagonal modal propagation-constant tensor.

The calculation uses the completed modal `YₘZₘ` product, so it is
independent of the decomposition formula's internal operator convention.
"""
function gamma(
        parameters::LineParameters{
        T, U, D}
) where {T <: Complex, U <: Real, D <: ModalDomain}
    impedance = parameters.Z.values
    admittance = parameters.Y.values
    n, _, nfrequencies = size(impedance)
    values = similar(impedance)
    product = Matrix{T}(undef, n, n)
    fill!(values, zero(T))
    @inbounds for frequency in 1:nfrequencies
        mul!(product,
            @view(admittance[:, :, frequency]),
            @view(impedance[:, :, frequency]))
        for mode in 1:n
            values[mode, mode, frequency] = sqrt(product[mode, mode])
        end
    end
    return values
end

"""
$(TYPEDSIGNATURES)

Return modal and phase-domain characteristic quantities.

The returned tuple contains modal `Z`, modal `Y`, modal characteristic
impedance, modal characteristic admittance, phase characteristic impedance,
and phase characteristic admittance. All values use the complete
frequency-dependent operators carried by the modal domain.
"""
function modal_quantities(
        parameters::LineParameters{
        T, U, D}
) where {T <: Complex, U <: Real, D <: ModalDomain}
    Zm = copy(parameters.Z.values)
    Ym = copy(parameters.Y.values)
    n, _, nfrequencies = size(Zm)
    Zc = zeros(T, n, n, nfrequencies)
    Yc = similar(Zc)
    Zp = similar(Zc)
    Yp = similar(Zc)
    maps = operators(parameters)

    @inbounds for frequency in 1:nfrequencies
        for mode in 1:n
            value = sqrt(Zm[mode, mode, frequency]) /
                    sqrt(Ym[mode, mode, frequency])
            Zc[mode, mode, frequency] = value
            Yc[mode, mode, frequency] = inv(value)
        end
        A = @view maps.voltage[:, :, frequency]
        B = @view maps.current[:, :, frequency]
        @views Zp[:, :, frequency] .= A \ (Zc[:, :, frequency] * B)
        @views Yp[:, :, frequency] .= B \ (Yc[:, :, frequency] * A)
    end
    return Zm, Ym, Zc, Yc, Zp, Yp
end
