"""
$(TYPEDSIGNATURES)

Return the diagonal modal propagation tensor, with dimensions
`(mode, mode, frequency)` and zero off-diagonal entries:

```math
\\gamma_{ii}(f)=\\sqrt{[Y_m(f)Z_m(f)]_{ii}}.
```

`parameters` contains modal series impedance ``Z_m`` and shunt admittance
``Y_m``. With a `:pul` basis (Ω/m and S/m), the result is in m⁻¹. With a
`:total` basis (Ω and S), it is the dimensionless propagation product
``γℓ`` for line length ``ℓ``. Each entry uses Julia's principal complex
square root; this function performs no additional branch tracking.
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

Return `(Zm, Ym, Zc, Yc, Zp, Yp)`: modal series impedance and shunt admittance,
followed by modal and phase-domain characteristic impedance and admittance.
Every array has dimensions `(conductor_or_mode, conductor_or_mode, frequency)`.

For each mode and frequency, calculate the diagonal characteristic quantities
using Julia's principal complex square roots:

```math
Z_{c,ii}=\\frac{\\sqrt{Z_{m,ii}}}{\\sqrt{Y_{m,ii}}},\\qquad
Y_{c,ii}=Z_{c,ii}^{-1}.
```

The voltage and current operators stored in `parameters` map phase quantities
to modal quantities, ``V_m=AV_p`` and ``I_m=BI_p``. Thus

```math
Z_p=A^{-1}Z_cB,\\qquad Y_p=B^{-1}Y_cA.
```

`Zm` and `Ym` retain the input basis: Ω/m and S/m for `:pul`, or Ω and S
for `:total`. Characteristic impedances `Zc` and `Zp` are in Ω;
characteristic admittances `Yc` and `Yp` are in S for either basis.
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
