function (F::AbstractTransformFormulation)(
        lp::LineParameters{
        Tc, U, ModalDomain, Basis},
) where {Tc <: Complex, U <: Real, Basis}
    throw(
        ErrorException(
        "Not yet implemented: inverse $(nameof(typeof(F)))( ::LineParameters{<:Complex,<:Real,ModalDomain} )",
    ),
    )
end
