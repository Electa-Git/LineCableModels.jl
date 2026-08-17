"""
    LineCableModels.Engine.Transforms

# Dependencies

$(IMPORTS)

"""
module Transforms

# Export public API
export Fortescue

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ...LineCableModels: description, PhaseDomain, ModalDomain, basis, nominal
import ..Engine:
                 AbstractTransformFormulation, LineParameters, SeriesImpedance,
                 ShuntAdmittance

#
using LinearAlgebra
# using GenericLinearAlgebra
using NLsolve

reciprocity(matrix) = (matrix + transpose(matrix)) / 2
reciprocity!(matrix) = (matrix .= (matrix .+ transpose(matrix)) ./ 2; matrix)

function ideal_transposition!(matrix::AbstractMatrix)
    n = LinearAlgebra.checksquare(matrix)
    coefficients = similar(diag(matrix))
    @inbounds for offset in 0:(n - 1)
        total = zero(eltype(matrix))
        for row in 1:n
            total += matrix[row, 1 + mod(row - 1 + offset, n)]
        end
        coefficients[offset + 1] = total / n
    end
    @inbounds for row in 1:n, column in 1:n

        matrix[row, column] = coefficients[mod1(column - row + 1, n)]
    end
    return matrix
end

function offdiagonal_ratio(matrix::AbstractMatrix)
    n = LinearAlgebra.checksquare(matrix)
    diagonal_max = maximum(abs(matrix[index, index]) for index in 1:n)
    offdiagonal_max = maximum(
        (abs(matrix[row, column]) for column in 1:n for row in 1:n if row != column);
        init = zero(diagonal_max)
    )
    return offdiagonal_max / max(diagonal_max, eps(float(real(diagonal_max))))
end

include("fortescue.jl")
include("eiglevenberg.jl")

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

end # module Transforms
