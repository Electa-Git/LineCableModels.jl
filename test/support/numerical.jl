@testmodule TestNumerics begin
    using LinearAlgebra

    precision_type(::Type{T}) where {T <: AbstractFloat} = T
    precision_type(::Type{Complex{T}}) where {T <: AbstractFloat} = T
    precision_type(::Type{<:Integer}) = Float64
    precision_type(::Type{<:Rational}) = Float64
    precision_type(value::Number) = precision_type(typeof(value))
    precision_type(values::AbstractArray) = precision_type(eltype(values))

    magnitude(value::Number) = abs(value)
    magnitude(values::AbstractArray) = norm(values)

    absolute_floor(::Type{T}; scale = one(T)) where {T <: AbstractFloat} = sqrt(eps(T)) *
                                                                           abs(scale)

    function tolerances(actual, expected; factor = 64, scale = nothing)
        T = promote_type(precision_type(actual), precision_type(expected))
        reference_scale = if isnothing(scale)
            max(magnitude(actual), magnitude(expected), floatmin(T))
        else
            convert(T, scale)
        end
        return (; rtol = factor * eps(T), atol = factor * eps(T) * reference_scale)
    end

    function isapprox_scaled(actual, expected; factor = 64, scale = nothing)
        tolerance = tolerances(actual, expected; factor, scale)
        return isapprox(actual, expected; tolerance...)
    end

    function normalized_residual(A, x, b)
        numerator = norm(A * x - b)
        denominator = max(norm(A) * norm(x) + norm(b), floatmin(Float64))
        return numerator / denominator
    end

    function normalized_eigen_residual(A, eigenvalue, eigenvector)
        numerator = norm(A * eigenvector - eigenvalue * eigenvector)
        denominator = max(
            norm(A) * norm(eigenvector) + abs(eigenvalue) * norm(eigenvector),
            floatmin(Float64)
        )
        return numerator / denominator
    end

    function is_positive_semidefinite(matrix; factor = 256)
        T = precision_type(matrix)
        scale = max(opnorm(matrix), floatmin(T))
        tolerance = factor * eps(T) * max(size(matrix)...) * scale
        return minimum(eigvals(Hermitian(matrix))) >= -tolerance
    end
end
