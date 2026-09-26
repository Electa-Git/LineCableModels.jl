@testmodule TestAssertions begin
    using LinearAlgebra

    function same_public_fields(left, right)
        typeof(left) === typeof(right) || return false
        return all(propertynames(left)) do name
            isequal(getproperty(left, name), getproperty(right, name))
        end
    end

    function is_symmetric_with_scale(matrix; factor = 64)
        T = typeof(float(real(zero(eltype(matrix)))))
        scale = max(norm(matrix), floatmin(T))
        return norm(matrix - transpose(matrix)) <= factor * eps(T) * scale
    end

    all_finite(values) = all(value -> isfinite(real(value)) && isfinite(imag(value)), values)
end
