module Beta
import ..Commons
function uses(m)
    Commons.shared(1) + Commons.undocumented(1) + Commons.untested(1) +
    Commons.ideal_transposition!(m)[1] + Commons.hidden_wrapper([1])[1]
end
c0() = 299792458
mu(T) = 4 * (one(T) * π) * (one(T) * 10)^(-7)
function average_cyclic!(values::AbstractMatrix)
    size = checksquare(values)
    sums = similar(diag(values))
    @inbounds for shift in 0:(size - 1)
        accumulator = zero(eltype(values))
        for i in 1:size
            accumulator += values[i, 1 + mod(i - 1 + shift, size)]
        end
        sums[shift + 1] = accumulator / size
    end
    @inbounds for i in 1:size, j in 1:size

        values[i, j] = sums[mod1(j - i + 1, size)]
    end
    return values
end
end
