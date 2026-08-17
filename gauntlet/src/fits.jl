"""Evaluate one PSCAD vector fit at positive frequencies in hertz."""
function evaluate(fit::Fit{T}, frequency::AbstractVector{<:Real}) where {T}
    f = T.(frequency)
    all(value -> isfinite(value) && value > zero(T), f) || throw(DomainError(
        f, "fit evaluation frequencies must be finite and positive"
    ))
    n = length(fit.columns)
    count = length(f)
    yc = Array{Complex{T}, 3}(undef, n, n, count)
    propagation = zeros(Complex{T}, n, n, count)

    @inbounds for index in eachindex(f)
        s = complex(zero(T), T(2π) * f[index])
        yc[:, :, index] .= fit.constant
        for column in 1:n
            model = fit.columns[column]
            for pole in eachindex(model.poles)
                yc[:, column, index] .+= model.residues[:, pole] ./
                                         (s - model.poles[pole])
            end
        end
        for group in fit.groups
            delay = exp(-s * group.delay)
            for pole in eachindex(group.poles)
                propagation[:, :, index] .+= group.residues[:, :, pole] .* delay ./
                                             (s - group.poles[pole])
            end
        end
    end
    return (; frequency = f, characteristic = yc, propagation)
end

evaluate(fit::Fit, frequency::Real) = evaluate(fit, [frequency])
