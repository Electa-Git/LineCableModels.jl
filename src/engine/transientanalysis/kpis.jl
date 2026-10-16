"""
$(TYPEDSIGNATURES)

Return `(value, time)` matrices indexed by end × terminal (2 × n) for the
largest absolute sample and its first occurrence in `response`.
`selector` is `voltage` or `current`. Values use \\[V\\] or \\[A\\], respectively.
Times use \\[s\\]. Ends are ordered sending, receiving.
"""
function peak(response::TransientResponse, selector::Union{typeof(voltage), typeof(current)})
    values = observe(response, selector)
    value = similar(values, size(values, 1), size(values, 2))
    time = similar(value)
    for terminal in axes(values, 2), line_end in axes(values, 1)
        trace = @view values[line_end, terminal, :]
        index = argmax(abs.(trace))
        value[line_end, terminal] = abs(trace[index])
        time[line_end, terminal] = response.t[index]
    end
    return (; value, time)
end

"""
$(TYPEDSIGNATURES)

Return the first sampled time \\[s\\] at which the absolute signal reaches
`fraction` of its absolute peak, indexed by end × terminal (2 × n).
`selector` is `voltage` or `current`. A zero peak gives `NaN`.
Signals that settle at zero can still have an arrival time.
`fraction` defaults to 0.1 and must be in (0, 1], otherwise throw `DomainError`.
"""
function arrival_time(response::TransientResponse,
        selector::Union{typeof(voltage), typeof(current)}; fraction=0.1)
    0 < fraction <= 1 || throw(DomainError(fraction, "fraction must be in (0, 1]"))
    values = observe(response, selector)
    time = fill!(similar(values, size(values, 1), size(values, 2)), NaN)
    for terminal in axes(values, 2), line_end in axes(values, 1)
        trace = @view values[line_end, terminal, :]
        amplitude = maximum(abs, trace)
        iszero(amplitude) && continue
        index = findfirst(x -> abs(x) >= fraction * amplitude, trace)
        time[line_end, terminal] = response.t[index]
    end
    return time
end

"""
$(TYPEDSIGNATURES)

Return logarithmic decrement \\[dimensionless\\], indexed by end × terminal
(2 × n), for `selector` (`voltage` or `current`).

```math
\\delta = \\ln(a_1/a_2)
```

Here ``a_1`` and ``a_2`` are the first two strict local maxima of the signal
minus its final value, strictly after [`arrival_time`](@ref) with fraction 0.1.
The final value is the mean of the last 10 % of samples, rounded up to at least
one sample. Fewer than two maxima gives `NaN`.
"""
function damping(response::TransientResponse, selector::Union{typeof(voltage), typeof(current)})
    values = observe(response, selector)
    arrival = arrival_time(response, selector)
    decrement = fill!(similar(arrival), NaN)
    tail = max(1, cld(length(response.t), 10))
    for terminal in axes(values, 2), line_end in axes(values, 1)
        trace = @view values[line_end, terminal, :]
        final = mean(@view trace[(end - tail + 1):end])
        first_maximum = zero(eltype(values))
        count = 0
        for index in 2:(length(trace) - 1)
            response.t[index] > arrival[line_end, terminal] || continue
            if trace[index - 1] < trace[index] > trace[index + 1]
                count += 1
                if count == 1
                    first_maximum = trace[index] - final
                else
                    decrement[line_end, terminal] = log(first_maximum / (trace[index] - final))
                    break
                end
            end
        end
    end
    return decrement
end

"""
$(TYPEDSIGNATURES)

Return final value minus `expected` (default 0), indexed by end × terminal
(2 × n). `selector` is `voltage` or `current`. Both `expected` and the returned
offset use \\[V\\] or \\[A\\], respectively. The final value is the mean of the
last 10 % of samples, rounded up to at least one sample.
"""
function dc_offset(response::TransientResponse,
        selector::Union{typeof(voltage), typeof(current)}; expected=0)
    values = observe(response, selector)
    tail = max(1, cld(length(response.t), 10))
    return dropdims(mean(@view(values[:, :, (end - tail + 1):end]); dims=3); dims=3) .- expected
end

"""
$(TYPEDSIGNATURES)

Compare terminal responses over their stored time samples. `selector` is
`voltage` or `current`. The first operand supplies the reference normalization.

```math
E_{\\mathrm{abs}} = \\sqrt{\\frac{1}{N_t}\\sum_k |x_k-r_k|^2},
\\qquad E_{\\mathrm{rel}} =
\\frac{E_{\\mathrm{abs}}}{\\sqrt{\\frac{1}{N_t}\\sum_k |r_k|^2}}.
```

Return [`RMSError`](@ref) with end × terminal (2 × n) matrices. Absolute errors
use \\[V\\] or \\[A\\]. Relative errors are fractions \\[dimensionless\\] and
are `missing` where the reference RMS is zero. Details retain `quantity`,
`normalization=:reference_rms`, `samples` and `time_window` \\[s\\].
Throw `ArgumentError` unless time axes match within relative tolerance 1e-12
(zero absolute tolerance) and coordinates match in order.
"""
function compare(reference::TransientResponse, result::TransientResponse,
        selector::Union{typeof(voltage), typeof(current)})
    length(reference.t) == length(result.t) &&
        all(isapprox.(reference.t, result.t; rtol=1e-12, atol=0)) ||
        throw(ArgumentError("reference and result time axes must match"))
    reference.coordinates == result.coordinates ||
        throw(ArgumentError("reference and result terminal coordinates must match"))
    left, right = observe(reference, selector), observe(result, selector)
    size(left) == size(right) ||
        throw(ArgumentError("reference and result response dimensions must match"))
    T = promote_type(eltype(left), eltype(right))
    absolute = Matrix{T}(undef, size(left, 1), size(left, 2))
    relative = Matrix{Union{Missing, T}}(missing, size(absolute))
    samples = length(reference.t)
    for terminal in axes(left, 2), line_end in axes(left, 1)
        a = @view left[line_end, terminal, :]
        b = @view right[line_end, terminal, :]
        error = sqrt(sum(abs2, b .- a) / samples)
        reference_rms = sqrt(sum(abs2, a) / samples)
        absolute[line_end, terminal] = error
        iszero(reference_rms) || (relative[line_end, terminal] = error / reference_rms)
    end
    return RMSError{T}(absolute, relative; details=ComputationDetails(;
        quantity=nameof(selector), normalization=:reference_rms, samples,
        time_window=(first(reference.t), last(reference.t))))
end
