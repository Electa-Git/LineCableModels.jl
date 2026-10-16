@testitem "TransientAnalysis / step and pulse KPIs" tags=[:unit, :engine] begin
    t = collect(0.0:0.1:1.0)
    values = repeat(reshape(Float64.(t .>= 0.3), 1, 1, :), 2, 2, 1)
    response = TransientResponse(t, values, copy(values),
        ["cable:1:core", "cable:2:core"], ComputationDetails())
    for selector in (voltage, current)
        extrema = peak(response, selector)
        @test size(extrema.value) == size(extrema.time) == (2, 2)
        @test extrema.value == ones(2, 2)
        @test extrema.time == fill(t[4], 2, 2)
        @test arrival_time(response, selector) == fill(t[4], 2, 2)
        @test arrival_time(response, selector; fraction=1) == fill(t[4], 2, 2)
        @test all(isnan, damping(response, selector))
        @test dc_offset(response, selector; expected=1) == zeros(2, 2)
        for fraction in (0, -0.1, 1.1, NaN, Inf)
            @test_throws DomainError arrival_time(response, selector; fraction)
        end
    end
    pulse = zeros(2, 2, length(t))
    pulse[1, 1, 4:5] .= -2
    pulse[2, 1, 7] = 3
    pulse[1, 2, end] = 2
    pulse_response = TransientResponse(t, pulse, copy(pulse), response.coordinates,
        ComputationDetails())
    for selector in (voltage, current)
        @test peak(pulse_response, selector).value == [2 2; 3 0]
        @test peak(pulse_response, selector).time == [t[4] t[end]; t[7] t[1]]
        @test arrival_time(pulse_response, selector)[1, 1] == t[4]
        @test arrival_time(pulse_response, selector)[2, 1] == t[7]
        @test isnan(arrival_time(pulse_response, selector)[2, 2])
        @test dc_offset(pulse_response, selector) == [0 1; 0 0]
        @test all(isnan, damping(pulse_response, selector))
    end
    single = TransientResponse([0.0], fill(2.0, 2, 1, 1), fill(3.0, 2, 1, 1),
        ["cable:1:core"], ComputationDetails())
    @test dc_offset(single, voltage) == fill(2.0, 2, 1)
    @test dc_offset(single, current) == fill(3.0, 2, 1)
    @test all(isnan, damping(single, voltage))
end

@testitem "TransientAnalysis / logarithmic decrement" tags=[:unit, :engine] begin
    t = collect(0.0:0.001:40.0)
    amplitude, decay, omega = 2.0, 0.5, 2π
    trace = amplitude .* exp.(-decay .* t) .* cos.(omega .* t)
    values = repeat(reshape(trace, 1, 1, :), 2, 2, 1)
    response = TransientResponse(t, values, 0.5values,
        ["cable:1:core", "cable:2:core"], ComputationDetails())
    for selector in (voltage, current)
        @test all(value -> isapprox(value, decay * 2π / omega; rtol=1e-3),
            damping(response, selector))
    end
end

@testitem "TransientAnalysis / response comparison" tags=[:unit, :engine] begin
    import LineCableModels.Engine: compare, RMSError

    t = collect(0.0:0.1:1.0)
    values = repeat(reshape(Float64.(t .>= 0.3), 1, 1, :), 2, 2, 1)
    coordinates = ["cable:1:core", "cable:2:core"]
    response = TransientResponse(t, values, -2values, coordinates, ComputationDetails())
    scaled = TransientResponse(t, 1.01values, -2.02values, coordinates, ComputationDetails())
    zero_response = TransientResponse(t, zero(values), zero(values), coordinates,
        ComputationDetails())
    for selector in (voltage, current)
        identical = compare(response, response, selector)
        @test identical isa RMSError
        @test identical.absolute == zeros(2, 2)
        @test identical.relative == zeros(2, 2)
        error = compare(response, scaled, selector)
        @test all(value -> isapprox(value, 0.01; rtol=1e-12), error.relative)
        reference_values = observe(response, selector)
        @test error.absolute[1, 1] ≈ 0.01sqrt(sum(abs2, reference_values[1, 1, :]) / length(t))
        @test details(error).data == (quantity=nameof(selector), normalization=:reference_rms,
            samples=length(t), time_window=(first(t), last(t)))
        @test all(ismissing, compare(zero_response, response, selector).relative)
        @test all(ismissing, compare(zero_response, zero_response, selector).relative)
        @test compare(zero_response, zero_response, selector).absolute == zeros(2, 2)
        for other_times in (t .+ 0.01, t[1:(end - 1)], reverse(t))
            other = TransientResponse(other_times, values, -2values, coordinates, ComputationDetails())
            @test_throws ArgumentError compare(response, other, selector)
        end
        reordered = TransientResponse(t, values, -2values, reverse(coordinates), ComputationDetails())
        @test_throws ArgumentError compare(response, reordered, selector)
        close = TransientResponse(t .* (1 + 5e-13), values, -2values, coordinates, ComputationDetails())
        @test compare(response, close, selector).absolute == zeros(2, 2)
    end
end
