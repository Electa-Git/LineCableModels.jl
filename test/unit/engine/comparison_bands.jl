@testitem "Engine / RMS bands slice stored frequencies and retain comparison meaning" tags=[:unit] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine: compare, RMSError, LineParametersBenchmark
    f = [0.1, 1.0, 40.0, 60.0, 100.0, 1000.0, 2000.0, 3000.0, 1e6]
    z = reshape(ComplexF64.(1:length(f)), 1, 1, :)
    y = (0.1+0.2im) .* z
    a = LineParameters(PhaseDomain, z, y, f)
    b = LineParameters(PhaseDomain, z .* 1.1, y .* 1.2, f)
    original = deepcopy((a.Z.values, a.Y.values, b.Z.values, b.Y.values))
    full = compare(a, b)
    @test full isa LineParametersBenchmark
    @test only(full.Z.relative) ≈ 0.1
    @test only(full.Y.relative) ≈ 0.2
    @test full.Y.details.data.sample_count == length(f)
    for (band, indices) in
        ((:dc, 1:5), (:harmonic, 3:7), (:narrow, 6:9), ((50.0, 2500.0), 3:7))
        result = compare(a, b; band)
        @test result.Z.details.data.indices == indices
        @test result.Z.details.data.actual_bounds == (f[first(indices)], f[last(indices)])
        @test result.Z.details.data.sample_count == length(indices)
        sliced_a = LineParameters(PhaseDomain, z[:, :, indices], y[:, :, indices], f[indices])
        sliced_b = LineParameters(PhaseDomain, 1.1z[:, :, indices], 1.2y[:, :, indices], f[indices])
        @test result.Z.absolute == compare(sliced_a, sliced_b).Z.absolute
        @test result.Y.relative == compare(sliced_a, sliced_b).Y.relative
    end
    @test compare(a, b; band = :harmonic, fundamental = 60, harmonics = 50).Z.details.data.indices ==
          4:8
    @test compare(a, b; band = (0.0, 1.0)).Z.details.data.indices == 1:2
    @test compare(a, b; band = (1000.0, 2e6)).Z.details.data.indices == 6:9
    @test compare(a, b; band = (55.0, 65.0)).Z.details.data.indices == 4:4
    @test compare(a, b; band = (55.0, 65.0)).Z.details.data.sample_count == 1
    for band in (:wide, (2e6, 3e6), (0.001, 0.01))
        result = compare(a, b; band)
        @test all(ismissing, result.Z.absolute)
        @test all(ismissing, result.Y.relative)
        @test all(==(:no_samples), result.Z.details.data.status)
        @test result.Z.details.data.sample_count == 0
        @test result.Z.details.data.actual_bounds === (missing, missing)
        @test result.Z.details.data.reason isa String
    end
    wider_f = [1e6, 1.1e6, 2e6]
    wider = LineParameters(PhaseDomain, ones(ComplexF64, 1, 1, 3), ones(ComplexF64, 1, 1, 3), wider_f)
    @test compare(wider, wider; band = :wide).Z.details.data.indices == 2:3
    @test_throws ArgumentError compare(a, b; band = :unknown)
    @test_throws ArgumentError compare(a, b; band = (100.0, 50.0))
    @test_throws ArgumentError compare(a, b; fundamental = 0)
    @test_throws ArgumentError compare(a, b; harmonics = 0)
    @test_throws ArgumentError compare(a, b; atol = -1)
    @test_throws ArgumentError compare(a, b; unsupported = (Y = true,))
    unavailable = compare(a, b; unsupported = (Y = "Backend does not provide this observable",))
    @test all(ismissing, unavailable.Y.relative)
    @test only(unavailable.Y.details.data.status) === :unsupported
    @test only(unavailable.Z.relative) ≈ 0.1
    declared = LineParameters(PhaseDomain, z, y, f;
        details = ComputationDetails(;comparison_unsupported = (G = "Dielectric conduction not represented",),))
    @test only(compare(declared, b, G).details.data.status) === :unsupported
    @test only(compare(declared, b).Y.details.data.status) === :compared
    # Numerical-zero policy is local to the selected band, not the full sweep.
    tiny = fill(1e-15+1e-15im, 1, 1, length(f))
    signal = copy(tiny)
    signal[1, 1, end] = 1.0
    quiet = LineParameters(PhaseDomain, tiny, tiny, f)
    noisy = LineParameters(PhaseDomain, signal, signal, f)
    @test ismissing(only(compare(quiet, noisy; band = :dc).Y.relative))
    @test only(compare(quiet, noisy; band = :dc).Y.details.data.status) === :reference_below_tolerance
    @test ismissing(only(compare(quiet, noisy).Y.relative))
    @test only(compare(quiet, noisy).Y.absolute) > 0.1
    @test only(compare(quiet, noisy; band = :dc, atol = 0.0).Y.details.data.status) === :compared
    zero = LineParameters(
        PhaseDomain, zeros(ComplexF64, 1, 1, length(f)), zeros(ComplexF64, 1, 1, length(f)), f)
    @test ismissing(only(compare(zero, noisy).Y.relative))
    tiny_capacitance = reshape(complex.(zeros(length(f)), 2π .* f .* 1e-17), 1, 1, :)
    negligible = LineParameters(PhaseDomain, z, tiny_capacitance, f)
    @test ismissing(only(compare(zero, negligible).Y.relative))
    @test only(compare(zero, negligible).Y.absolute) > 0
    @test ismissing(only(compare(zero, negligible, C).relative))
    @test only(compare(zero, negligible; atol = (C = 1e-18, G = 0.0)).Y.details.data.status) ===
          :reference_below_tolerance
    @test only(compare(quiet, noisy; band = :dc, atol = (Y = 1e-16,)).Y.details.data.status) ===
          :compared
    for quantity in (R, L, C, G)
        result = compare(a, b, quantity; band = :harmonic)
        @test result isa RMSError
        @test result.details.data.quantity === nameof(quantity)
        @test result.details.data.indices == 3:7
    end
    @test (a.Z.values, a.Y.values, b.Z.values, b.Y.values) == original
end

@testitem "Engine / RMS normalization retains ordered operands and every sample" tags=[:unit] begin
    using LineCableModels.Engine: compare
    function parameters(values)
        tensor = reshape(ComplexF64.(values), 1, 1, :)
        LineParameters(PhaseDomain, tensor, tensor, collect(1.0:length(values)))
    end
    a, b = parameters([1, 100]), parameters([2, 100])
    @test only(compare(a, b).Y.relative) ≈ sqrt(1/10001)
    @test only(compare(a, b; normalization = :pointwise).Y.relative) ≈ sqrt(1/2)
    @test only(compare(b, a).Y.relative) ≈ sqrt(1/10004)
    @test only(compare(b, a; normalization = :pointwise).Y.relative) ≈ sqrt(0.25/2)
    @test compare(a, b).Y.absolute == compare(a, b; normalization = :pointwise).Y.absolute
    @test compare(a, b).Y.details.data.normalization === :reference_rms
    @test compare(a, b; normalization = :pointwise).Y.details.data.normalization === :pointwise
    for normalization in (:reference_rms, :pointwise)
        @test only(compare(parameters([im, 2im]), parameters([-im, -2im]); normalization).Y.relative) ≈
              2
        @test ismissing(only(compare(parameters([0, 0]), parameters([0, 0]); normalization, atol = 0).Y.relative))
        @test ismissing(only(compare(parameters([0, 0]), parameters([0, 1]); normalization).Y.relative))
        @test ismissing(only(compare(a, b; normalization, band = :wide).Y.relative))
        @test ismissing(only(compare(a, b; normalization, unsupported = (Y = "unavailable",)).Y.relative))
    end
    # An unusable denominator does not change the selected sample population.
    @test ismissing(only(compare(parameters([0, 1]), parameters([0, 2]); normalization = :pointwise).Y.relative))
    @test ismissing(only(compare(parameters([0, 1]), parameters([1, 1]); normalization = :pointwise).Y.relative))
    @test ismissing(only(compare(parameters([0, 1]), parameters([1, 1])).Y.relative))
    @test_throws ArgumentError compare(a, b; normalization = :unknown)
end

@testitem "Engine / relative RMS requires two significant operands at every selected sample" tags=[:unit] begin
    using LineCableModels.Engine: compare
    tensor(values) = reshape(ComplexF64.(values), 1, 1, :)
    signal = tensor([2.0, -3.0, 4im])
    original = copy(signal)
    tolerance = [0.1, 0.2, 0.3]
    cases = (([0, 0, 0], :below_tolerance),
        ([0.01, -0.02, 0.03im], :below_tolerance),
        ([2, -0.2, 4im], :sample_below_tolerance),
        ([2, -3, 0.1im], :sample_below_tolerance))
    for normalization in (:reference_rms, :pointwise), (values, suffix) in cases
        quiet = tensor(values)
        for (reference, candidate, role) in ((signal, quiet, :candidate), (quiet, signal, :reference))
            error = compare(reference, candidate; normalization, atol=tolerance)
            @test ismissing(only(error.relative))
            @test only(error.absolute) ≈ sqrt(sum(abs2, signal .- quiet)/3)
            @test only(error.details.data.status) === Symbol(role, :_, suffix)
            @test occursin("no samples were omitted", only(error.details.data.normalization_reason))
            @test error.details.data.sample_count == 3
            @test quiet == tensor(values)
        end
    end
    for normalization in (:reference_rms, :pointwise)
        significant = tensor([0.1001, -0.2001, 0.3001im])
        @test only(compare(significant, significant; normalization, atol=tolerance).relative) == 0
        @test only(compare(signal, -signal; normalization, atol=tolerance).relative) ≈ 2
        @test ismissing(only(compare(tensor([0]), tensor([0]); normalization).relative))
        @test only(compare(signal, signal; normalization, atol=0).relative) == 0
    end
    @test signal == original

    # A quiet candidate invalidates only bands containing its quiet samples.
    f = [1.0, 50.0, 1e3, 1e7]
    z = ones(ComplexF64, 1, 1, length(f))
    y = tensor(2π .* f .* 1e-8im)
    tiny_y = tensor(2π .* f .* 1e-17im)
    reference = LineParameters(PhaseDomain, z, y, f)
    candidate = LineParameters(PhaseDomain, z, tiny_y, f)
    @test only(compare(reference, candidate, Y).details.data.status) === :candidate_below_tolerance
    @test ismissing(only(compare(reference, candidate, C).relative))
    @test only(compare(reference, candidate, Y; atol=(C=1e-18, G=0)).relative) ≈ 1-1e-9
    mixed_y = copy(y)
    mixed_y[1, 1, 1] = tiny_y[1, 1, 1]
    mixed = LineParameters(PhaseDomain, z, mixed_y, f)
    for normalization in (:reference_rms, :pointwise)
        for band in (:all, :dc)
            error = compare(reference, mixed, Y; band, normalization)
            @test ismissing(only(error.relative))
            @test only(error.details.data.status) === :candidate_sample_below_tolerance
            @test only(error.absolute) > 0
        end
        for band in (:harmonic, :narrow, :wide)
            error = compare(reference, mixed, Y; band, normalization)
            @test only(error.relative) == 0
            @test only(error.details.data.status) === :compared
        end
    end
    @test reference.Y.values == y
    @test candidate.Y.values == tiny_y
end

@testitem "Engine / lossless conductance normalization preserves measured differences" tags=[:unit] begin
    using LineCableModels.Engine: compare
    frequencies = [50.0, 500.0]
    z = ones(ComplexF64, 2, 2, 2)
    y = fill(1e-4im, 2, 2, 2)
    y[2, 2, :] .+= 1e-6
    reference = LineParameters(PhaseDomain, z, y, frequencies)
    candidate = LineParameters(PhaseDomain, z, y .+ 1e-14, frequencies)
    conductance = compare(reference, candidate, G)
    @test ismissing(conductance.relative[1, 1])
    @test conductance.absolute[1, 1] ≈ 1e-14
    @test conductance.details.data.status[1, 1] === :reference_below_tolerance
    @test conductance.details.data.normalization_reason[1, 1] isa String
    @test conductance.relative[2, 2] ≈ 1e-8
    @test conductance.details.data.normalization_reason[2, 2] === nothing
    @test all(!ismissing, compare(reference, candidate).Y.relative)
    @test reference.Y.values == y
end
