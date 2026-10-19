@testitem "PSCAD / fit report / native fitting and DC correction" tags=[:integration, :pscad] begin
    using LineCableModels.PSCAD: read_pscad_fit_report
    fixture = joinpath(@__DIR__, "..", "..", "fixtures", "pscad", "tline")

    report = read_pscad_fit_report(joinpath(fixture, "scan_violations.log"))
    @test report.admittance.before_dc == (maximum_error_percent = 0.155, rms_error_percent = nothing)
    @test report.admittance.after_dc === nothing
    @test report.admittance.poles == 13
    @test report.propagation.before_dc == (maximum_error_percent = 0.174, rms_error_percent = 0.0508)
    @test report.propagation.after_dc === nothing
    @test report.propagation.poles == [8, 20, 20, 20, 20]
    @test report.propagation.poles_source === :last_hmode
    @test report.maximum_residue_pole_ratio == 225520.39
    @test report.reports.before_dc === nothing
    @test report.reports.after_dc === nothing

    corrected = read_pscad_fit_report(joinpath(fixture, "dc_correction.log"))
    @test corrected.propagation.before_dc == (maximum_error_percent = 0.1443, rms_error_percent = 0.0357)
    @test corrected.propagation.after_dc == (maximum_error_percent = 0.14433, rms_error_percent = 0.03572)
    @test corrected.propagation.poles == [13, 20, 20, 20, 20]
    @test corrected.propagation.poles_source === :final_fit
    @test corrected.maximum_residue_pole_ratio == 360608.5
    @test corrected.admittance.after_dc === nothing
    @test occursin("Hphase:", corrected.reports.before_dc)
    @test occursin("Final Fit:", corrected.reports.after_dc)
    @test !occursin("Final Fit:", corrected.reports.before_dc)
end

@testitem "PSCAD / fit report / native passivity scan" tags=[:integration, :pscad] begin
    using LineCableModels.PSCAD: read_pscad_fit_report
    fixture = joinpath(@__DIR__, "..", "..", "fixtures", "pscad", "tline")
    report = read_pscad_fit_report(joinpath(fixture, "scan_violations.log"))
    scan = report.scan
    @test scan.band == (lower_Hz = 0.001, upper_Hz = 1e6)
    @test scan.ranges == [
        (lower_Hz = 0.001, upper_Hz = 1.97, smallest_eigenvalue = -0.349),
        (lower_Hz = 2.62, upper_Hz = 3.19, smallest_eigenvalue = -2.33e-7),
        (lower_Hz = 3.35, upper_Hz = 12.4, smallest_eigenvalue = -4.8e-6)]
    @test scan.smallest_eigenvalue == -0.349
    @test scan.outcome == "Passivity violations detected. See help for further details"
    @test startswith(scan.report, "Checking for Passivity Violations Between")
    @test !occursin('\r', scan.report)
    @test report.enforcement === nothing
    @test !haskey(first(scan.ranges), :relative_to_fitting_band)
end

@testitem "PSCAD / fit report / native enforcement iterations and outcomes" tags=[:integration, :pscad] begin
    using LineCableModels.PSCAD: read_pscad_fit_report
    fixture = joinpath(@__DIR__, "..", "..", "fixtures", "pscad", "tline")
    failed = read_pscad_fit_report(joinpath(fixture, "enforcement_failed.log"))
    record = failed.enforcement
    @test record.outcome == "Constrained solution failed"
    @test record.initial_smallest_eigenvalue == -0.0019472408
    @test record.final_smallest_eigenvalue == -0.14763037
    @test record.smallest_eigenvalue == -2.0589101
    @test record.maximum_abs_error_percent === nothing
    @test length(record.iterations) == 7
    @test first(record.iterations).main == 1
    @test first(record.iterations).inner == 1
    @test first(first(record.iterations).frequencies) == 0.0020540541
    @test first(first(record.iterations).eigenvalues) == -0.0024847484
    @test last(record.iterations) == (main = 3, inner = 1, frequencies = Float64[], eigenvalues = Float64[])
    @test minimum(record.iterations[end - 1].eigenvalues) == record.final_smallest_eigenvalue
    @test failed.scan === nothing
    @test failed.admittance === nothing
    @test failed.propagation === nothing

    successful = read_pscad_fit_report(joinpath(fixture, "enforcement_success.log")).enforcement
    # cspell:ignore sucessfully
    @test successful.outcome == "There are no passivity violations. Passivity is sucessfully enforced"
    @test successful.initial_smallest_eigenvalue == -6.5278033e-6
    @test successful.final_smallest_eigenvalue == -2.1426952e-9
    @test successful.maximum_abs_error_percent == 0.07203260719909321
    @test length(successful.iterations) == 10
    @test last(successful.iterations).main == 4
    @test last(successful.iterations).inner == 1
    @test isempty(last(successful.iterations).eigenvalues)
end

@testitem "PSCAD / fit report / partial and absent native output" tags=[:integration, :pscad] begin
    using LineCableModels.PSCAD: read_pscad_fit_report
    fixture = joinpath(@__DIR__, "..", "..", "fixtures", "pscad", "tline")
    partial = read_pscad_fit_report(joinpath(fixture, "cable_525kv_subsea_armoured_dc_bipole.log"))
    @test partial.admittance === nothing
    @test partial.propagation === nothing
    @test partial.maximum_residue_pole_ratio === nothing
    @test partial.scan === nothing
    @test partial.enforcement === nothing
    @test occursin("Error in INVCJD. Singularity occurs", partial.reports.fitting)
    @test partial.reports.before_dc === nothing
    @test partial.reports.after_dc === nothing
    mktemp() do path, io
        close(io)
        empty = read_pscad_fit_report(path)
        @test all(isnothing, (empty.admittance, empty.propagation, empty.maximum_residue_pole_ratio,
            empty.scan, empty.enforcement))
        @test all(isnothing, empty.reports)
        @test_throws SystemError read_pscad_fit_report(path * ".missing")
    end
    for name in ("scan_violations.log", "dc_correction.log", "enforcement_failed.log", "enforcement_success.log")
        bytes = read(joinpath(fixture, name), String)
        @test count("\r\n", bytes) > 0
        @test !occursin('\n', replace(bytes, "\r\n" => ""))
    end
end
