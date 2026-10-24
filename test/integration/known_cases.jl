@testitem "Engine / known cases" tags=[:integration, :engine] setup=[
    TestFixtures
] begin
    using LinearAlgebra: diag
    using TOML

    directory = joinpath(pkgdir(LineCableModels), "test", "integration", "known_cases")
    status_path = joinpath(directory, "status.toml")
    status = TOML.parsefile(status_path)
    frequencies = status["frequencies_hz"]
    cases = sort!(filter(!=("frequencies_hz"), collect(keys(status))))
    formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
        shunt_model = :equivalent,
        options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    for case in cases
        labels = String[]
        reference = Dict{Tuple{String, String, Float64}, Matrix{ComplexF64}}()
        for line in eachline(joinpath(directory, "$case.csv"))
            if startswith(line, "# Conductors in row/column order: ")
                labels = [String(last(split(entry, " = "))) for entry in
                          split(chopprefix(line, "# Conductors in row/column order: "), ", ")]
            elseif !startswith(line, "#") && !startswith(line, "source,")
                source, f, quantity, row, column, real_part, imaginary_part = split(line, ',')
                n = length(labels)
                matrix = get!(() -> zeros(ComplexF64, n, n), reference,
                    (String(source), String(quantity), parse(Float64, f)))
                matrix[parse(Int, row), parse(Int, column)] =
                    complex(parse(Float64, real_part), parse(Float64, imaginary_part))
            end
        end
        result = @inferred compute(TestFixtures.known_case_problem(Symbol(case); frequencies), formulation)
        @test result.terminals == labels
        for (k, f) in enumerate(frequencies), (quantity, selector) in (("Z", Z), ("Y", Y))
            computed = observe(result, selector)[:, :, k]
            @test isapprox(computed, reference[("analytical", quantity, f)]; rtol = 1e-12)
        end
        for source in filter(==("fem"), collect(keys(status[case]))), (quantity, selector) in (("Z", Z), ("Y", Y))
            computed = observe(result, selector)
            curve = [100 * maximum(abs, reference[(source, quantity, f)] - computed[:, :, k]) /
                     maximum(abs, diag(computed[:, :, k])) for (k, f) in enumerate(frequencies)]
            recorded = status[case][source][quantity]["e_percent"]
            @test all(curve .<= recorded .+ max.(1e-4 .* recorded, 1e-6))
        end
    end
end

@testitem "Engine / known cases / boundary shunt" tags=[:integration, :engine, :slow] setup=[
    TestFixtures
] begin
    using LinearAlgebra: diag
    using TOML

    directory = joinpath(pkgdir(LineCableModels), "test", "integration", "known_cases")
    status_path = joinpath(directory, "status.toml")
    status = TOML.parsefile(status_path)
    frequencies = status["frequencies_hz"]
    cases = sort!(filter(!=("frequencies_hz"), collect(keys(status))))
    formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
        shunt_model = :boundary,
        options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    for case in cases
        labels = String[]
        reference = Dict{Tuple{String, String, Float64}, Matrix{ComplexF64}}()
        for line in eachline(joinpath(directory, "$case.csv"))
            if startswith(line, "# Conductors in row/column order: ")
                labels = [String(last(split(entry, " = "))) for entry in
                          split(chopprefix(line, "# Conductors in row/column order: "), ", ")]
            elseif !startswith(line, "#") && !startswith(line, "source,")
                source, f, quantity, row, column, real_part, imaginary_part = split(line, ',')
                n = length(labels)
                matrix = get!(() -> zeros(ComplexF64, n, n), reference,
                    (String(source), String(quantity), parse(Float64, f)))
                matrix[parse(Int, row), parse(Int, column)] =
                    complex(parse(Float64, real_part), parse(Float64, imaginary_part))
            end
        end
        if any(key -> first(key) == "analytical_boundary", keys(reference))
            result = @inferred compute(TestFixtures.known_case_problem(Symbol(case); frequencies), formulation)
            for (k, f) in enumerate(frequencies), (quantity, selector) in (("Z", Z), ("Y", Y))
                computed = observe(result, selector)[:, :, k]
                @test isapprox(computed, reference[("analytical_boundary", quantity, f)]; rtol = 1e-12)
            end
            for source in filter(==("fem_boundary"), collect(keys(status[case]))), (quantity, selector) in (("Z", Z), ("Y", Y))
                computed = observe(result, selector)
                curve = [100 * maximum(abs, reference[("fem", quantity, f)] - computed[:, :, k]) /
                         maximum(abs, diag(computed[:, :, k])) for (k, f) in enumerate(frequencies)]
                recorded = status[case][source][quantity]["e_percent"]
                @test all(curve .<= recorded .+ max.(1e-4 .* recorded, 1e-6))
            end
        end
    end
end
