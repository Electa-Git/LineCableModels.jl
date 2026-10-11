# The first approved independent references (user, 2026-10-10). Each file in
# `fem_references/` gives the FEM and analytical Z and Y of one case of the FEM
# investigation at its 11 FEM frequencies. Its header states the source run folder.
@testitem "Engine / FEM references / three cheap cases" tags=[:integration, :engine] setup=[
    TestFixtures
] begin
    using LinearAlgebra: diag

    directory = joinpath(pkgdir(LineCableModels), "test", "integration", "fem_references")
    formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
        shunt_model = :equivalent,
        options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    for case in (:two_bare_wires, :two_insulated_wires, :solid_1000mm2_single)
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
        frequencies = sort!(unique(f for (_, _, f) in keys(reference)))
        @test length(frequencies) == 11
        result = compute(TestFixtures.fem_reference_problem(case; frequencies), formulation)
        @test details(result).data.coordinates == labels
        for (k, f) in enumerate(frequencies), (quantity, selector) in (("Z", Z), ("Y", Y))
            computed = observe(result, selector)[:, :, k]
            # The acceptance metric of the user: the largest entry difference from FEM over
            # the largest diagonal entry, at most 1 % at every frequency.
            fem = reference[("fem", quantity, f)]
            @test maximum(abs, computed - fem) / maximum(abs, diag(computed)) <= 0.01
            # Numerical preservation of the analytical values: `isapprox` on each matrix,
            # with the Frobenius norm, at rtol 1e-12.
            @test isapprox(computed, reference[("analytical", quantity, f)]; rtol = 1e-12)
        end
    end
end
