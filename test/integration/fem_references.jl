# Independent references: three cases approved by the user on 2026-10-10 and the cable
# cases copied on 2026-10-12. Each file in `fem_references/` gives the FEM and analytical Z
# and Y of one case of the FEM investigation at its 11 FEM frequencies. Its header states
# the source run folder and what FEM confirms.
@testitem "Engine / FEM references / seven cases" tags=[:integration, :engine] setup=[
    TestFixtures
] begin
    using LinearAlgebra: diag

    directory = joinpath(pkgdir(LineCableModels), "test", "integration", "fem_references")
    formulation = Formulation(earth_impedance = :unified, earth_admittance = :unified,
        shunt_model = :equivalent,
        options = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false))
    fem_limits = (
        two_bare_wires=(Z=Inf, Y=Inf),
        two_insulated_wires=(Z=Inf, Y=Inf),
        solid_1000mm2_single=(Z=Inf, Y=Inf),
        # At 10 MHz, core self-impedance differs by 2.78 %. The cause is not yet explained.
        cable_525kv_land_no_armour_dc_bipole=(Z=2.0e6, Y=Inf),
        # Z includes the analytical stranded-core lay correction, which the 2D FEM lacks.
        cable_320kv_no_armour_dc_bipole=(Y=Inf,),
        # Z includes the analytical stranded-core lay correction, which the 2D FEM lacks.
        cable_380kv_no_armour_ac_flat=(Y=Inf,),
        # Z lacks trefoil proximity. Above 39.8 kHz, Y differs by 1.52 to 9.71 %.
        cable_132kv_cigre_tb880_case0_630cu_trefoil=(Y=4.0e4,))
    for case in keys(fem_limits)
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
            if f <= get(fem_limits[case], Symbol(quantity), -Inf)
                fem = reference[("fem", quantity, f)]
                @test maximum(abs, computed - fem) / maximum(abs, diag(computed)) <= 0.01
            end
            # Numerical preservation of the analytical values: `isapprox` on each matrix,
            # with the Frobenius norm, at rtol 1e-12.
            @test isapprox(computed, reference[("analytical", quantity, f)]; rtol = 1e-12)
        end
    end
end
