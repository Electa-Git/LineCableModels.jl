@testitem "Engine / Pollaczek and ideal earth / armoured 320 kV line" tags=[:integration, :engine] setup=[TestFixtures] begin
    using LinearAlgebra: Hermitian, eigmin, opnorm
    frequencies = 10.0 .^ range(-3, log10(5e6); length=161)
    reference = TestFixtures.known_case_problem(:cable_320kv_armoured_dc_bipole; frequencies)
    selected = Formulation(earth_impedance=:pollaczek1926, earth_admittance=:ideal)
    minima = [Inf, Inf]
    shunt_roundoff = [0.0]
    for rho in (100.0, 1000.0, 10000.0)
        problem = LineParametersProblem(reference.system; reference.temperature,
            earth_props=homogeneous(; rho, eps_r=10.0), frequencies)
        result = @inferred compute(problem, selected)
        for (index, selector) in enumerate((Z, Y))
            matrices = observe(result, selector)
            @test all(isfinite, matrices)
            for matrix in eachslice(matrices; dims=3)
                value = eigmin(Hermitian((matrix+matrix')/2))
                minima[index] = min(minima[index], value)
                if selector === Z
                    @test value >= 0
                else
                    # Ideal earth and lossless dielectrics give G=0.
                    # Real P gives imaginary Y. Numerical inversion can leave B
                    # slightly asymmetric, with tiny negative Hermitian eigenvalues.
                    # Compare with the n*eps*||Y||_2 backward-error scale. Keep Y unchanged.
                    @test all(iszero, real.(matrix))
                    roundoff = size(matrix, 1)*eps(Float64)*opnorm(matrix)
                    shunt_roundoff[1] = max(shunt_roundoff[1], -value/roundoff)
                    @test value >= -roundoff
                end
                @test matrix ≈ transpose(matrix) rtol=1e-12
            end
        end
    end
    # The same complete public computation accepts the NLT's positive abscissa.
    shifted = LineParametersProblem(reference.system; reference.temperature,
        earth_props=reference.earth_props, frequencies=[1e-3, 1.0, 1e6, 5e6], sigma=80.0)
    result = compute(shifted, selected)
    @test all(isfinite, observe(result, Z))
    @test all(isfinite, observe(result, Y))
    @info "Pollaczek/ideal Hermitian minima, 1 mHz to 5 MHz" minimum_Z=minima[1] minimum_Y=minima[2] maximum_shunt_roundoff_ratio=shunt_roundoff[1]
end
