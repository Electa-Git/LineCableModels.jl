# Manual cross-backend check of the ordered mixed-layer earth matrices.
# julia --project=. test/gauntlet/unified_earth_references.jl
using LineCableModels, Test

const E=LineCableModels.Engine

function validate_references()
    geometry=E.EarthReturnGeometry([0.0, 1.0, 2.0], [1.2, -0.9, -1.4],
        [0.01, 0.025, 0.04])
    s=2pi*1e4im
    sigma=[0.0, 0.1]
    epsilon=8.8541878128e-12 .* [1.0, 8.0]
    mu=4pi*1e-7 .* [1.0, 3.0]
    state=(jω = s, Γ = 1e-4+2e-4im, sigma, epsilon, mu,
        gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
    @testset "Mixed layers, prescribed Γ and voltage references" begin
        for reference in (:deep, :interface, 3.0, :scalar)
            println((; reference, method = :quad))
            flush(stdout)
            quad=E.unified_earth!(E.EarthReturnWorkspace(geometry), state,
                E.computation_options(E.SpectralIntegral,
                    (method = :quad, options = (rtol = 1e-9,))); reference)
            for method in (:trapz, :cim)
                println((; reference, method))
                flush(stdout)
                actual=E.unified_earth!(E.EarthReturnWorkspace(geometry), state,
                    E.computation_options(E.SpectralIntegral,
                        (; method, options = (rtol = 1e-6,))); reference)
                for (kind, floor) in ((:Ze, 1e-14), (:Pe, 1e-14), (:Ye, 1e-10))
                    for (value, expected) in zip(getproperty(actual, kind),
                            getproperty(quad, kind)),
                        component in (real, imag)

                        @test isapprox(component(value), component(expected);
                            rtol = 1e-5, atol = floor)
                    end
                end
                println((; reference, method, passed = true))
                flush(stdout)
            end
        end
    end
end

validate_references()
