@testitem "BaseParams: calc_wirearray_coords unit tests" setup = [defaults] begin

    @testset "Basic Functionality" begin
        @testset "Standard 6-wire array at origin" begin
            let num_wires = 6, radius_wire = 0.001, radius_in = 0.01
                lay_radius = radius_in + radius_wire # 0.011
                coords = calc_wirearray_coords(num_wires, radius_wire, radius_in)

                @test length(coords) == num_wires
                @test coords isa Vector{Tuple{Float64,Float64}}

                # Expected coordinates for a 6-wire array (angle step = π/3)
                expected = [
                    (lay_radius, 0.0), # Angle 0
                    (lay_radius * cos(π / 3), lay_radius * sin(π / 3)),   # Angle π/3
                    (lay_radius * cos(2π / 3), lay_radius * sin(2π / 3)), # Angle 2π/3
                    (-lay_radius, 0.0), # Angle π
                    (lay_radius * cos(4π / 3), lay_radius * sin(4π / 3)), # Angle 4π/3
                    (lay_radius * cos(5π / 3), lay_radius * sin(5π / 3)), # Angle 5π/3
                ]

                @test length(coords) == length(expected)
                for (coord, exp_coord) in zip(coords, expected)
                    @test isapprox(coord[1], exp_coord[1]; atol=TEST_TOL)
                    @test isapprox(coord[2], exp_coord[2]; atol=TEST_TOL)
                end
            end
        end

        @testset "4-wire array with non-zero center" begin
            let num_wires = 4, radius_wire = 0.002, radius_in = 0.02, C = (0.1, -0.2)
                lay_radius = radius_in + radius_wire # 0.022
                coords = calc_wirearray_coords(num_wires, radius_wire, radius_in, C)

                @test length(coords) == num_wires

                # Expected coordinates for a 4-wire array (angle step = π/2)
                expected = [
                    (C[1] + lay_radius, C[2]),           # Angle 0
                    (C[1], C[2] + lay_radius),           # Angle π/2
                    (C[1] - lay_radius, C[2]),           # Angle π
                    (C[1], C[2] - lay_radius),           # Angle 3π/2
                ]
                @test length(coords) == length(expected)
                for (coord, exp_coord) in zip(coords, expected)
                    @test isapprox(coord[1], exp_coord[1]; atol=TEST_TOL)
                    @test isapprox(coord[2], exp_coord[2]; atol=TEST_TOL)
                end
            end
        end
    end

    @testset "Edge Cases" begin
        @testset "Single wire is always at the center" begin
            # A single wire's lay radius is defined as 0.
            coords = calc_wirearray_coords(1, 0.001, 0.01)
            @test coords == [(0.0, 0.0)]

            C = (10.0, -20.0)
            coords_C = calc_wirearray_coords(1, 0.001, 0.01, C)
            @test coords_C == [C]
        end

        @testset "Zero wires returns an empty vector" begin
            coords = calc_wirearray_coords(0, 0.001, 0.01)
            @test isempty(coords)
            @test coords isa Vector
        end

        @testset "Zero radii places all wires at the center" begin
            # If lay radius is zero, all wires should be at the center C.
            num_wires = 7
            coords = calc_wirearray_coords(num_wires, 0.0, 0.0)
            @test length(coords) == num_wires
            @test all(c -> c == (0.0, 0.0), coords)

            C = (1.0, 1.0)
            coords_C = calc_wirearray_coords(num_wires, 0.0, 0.0, C)
            @test length(coords_C) == num_wires
            @test all(c -> c == C, coords_C)
        end
    end

    @testset "Type Stability and Promotion" begin
        @testset "Base case: Float64 inputs" begin
            coords = calc_wirearray_coords(6, 0.001, 0.01)
            @test coords isa Vector{Tuple{Float64,Float64}}
            @test eltype(first(coords)) == Float64
        end

        @testset "Fully promoted: All inputs are Measurement" begin
            num_wires = 3
            rw = 0.001 ± 0.0001
            ri = 0.01 ± 0.0002
            C = (0.1 ± 0.01, -0.2 ± 0.02)
            coords = calc_wirearray_coords(num_wires, rw, ri, C)

            @test coords isa Vector{Tuple{Measurement{Float64},Measurement{Float64}}}
            @test eltype(first(coords)) == Measurement{Float64}

            # Check value and uncertainty propagation for the first wire (angle=0)
            lay_radius = rw + ri
            expected_x = C[1] + lay_radius
            expected_y = C[2] # sin(0) is 0, so lay_radius term is zero

            @test coords[1][1] ≈ expected_x
            @test coords[1][2] ≈ expected_y
        end

        @testset "Mixed types: radius_wire is Measurement" begin
            num_wires = 4
            rw = 0.001 ± 0.0001
            ri = 0.01 # Float64
            C = (0.1, -0.2) # Tuple{Float64, Float64}
            coords = calc_wirearray_coords(num_wires, rw, ri, C=C)

            @test coords isa Vector{Tuple{Measurement{Float64},Measurement{Float64}}}
            lay_radius_val = Measurements.value(rw) + ri

            # Wire 1 (angle 0)
            @test Measurements.value(coords[1][1]) ≈ C[1] + lay_radius_val atol = TEST_TOL
            @test Measurements.value(coords[1][2]) ≈ C[2] atol = TEST_TOL
            @test Measurements.uncertainty(coords[1][1]) > 0
            @test Measurements.uncertainty(coords[1][2]) == 0 # sin(0) = 0, no uncertainty propagation

            # Wire 2 (angle π/2)
            @test Measurements.value(coords[2][1]) ≈ C[1] atol = TEST_TOL
            @test Measurements.value(coords[2][2]) ≈ C[2] + lay_radius_val atol = TEST_TOL
            @test isapprox(Measurements.uncertainty(coords[2][1]), 0, atol=TEST_TOL) # cos(π/2) = 0, no uncertainty propagation
            @test Measurements.uncertainty(coords[2][2]) > 0
        end

        @testset "Mixed types: radius_in is Measurement" begin
            num_wires = 4
            rw = 0.001 # Float64
            ri = 0.01 ± 0.0002
            coords = calc_wirearray_coords(num_wires, rw, ri)

            @test coords isa Vector{Tuple{Measurement{Float64},Measurement{Float64}}}
            lay_radius_uncert = Measurements.uncertainty(ri)
            @test Measurements.uncertainty(coords[1][1]) ≈ lay_radius_uncert atol = TEST_TOL
        end

        @testset "Mixed types: Center C is Measurement" begin
            num_wires = 4
            rw = 0.001 # Float64
            ri = 0.01 # Float64
            C = (0.1 ± 0.01, -0.2 ± 0.02)
            # Use keyword argument version to test the helper method
            coords = calc_wirearray_coords(num_wires, rw, ri; C=C)

            @test coords isa Vector{Tuple{Measurement{Float64},Measurement{Float64}}}
            lay_radius = rw + ri

            # Check uncertainty propagation from center C
            @test Measurements.value(coords[1][1]) ≈ Measurements.value(C[1]) + lay_radius atol = TEST_TOL
            @test Measurements.value(coords[1][2]) ≈ Measurements.value(C[2]) atol = TEST_TOL
            @test Measurements.uncertainty(coords[1][1]) ≈ Measurements.uncertainty(C[1]) atol = TEST_TOL
            @test Measurements.uncertainty(coords[1][2]) ≈ Measurements.uncertainty(C[2]) atol = TEST_TOL
        end
    end
end