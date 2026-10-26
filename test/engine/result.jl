@testset "Engine results" begin
    empty = Array{ComplexF64, 3}(undef, 0, 0, 0)
    trace = (Zin = empty, Pin = empty, Zg = empty, Pg = empty)
    Z = fill(complex(1.0e-4, 1.0e-3), 2, 2, 3)
    Y = fill(complex(0.0, 1.0e-9), 2, 2, 3)
    frequencies = [50.0, 1.0e3, 1.0e4]
    labels = ["cable:1:core", "cable:2:core"]
    result = LineParameters{Float64}

    @testset "Z and Y are square, with one slice per frequency and one label per row" begin
        @test size(result(Z, Y, frequencies, labels, 1.0, FEM(), trace).Y) == (2, 2, 3)
        @test_throws DimensionMismatch result(
            Z[:, :, 1:2], Y, frequencies, labels, 1.0, FEM(), trace)
        @test_throws DimensionMismatch result(
            Z, Y[1:1, 1:1, :], frequencies, labels, 1.0, FEM(), trace)
        @test_throws DimensionMismatch result(
            Z, Y, frequencies, labels[1:1], 1.0, FEM(), trace)
        @test_throws DimensionMismatch result(
            Z, Y, frequencies, labels, 1.0, FEM(), (Z = Z[:, :, 1:2],))
    end

    @testset "Z and Y are finite" begin
        @test_throws DomainError result(
            fill(complex(NaN), 2, 2, 3), Y, frequencies, labels, 1.0, FEM(), trace)
    end
end
