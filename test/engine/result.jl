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

    @testset "a Z that is not passive warns once, with its frequencies" begin
        active = copy(Z)
        active[1, 1, 2:3] .= complex(-1.0e-4, 1.0e-3)
        @test_logs (:warn,
            "the series impedance Z is not passive at 2 of 3 frequencies, from 1000.0 Hz " *
            "to 10000.0 Hz; the results may be nonphysical") result(
            active, Y, frequencies, labels, 1.0, FEM(), trace)
        @test_logs result(
            Z, fill(complex(-1.0e-6, 1.0e-9), 2, 2, 3),
            frequencies, labels, 1.0, FEM(), trace)
    end

    @testset "the stored finite-element Z of $case is passive" for (case, _) in REFERENCES
        path = joinpath(@__DIR__, "..", "references", case, "matrices.csv")
        rows = [split(line, ',') for line in eachline(path) if startswith(line, "fem")]
        stored = sort(unique(parse(Float64, row[2]) for row in rows))
        n = maximum(parse(Int, row[4]) for row in rows)
        Zfem, Yfem = (zeros(ComplexF64, n, n, length(stored)) for _ in 1:2)
        for row in rows
            matrix = row[3] == "Z" ? Zfem : Yfem
            i, j = parse(Int, row[4]), parse(Int, row[5])
            k = searchsortedfirst(stored, parse(Float64, row[2]))
            matrix[i, j, k] = complex(parse(Float64, row[6]), parse(Float64, row[7]))
        end
        @test_logs result(
            Zfem, Yfem, stored, ["cable:$i:core" for i in 1:n], 1.0, FEM(), trace)
    end
end
