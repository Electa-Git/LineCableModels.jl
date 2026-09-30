# Arithmetic only: reads the four fixed executions, never launches a solver.
include("check_algebra.jl")
using LinearAlgebra, Printf

function petsc_matrix(path)
    open(path) do io
        header = Int.(ntoh.(read!(io, Vector{Int32}(undef, 4))))
        header[1] == 1211216 || error("Unexpected PETSc matrix header: $header")
        _, n, m, nnz = header
        counts = Int.(ntoh.(read!(io, Vector{Int32}(undef, n))))
        cols = Int.(ntoh.(read!(io, Vector{Int32}(undef, nnz)))) .+ 1
        words = ntoh.(read!(io, Vector{UInt64}(undef, 2nnz)))
        values = copy(reinterpret(ComplexF64, words))
        eof(io) || error("Trailing PETSc matrix bytes")
        return (; n, m, rowptr=vcat(1, 1 .+ cumsum(counts)), cols, values)
    end
end

function binary_vectors(path, n)
    vectors = Vector{ComplexF64}[]
    open(path) do io
        occursin("binary", readline(io)) || error("Expected native binary .res")
        readline(io) == "1.1 1" || error("Unexpected .res version")
        readline(io) == raw"$EndResFormat" || error("Invalid .res header")
        while !eof(io)
            line = readline(io)
            isempty(line) && continue
            startswith(line, raw"$Solution") || error("Unexpected .res section: $line")
            fields = split(readline(io))
            fields[1] == "0" || error("Unexpected field system")
            push!(vectors, read!(io, Vector{ComplexF64}(undef, n)))
            isempty(readline(io)) || error("Missing binary vector newline")
            readline(io) == raw"$EndSolution" || error("Invalid vector length")
        end
    end
    length(vectors) == 4 || error("Expected b1,x1,b2,x2; got $(length(vectors)) records")
    return vectors
end

function matrix_drift(A, B)
    (A.n, A.m, A.rowptr, A.cols) == (B.n, B.m, B.rowptr, B.cols) || error("Sparse patterns differ")
    d = abs.(A.values .- B.values)
    return Dict("changed_coefficients" => count(!iszero, d),
        "maximum_absolute_change" => maximum(d),
        "max_change_over_max_coefficient" => maximum(d)/maximum(abs, A.values))
end

function residual_diagnostics(A, b, x)
    # Sparse row sums only. Exact conversion of saved binary64 values to BigFloat.
    # Complex absolute value means sqrt(real(z)^2+imag(z)^2), not |re|+|im|.
    setprecision(BigFloat, 128) do
        xb = Complex{BigFloat}.(x)
        ax = abs.(xb)
        r = Vector{ComplexF64}(undef, A.n)
        eta = BigFloat(0)
        etarow = 0
        worst_numerator = BigFloat(0)
        worst_denominator = BigFloat(0)
        zeros_denominator = 0
        for i in 1:A.n
            bi = Complex{BigFloat}(b[i])
            ri = bi
            den = abs(bi)
            for k in A.rowptr[i]:A.rowptr[i+1]-1
                a = Complex{BigFloat}(A.values[k])
                j = A.cols[k]
                ri -= a * xb[j]
                den += abs(a) * ax[j]
            end
            # 0/0 -> 0; nonzero/0 -> Inf. Neither row is silently discarded.
            if iszero(den)
                zeros_denominator += 1
                ratio = iszero(ri) ? BigFloat(0) : BigFloat(Inf)
            else
                ratio = abs(ri)/den
            end
            if ratio > eta
                eta, etarow = ratio, i
                worst_numerator, worst_denominator = abs(ri), den
            end
            r[i] = ComplexF64(ri)
        end
        return Dict("residual_l2" => norm(r), "rhs_l2" => norm(b),
            "relative_residual_l2" => norm(r)/norm(b),
            "componentwise_backward_error" => Float64(eta), "worst_row" => etarow,
            "worst_row_residual_abs" => Float64(worst_numerator),
            "worst_row_denominator" => Float64(worst_denominator),
            "zero_denominator_rows" => zeros_denominator,
            "accumulation_precision_bits" => 128)
    end
end

function table_matrix(dir, name)
    path = joinpath(dir, "results/f0001-quasi-fw-b0000/matrices/$name.tsv")
    M = zeros(Complex{BigFloat}, 2, 2)
    for line in eachline(path)
        (isempty(line) || startswith(line, "#") || startswith(line, "row")) && continue
        parts = split(replace(line, raw"\t" => '\t'), '\t')
        i, j = parse.(Int, parts[1:2])
        M[i,j] = complex(parse(BigFloat, parts[end-1]), parse(BigFloat, parts[end]))
    end
    return M
end

function changes(dir0, dir1)
    P0, P1 = table_matrix(dir0, "P"), table_matrix(dir1, "P")
    Y0, Y1 = table_matrix(dir0, "Y"), table_matrix(dir1, "Y")
    C = [-Y0[1,a]*(P1[a,b]-P0[a,b])*Y1[b,2] for a in 1:2, b in 1:2]
    delta = Y1[1,2]-Y0[1,2]
    return Dict("delta_G12" => Float64(real(delta)),
        "identity_G12" => Float64(real(sum(C))),
        "identity_real_rounding_residual" => Float64(real(sum(C)-delta)),
        "P_entry_contributions_to_delta_G12" => vec(Float64.(real(C))),
        "contribution_order" => ["P11", "P21", "P12", "P22"],
        "delta_P_real" => vec(Float64.(real(P1-P0))),
        "delta_P_imag" => vec(Float64.(imag(P1-P0))))
end

function assess()
    setprecision(BigFloat, 256)
    report = Dict{String,Any}()
    for mesh in MESHES
        dirs = [joinpath(OUT, mesh, config) for config in CONFIGURATIONS]
        all(isfile(joinpath(d, "execution.toml")) for d in dirs) || error("Execution missing for $mesh")
        originals = [petsc_matrix(joinpath(d, "file_mat_before1.m.bin")) for d in dirs]
        vectors = [binary_vectors(joinpath(d, "study.res"), originals[k].n) for (k,d) in enumerate(dirs)]
        record(joinpath(OUT, mesh, "system-identity.toml"), Dict(
            "original_matrix_binary_equal" => digest(joinpath(dirs[1], "file_mat_before1.m.bin")) == digest(joinpath(dirs[2], "file_mat_before1.m.bin")),
            "preprocessing_equal" => digest(joinpath(dirs[1], "study.pre")) == digest(joinpath(dirs[2], "study.pre")),
            "rhs1_binary_equal" => reinterpret(UInt64, vectors[1][1]) == reinterpret(UInt64, vectors[2][1]),
            "rhs2_binary_equal" => reinterpret(UInt64, vectors[1][3]) == reinterpret(UInt64, vectors[2][3]),
            "second_column_working_matrix_binary_equal" => digest(joinpath(dirs[1], "file_mat_before2.m.bin")) == digest(joinpath(dirs[2], "file_mat_before2.m.bin"))))
        for (k, config) in enumerate(CONFIGURATIONS)
            dir, A, vs = dirs[k], originals[k], vectors[k]
            item = Dict{String,Any}("complex_dofs" => A.n, "stored_coefficients" => length(A.values))
            for tag in ("after1", "before2", "after2")
                item["matrix_"*tag*"_vs_original"] = matrix_drift(A, petsc_matrix(joinpath(dir, "file_mat_"*tag*".m.bin")))
            end
            for basis in 1:2
                say("ASSESS sparse original-system residual ", mesh, "/", config, " source=", basis)
                item["source$basis"] = residual_diagnostics(A, vs[2basis-1], vs[2basis])
                say("RESULT ", mesh, "/", config, " source=", basis, " ", item["source$basis"])
            end
            P, Y = table_matrix(dir, "P"), table_matrix(dir, "Y")
            item["P_real"] = vec(Float64.(real(P)))
            item["P_imag"] = vec(Float64.(imag(P)))
            item["Y_real"] = vec(Float64.(real(Y)))
            item["Y_imag"] = vec(Float64.(imag(Y)))
            item["matrix_entry_order"] = ["11", "21", "12", "22"]
            item["G12"] = Float64(real(Y[1,2]))
            item["saved_original_comparison"] = changes(joinpath(SAVED,mesh), dir)
            report[mesh*"/"*config] = item
            record(joinpath(dir,"assessment.toml"),item)
        end
        report[mesh*"/correction"] = changes(dirs...)
    end
    for config in CONFIGURATIONS
        report["corner_discrepancy/"*config] = changes(joinpath(OUT,"baseline",config), joinpath(OUT,"metric-1.35",config))
    end
    record(joinpath(OUT,"assessment.toml"),report)
    say("COMPLETE fixed-system arithmetic assessment. No further solves scheduled.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    assess()
end
