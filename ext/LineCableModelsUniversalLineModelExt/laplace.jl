"""
    LineCableModelsLine(problem, formulation)

Supply per-unit-length line matrices to UniversalLineModel on the Laplace plane.
Each exact step response computes its complete grid in one line-parameter scan.
"""
struct LineCableModelsLine{P <: LineParametersProblem, F <: LineParametersFormulation} <:
        UniversalLineModel.AbstractLineGeometry
    problem::P
    formulation::F
    scans::Dict{ComplexF64, Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}
end

function LineCableModelsLine(problem::LineParametersProblem,
        formulation::LineParametersFormulation)
    scans = Dict{ComplexF64, Tuple{Matrix{ComplexF64}, Matrix{ComplexF64}}}()
    return LineCableModelsLine(problem, formulation, scans)
end

function UniversalLineModel.nphases(g::LineCableModelsLine)
    if isempty(g.scans)
        point = complex(g.problem.sigma, 2π * first(g.problem.frequencies))
        return size(UniversalLineModel.zy(g, point)[1], 1)
    end
    return size(first(values(g.scans))[1], 1)
end

"""
    UniversalLineModel.zy(g::LineCableModelsLine, s::Number)

Return stored series impedance [Ω/m] and shunt admittance [S/m] at `s` [1/s].
A point without retained matrices takes the slow path: one line-parameter computation at
`frequency=imag(s)/(2π)` and `sigma=real(s)`, then retains the matrix pair.
"""
function UniversalLineModel.zy(g::LineCableModelsLine, s::Number)
    return get!(g.scans, s) do
        line = g.problem
        scan = compute(LineParametersProblem(line.system, line.temperature, line.earth_props,
            [imag(s) / (2π)], real(s)), g.formulation; options=(output_basis=:pul,))
        return (Matrix{ComplexF64}(scan.Z.values[:, :, 1]),
            Matrix{ComplexF64}(scan.Y.values[:, :, 1]))
    end
end

"""
    UniversalLineModel.zy(g::LineCableModelsLine, points::AbstractVector{<:Number})

Return matrix pairs in point order. Points share one Laplace abscissa. One
line-parameter scan evaluates the complete grid when any point is absent.

# Arguments

- `g`: physical line and its line-parameter formulation.
- `points`: nonempty Laplace grid \\[1/s\\] with positive imaginary parts.

# Returns

- Series impedance \\[Ω/m\\] and shunt admittance \\[S/m\\] pairs in point order.
"""
function UniversalLineModel.zy(g::LineCableModelsLine, points::AbstractVector{<:Number})
    isempty(points) && throw(ArgumentError("the Laplace grid must contain a point"))
    abscissa = real(first(points))
    all(point -> real(point) == abscissa, points) ||
        throw(ArgumentError("the Laplace grid must share one abscissa"))
    if !all(point -> haskey(g.scans, point), points)
        ordered = sort(points; by=imag)
        line = g.problem
        scan = compute(LineParametersProblem(line.system, line.temperature, line.earth_props,
            imag.(ordered) ./ (2π), abscissa), g.formulation; options=(output_basis=:pul,))
        for (k, point) in enumerate(ordered)
            g.scans[point] = (Matrix{ComplexF64}(scan.Z.values[:, :, k]),
                Matrix{ComplexF64}(scan.Y.values[:, :, k]))
        end
    end
    return [UniversalLineModel.zy(g, point) for point in points]
end

function wave_matrices(parameters::Tuple{<:AbstractMatrix, <:AbstractMatrix}, len::Real)
    Z, Y = parameters
    admittance = Z \ sqrt(Z * Y)
    propagation = exp(-sqrt(Y * Z) * len)
    return (admittance=(admittance + transpose(admittance)) / 2, propagation=propagation)
end

function dc_targets(g::LineCableModelsLine, len::Real)
    f = [1e-6, 1e-9, 1e-12, 1e-15, 1e-18]
    samples = UniversalLineModel.zy(g, 2π * im .* f)
    for (k, matrices) in enumerate(samples)
        all(matrix -> all(isfinite, matrix), matrices) ||
            throw(ArgumentError("DC line parameters are not finite at $(f[k]) Hz"))
    end
    pairs = ((1, 2), (2, 3), (3, 4), (4, 5))
    differences = map((1, 2)) do quantity
        return map(pairs) do (i, j)
            scale = maximum(abs, samples[j][quantity])
            difference = maximum(abs, samples[i][quantity] - samples[j][quantity])
            return iszero(scale) ? difference : difference / scale
        end
    end
    relative_differences = (impedance=differences[1], admittance=differences[2])
    accepted = findfirst(eachindex(pairs)) do k
        return all(values -> values[k] <= 1e-6, differences)
    end
    if isnothing(accepted)
        table = join(("$(f[i]) -> $(f[j]) Hz: Z=$(differences[1][k]), Y=$(differences[2][k])"
                      for (k, (i, j)) in enumerate(pairs)), "\n")
        throw(ArgumentError("DC line parameters do not converge within 1e-6 relative:\n$table"))
    end
    accepted_frequency = f[accepted + 1]
    matrices = wave_matrices(samples[accepted + 1], len)
    imaginary_relative = map(matrices) do matrix
        return maximum(abs, imag.(matrix)) / maximum(abs, real.(matrix))
    end
    all(ratio -> isfinite(ratio) && ratio < 1e-6, imaginary_relative) ||
        throw(ArgumentError("DC target imaginary parts at $(accepted_frequency) Hz exceed 1e-6 relative: $imaginary_relative"))
    targets = map(matrix -> real.(matrix), matrices)
    return targets, (; frequencies=f, pairs, accepted_frequency, relative_differences, imaginary_relative)
end

function dc_reference(g::LineCableModelsLine, data::LineFrequencyData)
    fmin = minimum(data.f)
    f = exp10.(range(log10(fmin) - 3, log10(fmin); length=31))
    f[1], f[end] = fmin / 1000, fmin
    points = 2π * im .* f
    samples = UniversalLineModel.zy(g, points)
    matrices = map(pair -> wave_matrices(pair, data.len), samples)
    return (admittance=(s=points, F=getproperty.(matrices, :admittance)),
        propagation=(s=points, F=getproperty.(matrices, :propagation)))
end

function UniversalLineModel.h_step_response_exact(g::LineCableModelsLine, len::Real,
        t::AbstractVector{<:Real}; N::Integer=4096, Tfactor::Real=1.25)
    T = Tfactor * maximum(t)
    Δω = π / T
    c = log(N^2) / T
    ω = [(2k - 1) * Δω for k in 1:N]
    points = [complex(c, value) for value in ω]
    if !all(point -> haskey(g.scans, point), points)
        measured = Base.@timed UniversalLineModel.zy(g, points)
        @info "ULM Laplace scan" N T wall_seconds=measured.time
    end
    return invoke(UniversalLineModel.h_step_response_exact,
        Tuple{UniversalLineModel.AbstractLineGeometry, Real, AbstractVector{<:Real}},
        g, len, t; N, Tfactor)
end
