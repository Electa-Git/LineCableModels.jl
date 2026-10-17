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

function UniversalLineModel.h_step_response_exact(g::LineCableModelsLine, len::Real,
        t::AbstractVector{<:Real}; N::Integer=4096, Tfactor::Real=1.25)
    T = Tfactor * maximum(t)
    Δω = π / T
    c = log(N^2) / T
    ω = [(2k - 1) * Δω for k in 1:N]
    points = [complex(c, value) for value in ω]
    if !all(point -> haskey(g.scans, point), points)
        line = g.problem
        measured = Base.@timed compute(LineParametersProblem(line.system, line.temperature,
            line.earth_props, ω ./ (2π), c), g.formulation; options=(output_basis=:pul,))
        scan = measured.value
        for (k, point) in enumerate(points)
            g.scans[point] = (Matrix{ComplexF64}(scan.Z.values[:, :, k]),
                Matrix{ComplexF64}(scan.Y.values[:, :, k]))
        end
        @info "ULM Laplace scan" N T wall_seconds=measured.time
    end
    return invoke(UniversalLineModel.h_step_response_exact,
        Tuple{UniversalLineModel.AbstractLineGeometry, Real, AbstractVector{<:Real}},
        g, len, t; N, Tfactor)
end
