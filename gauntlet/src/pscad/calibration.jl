"""Cross-check ordinary 60 Hz matrices against interpolated detailed evidence."""
function ordinary_consistency(name::Symbol, args...; kwargs...)
    ordinary_consistency(datasource(name), args...; kwargs...)
end

function ordinary_consistency(::PSCAD, raw::AbstractString)
    root = joinpath(abspath(raw), "cases")
    isdir(root) || throw(ArgumentError("raw dataset has no cases directory: $raw"))
    impedance = Float64[]
    admittance = Float64[]
    unavailable = 0
    for id in readdir(root)
        directory = joinpath(root, id)
        names = readdir(directory)
        any(endswith(lowercase(name), "_zm.out") for name in names) || continue
        files = Dict(name => joinpath(directory, name) for name in names)
        phase, _ = _phase_reference(files)
        ordinary = _ordinary(files)
        f = frequencies(phase)
        if ordinary === nothing || ordinary.Z === nothing ||
           !(first(f) < ordinary.frequency < last(f))
            unavailable += 1
            continue
        end
        z = _polar_log_interpolate(Array(Z(phase)), f, ordinary.frequency)
        y = _polar_log_interpolate(Array(Y(phase)), f, ordinary.frequency)
        push!(impedance, norm(z - ordinary.Z) / max(norm(ordinary.Z), eps()))
        push!(admittance, norm(y - ordinary.Y) / max(norm(ordinary.Y), eps()))
    end
    return (;
        method = "polar interpolation on log frequency",
        checked = length(impedance),
        unavailable,
        impedance_max = maximum(impedance; init = 0.0),
        impedance_median = isempty(impedance) ? 0.0 : median(impedance),
        admittance_max = maximum(admittance; init = 0.0),
        admittance_median = isempty(admittance) ? 0.0 : median(admittance)
    )
end
