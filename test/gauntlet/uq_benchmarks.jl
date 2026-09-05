const UQ_MONTE_CARLO_TRIALS = 512

function uq_inner_formulation()
    return Formulation(
        earth_impedance = :Pollaczek1926,
        earth_admittance = :IdealGround,
        insulation_admittance = formula(:Gustavsen2013),
        options = (
            kron_reduction = false,
            reduce_bundle = false,
            ideal_transposition = false
        )
    )
end

function uq_moment_tolerances(
        ; mean_relative::Real = 5.0e-2,
        std_relative::Real = 1.0e-1,
        minimum_speedup::Real = 2.0,
        timing_samples::Integer = 3,
        timing_seconds::Real = 20.0
)
    mean_relative > 0 || throw(ArgumentError("mean tolerance must be positive"))
    std_relative > 0 || throw(ArgumentError("std tolerance must be positive"))
    minimum_speedup > 1 || throw(ArgumentError(
        "the Monte Carlo/LEP timing ratio must be greater than one",
    ))
    timing_samples > 0 || throw(ArgumentError("timing samples must be positive"))
    timing_seconds > 0 || throw(ArgumentError("timing seconds must be positive"))
    absolute = (R = 1.0e-7, L = 1.0e-9, C = 1.0e-12, G = 1.0e-15)
    moment_limits(relative) = map(absolute) do floor
        (absolute = floor, relative = Float64(relative))
    end
    regression_limits = map(absolute) do floor
        (absolute = floor * 1.0e-3, relative = 1.0e-6)
    end
    return (
        reference = (
            mean = moment_limits(mean_relative),
            std = moment_limits(std_relative)
        ),
        regression = (
            mean = regression_limits,
            std = regression_limits
        ),
        performance = (
            minimum_speedup = Float64(minimum_speedup),
            samples = Int(timing_samples),
            seconds = Float64(timing_seconds)
        )
    )
end
